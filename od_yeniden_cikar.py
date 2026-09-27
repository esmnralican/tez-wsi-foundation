#!/usr/bin/env python3
"""
OD (optik yogunluk) tabanli maske ile secili slaytlari yeniden cikarir.

NEDEN
-----
Dort WSI'da (tumu TSS 86) kalem izinin yuksek doygunlugu, doygunluk tabanli
Otsu esigini yukari kaydirmis ve acik boyanmis GERCEK dokunun arka plan
olarak siniflandirilmasina yol acmistir:

    wsi              mevcut maske    OD maske    oran
    86-7711             %8.2          %70.2      8.5x
    86-7953            %14.7          %60.1      4.1x
    86-7701            %19.6          %46.3      2.4x
    86-8056            %27.3          %52.3      1.9x

Normal slaytlarda iki yontem ayni sonucu veriyor (1.0x), bu nedenle
degisiklik tum kohorta degil yalnizca bu slaytlara uygulanir.

OD maskesi renge ve doygunluga bagli degildir: doku, boyasi ne kadar acik
olursa olsun camdan koyudur.

NE YAPAR
--------
1. Secili slaytlarin eski patch'lerini ve manifest satirlarini kaldirir
2. OD maskesiyle yeniden tarar, kalem maskesini uygular
3. Patch'leri cikarir (blur QC + neg_frac olcumu ile birlikte)
4. final_patch_manifest.csv ve neg_frac.csv'yi gunceller

Kullanim:
    python od_yeniden_cikar.py                      # varsayilan 4 slayt
    python od_yeniden_cikar.py --wsi TCGA-86-7711   # tek slayt
    python od_yeniden_cikar.py --dry-run            # sadece olc, yazma
"""
import argparse
import csv
import shutil
import time
from pathlib import Path

import cv2
import numpy as np
import openslide
from PIL import Image

OUT = Path("~/tez/out/LUAD").expanduser()
SLIDES = Path("~/tcga_data/LUAD").expanduser()

HEDEF = [
    "TCGA-86-7711-01Z-00-DX1",
    "TCGA-86-7953-01Z-00-DX1",
    "TCGA-86-7701-01Z-00-DX1",
    "TCGA-86-8056-01Z-00-DX1",
]

CFG = {
    "target_mpp": 0.50,
    "patch_px": 224,
    "stride_um": 112.0,
    "mask_thumb_dim": 2048,
    "tissue_min": 0.80,
    "pen_max": 0.01,
    "cap": 2000,
    "blur_min": 100.0,
    "jpeg_quality": 95,
    "seed": 42,
    "od_min": 0.15,          # [olcum] normal slaytlarda Otsu ile ayni sonuc
}

_HE = np.array([[0.650, 0.704, 0.286], [0.072, 0.990, 0.105]])
_HE = _HE / np.linalg.norm(_HE, axis=1, keepdims=True)
_P = np.linalg.pinv(_HE.T)


# ------------------------------------------------------------------ maskeler

def od_maskeleri(slide, thumb_dim, od_min):
    """
    OD tabanli doku maskesi + H&E deconvolution tabanli kalem maskesi.

    Doku: ortalama optik yogunluk esigi. Renkten ve doygunluktan bagimsiz.
    Kalem: H&E ile aciklanamayan renk (kalinti > 0.55).
    """
    W, H = slide.dimensions
    sc = max(W, H) / thumb_dim
    thumb = np.array(slide.get_thumbnail((int(W / sc), int(H / sc)))
                     .convert("RGB"))

    od3 = -np.log10((thumb.astype(np.float64) + 1.0) / 256.0)
    od = od3.mean(axis=2)

    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    tis = (od > od_min).astype(np.uint8)
    tis = cv2.morphologyEx(tis, cv2.MORPH_CLOSE, k)
    tis = cv2.morphologyEx(tis, cv2.MORPH_OPEN, k) > 0

    # kalem: patches.py ile ayni
    hsv = cv2.cvtColor(thumb, cv2.COLOR_RGB2HSV)
    sat, val = hsv[..., 1], hsv[..., 2]
    flat = od3.reshape(-1, 3)
    dense = flat.sum(axis=1) > 0.15
    resid = np.zeros(flat.shape[0], dtype=np.float32)
    if dense.sum():
        sub = flat[dense]
        c = sub @ _P.T
        resid[dense] = (np.linalg.norm(sub - c @ _HE, axis=1)
                        / (np.linalg.norm(sub, axis=1) + 1e-8)).astype(np.float32)
    resid = resid.reshape(thumb.shape[:2])

    pen = (resid > 0.55) & (sat > 60) & (val > 30)
    k3 = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    k7 = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (7, 7))
    pen = cv2.morphologyEx(pen.astype(np.uint8), cv2.MORPH_OPEN, k3)
    pen = cv2.dilate(pen, k7) > 0

    return tis, pen, sc, thumb.shape[:2]


def pencere_ort(mask, xs, ys, w, h, mw, mh):
    ii = cv2.integral(mask.astype(np.uint8))
    x0 = np.clip(xs, 0, mw - 1); y0 = np.clip(ys, 0, mh - 1)
    x1 = np.clip(xs + w, 1, mw); y1 = np.clip(ys + h, 1, mh)
    area = np.maximum((x1 - x0) * (y1 - y0), 1)
    return (ii[y1, x1] - ii[y0, x1] - ii[y1, x0] + ii[y0, x0]) / area


def neg_frac(rgb):
    """H&E ile aciklanamayan renk orani -- negatif boya konsantrasyonu."""
    od = -np.log10((rgb.astype(np.float64) + 1.0) / 256.0).reshape(-1, 3)
    d = od.sum(axis=1) > 0.15
    if d.sum() < 50:
        return 0.0
    c = od[d] @ _P.T
    return float((c.min(axis=1) < -0.05).mean())


# -------------------------------------------------------------------- ana is

def isle(wid, qc, coh, dry=False):
    hits = list(SLIDES.rglob(wid + "*.svs"))
    if not hits:
        print(f"  [!] {wid}: svs bulunamadi")
        return None

    t0 = time.time()
    slide = openslide.OpenSlide(str(hits[0]))
    base_mpp = float(qc["base_mpp_x"])
    lvl = int(qc["chosen_level"])
    rsize = int(qc["read_size_px"])
    lvl_mpp = float(qc["chosen_mpp"])
    ds0 = slide.level_downsamples[lvl]
    size0 = int(round(rsize * ds0))
    step0 = int(round(CFG["stride_um"] / base_mpp))
    W, H = slide.dimensions

    tis, pen, sc, (mh, mw) = od_maskeleri(slide, CFG["mask_thumb_dim"],
                                          CFG["od_min"])

    xs0 = np.arange(0, max(1, W - size0 + 1), step0)
    ys0 = np.arange(0, max(1, H - size0 + 1), step0)
    gx, gy = np.meshgrid(xs0, ys0)
    gx, gy = gx.ravel(), gy.ravel()

    mx = (gx / sc).astype(np.int32)
    my = (gy / sc).astype(np.int32)
    ww = max(1, int(round(size0 / sc)))

    t_ratio = pencere_ort(tis, mx, my, ww, ww, mw, mh)
    p_ratio = pencere_ort(pen, mx, my, ww, ww, mw, mh)

    uygun = (t_ratio >= CFG["tissue_min"]) & (p_ratio <= CFG["pen_max"])
    n_aday = int(uygun.sum())
    print(f"  {wid}: doku %{tis.mean()*100:.1f} | grid {len(gx):,} | "
          f"aday {n_aday:,}")

    if dry or n_aday == 0:
        slide.close()
        return {"wsi_id": wid, "n_aday": n_aday, "n_yazilan": 0}

    cx, cy, ct = gx[uygun], gy[uygun], t_ratio[uygun]
    cp = p_ratio[uygun]

    # deterministik ornekleme -- blake2b (hash() surec basina degisir)
    import hashlib
    h = hashlib.blake2b(f"{CFG['seed']}:{wid}".encode(), digest_size=8).digest()
    rng = np.random.default_rng(int.from_bytes(h, "big") % (2 ** 32))
    idx = rng.permutation(n_aday)

    wdir = OUT / "patches" / coh["patient"] / wid
    if wdir.exists():
        shutil.rmtree(wdir)
    wdir.mkdir(parents=True, exist_ok=True)

    px = CFG["patch_px"]
    satirlar, negler = [], []
    n_blur = 0
    for i in idx:
        if len(satirlar) >= CFG["cap"]:
            break
        x0, y0 = int(cx[i]), int(cy[i])
        try:
            img = slide.read_region((x0, y0), lvl, (rsize, rsize))
        except Exception:
            continue
        arr = np.array(img.convert("RGB"))
        if rsize != px:
            interp = cv2.INTER_AREA if rsize > px else cv2.INTER_LINEAR
            arr = cv2.resize(arr, (px, px), interpolation=interp)

        gray = cv2.cvtColor(arr, cv2.COLOR_RGB2GRAY)
        blur = float(cv2.Laplacian(gray, cv2.CV_64F).var())
        if blur < CFG["blur_min"]:
            n_blur += 1
            continue

        nf = neg_frac(arr)
        fn = f"{wid}_x{x0}_y{y0}.jpg"
        Image.fromarray(arr).save(wdir / fn, quality=CFG["jpeg_quality"])

        satirlar.append({
            "wsi_id": wid, "patient": coh["patient"], "tss": coh["tss"],
            "label": coh["label"], "x0": x0, "y0": y0, "size_level0": size0,
            "read_level": lvl, "read_size_px": rsize,
            "base_mpp": round(base_mpp, 5), "level_mpp": round(lvl_mpp, 4),
            "target_mpp": CFG["target_mpp"],
            "physical_um": round(px * CFG["target_mpp"], 1),
            "tissue_thumb": round(float(ct[i]), 3), "texture_thumb": "",
            "pen_thumb": round(float(cp[i]), 3), "patch_file": fn,
            "tissue_patch": "", "ink_patch": "", "blur_var": round(blur, 1),
            "he_residual": "", "hema_frac": "", "written": 1,
        })
        negler.append({"patch_file": fn, "neg_frac": round(nf, 4)})

    slide.close()
    print(f"     -> {len(satirlar):,} patch yazildi "
          f"(blur eledi {n_blur}) {time.time()-t0:.0f}s")
    return {"wsi_id": wid, "n_aday": n_aday, "n_yazilan": len(satirlar),
            "satirlar": satirlar, "negler": negler}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wsi", action="append", default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    hedef = args.wsi or HEDEF

    qcs = {r["wsi_id"]: r for r in csv.DictReader(open(OUT / "wsi_qc.csv"))}
    cohs = {r["file_name"].split(".")[0]: r
            for r in csv.DictReader(open(OUT / "cohort.csv"))}

    print(f"OD maskesi ile yeniden cikarma | {len(hedef)} WSI")
    print(f"esikler: od>{CFG['od_min']} tissue>={CFG['tissue_min']} "
          f"pen<={CFG['pen_max']} blur>={CFG['blur_min']} cap={CFG['cap']}\n")

    sonuc = []
    for wid in hedef:
        if wid not in qcs or wid not in cohs:
            print(f"  [!] {wid}: kayitlarda yok")
            continue
        r = isle(wid, qcs[wid], cohs[wid], args.dry_run)
        if r:
            sonuc.append(r)

    if args.dry_run:
        print("\n(dry-run: hicbir sey yazilmadi)")
        return

    # --- manifestleri guncelle ---
    yeni_ids = {r["wsi_id"] for r in sonuc if r.get("satirlar")}
    if not yeni_ids:
        print("\nyazilacak satir yok")
        return

    man = OUT / "final_patch_manifest.csv"
    eski = list(csv.DictReader(open(man)))
    flds = list(eski[0].keys())
    tutulan = [r for r in eski if r["wsi_id"] not in yeni_ids]
    eklenen = [s for r in sonuc for s in r.get("satirlar", [])]

    shutil.copy2(man, str(man) + ".bak")
    with open(man, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=flds, extrasaction="ignore")
        w.writeheader(); w.writerows(tutulan); w.writerows(eklenen)

    nf = OUT / "neg_frac.csv"
    eski_nf = [r for r in csv.DictReader(open(nf))]
    eski_files = {s["patch_file"] for r in sonuc for s in r.get("satirlar", [])}
    silinecek = {r["patch_file"] for r in eski if r["wsi_id"] in yeni_ids}
    tut_nf = [r for r in eski_nf if r["patch_file"] not in silinecek]
    ek_nf = [s for r in sonuc for s in r.get("negler", [])]

    shutil.copy2(nf, str(nf) + ".bak")
    with open(nf, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["patch_file", "neg_frac"])
        w.writeheader(); w.writerows(tut_nf); w.writerows(ek_nf)

    print(f"\n{'wsi_id':<34}{'eski':>8}{'yeni':>8}{'kazanc':>9}")
    eski_say = {}
    for r in eski:
        eski_say[r["wsi_id"]] = eski_say.get(r["wsi_id"], 0) + 1
    for r in sonuc:
        e = eski_say.get(r["wsi_id"], 0); y = r["n_yazilan"]
        print(f"{r['wsi_id']:<34}{e:>8}{y:>8}{y-e:>+9}")

    print(f"\nmanifest: {len(tutulan)+len(eklenen):,} patch "
          f"(onceki {len(eski):,})")
    print(f"yedekler: {man}.bak  {nf}.bak")
    print("\nsiradaki adim: final_patch_manifest_clean.csv'yi yeniden uret")


if __name__ == "__main__":
    main()