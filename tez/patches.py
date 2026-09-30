"""
Patch aday tarama ve manifest uretimi.

UC ASAMA
--------
    scan      tavansiz tarama  -> candidates/<wsi_id>.csv + scan_summary.csv
    manifest  esik + tavan uygula -> canonical_patch_manifest.csv
    extract   goruntu yaz      -> patches/<case>/<wsi_id>/*.jpg

NEDEN UC ASAMA
--------------
scan hicbir esik/tavan uygulamaz, sadece OLCER. Onceki denemede tissue
threshold karsilastirmasi max_patches tavani yuzunden gecersiz olmustu
(469 x 100 = 46.900 tavana dayanmis, esik hicbir sey degistirmemisti).
Once dagilimi gorup sonra karar veriyoruz.

manifest, uc modelin (DINOv2/UNI/Virchow) ORTAK koordinat kaynagi.
Model degisirse manifest aynen kalir; patch parametresi degisirse
yeniden uretilir.

DIKKAT EDILEN NOKTALAR (docs/patch_kararlari.md)
------------------------------------------------
  1.1  fiziksel olcek (0.50 um/px), piksel boyutu degil
  1.2  level secimi asagi yonlu, yukari olcekleme yok
  1.3  blur QC HER ZAMAN resize sonrasi
  1.4  murekkep != hematoksilen -- renk tek basina kullanilmaz
  1.5  doku orani hem thumbnail hem patch uzerinde olculur
  1.6  ornekleme deterministik (seed + wsi_id)
  1.10 tss kolonu manifest'te tutulur (site-aware analizler icin)
"""
import csv
import hashlib
import json
import time
from collections import Counter
from pathlib import Path

import numpy as np

try:
    import cv2
except ImportError:
    cv2 = None
try:
    import openslide
except ImportError:
    openslide = None
try:
    from PIL import Image
except ImportError:
    Image = None


SCAN_FIELDS = [
    "wsi_id", "patient", "tss", "label",
    "base_mpp", "read_level", "read_size_px", "step_level0", "size_level0",
    "thumb_w", "thumb_h", "mask_scale",
    "tissue_area_mm2", "n_grid",
    "n_ge_30", "n_ge_50", "n_ge_60", "n_ge_70", "n_ge_80", "n_ge_90",
    "n_tex_ge_50", "n_both_70", "n_both_80",
    "pen_frac_slide", "pen_px_thumb", "n_pen_hit",
    "status", "note", "secs",
]

CAND_FIELDS = ["x0", "y0", "tissue", "texture", "pen"]

# scan asamasinda kaydedilen minimum doku orani.
# ESIK DEGIL -- sadece dosya boyutunu sinirlamak icin. Asil esik
# manifest asamasinda, dagilima bakilarak secilir.  [tahmin]
RECORD_MIN_TISSUE = 0.30


# --- Ruifrok & Johnston H&E boya vektorleri (optical density uzayinda) ---
_HE = np.array([[0.650, 0.704, 0.286],     # hematoksilen (cekirdek)
                [0.072, 0.990, 0.105]])    # eozin (sitoplazma/stroma)
_HE = _HE / np.linalg.norm(_HE, axis=1, keepdims=True)
_HE_PINV = np.linalg.pinv(_HE.T)           # (2,3)


# ------------------------------------------------------------- doku maskesi

def tissue_and_texture_mask(slide, thumb_dim=2048):
    """
    Thumbnail uzerinde iki maske:

      tissue   doygunluk tabanli (Otsu) -- beyaz cam / bos alan elenir
      texture  yerel standart sapma     -- DUZ renkli bolgeler elenir

    texture maskesi, hata 3.7'ye karsi: uniform turkuaz background
    doygunluk maskesini geciyordu ama hicbir yapi icermiyordu.
    Duz bolgede yerel std ~0'dir.

    Donus: (tissue_bool, texture_bool, pen_bool, scale, thumb_shape)
    """
    W, H = slide.dimensions
    scale = max(W, H) / thumb_dim
    size = (max(1, int(W / scale)), max(1, int(H / scale)))
    thumb = np.array(slide.get_thumbnail(size).convert("RGB"))

    hsv = cv2.cvtColor(thumb, cv2.COLOR_RGB2HSV)
    sat, val = hsv[..., 1], hsv[..., 2]

    _, tis = cv2.threshold(cv2.GaussianBlur(sat, (5, 5), 0), 0, 255,
                           cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    tis[val < 30] = 0          # cok koyu: coverslip kenari / artefakt
    tis[val > 235] = 0         # cok parlak: bos cam
    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    tis = cv2.morphologyEx(tis, cv2.MORPH_CLOSE, k)
    tis = cv2.morphologyEx(tis, cv2.MORPH_OPEN, k)

    # yerel std (7x7): E[x^2] - E[x]^2
    gray = cv2.cvtColor(thumb, cv2.COLOR_RGB2GRAY).astype(np.float32)
    mean = cv2.blur(gray, (7, 7))
    sq = cv2.blur(gray * gray, (7, 7))
    std = np.sqrt(np.maximum(sq - mean * mean, 0))
    tex = std > 4.0            # [tahmin] -- scan ciktisiyla kalibre edilecek

    # --- kalem izi tespiti: H&E DECONVOLUTION KALINTISI ---
    # Renk tonuna gore ayirmak calismiyor: hematoksilen (mor-mavi) ile mavi
    # kalem ayni hue araliginda. Bunun yerine "bu rengi H&E boyalariyla
    # aciklayabiliyor muyum?" sorusu soruluyor. Kalem H&E disi bir
    # pigmenttir -> kalinti yuksek cikar, boyama yogunlugundan bagimsiz.
    od = -np.log10((thumb.astype(np.float64) + 1.0) / 256.0)
    flat = od.reshape(-1, 3)
    dense = flat.sum(axis=1) > 0.15
    resid_map = np.zeros(flat.shape[0], dtype=np.float32)
    if dense.sum() > 0:
        sub = flat[dense]
        coeffs = sub @ _HE_PINV.T
        recon = coeffs @ _HE
        rr = np.linalg.norm(sub - recon, axis=1)
        dd = np.linalg.norm(sub, axis=1) + 1e-8
        resid_map[dense] = (rr / dd).astype(np.float32)
    resid_map = resid_map.reshape(thumb.shape[:2])

    # kalem: H&E ile aciklanamayan + doygun + yeterince koyu
        # [olcum] esik 0.55 -- 3 kalemli / 4 kalemsiz slayt karsilastirmasiyla
    # secildi. 0.55'te kalemsiz slaytlarda yakalanan doku orani <%0.2,
    # kalemli slaytlarda %5.6-63. 0.35'te iki grup ust uste biniyordu.
    pen = (resid_map > 0.55) & (sat > 60) & (val > 30)
    pen = cv2.morphologyEx(pen.astype(np.uint8), cv2.MORPH_OPEN,
                           cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3)))
    pen = cv2.dilate(pen, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (7, 7)))

    return tis > 0, tex, pen > 0, scale, thumb.shape[:2]


def _window_means(mask, xs, ys, wx, wy, mw, mh):
    """Integral image ile toplu pencere ortalamasi (hizli)."""
    ii = cv2.integral(mask.astype(np.uint8))       # (h+1, w+1)
    x0 = np.clip(xs, 0, mw - 1)
    y0 = np.clip(ys, 0, mh - 1)
    x1 = np.clip(xs + wx, 1, mw)
    y1 = np.clip(ys + wy, 1, mh)
    area = np.maximum((x1 - x0) * (y1 - y0), 1)
    s = (ii[y1, x1] - ii[y0, x1] - ii[y1, x0] + ii[y0, x0])
    return s / area


# -------------------------------------------------------------- ASAMA: scan

def scan_slide(svs_path, row, qc_row, cfg):
    """
    Bir WSI'i tarar. HICBIR esik/tavan uygulamaz.
    Donus: (summary_dict, candidate_rows)
    """
    rec = {f: "" for f in SCAN_FIELDS}
    rec.update(wsi_id=qc_row["wsi_id"], patient=row.get("patient", ""),
               tss=row.get("tss", ""), label=row.get("label", ""),
               status="ok", n_grid=0)

    t0 = time.time()
    try:
        slide = openslide.OpenSlide(str(svs_path))
    except Exception as e:
        rec.update(status="error", note=f"acilamadi: {e}")
        return rec, []

    base_mpp = float(qc_row["base_mpp_x"])
    level = int(qc_row["chosen_level"])
    read_size = int(qc_row["read_size_px"])
    ds0 = slide.level_downsamples[level]

    size0 = int(round(read_size * ds0))
    step0 = int(round(cfg["stride_um"] / base_mpp))
    W, H = slide.dimensions

    rec.update(base_mpp=round(base_mpp, 5), read_level=level,
               read_size_px=read_size, step_level0=step0, size_level0=size0)

    try:
        tis, tex, pen, mscale, (mh, mw) = tissue_and_texture_mask(
            slide, cfg["mask_thumb_dim"])
    except Exception as e:
        slide.close()
        rec.update(status="error", note=f"maske: {e}")
        return rec, []

    rec.update(thumb_w=mw, thumb_h=mh, mask_scale=round(mscale, 2))

    # slayt genelinde kalem izi: dokunun yuzde kaci isaretli
    tis_px = max(int(tis.sum()), 1)
    rec["pen_frac_slide"] = round(float((pen & tis).sum() / tis_px), 4)
    rec["pen_px_thumb"] = int(pen.sum())

    # thumbnail pikseli basina fiziksel alan -> toplam doku alani
    um_per_thumb_px = base_mpp * mscale
    mm2 = tis.sum() * (um_per_thumb_px ** 2) / 1e6
    rec["tissue_area_mm2"] = round(float(mm2), 1)

    # grid
    xs0 = np.arange(0, max(1, W - size0 + 1), step0)
    ys0 = np.arange(0, max(1, H - size0 + 1), step0)
    gx, gy = np.meshgrid(xs0, ys0)
    gx, gy = gx.ravel(), gy.ravel()
    rec["n_grid"] = int(gx.size)

    if gx.size == 0:
        slide.close()
        rec.update(status="empty", note="grid bos", secs=round(time.time()-t0, 1))
        return rec, []

    # thumbnail koordinatlarina donustur
    mx = (gx / mscale).astype(np.int32)
    my = (gy / mscale).astype(np.int32)
    wx = max(1, int(round(size0 / mscale)))
    wy = wx

    t_ratio = _window_means(tis, mx, my, wx, wy, mw, mh)
    x_ratio = _window_means(tex, mx, my, wx, wy, mw, mh)
    p_ratio = _window_means(pen, mx, my, wx, wy, mw, mh)

    for thr in (30, 50, 60, 70, 80, 90):
        rec[f"n_ge_{thr}"] = int((t_ratio >= thr / 100).sum())
    rec["n_tex_ge_50"] = int((x_ratio >= 0.50).sum())
    rec["n_both_70"] = int(((t_ratio >= 0.70) & (x_ratio >= 0.50)).sum())
    rec["n_both_80"] = int(((t_ratio >= 0.80) & (x_ratio >= 0.50)).sum())
    rec["n_pen_hit"] = int(((t_ratio >= 0.80) & (p_ratio > 0.01)).sum())

    keep = t_ratio >= RECORD_MIN_TISSUE
    cands = [
        {"x0": int(a), "y0": int(b), "tissue": round(float(c), 3),
         "texture": round(float(d), 3), "pen": round(float(e), 3)}
        for a, b, c, d, e in zip(gx[keep], gy[keep], t_ratio[keep],
                                 x_ratio[keep], p_ratio[keep])
    ]

    slide.close()
    rec["secs"] = round(time.time() - t0, 1)
    return rec, cands


def scan(cohort_csv, qc_csv, slides_root, out_dir, cfg, limit=0):
    """Tum kohortu tarar. Yeniden baslatilabilir."""
    rows = {r["file_name"].split(".")[0]: r
            for r in csv.DictReader(open(Path(cohort_csv).expanduser()))}
    qcs = list(csv.DictReader(open(Path(qc_csv).expanduser())))
    if limit:
        qcs = qcs[:limit]

    root = Path(slides_root).expanduser().resolve()
    out = Path(out_dir).expanduser()
    cand_dir = out / "candidates"
    cand_dir.mkdir(parents=True, exist_ok=True)

    sum_path = out / "scan_summary.csv"
    done = set()
    if sum_path.exists():
        done = {r["wsi_id"] for r in csv.DictReader(open(sum_path))}
        print(f"  {len(done)} WSI zaten taranmis, atlaniyor")

    fov = cfg["patch_px"] * cfg["target_mpp"]
    print(f"{len(qcs)} WSI taranacak | FOV {fov:.0f} um | "
          f"stride {cfg['stride_um']:.0f} um")
    print(f"TAVAN YOK, ESIK YOK -- sadece olcum. "
          f"tissue >= {RECORD_MIN_TISSUE} olanlar kaydediliyor.\n")

    results = []
    with open(sum_path, "a", newline="") as sf:
        sw = csv.DictWriter(sf, fieldnames=SCAN_FIELDS, extrasaction="ignore")
        if not done:
            sw.writeheader()

        for i, qc in enumerate(qcs, 1):
            wid = qc["wsi_id"]
            if wid in done:
                continue
            if qc["status"] != "ok":
                continue
            row = rows.get(wid, {})
            rel = row.get("rel_path") or (wid + ".svs")
            p = root / rel
            if not p.exists():
                hits = list(root.rglob(wid + "*.svs"))
                p = hits[0] if hits else p

            rec, cands = scan_slide(p, row, qc, cfg)
            sw.writerow(rec)
            sf.flush()
            results.append(rec)

            if cands:
                with open(cand_dir / f"{wid}.csv", "w", newline="") as cf:
                    cw = csv.DictWriter(cf, fieldnames=CAND_FIELDS)
                    cw.writeheader()
                    cw.writerows(cands)

            if i % 20 == 0 or i == len(qcs) or rec["status"] != "ok":
                print(f"[{i}/{len(qcs)}] {wid:32s} {rec['status']:7s} "
                      f"grid={rec['n_grid']:6} >=70:{rec['n_ge_70']:6} "
                      f">=80:{rec['n_ge_80']:6} doku={rec['tissue_area_mm2']:>7} mm2 "
                      f"{rec['secs']:5}s {rec['note']}")

    return results, sum_path


# ---------------------------------------------------------------- raporlama

def scan_report(results):
    ok = [r for r in results if r["status"] == "ok"]
    print(f"\n{'='*74}\nTARAMA OZETI\n{'='*74}")
    print(f"  toplam: {len(results)}")
    for s, n in Counter(r["status"] for r in results).most_common():
        print(f"      {n:>4}  {s}")
    if not ok:
        return

    def col(k):
        return np.array([float(r[k]) for r in ok])

    print(f"\n  doku alani (mm2):")
    a = col("tissue_area_mm2")
    for lbl, v in zip(("min", "p25", "ortanca", "p75", "max"),
                      np.percentile(a, [0, 25, 50, 75, 100])):
        print(f"      {lbl:>8}: {v:>8.0f}")

    print(f"\n  ESIGE GORE SLAYT BASINA ADAY PATCH SAYISI")
    print(f"  {'esik':<12}{'min':>8}{'p10':>8}{'ortanca':>9}{'p90':>8}{'max':>9}"
          f"{'toplam':>11}")
    for thr in (30, 50, 60, 70, 80, 90):
        v = col(f"n_ge_{thr}")
        q = np.percentile(v, [0, 10, 50, 90, 100])
        print(f"  tissue>={thr/100:<5.2f}{q[0]:>8.0f}{q[1]:>8.0f}{q[2]:>9.0f}"
              f"{q[3]:>8.0f}{q[4]:>9.0f}{v.sum():>11,.0f}")

    for name, key in (("+texture>=0.5", "n_both_70"), ("+texture>=0.5", "n_both_80")):
        v = col(key)
        q = np.percentile(v, [0, 10, 50, 90, 100])
        thr = key.split("_")[-1]
        print(f"  t>={int(thr)/100:.2f}{name:<6}{q[0]:>8.0f}{q[1]:>8.0f}{q[2]:>9.0f}"
              f"{q[3]:>8.0f}{q[4]:>9.0f}{v.sum():>11,.0f}")

    print(f"\n  TAVAN SECIMI: kac slayt hedefe ulasamiyor?")
    print(f"  {'esik':<12}{'<100':>8}{'<300':>8}{'<500':>8}{'<1000':>8}")
    for thr in (50, 70, 80):
        v = col(f"n_ge_{thr}")
        print(f"  tissue>={thr/100:<5.2f}"
              f"{(v < 100).sum():>8}{(v < 300).sum():>8}"
              f"{(v < 500).sum():>8}{(v < 1000).sum():>8}")

    print(f"\n  doku kapsami (patch basina 112x112 um = 0.0125 mm2):")
    med_area = np.median(a)
    for n in (100, 300, 500, 1000):
        cover = n * 0.0125 / med_area * 100
        print(f"      {n:>5} patch -> ortanca slaytin %{cover:.1f}'i")

    print(f"\n  KALEM IZI TARAMASI (H&E deconvolution kalintisi)")
    pf = col("pen_frac_slide")
    ph = col("n_pen_hit")
    print(f"      {(pf < 0.001).sum():>4} WSI  hic yok")
    print(f"      {((pf >= 0.001) & (pf < 0.01)).sum():>4} WSI  eser (<%1)")
    print(f"      {(pf >= 0.01).sum():>4} WSI  belirgin (>%1)")
    print(f"\n      kalem izli doku orani: ortanca %{np.median(pf)*100:.2f}, "
          f"p90 %{np.percentile(pf,90)*100:.2f}, max %{pf.max()*100:.1f}")
    print(f"      kalem degen aday patch: toplam {ph.sum():,.0f} "
          f"(ortanca {np.median(ph):.0f}/WSI)")

    worst = sorted(ok, key=lambda r: -float(r["pen_frac_slide"]))[:12]
    print(f"\n      en cok kalem izli 12 WSI:")
    for r in worst:
        print(f"        {r['wsi_id']:34s} tss={r['tss']:<4} "
              f"%{float(r['pen_frac_slide'])*100:5.2f}  "
              f"({r['n_pen_hit']} aday etkileniyor)")

    tssp = {}
    for r in ok:
        tssp.setdefault(r["tss"], []).append(float(r["pen_frac_slide"]))
    byt = sorted(((t, float(np.mean(v)), len(v)) for t, v in tssp.items()),
                 key=lambda x: -x[1])[:10]
    print(f"\n      merkez bazinda ortalama kalem orani (en yuksek 10):")
    for t, v, n in byt:
        print(f"        tss {t:<4} %{v*100:5.2f}  ({n} WSI)")

    secs = col("secs")
    print(f"\n  tarama suresi: ortanca {np.median(secs):.1f}s, "
          f"toplam {secs.sum()/60:.0f} dk")


# ---------------------------------------------------------- ASAMA: manifest

MANIFEST_FIELDS = [
    "wsi_id", "patient", "tss", "label",
    "x0", "y0", "size_level0", "read_level", "read_size_px",
    "base_mpp", "level_mpp", "target_mpp", "physical_um",
    "tissue_thumb", "texture_thumb", "pen_thumb", "patch_file",
]


def build_manifest(out_dir, cohort_csv, qc_csv, cfg,
                   tissue_min=0.80, cap=2000, oversample=1.5, pen_max=0.01):
    """
    candidates/ dosyalarindan esik + tavan uygulayarak manifest uretir.

    OVERSAMPLE: extract asamasinda blur/ink QC bazi patch'leri eleyecek.
    Hedef sayiya ulasabilmek icin cap * oversample kadar aday yazilir.
    Extract, cap'e ulasinca durur.

    DETERMINISTIK: seed, wsi_id ile birlestirilir. Ayni WSI her zaman ayni
    patch'leri verir -- islem sirasi degisse bile.

    Bu dosya UC MODELIN ORTAK koordinat kaynagidir.
    """
    out = Path(out_dir).expanduser()
    cand_dir = out / "candidates"

    coh = {r["file_name"].split(".")[0]: r
           for r in csv.DictReader(open(Path(cohort_csv).expanduser()))}
    qcs = {r["wsi_id"]: r
           for r in csv.DictReader(open(Path(qc_csv).expanduser()))}

    man_path = out / "canonical_patch_manifest.csv"
    n_total, n_short, per_slide, n_pen_drop = 0, [], [], 0

    with open(man_path, "w", newline="") as mf:
        mw = csv.DictWriter(mf, fieldnames=MANIFEST_FIELDS)
        mw.writeheader()

        for wid in sorted(qcs):
            qc = qcs[wid]
            if qc["status"] != "ok":
                continue
            cf = cand_dir / f"{wid}.csv"
            if not cf.exists():
                continue

            all_c = list(csv.DictReader(open(cf)))
            ok_t = [r for r in all_c if float(r["tissue"]) >= tissue_min]
            cands = [r for r in ok_t if float(r.get("pen", 0) or 0) <= pen_max]
            n_pen_drop += len(ok_t) - len(cands)
            if not cands:
                n_short.append((wid, 0))
                continue

            # deterministik karistirma
            # NOT: Python'da str hash'i surec basina rastgele tuzlanir
            # (PYTHONHASHSEED). hash() kullanmak "sabit seed ile
            # deterministik ornekleme" iddiasini gecersiz kilar.
            # blake2b surecten bagimsiz, her zaman ayni sonucu verir.
            h = hashlib.blake2b(f"{cfg['seed']}:{wid}".encode(),
                                digest_size=8).digest()
            rng = np.random.default_rng(int.from_bytes(h, "big") % (2 ** 32))
            idx = rng.permutation(len(cands))
            want = int(cap * oversample)
            picked = [cands[i] for i in idx[:want]]

            row = coh.get(wid, {})
            base_mpp = float(qc["base_mpp_x"])
            lvl = int(qc["chosen_level"])
            rsize = int(qc["read_size_px"])
            lvl_mpp = float(qc["chosen_mpp"])
            size0 = int(round(rsize * (lvl_mpp / base_mpp)))

            for c in picked:
                mw.writerow({
                    "wsi_id": wid,
                    "patient": row.get("patient", ""),
                    "tss": row.get("tss", ""),
                    "label": row.get("label", ""),
                    "x0": c["x0"], "y0": c["y0"], "size_level0": size0,
                    "read_level": lvl, "read_size_px": rsize,
                    "base_mpp": round(base_mpp, 5),
                    "level_mpp": round(lvl_mpp, 4),
                    "target_mpp": cfg["target_mpp"],
                    "physical_um": round(cfg["patch_px"] * cfg["target_mpp"], 1),
                    "tissue_thumb": c["tissue"], "texture_thumb": c["texture"],
                    "pen_thumb": c.get("pen", 0),
                    "patch_file": f"{wid}_x{c['x0']}_y{c['y0']}.jpg",
                })

            avail = len(cands)
            per_slide.append(min(avail, cap))
            n_total += len(picked)
            if avail < cap:
                n_short.append((wid, avail))

    print(f"\n{'='*74}\nMANIFEST\n{'='*74}")
    print(f"  tissue_min : {tissue_min}   cap: {cap}   oversample: {oversample}x")
    print(f"  pen_max    : {pen_max}  -> kalem degen {n_pen_drop:,} aday dislandi")
    print(f"  WSI        : {len(per_slide)}")
    print(f"  aday satir : {n_total:,}  (oversample dahil)")
    if per_slide:
        p = np.array(per_slide)
        print(f"  hedeflenen patch: {p.sum():,}  "
              f"(ortanca {np.median(p):.0f}/WSI)")
    if n_short:
        print(f"\n  [!] {len(n_short)} WSI cap'e ulasamiyor:")
        for wid, n in sorted(n_short, key=lambda x: x[1])[:10]:
            print(f"      {wid}  {n} aday")
    print(f"\n  manifest: {man_path}")
    return man_path


# ----------------------------------------------------------- ASAMA: extract

EXTRACT_FIELDS = MANIFEST_FIELDS + [
    "tissue_patch", "ink_patch", "blur_var",
    "he_residual", "hema_frac", "written"]

EXTRACT_QC_FIELDS = [
    "wsi_id", "patient", "tss", "label", "n_candidate", "n_written",
    "n_blur_rej", "n_ink_rej", "n_tissue_rej", "status", "note", "secs"]


def stain_metrics(rgb):
    """
    H&E stain deconvolution tabanli iki olcum.

    NEDEN: renk tonu (hue) murekkebi hematoksilenden ayiramiyor.
    Gorsel denetimde, cekirdek yogun gercek tumor dokusu (TCGA-83-5908)
    hue-tabanli ink maskesinde en yuksek skorlari aldi. Ton araligini
    daraltmak yetmedi -- boyama yogunlugu merkeze gore kaydigi icin
    sabit bir ton siniri her merkezde ayni yere dusmuyor.

    Deconvolution renge degil BOYA KIMYASINA bakar:
      OD = -log10(I/256)  ->  H ve E vektorlerinin gerdigi duzleme izdusum

      he_residual   patch'in rengi H&E ile ne kadar aciklanamiyor
                    gercek doku -> DUSUK (rengi H&E'den geliyor)
                    murekkep    -> YUKSEK (H&E disi pigment)
      hema_frac     hematoksilen konsantrasyonu yuksek piksel orani
                    = cekirdek/hucresellik olcusu
                    alveol bosluklu doku -> dusuk
                    tumor yatagi         -> yuksek

    Donus: (he_residual, hema_frac)
    """
    od = -np.log10((rgb.astype(np.float64) + 1.0) / 256.0)
    flat = od.reshape(-1, 3)

    # yalnizca boyali pikseller; bos beyaz alan kalintiyi yapay dusurur
    dense = flat.sum(axis=1) > 0.15
    if dense.sum() < 50:
        return 0.0, 0.0

    sub = flat[dense]
    coeffs = sub @ _HE_PINV.T                 # (N,2) H ve E konsantrasyonu
    recon = coeffs @ _HE
    resid = np.linalg.norm(sub - recon, axis=1)
    denom = np.linalg.norm(sub, axis=1) + 1e-8
    he_residual = float(np.median(resid / denom))

    # hematoksilen konsantrasyonu -> cekirdek benzeri piksel orani
    hema = coeffs[:, 0]
    hema_frac = float((hema > 0.35).sum() / flat.shape[0])

    return he_residual, hema_frac


def patch_metrics(rgb):
    """
    224x224 patch uzerinde bes olcum. HEPSI RESIZE SONRASI hesaplanir.

    tissue       doygunluk tabanli doku orani (patch'in KENDI pikselleri)
                 NOT: bu metrik boyama yogunluguna bagimli, merkeze gore
                 0.58-0.96 arasi degisiyor -- filtreleme icin KULLANILMAZ
    ink_hue      eski hue-tabanli olcum (karsilastirma icin tutuluyor)
    blur         Laplacian varyansi (odak disi bolgeler)
    he_residual  H&E ile aciklanamayan renk  -> murekkep gostergesi
    hema_frac    cekirdek yogunlugu          -> hucresellik gostergesi
    """
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    hue, sat, val = hsv[..., 0], hsv[..., 1], hsv[..., 2]

    tissue = ((sat > 20) & (val > 30) & (val < 235)).mean()

    ink_hue = ((sat > 90) & (val > 40) &
               (((hue >= 95) & (hue < 125)) |
                ((hue >= 35) & (hue < 95)))).mean()

    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    blur = float(cv2.Laplacian(gray, cv2.CV_64F).var())

    he_residual, hema_frac = stain_metrics(rgb)

    return float(tissue), float(ink_hue), blur, he_residual, hema_frac


def extract(out_dir, slides_root, cfg, cap=2000,
            blur_min=0.0, ink_max=1.0, tissue_min=0.0, limit=0,
            per_tss=0, res_max=1.0, hema_min=0.0):
    """
    Manifest koordinatlarindan gercek 224x224 JPEG yazar.

    VARSAYILAN: hicbir esik uygulanmaz (blur_min=0, ink_max=1, tissue_min=0).
    Tum olcumler kaydedilir. Esikler PILOT CIKTISINA BAKILARAK secilir ve
    sonradan final_patch_manifest.csv uzerinden filtrelenir -- yeniden
    cikarma gerekmez.
    """
    out = Path(out_dir).expanduser()
    root = Path(slides_root).expanduser().resolve()
    man = out / "canonical_patch_manifest.csv"
    if not man.exists():
        raise SystemExit(f"manifest yok: {man}")

    by_wsi = {}
    for r in csv.DictReader(open(man)):
        by_wsi.setdefault(r["wsi_id"], []).append(r)

    wids = sorted(by_wsi)

    if per_tss:
        # MERKEZ-KATMANLI PILOT: her TSS'ten per_tss adet WSI.
        # Alfabetik --limit tek merkezden ornek alir (TCGA-05-* gibi) ve
        # murekkep/boyama dagilimini temsil etmez.
        by_tss = {}
        for w in wids:
            t = by_wsi[w][0].get("tss", "??")
            by_tss.setdefault(t, []).append(w)
        wids = []
        for t in sorted(by_tss):
            wids.extend(by_tss[t][:per_tss])
        print(f"  merkez-katmanli pilot: {len(by_tss)} TSS x {per_tss} "
              f"= {len(wids)} WSI")
    elif limit:
        wids = wids[:limit]

    pdir = out / "patches"
    pdir.mkdir(parents=True, exist_ok=True)
    fin_path = out / "final_patch_manifest.csv"
    qc_path = out / "extract_qc.csv"

    done = set()
    if qc_path.exists():
        done = {r["wsi_id"] for r in csv.DictReader(open(qc_path))}
        print(f"  {len(done)} WSI zaten cikarilmis, atlaniyor")

    print(f"{len(wids)} WSI | cap {cap}/WSI | esikler: blur>={blur_min} "
          f"he_residual<={res_max} hema_frac>={hema_min}")
    if (blur_min == 0 and ink_max == 1 and tissue_min == 0
            and res_max == 1.0 and hema_min == 0.0):
        print("  -> OLCUM MODU: hicbir patch elenmiyor, metrikler kaydediliyor\n")

    ff = open(fin_path, "a", newline="")
    fw = csv.DictWriter(ff, fieldnames=EXTRACT_FIELDS, extrasaction="ignore")
    qf = open(qc_path, "a", newline="")
    qw = csv.DictWriter(qf, fieldnames=EXTRACT_QC_FIELDS, extrasaction="ignore")
    if not done:
        fw.writeheader()
        qw.writeheader()

    allm = []
    for i, wid in enumerate(wids, 1):
        if wid in done:
            continue
        rows = by_wsi[wid]
        t0 = time.time()
        qc = {"wsi_id": wid, "patient": rows[0]["patient"],
              "tss": rows[0]["tss"], "label": rows[0]["label"],
              "n_candidate": len(rows), "n_written": 0,
              "n_blur_rej": 0, "n_ink_rej": 0, "n_tissue_rej": 0,
              "status": "ok", "note": ""}

        hits = list(root.rglob(wid + "*.svs"))
        if not hits:
            qc.update(status="error", note="svs bulunamadi")
            qw.writerow(qc); qf.flush()
            continue
        try:
            slide = openslide.OpenSlide(str(hits[0]))
        except Exception as e:
            qc.update(status="error", note=str(e)[:60])
            qw.writerow(qc); qf.flush()
            continue

        wdir = pdir / rows[0]["patient"] / wid
        wdir.mkdir(parents=True, exist_ok=True)
        px = cfg["patch_px"]

        for r in rows:
            if qc["n_written"] >= cap:
                break
            rs = int(r["read_size_px"])
            try:
                img = slide.read_region((int(r["x0"]), int(r["y0"])),
                                        int(r["read_level"]), (rs, rs))
            except Exception:
                continue
            arr = np.array(img.convert("RGB"))

            # --- ONCE resize, SONRA QC (hata 3.3 / kural 1.3) ---
            if rs != px:
                interp = cv2.INTER_AREA if rs > px else cv2.INTER_LINEAR
                arr = cv2.resize(arr, (px, px), interpolation=interp)

            tis, ink, blur, heres, hema = patch_metrics(arr)
            if blur < blur_min:
                qc["n_blur_rej"] += 1
                continue
            if heres > res_max:
                qc["n_ink_rej"] += 1
                continue
            if hema < hema_min:
                qc["n_tissue_rej"] += 1
                continue
            if ink_max < 1.0 and ink > ink_max:
                qc["n_ink_rej"] += 1
                continue
            if tis < tissue_min:
                qc["n_tissue_rej"] += 1
                continue

            Image.fromarray(arr).save(wdir / r["patch_file"],
                                      quality=cfg["jpeg_quality"])
            out_row = dict(r)
            out_row.update(tissue_patch=round(tis, 4), ink_patch=round(ink, 4),
                           blur_var=round(blur, 1),
                           he_residual=round(heres, 4),
                           hema_frac=round(hema, 4), written=1)
            fw.writerow(out_row)
            allm.append((tis, ink, blur, heres, hema, qc["tss"]))
            qc["n_written"] += 1

        slide.close()
        qc["secs"] = round(time.time() - t0, 1)
        if qc["n_written"] < cap * 0.5:
            qc["status"] = "low"
        qw.writerow(qc)
        ff.flush(); qf.flush()

        if i % 10 == 0 or i == len(wids) or qc["status"] != "ok":
            print(f"[{i}/{len(wids)}] {wid:32s} {qc['status']:6s} "
                  f"yazilan={qc['n_written']:5} blur_rej={qc['n_blur_rej']:4} "
                  f"{qc['secs']:6.1f}s {qc['note']}")

    ff.close(); qf.close()
    return allm, fin_path


def extract_report(metrics):
    """Pilot ciktisi: esikleri BU dagilima bakarak sec."""
    if not metrics:
        return
    tss_col = [r[-1] for r in metrics]
    m = np.array([r[:5] for r in metrics], dtype=float)
    print(f"\n{'='*74}\nPATCH METRIK DAGILIMI  (n={len(m):,})\n{'='*74}")
    print(f"  {'metrik':<14}{'min':>10}{'p1':>10}{'p5':>10}{'p25':>10}"
          f"{'ortanca':>10}{'p95':>10}{'max':>10}")
    for j, name in enumerate(("tissue_patch", "ink_hue", "blur_var",
                              "he_residual", "hema_frac")):
        q = np.percentile(m[:, j], [0, 1, 5, 25, 50, 95, 100])
        print(f"  {name:<14}" + "".join(f"{v:>10.2f}" for v in q))

    print(f"\n  ESIK SECIMI -- kac patch elenir?")
    print(f"  {'blur_min':<12}{'elenen':>10}{'oran':>8}")
    for t in (10, 25, 50, 100, 200):
        n = (m[:, 2] < t).sum()
        print(f"  {t:<12}{n:>10,}{n/len(m)*100:>7.1f}%")
    print(f"\n  {'ink_max':<12}{'elenen':>10}{'oran':>8}")
    for t in (0.02, 0.05, 0.10, 0.20):
        n = (m[:, 1] > t).sum()
        print(f"  {t:<12}{n:>10,}{n/len(m)*100:>7.1f}%")
    print(f"\n  {'he_res_max':<12}{'elenen':>10}{'oran':>8}   (murekkep)")
    for t in (0.10, 0.15, 0.20, 0.25, 0.30):
        n = (m[:, 3] > t).sum()
        print(f"  {t:<12}{n:>10,}{n/len(m)*100:>7.1f}%")
    print(f"\n  {'hema_min':<12}{'elenen':>10}{'oran':>8}   (hucresellik)")
    for t in (0.02, 0.05, 0.10, 0.20):
        n = (m[:, 4] < t).sum()
        print(f"  {t:<12}{n:>10,}{n/len(m)*100:>7.1f}%")

    # --- merkez kirilimi: boyama/murekkep merkeze bagli mi? ---
    uniq = sorted(set(tss_col))
    if len(uniq) > 1:
        print(f"\n  MERKEZ (TSS) KIRILIMI -- {len(uniq)} merkez")
        print(f"  {'tss':<6}{'n':>8}{'tissue_ort':>12}{'hema_ort':>10}"
              f"{'res_ort':>10}{'res_p95':>10}{'blur_ort':>10}")
        arr = np.array(tss_col)
        rows = []
        for t in uniq:
            sel = m[arr == t]
            if len(sel) < 5:
                continue
            rows.append((t, len(sel), np.median(sel[:, 0]),
                         np.median(sel[:, 4]), np.median(sel[:, 3]),
                         np.percentile(sel[:, 3], 95), np.median(sel[:, 2])))
        for r in rows:
            print(f"  {r[0]:<6}{r[1]:>8}{r[2]:>12.2f}{r[3]:>10.3f}"
                  f"{r[4]:>10.3f}{r[5]:>10.3f}{r[6]:>10.0f}")
        if rows:
            tv = [r[2] for r in rows]
            iv = [r[5] for r in rows]
            hv = [r[3] for r in rows]
            print(f"\n  tissue_patch ortancasi merkezler arasi: "
                  f"{min(tv):.2f} - {max(tv):.2f}")
            print(f"  hema_frac ortancasi merkezler arasi: "
                  f"{min(hv):.3f} - {max(hv):.3f}")
            print(f"  he_residual p95 merkezler arasi: "
                  f"{min(iv):.3f} - {max(iv):.3f}")
            if max(tv) - min(tv) > 0.15:
                print("  [!] tissue_patch merkeze gore belirgin degisiyor "
                      "-> bu metrikle filtreleme merkez-bagimli secilim yaratir")