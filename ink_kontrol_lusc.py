#!/usr/bin/env python3
"""
Kalem izi tespitini manuel dogrulama araci.

Ne yapar:
  1. scan_summary.csv'den kalem orani listesini basar
  2. Her slayt icin bir onizleme PNG'i uretir:
       sol   = orijinal thumbnail
       sag   = kalem maskesi kirmizi ile isaretlenmis
  3. Kontrol listesi CSV'si yazar (sen doldurursun)

Boylece ImageScope acmadan, tek klasorde hepsine bakabilirsin.

Kullanim:
    # en cok kalem izli 25 slayt (varsayilan)
    python ink_kontrol.py

    # esigin hemen altindakiler -- kacan var mi kontrolu
    python ink_kontrol.py --mod sinir

    # rastgele temiz slaytlar -- yanlis pozitif kontrolu
    python ink_kontrol.py --mod temiz

    # tek slayt
    python ink_kontrol.py --wsi TCGA-86-7953
"""
import argparse
import csv
from pathlib import Path

import cv2
import numpy as np
import openslide

OUT = Path("~/tez/out/LUSC").expanduser()
SLIDES = Path("~/tcga_data/LUSC").expanduser()
THUMB = 1400

# patches.py ile AYNI olmali
PEN_RESIDUAL = 0.55
PEN_SAT = 60
PEN_VAL = 30

_HE = np.array([[0.650, 0.704, 0.286], [0.072, 0.990, 0.105]])
_HE = _HE / np.linalg.norm(_HE, axis=1, keepdims=True)
_P = np.linalg.pinv(_HE.T)


def maskeler(thumb):
    """Doku ve kalem maskeleri -- patches.py mantiginin aynisi."""
    hsv = cv2.cvtColor(thumb, cv2.COLOR_RGB2HSV)
    sat, val = hsv[..., 1], hsv[..., 2]

    _, tis = cv2.threshold(cv2.GaussianBlur(sat, (5, 5), 0), 0, 255,
                           cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    tis[val < 30] = 0
    tis[val > 235] = 0
    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    tis = cv2.morphologyEx(tis, cv2.MORPH_CLOSE, k)
    tis = cv2.morphologyEx(tis, cv2.MORPH_OPEN, k)

    od = -np.log10((thumb.astype(np.float64) + 1.0) / 256.0).reshape(-1, 3)
    dense = od.sum(axis=1) > 0.15
    resid = np.zeros(od.shape[0], dtype=np.float32)
    if dense.sum():
        sub = od[dense]
        c = sub @ _P.T
        resid[dense] = (np.linalg.norm(sub - c @ _HE, axis=1)
                        / (np.linalg.norm(sub, axis=1) + 1e-8)).astype(np.float32)
    resid = resid.reshape(thumb.shape[:2])

    pen = (resid > PEN_RESIDUAL) & (sat > PEN_SAT) & (val > PEN_VAL)
    k3 = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    k7 = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (7, 7))
    pen = cv2.morphologyEx(pen.astype(np.uint8), cv2.MORPH_OPEN, k3)
    pen = cv2.dilate(pen, k7) > 0

    return tis > 0, pen


def onizleme(wsi_id, oran, hedef):
    hits = list(SLIDES.rglob(wsi_id + "*.svs"))
    if not hits:
        print(f"  [!] {wsi_id}: svs bulunamadi")
        return False

    s = openslide.OpenSlide(str(hits[0]))
    w, h = s.dimensions
    sc = max(w, h) / THUMB
    thumb = np.array(s.get_thumbnail((int(w / sc), int(h / sc))).convert("RGB"))
    s.close()

    tis, pen = maskeler(thumb)
    sag = thumb.copy()
    sag[pen] = (0.35 * sag[pen] + 0.65 * np.array([255, 40, 40])).astype(np.uint8)

    # ayirici cizgi ile yan yana
    bosluk = np.full((thumb.shape[0], 8, 3), 40, dtype=np.uint8)
    birlesik = np.hstack([thumb, bosluk, sag])

    # ust serit: bilgi
    serit = np.full((34, birlesik.shape[1], 3), 25, dtype=np.uint8)
    cv2.putText(serit, f"{wsi_id}   kalem: %{oran*100:.2f}",
                (10, 23), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (240, 240, 240), 1,
                cv2.LINE_AA)
    cv2.putText(serit, "sol: orijinal   |   sag: kirmizi = kalem maskesi",
                (birlesik.shape[1] - 430, 23), cv2.FONT_HERSHEY_SIMPLEX,
                0.45, (150, 150, 150), 1, cv2.LINE_AA)
    birlesik = np.vstack([serit, birlesik])

    ad = f"{oran*100:06.2f}_{wsi_id}.png"
    cv2.imwrite(str(hedef / ad), cv2.cvtColor(birlesik, cv2.COLOR_RGB2BGR))
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mod", choices=["kalemli", "sinir", "temiz"],
                    default="kalemli",
                    help="kalemli: en yuksek | sinir: esik civari | temiz: %0 olanlar")
    ap.add_argument("--n", type=int, default=25)
    ap.add_argument("--wsi", default=None, help="tek slayt")
    args = ap.parse_args()

    rows = [r for r in csv.DictReader(open(OUT / "scan_summary.csv"))
            if r["status"] == "ok"]
    for r in rows:
        r["_pen"] = float(r["pen_frac_slide"])

    if args.wsi:
        sec = [r for r in rows if r["wsi_id"].startswith(args.wsi)]
        klasor = "tek"
    elif args.mod == "kalemli":
        sec = sorted(rows, key=lambda r: -r["_pen"])[:args.n]
        klasor = "kalemli"
    elif args.mod == "sinir":
        # esigin hemen altinda kalanlar: %0.1 - %2 arasi
        ara = [r for r in rows if 0.001 <= r["_pen"] < 0.02]
        sec = sorted(ara, key=lambda r: -r["_pen"])[:args.n]
        klasor = "sinir"
    else:
        temiz = [r for r in rows if r["_pen"] < 0.0005]
        rng = np.random.default_rng(42)
        idx = rng.permutation(len(temiz))[:args.n]
        sec = [temiz[i] for i in idx]
        klasor = "temiz"

    if not sec:
        print("secilen slayt yok")
        return

    hedef = OUT / "ink_kontrol" / klasor
    hedef.mkdir(parents=True, exist_ok=True)

    print(f"\n{'wsi_id':<34}{'tss':>5}{'kalem':>9}{'etkilenen aday':>16}")
    print("-" * 66)
    kayit = []
    for r in sec:
        print(f"{r['wsi_id']:<34}{r['tss']:>5}"
              f"{r['_pen']*100:>8.2f}%{int(r['n_pen_hit']):>16,}")
        if onizleme(r["wsi_id"], r["_pen"], hedef):
            kayit.append({"wsi_id": r["wsi_id"], "tss": r["tss"],
                          "kalem_orani": round(r["_pen"], 4),
                          "etkilenen_aday": r["n_pen_hit"],
                          "GERCEKTEN_KALEM_VAR_MI": ""})

    # kontrol listesi -- son kolonu sen doldur
    kp = hedef / "kontrol_listesi.csv"
    with open(kp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(kayit[0].keys()))
        w.writeheader()
        w.writerows(kayit)

    print(f"\n{len(kayit)} onizleme: {hedef}")
    print(f"kontrol listesi: {kp}")
    print("\nWindows'tan:")
    print(f"  \\\\wsl.localhost\\Ubuntu-22.04{str(hedef).replace('/', chr(92))}")
    print("\nDosya adlari kalem oraniyla basliyor -> klasorde sirali gorunur.")
    print("Her goruntude SOL orijinal, SAG kalem maskesi kirmizi.")
    print("kontrol_listesi.csv'nin son kolonuna evet/hayir yaz.")


if __name__ == "__main__":
    main()