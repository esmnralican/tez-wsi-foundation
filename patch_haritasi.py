#!/usr/bin/env python3
"""
Patch konum haritasi.

Secilen patch koordinatlarini slaydin onizlemesi uzerine cizer.
"Patch'ler nereden aliniyor?" sorusunu tahminle degil gorerek cevaplar.

Renkler:
  yesil nokta   secilen ve yazilan patch
  kirmizi alan  kalem izi maskesi (dislanmis bolge)

Kullanim:
    python patch_haritasi.py TCGA-86-8074
    python patch_haritasi.py TCGA-86-8074 TCGA-05-4244 TCGA-50-5044
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import cv2
import numpy as np
import openslide

OUT_DIR = Path("~/tez/out/LUAD").expanduser()
SLIDES = Path("~/tcga_data/LUAD").expanduser()
THUMB = 1600

_HE = np.array([[0.650, 0.704, 0.286], [0.072, 0.990, 0.105]])
_HE = _HE / np.linalg.norm(_HE, axis=1, keepdims=True)
_P = np.linalg.pinv(_HE.T)


def pen_mask(thumb):
    hsv = cv2.cvtColor(thumb, cv2.COLOR_RGB2HSV)
    sat, val = hsv[..., 1], hsv[..., 2]
    od = -np.log10((thumb.astype(np.float64) + 1.0) / 256.0).reshape(-1, 3)
    d = od.sum(axis=1) > 0.15
    r = np.zeros(od.shape[0], dtype=np.float32)
    if d.sum():
        sub = od[d]
        c = sub @ _P.T
        r[d] = (np.linalg.norm(sub - c @ _HE, axis=1)
                / (np.linalg.norm(sub, axis=1) + 1e-8)).astype(np.float32)
    r = r.reshape(thumb.shape[:2])
    pen = (r > 0.55) & (sat > 60) & (val > 30)
    k3 = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    k7 = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (7, 7))
    pen = cv2.morphologyEx(pen.astype(np.uint8), cv2.MORPH_OPEN, k3)
    return cv2.dilate(pen, k7) > 0


def harita(prefix):
    # yazilan patch'lerin koordinatlari
    coords = []
    man = OUT_DIR / "final_patch_manifest.csv"
    for r in csv.DictReader(open(man)):
        if r["wsi_id"].startswith(prefix):
            coords.append((int(r["x0"]), int(r["y0"]), int(r["size_level0"])))
    if not coords:
        print(f"[!] {prefix}: manifest'te kayit yok")
        return

    hits = list(SLIDES.rglob(prefix + "*.svs"))
    if not hits:
        print(f"[!] {prefix}: svs bulunamadi")
        return

    s = openslide.OpenSlide(str(hits[0]))
    W, H = s.dimensions
    sc = max(W, H) / THUMB
    thumb = np.array(s.get_thumbnail((int(W / sc), int(H / sc))).convert("RGB"))
    s.close()

    pen = pen_mask(thumb)
    vis = thumb.copy()
    vis[pen] = (0.55 * vis[pen] + 0.45 * np.array([255, 60, 60])).astype(np.uint8)

    xs, ys = [], []
    for x0, y0, sz in coords:
        cx, cy = int((x0 + sz / 2) / sc), int((y0 + sz / 2) / sc)
        if 0 <= cx < vis.shape[1] and 0 <= cy < vis.shape[0]:
            cv2.circle(vis, (cx, cy), 2, (0, 200, 0), -1)
            xs.append(x0 + sz / 2)
            ys.append(y0 + sz / 2)

    out = OUT_DIR / "haritalar"
    out.mkdir(exist_ok=True)
    p = out / f"{prefix}.png"
    cv2.imwrite(str(p), cv2.cvtColor(vis, cv2.COLOR_RGB2BGR))

    # kenar yogunlasmasi var mi: merkeze uzaklik dagilimi
    xs, ys = np.array(xs), np.array(ys)
    cx, cy = W / 2, H / 2
    d = np.sqrt(((xs - cx) / (W / 2)) ** 2 + ((ys - cy) / (H / 2)) ** 2)
    print(f"{prefix}: {len(coords)} patch  ->  {p.name}")
    print(f"   x araligi: {xs.min():.0f} - {xs.max():.0f}  (slayt 0 - {W})")
    print(f"   y araligi: {ys.min():.0f} - {ys.max():.0f}  (slayt 0 - {H})")
    print(f"   merkeze normalize uzaklik: ortanca {np.median(d):.2f}, "
          f"p90 {np.percentile(d, 90):.2f}")
    print(f"   (1.0'a yakin = kenarlarda, 0.5 civari = homojen dagilmis)")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("kullanim: python patch_haritasi.py TCGA-86-8074 [...]")
    for p in sys.argv[1:]:
        harita(p)
        print()