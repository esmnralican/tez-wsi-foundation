#!/usr/bin/env python3
"""
Merkez (TSS) bazinda hata orani, tohumlar arasi toplanmis.

Her tohum her slayda tam bir tahmin verir (slayt bes fold'un yalnizca
birinde testtedir). Dolayisiyla her (merkez, tohum) cifti icin bir hata
orani hesaplanir; tablo bu oranlarin tohumlar arasi ortalamasini, en
kucuk-en buyuk araligini ve slayt sayisina dayali Wilson guven araligini
verir.

Tohum tekrari TOHUM gurultusunu azaltir, BINOM gurultusunu azaltmaz:
38 slaytlik bir merkezde %5 hata iki slayt demektir. Wilson araligi bu
ikinci kaynagi gosterir ve merkez basina slayt sayisi kucuk oldugunda
genis kalir -- bu bir olcum hatasi degil, verinin siniri.

Kullanim:
    python merkez_hata.py --model dinov2
    python merkez_hata.py --model uni --en-az 20
    python merkez_hata.py --model abl_n100     # ablasyon kosulari da okunur
"""

import argparse
import csv
import glob
import re
from collections import defaultdict

import numpy as np


def wilson(k, n, z=1.96):
    """Oran icin Wilson skor araligi. Kucuk n'de normal yaklasimdan iyi."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    orta = (p + z * z / (2 * n)) / d
    yari = z / d * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, orta - yari), min(1.0, orta + yari))


def oku(model, mod):
    hata = defaultdict(lambda: [0, 0])   # (tohum, tss) -> [yanlis, toplam]
    sinif = {}
    slaytlar = defaultdict(set)
    tohumlar = set()

    for f in sorted(glob.glob("out/mil/*/predictions.csv")):
        tag = f.split("/")[2]
        m = re.match(rf"^{re.escape(model)}_{mod}(_k5|_s(\d+))$", tag)
        if not m:
            continue
        tohum = int(m.group(2)) if m.group(2) else 42
        tohumlar.add(tohum)
        for r in csv.DictReader(open(f)):
            tss = r["tss"]
            hata[(tohum, tss)][1] += 1
            if r["tahmin"] != r["gercek"]:
                hata[(tohum, tss)][0] += 1
            sinif[tss] = r["gercek"]
            slaytlar[tss].add(r["wsi_id"])

    return hata, sinif, {k: len(v) for k, v in slaytlar.items()}, sorted(tohumlar)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True,
                    help="kosu adi: dinov2, uni, virchow, abl_n100 ...")
    ap.add_argument("--en-az", type=int, default=20,
                    help="en az bu kadar slayti olan merkezler listelenir")
    a = ap.parse_args()

    for mod in ("patient", "site"):
        hata, sinif, n_slayt, tohumlar = oku(a.model, mod)
        if not hata:
            print(f"\n[{a.model} / {mod}] kosu bulunamadi")
            continue

        print(f"\n=== {a.model.upper()} / {mod.upper()} "
              f"- MERKEZ BAZINDA HATA ===")
        print(f"tohum: {tohumlar}  (n = {len(tohumlar)})   "
              f"merkez: {len(sinif)}   "
              f"listelenen: slayt >= {a.en_az}")
        print(f"\n{'tss':>4} {'sinif':<6} {'slayt':>6} {'hata':>8} "
              f"{'tohum araligi':>17} {'%95 Wilson':>16}")
        print("-" * 64)

        satirlar = []
        for tss in sinif:
            oranlar = [hata[(t, tss)][0] / hata[(t, tss)][1]
                       for t in tohumlar if hata[(t, tss)][1] > 0]
            if not oranlar or n_slayt[tss] < a.en_az:
                continue
            ort = float(np.mean(oranlar))
            lo, hi = wilson(int(round(ort * n_slayt[tss])), n_slayt[tss])
            satirlar.append((n_slayt[tss], tss, sinif[tss], ort,
                             min(oranlar), max(oranlar), lo, hi))

        for n, tss, sn, ort, mn, mx, lo, hi in sorted(satirlar, reverse=True):
            print(f"{tss:>4} {sn:<6} {n:>6} {ort*100:>7.1f}% "
                  f"{mn*100:>7.1f} - {mx*100:<7.1f} "
                  f"[{lo*100:>5.1f}, {hi*100:>5.1f}]")

        y = sum(hata[(t, s)][0] for t in tohumlar for s in sinif)
        tp = sum(hata[(t, s)][1] for t in tohumlar for s in sinif)
        print("-" * 64)
        print(f"{'TUM':>4} {'':<6} {sum(n_slayt.values()):>6} "
              f"{y / tp * 100:>7.1f}%   (tum merkezler, tum tohumlar)")

    print("\nTohum araligi dar + Wilson araligi genis = tohum degil, "
          "slayt sayisi sinirliyor.")


if __name__ == "__main__":
    main()
