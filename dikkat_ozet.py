#!/usr/bin/env python3
"""
Dikkat yogunlasmasi ozeti: out/dikkat/ altindaki tum kosulari toplar.

Her kosu bir (model, sema, fold) ucgenidir. Slayt seviyesindeki olculerin
fold ici ORTANCASI alinir (slayt dagilimi carpik oldugu icin ortalama degil),
sonra fold'lar arasi ortalama ve en kucuk-en buyuk aralik raporlanir.

Kullanim: python dikkat_ozet.py
"""

import csv
import glob
import re
from collections import defaultdict

import numpy as np


def main():
    g = defaultdict(list)   # (model, sema) -> [(fold, ust1, etkin, dogruluk, n)]
    for f in sorted(glob.glob("out/dikkat/*/dikkat_ozet.csv")):
        ad = f.split("/")[2]
        m = re.match(r"^(.+)_(patient|site)_k5_fold(\d+)(?:_s(\d+))?$", ad)
        if not m:
            print(f"[atlandi] {ad}")
            continue
        model, sema, fold = m.group(1), m.group(2), int(m.group(3))
        tohum = int(m.group(4)) if m.group(4) else 42
        fold = (tohum, fold)
        r = list(csv.DictReader(open(f)))
        if not r:
            print(f"[bos] {ad}")
            continue
        ust1 = float(np.median([float(x["ust1_pct_dikkat"]) for x in r]))
        etkin = float(np.median([float(x["etkin_patch"]) for x in r]))
        dog = float(np.mean([int(x["dogru"]) for x in r]))
        g[(model, sema)].append((fold, ust1, etkin, dog, len(r)))

    if not g:
        raise SystemExit("out/dikkat/ altinda ozet bulunamadi")

    print(f"\n{'model':<12} {'sema':<9} {'fold':>6} {'ust%1 dikkat':>13} "
          f"{'etkin patch':>12} {'dogruluk':>9} {'slayt':>6}")
    print("=" * 72)
    for anahtar in sorted(g):
        model, sema = anahtar
        v = sorted(g[anahtar])
        for fold, u, e, d, n in v:
            print(f"{model:<12} {sema:<9} {str(fold):>10} {u*100:>12.1f}% "
                  f"{e:>12.0f} {d*100:>8.1f}% {n:>6}")
        u = np.array([x[1] for x in v])
        e = np.array([x[2] for x in v])
        d = np.array([x[3] for x in v])
        print(f"{model:<12} {sema:<9} {'ORT':>6} {u.mean()*100:>12.1f}% "
              f"{e.mean():>12.0f} {d.mean()*100:>8.1f}%")
        print(f"{'':<12} {'':<9} {'aralik':>6} "
              f"{u.min()*100:>6.1f}-{u.max()*100:<5.1f} "
              f"{e.min():>7.0f}-{e.max():<5.0f}")
        print("-" * 72)

    print("\nDuzgun dagilimda: ust %1 payi = %1, etkin patch = 2000")


if __name__ == "__main__":
    main()
