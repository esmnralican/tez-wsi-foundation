#!/usr/bin/env python3
"""
Dikkat yogunlasmasinin modeller arasi eslestirilmis testi.

Veri iki makineye dagilmis durumda (DINOv2 bir yerde, UNI/Virchow
digerinde). Bu yuzden iki asamali calisir:

  1) Her makinede  --etiket <ad>  ile calistir: yerel out/dikkat/
     altindaki her kosu icin tek satirlik ozet uretir ve
     out/dikkat_kompakt_<ad>.csv dosyasina yazar.
  2) Kompakt dosyalari tek makinede topla (git ile) ve etiketsiz
     calistir: out/dikkat_kompakt_*.csv dosyalarinin tamamini okur,
     eslestirilmis testleri uygular.

Kosular ayni fold'larda yapildigindan karsilastirma eslestirilmistir;
fold'larin egitim kumeleri ortustugunden varyans Nadeau-Bengio ile
duzeltilir (bkz. karsilastir.py).

Kullanim:
    python dikkat_test.py --etiket monster     # 1. asama
    python dikkat_test.py                      # 2. asama
"""

import argparse
import csv
import glob
import re
from collections import defaultdict
from itertools import combinations

import numpy as np
from scipy import stats

N_TEST = 189
N_TRAIN = 755

OLCULER = {
    "ust1_pct_dikkat": ("ust %1 dikkat payi (puan)", 100.0),
    "etkin_patch": ("etkin patch sayisi", 1.0),
}


def ozet_uret(etiket):
    """out/dikkat/ altindaki her kosuyu tek satira indirger."""
    satirlar = []
    for f in sorted(glob.glob("out/dikkat/*/dikkat_ozet.csv")):
        ad = f.split("/")[2]
        m = re.match(r"^(.+)_(patient|site)_k5_fold(\d+)(?:_s(\d+))?$", ad)
        if not m:
            print(f"[atlandi] {ad}")
            continue
        model, sema, fold = m.group(1), m.group(2), int(m.group(3))
        tohum = int(m.group(4)) if m.group(4) else 42
        r = list(csv.DictReader(open(f)))
        if not r:
            continue
        # slayt dagilimi carpik -> ortalama degil ORTANCA
        satirlar.append({
            "model": model, "sema": sema, "tohum": tohum, "fold": fold,
            "ust1_pct_dikkat": np.median(
                [float(x["ust1_pct_dikkat"]) for x in r]),
            "etkin_patch": np.median([float(x["etkin_patch"]) for x in r]),
            "dogruluk": np.mean([int(x["dogru"]) for x in r]),
            "n_slayt": len(r),
        })

    yol = f"out/dikkat_kompakt_{etiket}.csv"
    with open(yol, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(satirlar[0]))
        w.writeheader()
        w.writerows(satirlar)
    print(f"{len(satirlar)} kosu -> {yol}")


def kompakt_oku():
    d = defaultdict(dict)       # (model, sema) -> {(tohum, fold): satir}
    dosyalar = sorted(glob.glob("out/dikkat_kompakt_*.csv"))
    if not dosyalar:
        raise SystemExit("out/dikkat_kompakt_*.csv bulunamadi "
                         "-- once --etiket ile calistir")
    for yol in dosyalar:
        for r in csv.DictReader(open(yol)):
            anahtar = (int(r["tohum"]), int(r["fold"]))
            d[(r["model"], r["sema"])][anahtar] = r
    print("okunan dosyalar: " + ", ".join(dosyalar))
    return d


def test(a, b):
    """Eslestirilmis fark, Nadeau-Bengio duzeltmeli."""
    ortak = sorted(set(a) & set(b))
    if len(ortak) < 3:
        return None
    f = np.array([a[k] - b[k] for k in ortak])
    n = len(f)
    se = np.sqrt(f.var(ddof=1) * (1 / n + N_TEST / N_TRAIN))
    t = f.mean() / se if se > 0 else 0.0
    p = 2 * stats.t.sf(abs(t), df=n - 1)
    k = stats.t.ppf(0.975, df=n - 1)
    return {"n": n, "fark": f.mean(), "alt": f.mean() - k * se,
            "ust": f.mean() + k * se, "p": p}


def holm(pler, alfa=0.05):
    """Dondurulen: p -> Holm sonrasi anlamli mi."""
    m = len(pler)
    sonuc, esikler = {}, {}
    for i, p in enumerate(sorted(pler)):
        esikler[p] = alfa / (m - i)
    gecti = True
    for p in sorted(pler):
        gecti = gecti and (p < esikler[p])
        sonuc[p] = gecti
    return sonuc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--etiket", default=None,
                    help="1. asama: bu makinenin ozetini uret")
    a = ap.parse_args()

    if a.etiket:
        ozet_uret(a.etiket)
        return

    d = kompakt_oku()
    for sema in ("site", "patient"):
        modeller = sorted({m for (m, s) in d if s == sema})
        if len(modeller) < 2:
            continue

        print(f"\n{'='*66}\n{sema.upper()} BOLME\n{'='*66}")
        for olcu, (baslik, _) in OLCULER.items():
            print(f"\n--- {baslik} ---")
            # ust1 kaynak dosyada ORAN (0.52) olarak duruyor -> puana cevir
            olcek = OLCULER[olcu][1]
            deger = {m: {k: float(v[olcu]) * olcek
                         for k, v in d[(m, sema)].items()} for m in modeller}

            for m in sorted(modeller, key=lambda x: -np.mean(
                    list(deger[x].values()))):
                v = np.array(list(deger[m].values()))
                birim = "%" if olcu.endswith("pct_dikkat") else ""
                print(f"  {m:<9} n={len(v):<3} ort {v.mean():>7.1f}{birim}"
                      f"   aralik {v.min():>6.1f} - {v.max():<6.1f}")

            ciftler = list(combinations(
                sorted(modeller, key=lambda x: -np.mean(
                    list(deger[x].values()))), 2))
            sonuclar = [(A, B, test(deger[A], deger[B])) for A, B in ciftler]
            sonuclar = [s for s in sonuclar if s[2]]
            if not sonuclar:
                continue
            hsonuc = holm([s[2]["p"] for s in sonuclar])

            print()
            for A, B, r in sonuclar:
                etiket = "FARKLI" if r["p"] < 0.05 else "ayirt edilemez"
                hm = "  (Holm sonrasi da)" if hsonuc[r["p"]] else ""
                print(f"  {A:>8} - {B:<8} {r['fark']:>8.1f}  "
                      f"[{r['alt']:>7.1f}, {r['ust']:>7.1f}]  "
                      f"p={r['p']:.4f}  n={r['n']}  {etiket}{hm}")

    print("\nHolm: sema basina olcu basina 3 karsilastirma uzerinden.")


if __name__ == "__main__":
    main()
