#!/usr/bin/env python3
"""
Kosullar arasi eslestirilmis karsilastirma (Nadeau-Bengio duzeltmeli).

Modeller ayni fold'larda degerlendirildiginden skorlar bagimlidir.
Bagimsiz varsayan karsilastirma yanlis sonuc verir. Bu script her
(tohum, fold) ciftinde farki alir ve CV fold'larinin ortusen egitim
kumelerinden gelen bagimliligi Nadeau-Bengio duzeltmesiyle hesaba katar:

    duzeltilmis_var = var(farklar) * (1/n + n_test/n_egitim)

Kaynak: Nadeau & Bengio 2003; scikit-learn model_selection dokumantasyonu.

Kullanim:
    python karsilastir.py                      # tum kosullar, site bolmesi
    python karsilastir.py --mod patient
    python karsilastir.py --abl uni --mod site # ornekleme yogunlugu ablasyonu
"""

import argparse
import csv
import glob
import re
from collections import defaultdict

import numpy as np
from scipy import stats

N_TEST = 189      # fold basina ortalama test slaydi
N_TRAIN = 755     # fold basina egitim + ic dogrulama


def kosullari_oku(mod, abl=None, secili=None):
    """
    {kosul: {(tohum, fold): auc}}

    abl verilirse (model adi, orn. 'uni') ablasyon kosullari (abl_n*)
    ve o modelin tam kosusu okunur; normal kosullar atlanir. Ablasyon
    noktalari ayni fold/tohumda olctugu icin eslestirilmis karsilastirma
    ayni sekilde gecerlidir -- 2.000 patch noktasi modelin kendi kosusudur.
    """
    d = defaultdict(dict)
    for f in glob.glob("out/mil/*/fold_metrics.csv"):
        tag = f.split("/")[2]
        m = re.match(rf"^(.+?)_{mod}(_k5|_s(\d+))$", tag)
        if not m:
            continue
        kosul = m.group(1)
        if secili is not None:
            if kosul not in secili:
                continue
        elif abl:
            if not (kosul.startswith("abl_n") or kosul == abl):
                continue
        elif kosul.startswith(("abl_", "renk_")) or "_n100" in kosul or "_cls" in kosul:
            continue
        tohum = int(m.group(3)) if m.group(3) else 42
        for r in csv.DictReader(open(f)):
            d[kosul][(tohum, int(r["fold"]))] = float(r["auc"])
    return d


def karsilastir(a, b):
    """
    Eslestirilmis fark, duzeltilmis guven araligi ve p degeri.
    Ortak (tohum, fold) anahtarlari kullanilir.
    """
    ortak = sorted(set(a) & set(b))
    if len(ortak) < 3:
        return None
    farklar = np.array([a[k] - b[k] for k in ortak])
    n = len(farklar)

    duz_var = farklar.var(ddof=1) * (1 / n + N_TEST / N_TRAIN)
    duz_se = np.sqrt(duz_var)
    t_ist = farklar.mean() / duz_se if duz_se > 0 else 0.0
    p = 2 * stats.t.sf(abs(t_ist), df=n - 1)
    kritik = stats.t.ppf(0.975, df=n - 1)

    return {
        "n": n,
        "fark": farklar.mean(),
        "se": duz_se,
        "alt": farklar.mean() - kritik * duz_se,
        "ust": farklar.mean() + kritik * duz_se,
        "p": p,
        "ham_se": farklar.std(ddof=1) / np.sqrt(n),   # duzeltmesiz, kiyas icin
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mod", default="site", choices=["site", "patient"])
    ap.add_argument("--abl", metavar="MODEL", default=None,
                    help="ornekleme yogunlugu ablasyonu; MODEL = 2.000 patch "
                         "referansi (orn. uni, dinov2)")
    ap.add_argument("--kosul", metavar="A,B", default=None,
                    help="yalnizca bu kosullari karsilastir "
                         "(diger filtreleri atlar)")
    a = ap.parse_args()

    secili = set(a.kosul.split(",")) if a.kosul else None
    d = kosullari_oku(a.mod, a.abl, secili)
    if not d:
        raise SystemExit(f"out/mil/ altinda '{a.mod}' kosulu bulunamadi")
    if a.abl:
        print(f"[ablasyon modu] '{a.abl}' satiri = 2.000 patch (tam kosu)\n")

    print(f"=== {a.mod.upper()} BOLME - KOSUL OZETI ===")
    print(f"{'kosul':<16} {'n':>4} {'ortalama AUC':>14}")
    print("-" * 36)
    if a.abl:
        # egri okunabilsin diye patch sayisina gore sirala (model = 2.000)
        sirali = sorted(d, key=lambda k: 2000 if k == a.abl
                        else int(re.match(r"abl_n(\d+)", k).group(1)))
    else:
        sirali = sorted(d, key=lambda k: -np.mean(list(d[k].values())))
    for k in sirali:
        v = np.array(list(d[k].values()))
        print(f"{k:<16} {len(v):>4} {v.mean():>14.4f}")

    print(f"\n=== ESLESTIRILMIS KARSILASTIRMALAR ===")
    print(f"{'A':<14} {'B':<14} {'fark':>8} {'%95 guven araligi':>22} "
          f"{'p':>8}  sonuc")
    print("-" * 82)

    for i, ka in enumerate(sirali):
        for kb in sirali[i + 1:]:
            r = karsilastir(d[ka], d[kb])
            if not r:
                continue
            anlamli = r["p"] < 0.05
            sonuc = "FARKLI" if anlamli else "ayirt edilemez"
            print(f"{ka:<14} {kb:<14} {r['fark']:>8.4f} "
                  f"[{r['alt']:>8.4f}, {r['ust']:>8.4f}] "
                  f"{r['p']:>8.4f}  {sonuc}")

    print("\nGuven araligi 0'i iceriyorsa fark gosterilemedi demektir.")
    print("Aralik dar ve 0'i iceriyorsa 'fark varsa bile su sinirin altinda'")
    print("seklinde raporlanabilir.")

    # duzeltmenin etkisi
    ornek = karsilastir(d[sirali[0]], d[sirali[-1]])
    if ornek:
        print(f"\nDuzeltmenin etkisi ({sirali[0]} vs {sirali[-1]}):")
        print(f"  duzeltmesiz standart hata : {ornek['ham_se']:.4f}")
        print(f"  Nadeau-Bengio duzeltmeli  : {ornek['se']:.4f}  "
              f"({ornek['se']/ornek['ham_se']:.1f}x)")


if __name__ == "__main__":
    main()
