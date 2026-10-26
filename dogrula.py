#!/usr/bin/env python3
"""
Istatistik dogrulama -- bagimsiz yeniden hesaplama.

Bu script karsilastir.py'yi IMPORT ETMEZ. Ham fold_metrics.csv
dosyalarindan baslayarak her seyi yeniden hesaplar ve bildiride
raporlanan degerlerle karsilastirir. Amac, formulun aciklanmasi degil,
kodun o formulu gercekten uyguladiginin gosterilmesidir.

Bes kontrol:

  1. Eslestirilmis test mekanizmasi: duzeltme kapaliyken (carpan = 1/n)
     sonucumuz scipy.stats.ttest_rel ile BIREBIR ayni olmali. Ayni
     cikiyorsa eslestirme ve t hesabi dogru demektir; geriye yalnizca
     duzeltme carpani kalir.

  2. Nadeau-Bengio adim adim: ara buyuklukler tek tek yazdirilir
     (var(d), 1/n, n_test/n_train, carpan, duzeltilmis var, SE, t, p, GA)
     ve bildirideki degerle karsilastirilir.

  3. n_train duyarliligi: 755 (test disi her sey) ve 642 (yalnizca
     agirlik egitimi) ile ayni hesap tekrarlanir, SE orani raporlanir.

  4. Holm adim adim: her sirada p degeri ve esigi yazdirilir, nerede
     durdugu gorunur.

  5. Holm, her iki n_train ile: sonucun degisip degismedigi.

Kullanim:
    python dogrula.py                 # merkez-korumali
    python dogrula.py --mod patient
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
N_TRAIN_ANA = 755      # test disi her sey (egitim + ic dogrulama)
N_TRAIN_ALT = 642      # yalnizca agirlik egitimi

# Bildiride raporlanan degerler -- script bunlari DOGRULAR, kullanmaz.
IDDIA = {
    "site":    {"se": 0.0228, "ham_se": 0.0093, "carpan": 2.5},
    "patient": {"se": 0.0067, "ham_se": 0.0027, "carpan": 2.5},
}
IDDIA_HOLM = {"site": 0, "patient": 18}      # Holm sonrasi anlamli sayisi
IDDIA_ORAN = 3.4                              # site SE / patient SE


def oku(mod):
    """Ham fold_metrics.csv -> {kosul: {(tohum, fold): auc}}"""
    d = defaultdict(dict)
    for f in glob.glob("out/mil/*/fold_metrics.csv"):
        tag = f.split("/")[2]
        m = re.match(rf"^(.+?)_{mod}(_k5|_s(\d+))$", tag)
        if not m:
            continue
        kosul = m.group(1)
        if (kosul.startswith(("abl_", "renk_")) or "_n100" in kosul
                or "_cls" in kosul):
            continue
        tohum = int(m.group(3)) if m.group(3) else 42
        for r in csv.DictReader(open(f)):
            d[kosul][(tohum, int(r["fold"]))] = float(r["auc"])
    return d


def nb(a, b, n_train, ayrinti=False):
    """Nadeau-Bengio duzeltmeli eslestirilmis t testi."""
    ortak = sorted(set(a) & set(b))
    f = np.array([a[k] - b[k] for k in ortak])
    n = len(f)
    var = f.var(ddof=1)
    t1, t2 = 1 / n, N_TEST / n_train
    carpan = t1 + t2
    duz_var = var * carpan
    se = np.sqrt(duz_var)
    t = f.mean() / se
    p = 2 * stats.t.sf(abs(t), df=n - 1)
    k = stats.t.ppf(0.975, df=n - 1)

    if ayrinti:
        print(f"    n                        = {n}")
        print(f"    ortalama fark            = {f.mean():.6f}")
        print(f"    var(farklar), ddof=1     = {var:.8f}")
        print(f"    1/n                      = {t1:.6f}")
        print(f"    n_test/n_train           = {N_TEST}/{n_train}"
              f" = {t2:.6f}")
        print(f"    carpan                   = {carpan:.6f}")
        print(f"    duzeltilmis var          = {duz_var:.8f}")
        print(f"    duzeltilmis SE           = {se:.6f}")
        print(f"    duzeltmesiz SE (sd/sqrt(n)) = "
              f"{np.sqrt(var / n):.6f}")
        print(f"    t = {t:.4f}   df = {n-1}   p = {p:.6f}")
        print(f"    %95 GA = [{f.mean()-k*se:.4f}, {f.mean()+k*se:.4f}]")
    return {"n": n, "fark": f.mean(), "se": se, "p": p,
            "ham_se": np.sqrt(var / n), "farklar": f}


def holm(pler, alfa=0.05, ayrinti=False):
    """Holm-Bonferroni. Donen: anlamli sayisi. Sinir civari yazdirilir."""
    m = len(pler)
    sirali = sorted(pler)
    dur = m
    for i, p in enumerate(sirali):
        if p >= alfa / (m - i):
            dur = i
            break
    if ayrinti:
        print(f"    {'sira':>5} {'p':>10} {'esik':>22} {'gecti':>7}")
        gosterilecek = sorted({0, 1, max(0, dur - 2), max(0, dur - 1),
                               dur, min(dur + 1, m - 1)})
        onceki = -1
        for i in gosterilecek:
            if i >= m:
                continue
            if i > onceki + 1:
                print(f"    {'...':>5}")
            esik = alfa / (m - i)
            print(f"    {i+1:>5} {sirali[i]:>10.6f} "
                  f"{f'0.05/{m-i} = {esik:.6f}':>22} "
                  f"{'evet' if sirali[i] < esik else 'HAYIR':>7}")
            onceki = i
        print(f"    -> {dur}. sirada duruldu: {dur} karsilastirma anlamli")
    return dur


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mod", default="site", choices=["site", "patient"])
    a = ap.parse_args()
    d = oku(a.mod)
    if len(d) < 9:
        raise SystemExit(f"9 kosul bekleniyordu, {len(d)} bulundu: {sorted(d)}")

    sirali = sorted(d, key=lambda k: -np.mean(list(d[k].values())))
    A, B = sirali[0], sirali[-1]        # en buyuk kontrast
    iddia = IDDIA[a.mod]

    print(f"\n{'='*64}\n{a.mod.upper()} BOLME -- {len(d)} kosul, "
          f"ham fold_metrics.csv dosyalarindan\n{'='*64}")

    # ---- 1
    print(f"\n[1] ESLESTIRILMIS TEST MEKANIZMASI")
    print(f"    Duzeltme kapali (carpan = 1/n) -> scipy.ttest_rel ile "
          f"ayni cikmali.")
    ortak = sorted(set(d[A]) & set(d[B]))
    x = np.array([d[A][k] for k in ortak])
    y = np.array([d[B][k] for k in ortak])
    f = x - y
    n = len(f)
    se0 = np.sqrt(f.var(ddof=1) / n)
    t0 = f.mean() / se0
    p0 = 2 * stats.t.sf(abs(t0), df=n - 1)
    sc = stats.ttest_rel(x, y)
    print(f"    kendi  t = {t0:.10f}   p = {p0:.10f}")
    print(f"    scipy  t = {sc.statistic:.10f}   p = {sc.pvalue:.10f}")
    fark_t = abs(t0 - sc.statistic)
    print(f"    fark   t = {fark_t:.2e}   -> "
          f"{'DOGRULANDI' if fark_t < 1e-9 else 'UYUSMUYOR'}")

    # ---- 2
    print(f"\n[2] NADEAU-BENGIO ADIM ADIM  ({A} vs {B})")
    r = nb(d[A], d[B], N_TRAIN_ANA, ayrinti=True)
    print(f"    bildiride raporlanan SE  = {iddia['se']}")
    print(f"    hesaplanan SE            = {r['se']:.4f}   -> "
          f"{'ESLESIYOR' if abs(r['se']-iddia['se']) < 5e-5 else 'UYUSMUYOR'}")
    oran = r["se"] / r["ham_se"]
    print(f"    duzeltme carpani = {oran:.2f}x  (bildiride "
          f"{iddia['carpan']}x, teorik sqrt(1+n*{N_TEST}/{N_TRAIN_ANA}) = "
          f"{np.sqrt(1+n*N_TEST/N_TRAIN_ANA):.2f}x)")

    # ---- 3
    print(f"\n[3] n_train DUYARLILIGI")
    r2 = nb(d[A], d[B], N_TRAIN_ALT)
    print(f"    n_train=755 : SE {r['se']:.6f}   p {r['p']:.6f}")
    print(f"    n_train=642 : SE {r2['se']:.6f}   p {r2['p']:.6f}")
    print(f"    SE orani    : {r2['se']/r['se']:.4f}  "
          f"-> %{100*(r2['se']/r['se']-1):.1f}  "
          f"(bildiride %7.1)")

    # ---- 4 ve 5
    for ad, ntr in (("755 (bildiride kullanilan)", N_TRAIN_ANA),
                    ("642 (muhafazakar)", N_TRAIN_ALT)):
        pler = [nb(d[x], d[y], ntr)["p"] for x, y in combinations(sirali, 2)]
        print(f"\n[{'4' if ntr == N_TRAIN_ANA else '5'}] HOLM, n_train = {ad}"
              f"   ({len(pler)} karsilastirma)")
        kalan = holm(pler, ayrinti=(ntr == N_TRAIN_ANA))
        ham = sum(1 for p in pler if p < 0.05)
        print(f"    ham p<0.05       : {ham}/{len(pler)}")
        print(f"    Holm sonrasi     : {kalan}/{len(pler)}")
        if ntr == N_TRAIN_ANA:
            bek = IDDIA_HOLM[a.mod]
            print(f"    bildiride        : {bek}  -> "
                  f"{'ESLESIYOR' if kalan == bek else 'UYUSMUYOR'}")

    print(f"\n{'='*64}")


if __name__ == "__main__":
    main()