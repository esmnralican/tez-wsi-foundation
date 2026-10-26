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

  7. Tablo 1: dokuz kosulun ortalama AUC degerleri ham veriden
     yeniden hesaplanir ve bildiride basili degerlerle karsilastirilir
     (18 sayi). Teskkur bolumundeki "deney loglarina karsi
     kontrol edildi" taahhudunun en kalabalik parcasi budur.

  6. Fold boyutlari: N_TEST ve N_TRAIN_ANA sabitleri splits
     dosyalarindan olculen GERCEK fold boyutlariyla karsilastirilir.
     Ayrica fold atamasinin tum kosullarda ozdes oldugu gosterilir
     (eslestirilmis testin on kosulu). Sabitler [2], [3] ve [5]
     kontrollerinde kullanildigi icin bu onlarin dayanagidir.

Kullanim:
    python dogrula.py                 # merkez-korumali
    python dogrula.py --mod patient
"""

import argparse
import csv
import glob
import hashlib
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
    "site":    {"se": 0.0228, "ham_se": 0.0093, "carpan": 2.45},
    "patient": {"se": 0.0067, "ham_se": 0.0027, "carpan": 2.45},
}
IDDIA_HOLM = {"site": 0, "patient": 18}      # Holm sonrasi anlamli sayisi

# Bildiri Tablo 1: kosul -> (hasta-katmanli, merkez-korumali) ortalama AUC.
# Elle girilmistir; script bunlari ham veriden yeniden hesaplayip
# KARSILASTIRIR. Bildiri 4 haneye yuvarladigi icin tolerans 5e-5.
IDDIA_TABLO1 = {
    "uni":         (0.9826, 0.9535),
    "uni_gri":     (0.9826, 0.9543),
    "uni_mac":     (0.9795, 0.9538),
    "virchow":     (0.9776, 0.9484),
    "virchow_gri": (0.9777, 0.9500),
    "virchow_mac": (0.9774, 0.9526),
    "dinov2":      (0.9411, 0.9017),
    "dinov2_gri":  (0.9397, 0.8993),
    "dinov2_mac":  (0.9453, 0.9111),
}
IDDIA_ORAN = 3.4                              # site SE / patient SE

# Analize girmesi beklenen TAM kosul kumesi. Kara liste degil BEYAZ LISTE:
# beklenmeyen bir etiket sessizce analize giremez.
BEKLENEN = {"uni", "uni_gri", "uni_mac",
            "virchow", "virchow_gri", "virchow_mac",
            "dinov2", "dinov2_gri", "dinov2_mac"}
BEKLENEN_TOHUM = {42, 7, 123, 2024}
BEKLENEN_FOLD = set(range(5))


def oku(mod):
    """Ham fold_metrics.csv -> (veri, yinelenenler, analiz disi etiketler)."""
    d = defaultdict(dict)
    yinelenen, atlanan = [], []
    for f in sorted(glob.glob("out/mil/*/fold_metrics.csv")):
        tag = f.split("/")[2]
        m = re.match(rf"^(.+?)_{mod}(_k5|_s(\d+))$", tag)
        if not m:
            continue
        kosul = m.group(1)
        if kosul not in BEKLENEN:
            atlanan.append(tag)
            continue
        tohum = int(m.group(3)) if m.group(3) else 42
        for r in csv.DictReader(open(f)):
            anahtar = (tohum, int(r["fold"]))
            if anahtar in d[kosul]:
                yinelenen.append((kosul, anahtar, tag))
            d[kosul][anahtar] = float(r["auc"])
    return d, yinelenen, atlanan


def butunluk(d, yinelenen, atlanan):
    """Analizden ONCE veri butunlugu. Herhangi bir ihlalde durur."""
    print("\n[0] VERI BUTUNLUGU")
    hata = []

    eksik, fazla = BEKLENEN - set(d), set(d) - BEKLENEN
    print(f"    beklenen / bulunan kosul : {len(BEKLENEN)} / {len(d)}")
    if eksik:
        hata.append(f"eksik kosul: {sorted(eksik)}")
    if fazla:
        hata.append(f"fazla kosul: {sorted(fazla)}")
    print(f"    eksik veya fazla         : "
          f"{'YOK' if not (eksik or fazla) else 'VAR'}")

    boyut = {k: len(v) for k, v in d.items()}
    print(f"    kosul basina olcum       : {sorted(set(boyut.values()))}")
    if set(boyut.values()) != {20}:
        hata.append(f"her kosulda 20 olcum beklenir: {boyut}")

    kumeler = {k: frozenset(v) for k, v in d.items()}
    ayni = len(set(kumeler.values())) == 1
    print(f"    (tohum, fold) kumeleri   : "
          f"{'hepsi ozdes' if ayni else 'FARKLI'}")
    if not ayni:
        hata.append("kosullarin anahtar kumeleri ozdes degil")

    if d:
        ornek = next(iter(kumeler.values()))
        t = {a for a, _ in ornek}
        fl = {b for _, b in ornek}
        print(f"    tohumlar                 : {sorted(t)}")
        print(f"    fold'lar                 : {sorted(fl)}")
        if t != BEKLENEN_TOHUM:
            hata.append(f"tohum kumesi farkli: {sorted(t)}")
        if fl != BEKLENEN_FOLD:
            hata.append(f"fold kumesi farkli: {sorted(fl)}")

    print(f"    yinelenen kayit          : {len(yinelenen)}")
    if yinelenen:
        hata.append(f"yinelenen (tohum, fold): {yinelenen[:5]}")
    if atlanan:
        print(f"    beyaz liste disi etiket  : {len(atlanan)} "
              f"(orn. {atlanan[:3]})")

    if hata:
        print("\n    BUTUNLUK HATASI -- analiz durduruldu:")
        for h in hata:
            print(f"      - {h}")
        raise SystemExit(1)
    print("    -> butunluk DOGRULANDI")


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


def fold_boyut(mod):
    """[6] Gercek fold boyutlari ve fold atamasinin ozdesligi."""
    print("\n[6] FOLD BOYUTLARI -- splits dosyalarindan")
    dosyalar = sorted(glob.glob(f"out/splits_*_{mod}_k5.csv"))
    if not dosyalar:
        print(f"    out/splits_*_{mod}_k5.csv bulunamadi -- ATLANDI")
        return

    # Fold atamasi: YALNIZCA wsi_id -> fold. npy sutunu kodlayiciya ozgu
    # yol icerdiginden dosyalarin bayt duzeyinde karsilastirilmasi
    # yanlistir: atama ozdes olsa bile hash farkli cikar.
    atama, imza = {}, {}
    for yol in dosyalar:
        r = list(csv.DictReader(open(yol)))
        a = {x["wsi_id"]: int(x["fold"]) for x in r}
        if len(a) != len(r):
            print(f"    UYARI {yol}: yinelenen wsi_id "
                  f"({len(r)} satir, {len(a)} benzersiz)")
        atama[yol] = a
        imza[yol] = hashlib.md5(
            repr(sorted(a.items())).encode()).hexdigest()[:12]

    kosul = [yol.split("splits_")[1].rsplit(f"_{mod}_k5", 1)[0]
             for yol in dosyalar]
    tek = sorted(set(imza.values()))
    print(f"    splits dosyasi           : {len(dosyalar)}")
    print(f"    kosullar                 : {', '.join(sorted(kosul))}")
    print(f"    fold atama imzasi        : {tek[0]}"
          f"{'  (hepsi ozdes)' if len(tek) == 1 else ''}")
    if len(tek) > 1:
        print(f"    IMZA SAYISI = {len(tek)} -- FOLD ATAMASI OZDES DEGIL")
        for yol in dosyalar:
            print(f"      {imza[yol]}  {yol}")
        print("    eslestirilmis test varsayimi gecersiz -- durduruldu")
        raise SystemExit(1)

    a = atama[dosyalar[0]]
    n = len(a)
    boy = [sum(1 for v in a.values() if v == k)
           for k in sorted(set(a.values()))]
    print(f"    kohort                   : {n} slayt"
          f"{'' if n == 944 else '  <-- 944 BEKLENIYORDU'}")
    print(f"    fold sayisi              : {len(boy)}")
    print("    fold  n_test  n_train")
    for k, b in enumerate(boy):
        print(f"      {k}    {b:4d}    {n - b:4d}")

    ort_test = sum(boy) / len(boy)
    ort_tr = n - ort_test
    print(f"    ortalama n_test          : {ort_test:.1f}"
          f"   (sabit {N_TEST})")
    print(f"    ortalama n_train         : {ort_tr:.1f}"
          f"   (sabit {N_TRAIN_ANA})")

    nn = 20
    def carpan(nt, ntr):
        return 1 / nn + nt / ntr
    c_sabit = carpan(N_TEST, N_TRAIN_ANA)
    c_ger = carpan(ort_test, ort_tr)
    oran = np.sqrt(c_ger / c_sabit)
    print(f"    carpan sabit / gercek    : {c_sabit:.6f} / {c_ger:.6f}")
    print(f"    SE orani (gercek/sabit)  : {oran:.4f}")

    en_k, en_b = min(boy), max(boy)
    o1 = np.sqrt(carpan(en_k, n - en_k) / c_sabit)
    o2 = np.sqrt(carpan(en_b, n - en_b) / c_sabit)
    print(f"    fold yayilimi            : {en_k} - {en_b} slayt"
          f"   -> SE orani {min(o1, o2):.4f} - {max(o1, o2):.4f}")
    print(f"    -> sabitler {'DOGRULANDI' if abs(oran-1) < 0.02 else 'SAPIYOR'}"
          f"  (ortalamada %{100*abs(oran-1):.1f} fark)")


def tablo1(d, mod):
    """[7] Bildiri Tablo 1'in ilgili kolonu ham veriden yeniden."""
    print("\n[7] TABLO 1 -- ortalama AUC (n=20)")
    sut = 1 if mod == "site" else 0
    print("    kosul          hesaplanan  bildiride   fark")
    kotu = []
    for k in sorted(IDDIA_TABLO1, key=lambda x: -IDDIA_TABLO1[x][sut]):
        if k not in d:
            print(f"    {k:<14} --- veri yok ---")
            kotu.append(k)
            continue
        h = float(np.mean(list(d[k].values())))
        b = IDDIA_TABLO1[k][sut]
        fark = h - b
        im = "ESLESIYOR" if abs(fark) < 5e-5 else "UYUSMUYOR"
        if im == "UYUSMUYOR":
            kotu.append(k)
        print(f"    {k:<14} {h:.6f}    {b:.4f}   "
              f"{fark:+.6f}  {im}")
    print(f"    -> 9 degerin {9-len(kotu)}'u esliyor"
          f"{'' if not kotu else '   SAPAN: ' + ', '.join(kotu)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mod", default="site", choices=["site", "patient"])
    ap.add_argument("--yalniz-fold", action="store_true",
                    help="yalnizca [6] FOLD BOYUTLARI: fold_metrics.csv "
                         "dosyalari eksik olan makinelerde kullanilir")
    a = ap.parse_args()

    if a.yalniz_fold:
        print(f"\n{'='*64}\n{a.mod.upper()} BOLME -- yalnizca fold "
              f"boyutlari\n{'='*64}")
        fold_boyut(a.mod)
        print(f"\n{'='*64}")
        return

    d, yinelenen, atlanan = oku(a.mod)
    iddia = IDDIA[a.mod]

    print(f"\n{'='*64}\n{a.mod.upper()} BOLME -- "
          f"ham fold_metrics.csv dosyalarindan\n{'='*64}")
    butunluk(d, yinelenen, atlanan)

    sirali = sorted(d, key=lambda k: -np.mean(list(d[k].values())))
    A, B = sirali[0], sirali[-1]        # en buyuk kontrast

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

    fold_boyut(a.mod)
    tablo1(d, a.mod)

    print(f"\n{'='*64}")


if __name__ == "__main__":
    main()