"""
KOD NEDEN VAR
Tek bir soruyu doğrudan sormak: embedding'den hangi hastaneden geldiği anlaşılıyor mu?

Şimdiye kadar bunu dolaylı ölçtük — iki bölme arasındaki farka bakarak. Bu script doğrudan soruyor: aynı vektörlerden sınıf yerine TSS tahmin et, doğruluğa bak.

Sonuç	Anlamı
Şans seviyesinde	İmza yok, endişe yersizdi
%50-70	İmza var ama zayıf
%80+	Merkez bilgisi embedding'lere kodlanmış

Literatürde foundation model'lerde %80-99 çıkmış.

Üç faydası:

1. Tek sayı. Tezde "merkez imzası embedding'lerde %X doğrulukla tespit edilebilmektedir" diye yazarsın. Dolaylı çıkarım değil, doğrudan ölçüm.

2. Model karşılaştırması. İleride UNI ve Virchow'da da çalıştırırsın. Literatür "daha iyi model = daha güçlü imza" diyor; kendi verinde test etmiş olursun. Bu tezin özgün tarafı olabilir.

3. Ablasyonların ölçüsü. Gri tonlamada bu sayı düşüyor mu? Macenko'da? "İşe yaradı mı" sorusunun cevabı bu sayının değişimi olur.

Çözüm değil, teşhis. Ama olmadan diğer müdahalelerin etkisini ölçemezsin

"""#!/usr/bin/env python3


"""
Site prediction baseline.

Ayni embedding'lerden sinif yerine kaynak merkez (TSS) tahmin edilir.
Yuksek dogruluk, merkez bilgisinin embedding'lerde kodlandigini
dogrudan gosterir (Howard et al. 2021; arXiv:2411.05489).

Iki seviye:
  patch  : her patch bir ornek (literaturdeki standart)
  slide  : slayt ortalamasi bir ornek

Bolme hasta seviyesindedir -- ayni slaytin patch'leri hem egitimde
hem testte bulunmaz.

Kullanim:
  python site_predict.py --emb dinov2_n100
  python site_predict.py --emb dinov2_n100 --min-slayt 10 --patch-basina 20
"""

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

ROOT = Path(__file__).resolve().parent

def yukle_duzgun(emb, patch_basina, seed):
    rng = np.random.default_rng(seed)
    P = {"X": [], "tss": [], "sinif": [], "hasta": []}
    S = {"X": [], "tss": [], "sinif": [], "hasta": []}

    for lab in ("LUAD", "LUSC"):
        d = ROOT / "out" / lab / "embeddings" / emb
        if not d.exists():
            raise SystemExit(f"yok: {d}")
        for npy in sorted(d.glob("*.npy")):
            wsi = npy.stem
            t = wsi.split("-")[1]
            h = "-".join(wsi.split("-")[:3])
            a = np.load(npy).astype(np.float32)

            S["X"].append(a.mean(axis=0)); S["tss"].append(t)
            S["sinif"].append(lab); S["hasta"].append(h)

            n = min(patch_basina, len(a))
            for i in rng.choice(len(a), n, replace=False):
                P["X"].append(a[i]); P["tss"].append(t)
                P["sinif"].append(lab); P["hasta"].append(h)

    for D in (P, S):
        for k in D:
            D[k] = np.array(D[k])
    return P, S


def sans(y):
    """En sik sinifin orani -- taban cizgisi."""
    return Counter(y).most_common(1)[0][1] / len(y)


def olc(X, y, gruplar, ad, k=5):
    if len(set(y)) < 2:
        print(f"  {ad:<34} tek sinif, atlandi")
        return None
    kf = GroupKFold(n_splits=min(k, len(set(gruplar))))
    model = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=1000, n_jobs=-1))
    tah = cross_val_predict(model, X, y, cv=kf, groups=gruplar, n_jobs=1)
    dog = (tah == y).mean()
    s = sans(y)
    print(f"  {ad:<34} {dog:>6.3f}   (sans {s:.3f})   "
          f"{len(set(y))} sinif, n={len(y):,}")
    return dog


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb", required=True)
    ap.add_argument("--patch-basina", type=int, default=20)
    ap.add_argument("--min-slayt", type=int, default=10,
                    help="bu sayidan az slaydi olan merkezler elenir")
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()

    P, S = yukle_duzgun(a.emb, a.patch_basina, a.seed)
    print(f"embedding: {a.emb}")
    print(f"slayt {len(S['X'])} | patch {len(P['X']):,} "
          f"(slayt basina {a.patch_basina}) | boyut {P['X'].shape[1]}")

    # kucuk merkezleri ele
    say = Counter(S["tss"])
    gecerli = {t for t, n in say.items() if n >= a.min_slayt}
    print(f"merkez: {len(say)} toplam, {len(gecerli)} tanesi "
          f">= {a.min_slayt} slayt\n")

    mp = np.isin(P["tss"], list(gecerli))
    ms = np.isin(S["tss"], list(gecerli))

    print("MERKEZ (TSS) TAHMINI  -- yuksek deger = imza var")
    olc(P["X"][mp], P["tss"][mp], P["hasta"][mp], "patch seviyesi")
    olc(S["X"][ms], S["tss"][ms], S["hasta"][ms], "slayt seviyesi (ortalama)")

    print("\nSINIF TAHMINI  -- karsilastirma icin")
    olc(P["X"], P["sinif"], P["hasta"], "patch seviyesi")
    olc(S["X"], S["sinif"], S["hasta"], "slayt seviyesi (ortalama)")


if __name__ == "__main__":
    main()