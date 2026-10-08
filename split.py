#!/usr/bin/env python3
"""
Hasta seviyesinde k-fold bolme.

Iki mod:
  patient : sinifa gore katmanli, hasta seviyesinde (literatur standardi)
  site    : merkez (TSS) seviyesinde -- hicbir merkez hem egitimde hem
            testte bulunmaz. Merkezler tek sinifli oldugu icin her fold'a
            her siniftan merkez dagitilir (leave-sites-out mantigi).

Tum bilgi barkoddan turetiliyor:
  TCGA-05-4244-01Z-00-DX1  ->  tss=05, hasta=TCGA-05-4244

Kullanim:
  python split.py --emb dinov2_n100 --k 5 --mod patient
  python split.py --emb dinov2_n100 --k 5 --mod site
"""

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
SINIF = {"LUAD": 0, "LUSC": 1}


def kayitlari_topla(emb):
    """Embedding klasorlerini tarar, barkodu ayristirir."""
    kayit = []
    for lab in SINIF:
        d = ROOT / "out" / lab / "embeddings" / emb
        if not d.exists():
            raise SystemExit(f"embedding klasoru yok: {d}")
        for npy in sorted(d.glob("*.npy")):
            wsi = npy.stem
            p = wsi.split("-")
            kayit.append({
                "wsi_id": wsi,
                "patient": "-".join(p[:3]),      # TCGA-05-4244
                "tss": p[1],                     # 05
                "label": lab,
                "y": SINIF[lab],
                "npy": str(npy.relative_to(ROOT)),
                "n_patch": int(np.load(npy, mmap_mode="r").shape[0]),
            })
    return kayit


def bol_patient(kayit, k, seed):
    """Hasta seviyesinde, sinifa gore katmanli."""
    from sklearn.model_selection import StratifiedKFold

    hastalar = sorted({r["patient"]: r["y"] for r in kayit}.items())
    isim = np.array([h for h, _ in hastalar])
    y = np.array([v for _, v in hastalar])

    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    fold = {}
    for f, (_, test_i) in enumerate(skf.split(isim, y)):
        for i in test_i:
            fold[isim[i]] = f
    for r in kayit:
        r["fold"] = fold[r["patient"]]
    return kayit


def bol_site(kayit, k, seed):
    """
    Merkez seviyesinde. Merkezler tek sinifli oldugu icin klasik LOSO
    uygulanamaz: her fold'a her siniftan merkez atanir.

    Acgozlu denge: her sinif icin merkezler buyukten kucuge siralanir,
    siradaki merkez o sinifta en az WSI'ya sahip fold'a verilir.
    """
    rng = np.random.default_rng(seed)
    merkez = defaultdict(lambda: {"n": 0, "y": None})
    for r in kayit:
        merkez[r["tss"]]["n"] += 1
        merkez[r["tss"]]["y"] = r["y"]

    atama = {}
    for sinif in (0, 1):
        liste = [(m, v["n"]) for m, v in merkez.items() if v["y"] == sinif]
        rng.shuffle(liste)                       # esit buyuklukte baglari kir
        liste.sort(key=lambda t: -t[1])
        yuk = [0] * k
        for m, n in liste:
            f = int(np.argmin(yuk))
            atama[m] = f
            yuk[f] += n

    for r in kayit:
        r["fold"] = atama[r["tss"]]
    return kayit


def ozet(kayit, k, mod):
    print(f"\n{'fold':>5} {'WSI':>6} {'LUAD':>6} {'LUSC':>6} {'merkez':>8} "
          f"{'patch':>10}")
    for f in range(k):
        g = [r for r in kayit if r["fold"] == f]
        c = Counter(r["label"] for r in g)
        print(f"{f:>5} {len(g):>6} {c['LUAD']:>6} {c['LUSC']:>6} "
              f"{len({r['tss'] for r in g}):>8} "
              f"{sum(r['n_patch'] for r in g):>10,}")

    # sizinti kontrolu
    hata = 0
    for alan in ("patient", "tss"):
        if alan == "tss" and mod != "site":
            continue
        d = defaultdict(set)
        for r in kayit:
            d[r[alan]].add(r["fold"])
        tasan = {a: v for a, v in d.items() if len(v) > 1}
        if tasan:
            print(f"\n  [!] {len(tasan)} {alan} birden fazla fold'da: "
                  f"{list(tasan)[:5]}")
            hata += 1
        else:
            print(f"  [OK] her {alan} tek fold'da")
    return hata == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb", required=True, help="orn. dinov2_n100")
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--mod", default="patient", choices=["patient", "site"])
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()

    kayit = kayitlari_topla(a.emb)
    c = Counter(r["label"] for r in kayit)
    print(f"{len(kayit)} WSI | LUAD {c['LUAD']} / LUSC {c['LUSC']} | "
          f"{len({r['patient'] for r in kayit})} hasta | "
          f"{len({r['tss'] for r in kayit})} merkez")

    kayit = (bol_site if a.mod == "site" else bol_patient)(kayit, a.k, a.seed)
    temiz = ozet(kayit, a.k, a.mod)

    out = ROOT / "out" / f"splits_{a.emb}_{a.mod}_k{a.k}.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(kayit[0].keys()))
        w.writeheader()
        w.writerows(kayit)
    print(f"\ncikti: {out}")
    if not temiz:
        raise SystemExit("[!] sizinti var, egitime gecme")


if __name__ == "__main__":
    main()