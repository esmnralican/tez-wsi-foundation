#!/usr/bin/env python3
"""
ABMIL dikkat agirliklarinin gorsellestirilmesi.

Bir fold uzerinde model egitir, o fold'un test slaytlarinda dikkat
agirliklarini cikarir ve en yuksek / en dusuk agirlikli patch'leri yan
yana koyarak kontak sayfasi uretir.

Amac: modelin karara katkisi en yuksek bolgelerin morfolojik yapi mi
(glandüler organizasyon, keratinizasyon, hucreler arasi kopru) yoksa
edinim artefakti mi (murekkep kenari, arka plan, kesit hatasi)
gosterdigini nitel olarak degerlendirmek.

Ayrica dikkatin ne kadar yogunlastigini olcer: toplam dikkatin yuzde
kaci en yuksek %1'lik patch dilimindedir.

Kullanim:
  python dikkat.py --splits out/splits_dinov2_site_k5.csv --fold 0
  python dikkat.py --splits out/splits_dinov2_site_k5.csv --fold 0 \
                   --n-slayt 8 --topk 12

Cikti:
  out/dikkat/<tag>_fold<k>/
    <wsi_id>_dogru.jpg      kontak sayfasi (ust: en yuksek, alt: en dusuk)
    dikkat_ozet.csv         slayt basina yogunlasma olculeri
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import torch
from PIL import Image, ImageDraw

from train_mil import ABMIL, bag_yukle, ic_dogrulama_ayir, egit, esik_sec

ROOT = Path(__file__).resolve().parent
PATCH_PX = 224
BOSLUK = 6
ETIKET_Y = 26


# ---------------------------------------------------------------- yardimci

def json_oku(npy_yolu):
    """Embedding'in yanindaki .json: patch sirasi ve hasta kimligi."""
    j = Path(npy_yolu).with_suffix(".json")
    with open(j) as f:
        return json.load(f)


def patch_yolu(label, patient, wsi_id, patch_file):
    return ROOT / "out" / label / "patches" / patient / wsi_id / patch_file


def kontak_sayfasi(ust, alt, baslik, cikti, sutun):
    """
    Iki satir grubundan kontak sayfasi kurar.
    ust/alt: [(PIL.Image, "etiket"), ...]
    """
    n = max(len(ust), len(alt))
    w = n * PATCH_PX + (n + 1) * BOSLUK
    satir_h = PATCH_PX + ETIKET_Y + BOSLUK
    h = 2 * satir_h + 3 * BOSLUK + 34

    tuval = Image.new("RGB", (w, h), (255, 255, 255))
    ciz = ImageDraw.Draw(tuval)
    ciz.text((BOSLUK, 6), baslik, fill=(0, 0, 0))

    for sira, (grup, ad) in enumerate(((ust, "EN YUKSEK DIKKAT"),
                                       (alt, "EN DUSUK DIKKAT"))):
        y0 = 34 + sira * (satir_h + BOSLUK)
        ciz.text((BOSLUK, y0 - 2), ad, fill=(60, 60, 60))
        for i, (img, et) in enumerate(grup):
            x = BOSLUK + i * (PATCH_PX + BOSLUK)
            tuval.paste(img, (x, y0 + 14))
            ciz.text((x + 2, y0 + 14 + PATCH_PX + 4), et, fill=(0, 0, 0))

    tuval.save(cikti, quality=92)


def yogunlasma(a):
    """
    Dikkat dagiliminin yogunlasma olculeri.
    a: (n,) softmax ciktisi, toplami 1.
    """
    n = len(a)
    s = np.sort(a)[::-1]
    ust1 = s[:max(1, n // 100)].sum()
    ust5 = s[:max(1, n // 20)].sum()
    # etkin patch sayisi (perplexity): duzgun dagilimda n, tek patch'te 1
    p = a[a > 0]
    etkin = float(np.exp(-(p * np.log(p)).sum()))
    return ust1, ust5, etkin


# ---------------------------------------------------------------- ana

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", required=True)
    ap.add_argument("--fold", type=int, default=0)
    ap.add_argument("--n-slayt", type=int, default=6,
                    help="gorsellestirilecek test slaydi sayisi")
    ap.add_argument("--topk", type=int, default=10,
                    help="satir basina patch sayisi")
    ap.add_argument("--yalnizca-dogru", action="store_true",
                    help="sadece dogru siniflandirilan slaytlari goster")
    ap.add_argument("--epoch", type=int, default=30)
    ap.add_argument("--lr", type=float, default=2e-4)
    ap.add_argument("--wd", type=float, default=1e-5)
    ap.add_argument("--accum", type=int, default=8)
    ap.add_argument("--sabir", type=int, default=7)
    ap.add_argument("--val-oran", type=float, default=0.15)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    satirlar = list(csv.DictReader(open(a.splits)))
    for r in satirlar:
        r["y"], r["fold"] = int(r["y"]), int(r["fold"])

    ad = Path(a.splits).stem.replace("splits_", "")
    ic_mod = "site" if "_site_" in Path(a.splits).name else "patient"
    out = ROOT / "out" / "dikkat" / f"{ad}_fold{a.fold}"
    out.mkdir(parents=True, exist_ok=True)

    print(f"splits : {a.splits}")
    print(f"fold   : {a.fold} | ic dogrulama: {ic_mod} | cihaz {device}")
    print(f"cikti  : {out}\n")

    tum = bag_yukle(satirlar, None)
    test_i = [i for i, r in enumerate(satirlar) if r["fold"] == a.fold]
    kalan = [tum[i] for i, r in enumerate(satirlar) if r["fold"] != a.fold]
    test = [tum[i] for i in test_i]

    rng = np.random.default_rng(a.seed)
    tren, val = ic_dogrulama_ayir(kalan, ic_mod, a.val_oran, rng)
    print(f"tren {len(tren)}  val {len(val)}  test {len(test)}")

    print("\nmodel egitiliyor...")
    model, en_iyi = egit(tren, val, device, a.epoch, a.lr, a.wd,
                         a.accum, a.sabir, a.seed + a.fold, sessiz=True)
    esik = esik_sec(model, val, device)
    print(f"  en iyi val_auc {en_iyi:.4f} | esik {esik:.3f}\n")

    # --------------------------------------------------- dikkat cikarimi
    model.eval()
    kayit = []
    secilen = []

    with torch.inference_mode():
        for b, i in zip(test, test_i):
            logit, att = model(b["x"].to(device))
            skor = torch.sigmoid(logit).item()
            att = att.cpu().numpy()
            tahmin = "LUSC" if skor >= esik else "LUAD"
            dogru = tahmin == b["label"]

            u1, u5, etkin = yogunlasma(att)
            kayit.append({
                "wsi_id": b["wsi_id"], "tss": b["tss"], "gercek": b["label"],
                "tahmin": tahmin, "dogru": int(dogru),
                "skor": round(skor, 4), "n_patch": len(att),
                "ust1_pct_dikkat": round(float(u1), 4),
                "ust5_pct_dikkat": round(float(u5), 4),
                "etkin_patch": round(etkin, 1),
            })
            if (not a.yalnizca_dogru) or dogru:
                secilen.append((b, satirlar[i], att, skor, tahmin, dogru))

    with open(out / "dikkat_ozet.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(kayit[0].keys()))
        w.writeheader(); w.writerows(kayit)

    # ozet istatistik
    u1 = np.array([k["ust1_pct_dikkat"] for k in kayit])
    etk = np.array([k["etkin_patch"] for k in kayit])
    npatch = np.array([k["n_patch"] for k in kayit])
    print("DIKKAT YOGUNLASMASI (test slaytlari)")
    print(f"  en yuksek %1'lik patch dilimi toplam dikkatin: "
          f"ortanca %{np.median(u1)*100:.1f}  (duzgun dagilimda %1)")
    print(f"  etkin patch sayisi: ortanca {np.median(etk):.0f} / "
          f"{int(np.median(npatch))}")
    print(f"  dogru siniflandirilan: {sum(k['dogru'] for k in kayit)}"
          f"/{len(kayit)}\n")

    # --------------------------------------------------- gorsellestirme
    # once en yuksek dikkat yogunlasmasina sahip slaytlari sec
    secilen.sort(key=lambda t: -yogunlasma(t[2])[0])
    secilen = secilen[:a.n_slayt]

    print("kontak sayfalari uretiliyor...")
    for b, satir, att, skor, tahmin, dogru in secilen:
        meta = json_oku(ROOT / satir["npy"])
        dosyalar = meta["patch_files"]
        if len(dosyalar) != len(att):
            print(f"  [!] {b['wsi_id']}: json {len(dosyalar)} / "
                  f"dikkat {len(att)} uyusmuyor, atlandi")
            continue

        sira = np.argsort(att)[::-1]
        ust_i, alt_i = sira[:a.topk], sira[-a.topk:]

        def grup(idx):
            g = []
            for j in idx:
                p = patch_yolu(b["label"], b["patient"], b["wsi_id"],
                               dosyalar[j])
                if not p.exists():
                    continue
                img = Image.open(p).convert("RGB")
                if img.size != (PATCH_PX, PATCH_PX):
                    img = img.resize((PATCH_PX, PATCH_PX))
                g.append((img, f"{att[j]*len(att):.1f}x"))
            return g

        u1, u5, etkin = yogunlasma(att)
        baslik = (f"{b['wsi_id']}  TSS {b['tss']}  gercek {b['label']}  "
                  f"tahmin {tahmin} ({skor:.3f})  |  n={len(att)}  "
                  f"ust%1 dikkat %{u1*100:.1f}  etkin patch {etkin:.0f}")
        dosya = out / f"{b['wsi_id']}_{'dogru' if dogru else 'yanlis'}.jpg"
        kontak_sayfasi(grup(ust_i), grup(alt_i), baslik, dosya, a.topk)
        print(f"  {dosya.name}")

    print(f"\nEtiketler ortalama dikkate gore normalize edilmistir: "
          f"'3.0x' o patch'in duzgun dagilimdan 3 kat fazla agirlik "
          f"aldigini gosterir.")
    print(f"cikti: {out}")


if __name__ == "__main__":
    main()