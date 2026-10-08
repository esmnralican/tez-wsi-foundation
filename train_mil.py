#!/usr/bin/env python3
"""
ABMIL (Ilse et al. 2018, gated attention) ile slayt seviyesinde
LUAD / LUSC siniflandirmasi.

Her slayt bir bag, her patch embedding'i bir instance.
Etiket slayt (= hasta) seviyesindedir.

Kullanim:
  python train_mil.py --splits out/splits_dinov2_n100_patient_k5.csv
  python train_mil.py --splits out/splits_dinov2_n100_site_k5.csv

Cikti:
  out/mil/<tag>/fold_metrics.csv
  out/mil/<tag>/predictions.csv
  out/mil/<tag>/config.json
"""

#!/usr/bin/env python3
"""
ABMIL (Ilse et al. 2018, gated attention) ile slayt seviyesinde
LUAD / LUSC siniflandirmasi.

Ic dogrulama kumesi, dis bolmeyle ayni seviyede ayrilir:
  patient bolme  -> ic dogrulama hasta seviyesinde
  site bolme     -> ic dogrulama MERKEZ seviyesinde

Aksi halde erken durdurma, egitim merkezlerinde iyi olan (yani merkez
kestirmesini en iyi kullanan) checkpoint'i secer ve yabanci merkezlerde
test edilir -- sistematik olarak yanli model secimi.

Kullanim:
  python train_mil.py --splits out/splits_dinov2_patient_k5.csv
  python train_mil.py --splits out/splits_dinov2_site_k5.csv
  python train_mil.py --splits ... --max-patch 500 --tag dinov2_site_n500
"""

import argparse
import csv
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                             f1_score, roc_auc_score)

ROOT = Path(__file__).resolve().parent


# ---------------------------------------------------------------- model

class ABMIL(nn.Module):
    """Gated attention MIL. Ilse, Tomczak & Welling, ICML 2018."""

    def __init__(self, in_dim=1024, hid=512, att=384, dropout=0.25):
        super().__init__()
        self.govde = nn.Sequential(
            nn.Linear(in_dim, hid), nn.ReLU(), nn.Dropout(dropout))
        self.att_V = nn.Sequential(nn.Linear(hid, att), nn.Tanh())
        self.att_U = nn.Sequential(nn.Linear(hid, att), nn.Sigmoid())
        self.att_w = nn.Linear(att, 1)
        self.bas = nn.Linear(hid, 1)

    def forward(self, x):                      # x: (n_patch, in_dim)
        h = self.govde(x)
        a = self.att_w(self.att_V(h) * self.att_U(h))
        a = torch.softmax(a, dim=0)
        m = (a * h).sum(dim=0)
        return self.bas(m).squeeze(), a.squeeze(1)


# ---------------------------------------------------------------- veri

def bag_yukle(satirlar, max_patch=None):
    torbalar = []
    for r in satirlar:
        x = np.load(ROOT / r["npy"]).astype(np.float32)
        if max_patch:
            x = x[:max_patch]
        torbalar.append({
            "wsi_id": r["wsi_id"], "patient": r["patient"], "tss": r["tss"],
            "label": r["label"], "y": float(r["y"]),
            "x": torch.from_numpy(x),
        })
    return torbalar


def ic_dogrulama_ayir(kalan, mod, oran, rng):
    """
    Egitim kumesinden ic dogrulama ayirir.

    mod='patient' : hasta seviyesinde
    mod='site'    : MERKEZ seviyesinde -- ic dogrulama merkezleri egitimde
                    bulunmaz, boylece erken durdurma merkezler arasi
                    genellemeye gore secim yapar
    """
    if mod == "patient":
        birim = sorted({b["patient"] for b in kalan})
        rng.shuffle(birim)
        n = max(1, int(len(birim) * oran))
        secili = set(birim[:n])
        anahtar = "patient"
    else:
        gruplar = defaultdict(list)
        for b in kalan:
            gruplar[b["tss"]].append(b)
        secili = set()
        for sinif in (0.0, 1.0):
            ms = [(t, len(v)) for t, v in gruplar.items() if v[0]["y"] == sinif]
            rng.shuffle(ms)
            hedef = oran * sum(n for _, n in ms)
            top = 0
            for t, n in ms:
                if top >= hedef:
                    break
                secili.add(t)
                top += n
        anahtar = "tss"

    val = [b for b in kalan if b[anahtar] in secili]
    tren = [b for b in kalan if b[anahtar] not in secili]

    if not val or len({b["y"] for b in val}) < 2 or not tren:
        raise SystemExit("[!] ic dogrulama kumesi kurulamadi "
                         "(bos ya da tek sinifli) -- --val-oran degerini artir")
    return tren, val


def degerlendir(model, torbalar, device):
    model.eval()
    skor, gercek = [], []
    with torch.inference_mode():
        for b in torbalar:
            logit, _ = model(b["x"].to(device))
            skor.append(torch.sigmoid(logit).item())
            gercek.append(b["y"])
    skor, gercek = np.array(skor), np.array(gercek)
    tahmin = (skor >= 0.5).astype(int)
    return {
        "auc": roc_auc_score(gercek, skor) if len(set(gercek)) > 1 else float("nan"),
        "acc": accuracy_score(gercek, tahmin),
        "bal_acc": balanced_accuracy_score(gercek, tahmin),
        "f1": f1_score(gercek, tahmin, zero_division=0),
    }, skor


def egit(tren, val, device, epoch, lr, wd, accum, sabir, tohum, sessiz):
    torch.manual_seed(tohum)
    model = ABMIL(in_dim=tren[0]["x"].shape[1]).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=wd)
    kayip_fn = nn.BCEWithLogitsLoss()
    rng = np.random.default_rng(tohum)

    en_iyi, en_iyi_durum, bekleme = -1.0, None, 0
    for ep in range(1, epoch + 1):
        model.train()
        sira = rng.permutation(len(tren))
        opt.zero_grad()
        toplam = 0.0
        for adim, i in enumerate(sira, 1):
            b = tren[i]
            logit, _ = model(b["x"].to(device))
            kayip = kayip_fn(logit, torch.tensor(b["y"], device=device))
            (kayip / accum).backward()
            toplam += kayip.item()
            if adim % accum == 0 or adim == len(sira):
                opt.step()
                opt.zero_grad()

        m, _ = degerlendir(model, val, device)
        if m["auc"] > en_iyi:
            en_iyi = m["auc"]
            en_iyi_durum = {k: v.detach().clone()
                            for k, v in model.state_dict().items()}
            bekleme = 0
        else:
            bekleme += 1

        if not sessiz:
            print(f"      ep {ep:>3}  kayip {toplam/len(sira):.4f}  "
                  f"val_auc {m['auc']:.4f}" + ("  *" if bekleme == 0 else ""))
        if bekleme >= sabir:
            if not sessiz:
                print(f"      erken durdurma (en iyi val_auc {en_iyi:.4f})")
            break

    model.load_state_dict(en_iyi_durum)
    return model, en_iyi


# ---------------------------------------------------------------- ana

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", required=True)
    ap.add_argument("--max-patch", type=int, default=None,
                    help="bag basina ilk N patch (ablasyon icin)")
    ap.add_argument("--ic-mod", default="auto",
                    choices=["auto", "patient", "site"],
                    help="ic dogrulama hangi seviyede ayrilsin")
    ap.add_argument("--epoch", type=int, default=30)
    ap.add_argument("--lr", type=float, default=2e-4)
    ap.add_argument("--wd", type=float, default=1e-5)
    ap.add_argument("--accum", type=int, default=8)
    ap.add_argument("--sabir", type=int, default=7)
    ap.add_argument("--val-oran", type=float, default=0.15)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--sessiz", action="store_true",
                    help="epoch satirlarini basma (ablasyon icin)")
    a = ap.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    satirlar = list(csv.DictReader(open(a.splits)))
    for r in satirlar:
        r["y"], r["fold"] = int(r["y"]), int(r["fold"])
    k = max(r["fold"] for r in satirlar) + 1

    ic_mod = a.ic_mod
    if ic_mod == "auto":
        ic_mod = "site" if "_site_" in Path(a.splits).name else "patient"

    tag = a.tag or Path(a.splits).stem.replace("splits_", "")
    if a.max_patch and not a.tag:
        tag += f"_n{a.max_patch}"
    out = ROOT / "out" / "mil" / tag
    out.mkdir(parents=True, exist_ok=True)

    print(f"splits  : {a.splits}")
    print(f"ic dogr.: {ic_mod} seviyesinde ayriliyor")
    print(f"veri    : {len(satirlar)} WSI | {k} fold | cihaz {device}"
          + (f" | bag tavani {a.max_patch}" if a.max_patch else ""))
    print(f"cikti   : {out}\n")

    tum = bag_yukle(satirlar, a.max_patch)
    rng = np.random.default_rng(a.seed)
    metrikler, tahminler = [], []
    t0 = time.time()

    for f in range(k):
        test = [b for b, r in zip(tum, satirlar) if r["fold"] == f]
        kalan = [b for b, r in zip(tum, satirlar) if r["fold"] != f]
        tren, val = ic_dogrulama_ayir(kalan, ic_mod, a.val_oran, rng)

        print(f"  fold {f}:  tren {len(tren)} ({len({b['tss'] for b in tren})} merkez)"
              f"  val {len(val)} ({len({b['tss'] for b in val})} merkez)"
              f"  test {len(test)} ({len({b['tss'] for b in test})} merkez)")

        model, en_iyi = egit(tren, val, device, a.epoch, a.lr, a.wd,
                             a.accum, a.sabir, a.seed + f, a.sessiz)

        m, skor = degerlendir(model, test, device)
        m["fold"], m["val_auc"] = f, en_iyi
        metrikler.append(m)
        for b, s in zip(test, skor):
            tahminler.append({"fold": f, "wsi_id": b["wsi_id"],
                              "patient": b["patient"], "tss": b["tss"],
                              "gercek": b["label"],
                              "skor_LUSC": round(float(s), 4),
                              "tahmin": "LUSC" if s >= 0.5 else "LUAD"})
        print(f"      TEST  auc {m['auc']:.4f}  acc {m['acc']:.4f}  "
              f"bal_acc {m['bal_acc']:.4f}  f1 {m['f1']:.4f}\n")

    print("=" * 62)
    print(f"{'fold':>5} {'auc':>8} {'acc':>8} {'bal_acc':>9} {'f1':>8}")
    for m in metrikler:
        print(f"{m['fold']:>5} {m['auc']:>8.4f} {m['acc']:>8.4f} "
              f"{m['bal_acc']:>9.4f} {m['f1']:>8.4f}")
    alan = ("auc", "acc", "bal_acc", "f1")
    ort = {x: float(np.mean([m[x] for m in metrikler])) for x in alan}
    std = {x: float(np.std([m[x] for m in metrikler])) for x in alan}
    print("-" * 62)
    print(f"{'ort':>5} {ort['auc']:>8.4f} {ort['acc']:>8.4f} "
          f"{ort['bal_acc']:>9.4f} {ort['f1']:>8.4f}")
    print(f"{'std':>5} {std['auc']:>8.4f} {std['acc']:>8.4f} "
          f"{std['bal_acc']:>9.4f} {std['f1']:>8.4f}")
    print("=" * 62)
    print(f"sure: {(time.time()-t0)/60:.1f} dk")

    with open(out / "fold_metrics.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["fold", "auc", "acc", "bal_acc",
                                           "f1", "val_auc"])
        w.writeheader(); w.writerows(metrikler)
    with open(out / "predictions.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(tahminler[0].keys()))
        w.writeheader(); w.writerows(tahminler)
    with open(out / "config.json", "w") as fh:
        json.dump({**vars(a), "ic_mod": ic_mod,
                   "ortalama": ort, "std": std}, fh, indent=2)


if __name__ == "__main__":
    main()