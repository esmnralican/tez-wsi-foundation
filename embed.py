#!/usr/bin/env python3
"""
Foundation model embedding cikarma.

Slayt basina bir .npy dosyasi uretir: (n_patch, embed_dim), fp16.
Patch sirasi manifest sirasidir ve .json yan dosyasinda kaydedilir --
boylece alt kume ablasyonu embedding'i yeniden cikarmadan yapilabilir
(ilk N satir alinir).

Kullanim:
    python embed.py --model dinov2 --label LUAD --limit 100
    python embed.py --model dinov2 --label LUSC --limit 100
    python embed.py --model uni    --label LUAD            # tum patch'ler
    python embed.py --model virchow --label LUAD --batch 16

Cikti:
    out/<LABEL>/embeddings/<model>/<wsi_id>.npy
    out/<LABEL>/embeddings/<model>/<wsi_id>.json
    out/<LABEL>/embeddings/<model>/_log.csv
"""

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
from PIL import Image
from torch.utils.data import DataLoader, Dataset

# ---------------------------------------------------------------- ayarlar

ROOT = Path(__file__).resolve().parent
IMAGENET_MEAN = (0.485, 0.456, 0.406)
IMAGENET_STD = (0.229, 0.224, 0.225)

MODELLER = {
    # ad       : (embed_dim, aciklama)
    "dinov2": (1024, "DINOv2 ViT-L/14 -- dogal goruntu, genel amacli"),
    "uni": (1024, "UNI ViT-L/16 -- patoloji, 100K WSI"),
    "virchow": (2560, "Virchow ViT-H/14 -- patoloji, 1.5M WSI (CLS + mean)"),
}


# ---------------------------------------------------------------- veri

class PatchSeti(Dataset):
    """Bir slaydin patch'leri. Sira manifest sirasidir."""

    def __init__(self, yollar):
        self.yollar = yollar
        import torchvision.transforms as T
        self.tf = T.Compose([
            T.Resize((224, 224)),          # zaten 224, guvenlik icin
            T.ToTensor(),
            T.Normalize(IMAGENET_MEAN, IMAGENET_STD),
        ])

    def __len__(self):
        return len(self.yollar)

    def __getitem__(self, i):
        img = Image.open(self.yollar[i]).convert("RGB")
        return self.tf(img)


def manifest_oku(out_dir, limit=None, seed=42):
    """
    Manifest'i slayt bazinda gruplar.

    limit verilirse her slayttan ilk N patch alinir. Manifest zaten
    deterministik siradadir (cikarma asamasinda blake2b seed=42 ile
    secilmisti), bu yuzden burada yeniden rastgelelik eklenmez --
    ayni limit her zaman ayni patch'leri verir.
    """
    man = out_dir / "final_patch_manifest_clean.csv"
    if not man.exists():
        sys.exit(f"manifest bulunamadi: {man}")

    slaytlar = {}
    with open(man) as f:
        for r in csv.DictReader(f):
            slaytlar.setdefault(r["wsi_id"], []).append(r)

    if limit:
        for w in slaytlar:
            slaytlar[w] = slaytlar[w][:limit]

    return slaytlar


def patch_yolu(out_dir, satir):
    """patches/<patient>/<wsi_id>/<patch_file>"""
    return out_dir / "patches" / satir["patient"] / satir["wsi_id"] / satir["patch_file"]


# ---------------------------------------------------------------- model

def model_yukle(ad, device):
    if ad == "dinov2":
        m = torch.hub.load("facebookresearch/dinov2", "dinov2_vitl14")
        cikar = lambda x: m(x)                      # CLS token, 1024

    elif ad == "uni":
        import timm
        m = timm.create_model("hf-hub:MahmoodLab/UNI", pretrained=True,
                              init_values=1e-5, dynamic_img_size=True,
                              num_classes=0)
        cikar = lambda x: m(x)                      # pooled CLS, 1024

    elif ad == "virchow":
        import timm
        from timm.layers import SwiGLUPacked
        m = timm.create_model("hf-hub:paige-ai/Virchow", pretrained=True,
                              mlp_layer=SwiGLUPacked, act_layer=torch.nn.SiLU)

        def cikar(x):
            o = m.forward_features(x)               # (B, 257, 1280)
            cls = o[:, 0]                           # (B, 1280)
            patch_tok = o[:, 1:].mean(dim=1)        # (B, 1280)
            return torch.cat([cls, patch_tok], dim=-1)   # (B, 2560)

    else:
        sys.exit(f"bilinmeyen model: {ad}")

    m.eval().to(device)
    for p in m.parameters():
        p.requires_grad_(False)
    return m, cikar


# ---------------------------------------------------------------- ana

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, choices=list(MODELLER))
    ap.add_argument("--label", required=True, help="LUAD veya LUSC")
    ap.add_argument("--out", default=None, help="varsayilan: out/<LABEL>")
    ap.add_argument("--limit", type=int, default=None,
                    help="slayt basina en fazla N patch (pilot icin 100)")
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--fp32", action="store_true", help="fp16 yerine fp32")
    ap.add_argument("--yeniden", action="store_true",
                    help="var olan .npy dosyalarini da yeniden hesapla")
    a = ap.parse_args()

    out_dir = Path(a.out) if a.out else ROOT / "out" / a.label
    emb_dir = out_dir / "embeddings" / a.model
    if a.limit:
        emb_dir = out_dir / "embeddings" / f"{a.model}_n{a.limit}"
    emb_dir.mkdir(parents=True, exist_ok=True)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cpu":
        print("[!] GPU bulunamadi -- CPU ile cok yavas olacak")
    else:
        g = torch.cuda.get_device_properties(0)
        print(f"GPU  : {g.name}  {g.total_memory/1e9:.1f} GB")

    slaytlar = manifest_oku(out_dir, a.limit)
    toplam = sum(len(v) for v in slaytlar.values())
    bekl_dim = MODELLER[a.model][0]

    print(f"model: {a.model}  ({MODELLER[a.model][1]})")
    print(f"veri : {len(slaytlar)} WSI / {toplam:,} patch"
          + (f"  (limit {a.limit}/slayt)" if a.limit else ""))
    print(f"cikti: {emb_dir}")

    model, cikar = model_yukle(a.model, device)
    dtype = torch.float32 if a.fp32 else torch.float16
    kayit_dtype = np.float32 if a.fp32 else np.float16

    log_yolu = emb_dir / "_log.csv"
    yeni_log = not log_yolu.exists()
    log = open(log_yolu, "a", newline="")
    logw = csv.writer(log)
    if yeni_log:
        logw.writerow(["wsi_id", "n_patch", "embed_dim", "saniye", "durum"])

    t_bas = time.time()
    islenen = atlanan = hatali = 0

    for i, (wsi, satirlar) in enumerate(sorted(slaytlar.items()), 1):
        npy = emb_dir / f"{wsi}.npy"
        if npy.exists() and not a.yeniden:
            atlanan += 1
            continue

        yollar = [patch_yolu(out_dir, r) for r in satirlar]
        eksik = [p for p in yollar if not p.exists()]
        if eksik:
            print(f"  [!] {wsi}: {len(eksik)} patch diskte yok, atlandi")
            logw.writerow([wsi, len(satirlar), 0, 0, f"eksik_{len(eksik)}"])
            hatali += 1
            continue

        t0 = time.time()
        dl = DataLoader(PatchSeti(yollar), batch_size=a.batch, shuffle=False,
                        num_workers=a.workers, pin_memory=(device == "cuda"))

        parcalar = []
        try:
            with torch.inference_mode():
                for x in dl:
                    x = x.to(device, non_blocking=True)
                    with torch.autocast(device_type=device, dtype=dtype,
                                        enabled=(device == "cuda" and not a.fp32)):
                        z = cikar(x)
                    parcalar.append(z.float().cpu().numpy())
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            sys.exit(f"\n[!] VRAM yetmedi. --batch {max(4, a.batch//2)} ile dene.")

        emb = np.concatenate(parcalar, axis=0).astype(kayit_dtype)

        if emb.shape != (len(yollar), bekl_dim):
            print(f"  [!] {wsi}: beklenmeyen boyut {emb.shape}, "
                  f"beklenen ({len(yollar)}, {bekl_dim})")
            hatali += 1
            continue

        np.save(npy, emb)
        with open(emb_dir / f"{wsi}.json", "w") as f:
            json.dump({
                "wsi_id": wsi,
                "patient": satirlar[0]["patient"],
                "model": a.model,
                "embed_dim": bekl_dim,
                "n_patch": len(yollar),
                "dtype": str(kayit_dtype.__name__),
                "patch_files": [r["patch_file"] for r in satirlar],
            }, f)

        sn = time.time() - t0
        logw.writerow([wsi, len(yollar), bekl_dim, round(sn, 1), "ok"])
        log.flush()
        islenen += 1

        hiz = len(yollar) / sn
        kalan = len(slaytlar) - i
        print(f"  [{i:>3}/{len(slaytlar)}] {wsi}  {len(yollar):>5} patch  "
              f"{sn:>5.1f}s  {hiz:>5.1f}/s  ~{kalan*sn/60:.0f} dk kaldi")

    log.close()
    dk = (time.time() - t_bas) / 60
    print(f"\nislenen {islenen} | atlanan {atlanan} | hatali {hatali} | "
          f"toplam {dk:.1f} dk")
    print(f"cikti: {emb_dir}")


if __name__ == "__main__":
    main()