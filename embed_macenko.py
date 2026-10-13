#!/usr/bin/env python3
"""
Macenko boya normalizasyonu ile embedding cikarma.

embed.py'yi degistirmez; ondan yardimci fonksiyonlari ice aktarir ve
yalnizca patch okuma asamasina Macenko normalizasyonu ekler.

Referans goruntu TEK ve SINIF-BAGIMSIZDIR. Sinifa ozgu referans
kullanmak (ornegin LUAD icin lung_aca, LUSC icin lung_scc) etiket
sizintisidir ve bu calismada yapilmamistir.

Kullanim:

  # 1) referans goruntuyu sec (bir kez)
  python embed_macenko.py --ref-sec

  # 2) embedding cikar
  python embed_macenko.py --model dinov2 --label LUAD --limit 100 \
                          --ref macenko_ref.png --batch 64 --workers 8

Cikti:
  out/<LABEL>/embeddings/<model>_mac[_nN]/<wsi_id>.npy
  out/<LABEL>/embeddings/<model>_mac[_nN]/<wsi_id>.json
  out/<LABEL>/embeddings/<model>_mac[_nN]/_log.csv
"""

import argparse
import csv
import json
import random
import sys
import time
from pathlib import Path

import numpy as np
import torch
from PIL import Image
from torch.utils.data import DataLoader, Dataset

from embed import (ROOT, IMAGENET_MEAN, IMAGENET_STD, MODELLER,
                   manifest_oku, patch_yolu, model_yukle)


# ---------------------------------------------------------------- referans

def referans_sec(cikti="macenko_ref.png", seed=42, aday=300):
    """
    Kohorttan tek bir referans patch secer.

    Olcut: ortalama parlaklik 140-190 araliginda (ne cok koyu ne cok
    acik) ve bu arada kontrasti (std) en yuksek olan. Her iki siniftan
    esit sayida aday incelenir; secim deterministiktir.
    """
    rng = random.Random(seed)
    adaylar = []
    for lab in ("LUAD", "LUSC"):
        man = ROOT / "out" / lab / "final_patch_manifest_clean.csv"
        if not man.exists():
            sys.exit(f"manifest yok: {man}")
        rows = list(csv.DictReader(open(man)))
        for r in rng.sample(rows, min(aday, len(rows))):
            adaylar.append((lab, r))

    en_iyi, en_iyi_std = None, -1.0
    for lab, r in adaylar:
        p = (ROOT / "out" / lab / "patches" / r["patient"]
             / r["wsi_id"] / r["patch_file"])
        if not p.exists():
            continue
        a = np.asarray(Image.open(p).convert("RGB"), np.float32)
        if not (140 < a.mean() < 190):
            continue
        s = float(a.std())
        if s > en_iyi_std:
            en_iyi_std, en_iyi = s, (lab, r, p)

    if en_iyi is None:
        sys.exit("uygun referans bulunamadi -- parlaklik araligini genislet")

    lab, r, p = en_iyi
    Image.open(p).convert("RGB").save(cikti)

    bilgi = {
        "referans_dosya": cikti,
        "kaynak_kohort": lab,
        "wsi_id": r["wsi_id"],
        "patch_file": r["patch_file"],
        "kontrast_std": round(en_iyi_std, 2),
        "secim_tohumu": seed,
        "aday_sayisi_per_kohort": aday,
    }
    with open(Path(cikti).with_suffix(".json"), "w") as f:
        json.dump(bilgi, f, indent=2)

    print("referans secildi")
    for k, v in bilgi.items():
        print(f"  {k:<24} {v}")
    return cikti


# ---------------------------------------------------------------- veri

class MacenkoPatchSeti(Dataset):
    """
    Bir slaydin patch'leri, Macenko normalizasyonu uygulanmis.

    Normalizer her DataLoader worker'inda tembel kurulur; torchstain
    nesnesi worker'lar arasi pickle edilmez.
    """

    def __init__(self, yollar, ref_yolu):
        self.yollar = yollar
        self.ref_yolu = ref_yolu
        self._norm = None
        import torchvision.transforms as T
        self.tf = T.Compose([
            T.Resize((224, 224)),
            T.ToTensor(),
            T.Normalize(IMAGENET_MEAN, IMAGENET_STD),
        ])

    def _normalizer(self):
        if self._norm is None:
            import torchstain
            ref = np.asarray(Image.open(self.ref_yolu).convert("RGB"))
            n = torchstain.normalizers.MacenkoNormalizer(backend="numpy")
            n.fit(ref)
            self._norm = n
        return self._norm

    def __len__(self):
        return len(self.yollar)

    def __getitem__(self, i):
        img = Image.open(self.yollar[i]).convert("RGB")
        try:
            cikis = self._normalizer().normalize(I=np.asarray(img),
                                                 stains=False)
            arr = cikis[0] if isinstance(cikis, tuple) else cikis
            arr = np.asarray(arr)
            if arr.dtype != np.uint8:
                arr = np.clip(arr, 0, 255).astype(np.uint8)
            img = Image.fromarray(arr)
        except Exception:
            # SVD yakinsamazsa veya patch neredeyse tek renkse orijinali kullan
            pass
        return self.tf(img)


# ---------------------------------------------------------------- ana

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-sec", action="store_true",
                    help="referans goruntuyu sec ve cik")
    ap.add_argument("--ref", default="macenko_ref.png",
                    help="Macenko referans goruntusu")
    ap.add_argument("--model", choices=list(MODELLER))
    ap.add_argument("--label", help="LUAD veya LUSC")
    ap.add_argument("--out", default=None)
    ap.add_argument("--limit", type=int, default=None,
                    help="slayt basina en fazla N patch")
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--fp32", action="store_true")
    ap.add_argument("--yeniden", action="store_true")
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()

    if a.ref_sec:
        referans_sec(a.ref, a.seed)
        return

    if not (a.model and a.label):
        sys.exit("--model ve --label gerekli (ya da --ref-sec)")
    if not Path(a.ref).exists():
        sys.exit(f"referans bulunamadi: {a.ref}\n"
                 f"once calistir: python embed_macenko.py --ref-sec")

    try:
        import torchstain  # noqa: F401
    except ImportError:
        sys.exit("torchstain kurulu degil:  pip install torchstain")

    out_dir = Path(a.out) if a.out else ROOT / "out" / a.label
    ad = f"{a.model}_mac"
    emb_dir = out_dir / "embeddings" / (f"{ad}_n{a.limit}" if a.limit else ad)
    emb_dir.mkdir(parents=True, exist_ok=True)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda":
        g = torch.cuda.get_device_properties(0)
        print(f"GPU  : {g.name}  {g.total_memory/1e9:.1f} GB")
    else:
        print("[!] GPU bulunamadi -- CPU ile cok yavas olacak")

    slaytlar = manifest_oku(out_dir, a.limit)
    toplam = sum(len(v) for v in slaytlar.values())
    bekl_dim = MODELLER[a.model][0]

    ref_bilgi = Path(a.ref).with_suffix(".json")
    print(f"model    : {a.model}  ({MODELLER[a.model][1]})")
    print(f"normaliz.: Macenko, referans {a.ref}"
          + (f"  ({json.load(open(ref_bilgi))['wsi_id']})"
             if ref_bilgi.exists() else ""))
    print(f"veri     : {len(slaytlar)} WSI / {toplam:,} patch"
          + (f"  (limit {a.limit}/slayt)" if a.limit else ""))
    print(f"cikti    : {emb_dir}")

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
        dl = DataLoader(MacenkoPatchSeti(yollar, a.ref),
                        batch_size=a.batch, shuffle=False,
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
            print(f"  [!] {wsi}: beklenmeyen boyut {emb.shape}")
            hatali += 1
            continue

        np.save(npy, emb)
        with open(emb_dir / f"{wsi}.json", "w") as f:
            json.dump({
                "wsi_id": wsi,
                "patient": satirlar[0]["patient"],
                "model": a.model,
                "normalizasyon": "macenko",
                "referans": a.ref,
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