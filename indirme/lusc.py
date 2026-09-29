#!/usr/bin/env python3
"""
TCGA-LUSC tanisal WSI indirici.

Ikizi: luad.py (ayni mantik, TCGA-LUSC icin). Birinde degisiklik
yaparsan digerine de uygula.

Iki asamali:
  1. GDC API'den tanisal slayt kayitlarini sorgular -> CSV
  2. Barkod kurallarina uyanlari indirir, MD5 ile dogrular

Eski indirme kodundan farklari:
  - CSV'yi kendisi uretir (elle hazirlamak gerekmiyor)
  - Barkod filtresi indirmeden ONCE uygulanir: yalnizca sample_type=01 ve
    DX1 olanlar iner. 02Z (nuks tumor) ve DX10 gibi dosyalar CSV'ye
    kaydedilir ama indirilmez -- diskte ~10 GB tasarruf, kayit yine tam
  - Kesintiden sonra kaldigi yerden devam eder
  - AWS S3 (hizli) once denenir, basarisizsa GDC API'ye duser

Kullanim:
    python lusc.py                       # indir
    python lusc.py --dry-run             # sadece listele, kac GB
    python lusc.py --verify-only         # mevcut dosyalari MD5 ile dogrula
    python lusc.py --out /mnt/hdd/LUSC   # farkli klasor
"""
import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time

import requests

GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data"
GDC_LEGACY = "https://api.gdc.cancer.gov/legacy/data"
AWS_BASE = "https://tcga-2-open.s3.us-east-1.amazonaws.com"

CHUNK = 1024 * 1024
TIMEOUT = (30, 600)
MAX_RETRY = 12

BARCODE = re.compile(
    r"^(?P<patient>TCGA-(?P<tss>[0-9A-Z]{2})-[0-9A-Z]{4})"
    r"-(?P<sample>\d{2})[A-Z]-\d{2}-(?P<kind>[A-Z]{2})(?P<no>\d+)\."
)


# ------------------------------------------------------------------ sorgu

def gdc_sorgula(project):
    """Projedeki tum tanisal (diagnostic) slaytlari listeler."""
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id",
                                     "value": [project]}},
            {"op": "in", "content": {"field": "data_format",
                                     "value": ["SVS"]}},
            {"op": "in", "content": {"field": "experimental_strategy",
                                     "value": ["Diagnostic Slide"]}},
        ],
    }
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,md5sum,file_size,cases.submitter_id",
        "format": "JSON",
        "size": "5000",
    }
    print(f"GDC sorgulaniyor: {project} ...")
    r = requests.get(GDC_FILES, params=params, timeout=TIMEOUT)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    print(f"  {len(hits)} tanisal slayt kaydi bulundu")

    rows = []
    for h in hits:
        fn = h["file_name"]
        m = BARCODE.match(fn)
        case = (h.get("cases") or [{}])[0].get("submitter_id", "")
        if m:
            sample, kind, no = m["sample"], m["kind"], int(m["no"])
            patient, tss = m["patient"], m["tss"]
        else:
            sample = kind = tss = ""
            no = -1
            patient = case

        # indirme kurali: primary tumor + birinci tanisal slayt
        indir = (sample == "01" and kind == "DX" and no == 1)
        if not indir:
            if not m:
                sebep = "barcode_parse_failed"
            elif kind != "DX":
                sebep = f"not_diagnostic_slide({kind})"
            elif no != 1:
                sebep = f"secondary_diagnostic_slide({kind}{no})"
            else:
                sebep = f"non_primary_sample({sample})"
        else:
            sebep = ""

        rows.append({
            "file_id": h["file_id"], "file_name": fn,
            "case_barcode": case or patient, "patient": patient, "tss": tss,
            "sample": sample, "slide": f"{kind}{no}" if m else "",
            "file_size": int(h["file_size"]), "md5sum": h["md5sum"],
            "aws_s3_url": f"{AWS_BASE}/{h['file_id']}/{fn}",
            "gdc_data_url": f"{GDC_DATA}/{h['file_id']}?related_files=false",
            "download": "1" if indir else "0", "skip_reason": sebep,
        })

    rows.sort(key=lambda r: r["file_name"])
    return rows


# ---------------------------------------------------------------- indirme

def md5_hesapla(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            c = f.read(CHUNK)
            if not c:
                break
            h.update(c)
    return h.hexdigest()


def indir(urls, dest, boyut, md5):
    """
    Devam edebilir indirme.

    Buyuk WSI dosyalari (0.5-2 GB) indirilirken baglanti kopabiliyor.
    Her denemede bastan baslamak yerine HTTP Range ile kalinan yerden
    devam edilir; .part dosyasi denemeler arasinda korunur.

    Kaynak sirasi: AWS S3 (hizli, Registry of Open Data) -> GDC API.
    """
    tmp = dest + ".part"

    for deneme in range(1, MAX_RETRY + 1):
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        url = urls[(deneme - 1) % len(urls)]
        kaynak = "AWS" if "amazonaws" in url else "GDC"

        # kalinan yerden devam
        mevcut = os.path.getsize(tmp) if os.path.exists(tmp) else 0
        if mevcut > boyut:
            os.remove(tmp); mevcut = 0

        headers = {}
        mod = "wb"
        if mevcut:
            headers["Range"] = f"bytes={mevcut}-"
            mod = "ab"

        try:
            t0 = time.time()
            with requests.get(url, stream=True, timeout=TIMEOUT,
                              headers=headers) as r:
                if mevcut and r.status_code == 200:
                    # sunucu Range'i yok saydi -> bastan yaz
                    mevcut, mod = 0, "wb"
                elif mevcut and r.status_code != 206:
                    r.raise_for_status()
                else:
                    r.raise_for_status()

                inen = mevcut
                son = t0
                with open(tmp, mod) as f:
                    for c in r.iter_content(chunk_size=CHUNK):
                        if not c:
                            continue
                        f.write(c)
                        inen += len(c)
                        if time.time() - son > 2:
                            son = time.time()
                            pc = inen / boyut * 100
                            hz = (inen - mevcut) / (son - t0) / 1e6
                            kalan = (boyut - inen) / max(hz * 1e6, 1) / 60
                            print(f"\r    {kaynak} %{pc:5.1f}  {hz:5.1f} MB/s"
                                  f"  ~{kalan:.0f} dk kaldi   ",
                                  end="", flush=True)
            print("\r" + " " * 56 + "\r", end="")

            gercek = os.path.getsize(tmp)
            if gercek != boyut:
                eksik = boyut - gercek
                print(f"    [{deneme}/{MAX_RETRY}] {kaynak}: baglanti koptu, "
                      f"{gercek/1e6:.0f}/{boyut/1e6:.0f} MB "
                      f"({eksik/1e6:.0f} MB eksik) -- devam edilecek")
                time.sleep(3)
                continue          # .part SILINMIYOR, kalinan yerden devam

            if md5_hesapla(tmp).lower() != md5.lower():
                print(f"    [{deneme}/{MAX_RETRY}] {kaynak}: MD5 uyusmadi, "
                      f"bastan indirilecek")
                os.remove(tmp); time.sleep(3); continue

            os.replace(tmp, dest)
            print(f"    OK ({kaynak}, {boyut/1e6:.0f} MB, "
                  f"{time.time()-t0:.0f}s)")
            return True

        except Exception as e:
            print(f"    [{deneme}/{MAX_RETRY}] {kaynak}: "
                  f"{type(e).__name__}: {str(e)[:60]}")
            time.sleep(5)         # .part korunuyor

    return False


# -------------------------------------------------------------------- main

PROJE = "TCGA-LUSC"
VARSAYILAN_KLASOR = "~/tcga_data/LUSC"


def main():
    ap = argparse.ArgumentParser(description=f"{PROJE} WSI indirme")
    ap.add_argument("--project", default=PROJE)
    ap.add_argument("--out", default=VARSAYILAN_KLASOR)
    ap.add_argument("--refresh", action="store_true",
                    help="CSV'yi GDC'den yeniden sorgula")
    ap.add_argument("--dry-run", action="store_true",
                    help="sadece listele, indirme")
    ap.add_argument("--verify-only", action="store_true",
                    help="mevcut dosyalari dogrula, eksikleri indirme")
    args = ap.parse_args()

    out = os.path.expanduser(args.out)
    label = args.project.split("-")[1]
    os.makedirs(out, exist_ok=True)
    csv_path = os.path.join(out, f"tcga_{label.lower()}_diagnostic_slides.csv")
    dl_dir = os.path.join(out, f"tcga_{label.lower()}_wsi_dx1_only_files")

    # --- CSV ---
    if args.refresh or not os.path.exists(csv_path):
        rows = gdc_sorgula(args.project)
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"  CSV: {csv_path}")
    else:
        rows = list(csv.DictReader(open(csv_path)))
        print(f"CSV okundu: {csv_path} ({len(rows)} kayit)")
        print("  (yeniden sorgulamak icin --refresh)")

    inecek = [r for r in rows if r["download"] == "1"]
    atlanan = [r for r in rows if r["download"] != "1"]
    toplam = sum(int(r["file_size"]) for r in inecek)
    hasta = len({r["patient"] for r in inecek})

    print(f"\n{'='*66}\n{args.project}\n{'='*66}")
    print(f"  toplam tanisal slayt : {len(rows)}")
    print(f"  indirilecek (01+DX1) : {len(inecek)}  /  {hasta} hasta")
    print(f"  toplam boyut         : {toplam/1e9:.1f} GB")
    if atlanan:
        from collections import Counter
        print(f"\n  indirilmeyecek ({len(atlanan)}):")
        for s, n in Counter(r["skip_reason"] for r in atlanan).most_common():
            print(f"      {n:>4}  {s}")

    # hasta basina birden fazla DX1 var mi
    from collections import Counter
    cok = [p for p, n in Counter(r["patient"] for r in inecek).items() if n > 1]
    if cok:
        print(f"\n  [!] {len(cok)} hastada birden fazla DX1: {cok[:5]}")

    if args.dry_run:
        print("\n(dry-run: indirme yapilmadi)")
        return

    # --- disk kontrolu ---
    st = os.statvfs(out)
    bos = st.f_bavail * st.f_frsize
    print(f"\n  bos disk alani: {bos/1e9:.0f} GB")
    if bos < toplam * 1.05:
        sys.exit(f"[DUR] yetersiz alan: {toplam/1e9:.0f} GB gerekiyor")

    # --- indirme ---
    print(f"\n{'='*66}\nINDIRME\n{'='*66}")
    ok = eksik = basarisiz = bozuk = 0
    t0 = time.time()

    for i, r in enumerate(inecek, 1):
        dest = os.path.join(dl_dir, r["case_barcode"], r["file_name"])
        boyut = int(r["file_size"])
        print(f"[{i}/{len(inecek)}] {r['file_name'][:52]}")

        if os.path.exists(dest):
            if os.path.getsize(dest) == boyut:
                if md5_hesapla(dest).lower() == r["md5sum"].lower():
                    print("    zaten var, dogrulandi")
                    ok += 1
                    continue
                print("    MD5 hatali, siliniyor")
                bozuk += 1
            else:
                print("    eksik/bozuk, siliniyor")
                bozuk += 1
            os.remove(dest)

        if args.verify_only:
            eksik += 1
            continue

        if indir([r["aws_s3_url"], r["gdc_data_url"]], dest, boyut,
                 r["md5sum"]):
            ok += 1
        else:
            basarisiz += 1
            print(f"    [!] BASARISIZ: {r['file_name']}")

        if i % 25 == 0:
            gecen = time.time() - t0
            kalan = (len(inecek) - i) * gecen / i / 3600
            print(f"\n  --- {i}/{len(inecek)} | gecen {gecen/3600:.1f} sa | "
                  f"tahmini kalan {kalan:.1f} sa ---\n")

    print(f"\n{'='*66}\nOZET\n{'='*66}")
    print(f"  dogrulanan  : {ok}")
    print(f"  bozuk bulup yeniden inen : {bozuk}")
    if args.verify_only:
        print(f"  eksik       : {eksik}")
    print(f"  basarisiz   : {basarisiz}")
    print(f"  sure        : {(time.time()-t0)/3600:.1f} saat")

    if basarisiz == 0 and not args.verify_only:
        print(f"\n  Tum dosyalar indirildi ve dogrulandi.")
        print(f"\n  siradaki adim:")
        print(f"    python run.py cohort --data {out} --label {label} \\")
        print(f"        --primary tcga_{label.lower()}_wsi_dx1_only_files \\")
        print(f"        --out ~/tez/out/{label}")
    elif basarisiz:
        print(f"\n  [!] {basarisiz} dosya inmedi. Script'i tekrar calistir.")


if __name__ == "__main__":
    main()