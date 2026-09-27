#!/usr/bin/env python3
"""
Tez pipeline giris noktasi.

Her asama bir alt komut. Yeni asama eklendiginde yeni dosya degil,
tez/ altinda ilgili modul buyur ve buraya bir alt komut eklenir.

    python run.py cohort  --data ~/tcga_data/LUAD --label LUAD \
                          --primary tcga_luad_wsi_dx1_only_files \
                          --out ~/tez/out/LUAD

    python run.py check   --out ~/tez/out/LUAD --expect-wsi 468

    python run.py sites   --cohort ~/tez/out/LUAD/cohort.csv \
                          --cohort ~/tez/out/LUSC/cohort.csv

Kural: bir asama tamamlanmadan ve check gecmeden sonrakine gecilmez.
"""
import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from tez import barcode, cohort as cohort_mod, wsi as wsi_mod, patches as patch_mod


# ------------------------------------------------------------------ cohort

def cmd_cohort(args):
    print("barcode.py kendi testi:")
    fails = barcode.self_test(verbose=True)
    if fails:
        sys.exit("\n[DUR] barcode self-test basarisiz. Kohort kurulmadi.")
    print("  -> gecti\n")

    rows, excl = cohort_mod.scan(args.data, args.label, args.primary)
    rows, excl = cohort_mod.dedupe_patients(rows, excl)
    cohort_mod.report(rows, excl, args.label)

    cp, ep = cohort_mod.write(rows, excl, args.out)
    print(f"\n  cohort.csv   : {cp}")
    print(f"  excluded.csv : {ep}")
    print(f"\n  siradaki adim: python run.py check --out {args.out} "
          f"--expect-wsi {len(rows)}")


# ------------------------------------------------------------------ wsi-qc

def cmd_wsi_qc(args):
    """Asama 2: fiziksel olcek analizi."""
    out = Path(args.out).expanduser()
    cohort_csv = out / "cohort.csv"
    if not cohort_csv.exists():
        sys.exit(f"[DUR] {cohort_csv} yok -- once 'run.py cohort' calistir")

    results, qc_path = wsi_mod.run(
        cohort_csv, args.slides, out,
        target_mpp=args.target_mpp, patch_px=args.patch_px,
        limit=args.limit, read_test_on=not args.no_read_test)

    wsi_mod.report(results, args.target_mpp, args.patch_px)
    print(f"\n  wsi_qc.csv : {qc_path}")
    print(f"\n  siradaki adim: python run.py check --out {out}")


# -------------------------------------------------------------- patch-scan

PATCH_CFG = {
    "target_mpp": 0.50,      # um/px  -- canonical fiziksel olcek
    "patch_px": 224,         # px     -- model input (uc model de 224)
    "stride_um": 112.0,      # um     -- 224 * 0.50, overlap yok
    "mask_thumb_dim": 2048,
    "seed": 42,              # wsi_id ile birlestirilir -> deterministik
    "jpeg_quality": 95,
}


def cmd_patch_scan(args):
    """Asama 3a: tavansiz aday tarama. Esik/tavan UYGULAMAZ, olcer."""
    out = Path(args.out).expanduser()
    cohort_csv = out / "cohort.csv"
    qc_csv = out / "wsi_qc.csv"
    for p in (cohort_csv, qc_csv):
        if not p.exists():
            sys.exit(f"[DUR] {p} yok -- onceki asamayi calistir")

    cfg = dict(PATCH_CFG)
    cfg["target_mpp"] = args.target_mpp
    cfg["patch_px"] = args.patch_px
    cfg["stride_um"] = args.stride_um or args.patch_px * args.target_mpp

    results, sp = patch_mod.scan(cohort_csv, qc_csv, args.slides, out,
                                 cfg, limit=args.limit)
    if results:
        patch_mod.scan_report(results)
    print(f"\n  scan_summary.csv : {sp}")
    print(f"  candidates/      : {out / 'candidates'}")
    print(f"\n  Bu ciktiya bakarak tissue_threshold ve max_patches secilecek.")


# ---------------------------------------------------- patch-manifest / extract

def cmd_patch_manifest(args):
    """Asama 3b: esik + tavan uygulayip ortak koordinat kaynagini uret."""
    out = Path(args.out).expanduser()
    for n in ("cohort.csv", "wsi_qc.csv"):
        if not (out / n).exists():
            sys.exit(f"[DUR] {out/n} yok")
    if not (out / "candidates").is_dir():
        sys.exit(f"[DUR] {out/'candidates'} yok -- once patch-scan calistir")

    cfg = dict(PATCH_CFG)
    cfg["target_mpp"] = args.target_mpp
    cfg["patch_px"] = args.patch_px

    patch_mod.build_manifest(out, out / "cohort.csv", out / "wsi_qc.csv", cfg,
                             tissue_min=args.tissue_min, cap=args.cap,
                             oversample=args.oversample,
                             pen_max=args.pen_max)
    print(f"\n  siradaki adim: python run.py patch-extract "
          f"--slides <yol> --out {out} --limit 10")


def cmd_patch_extract(args):
    """Asama 3c: gercek JPEG patch'leri yaz."""
    out = Path(args.out).expanduser()
    cfg = dict(PATCH_CFG)
    cfg["target_mpp"] = args.target_mpp
    cfg["patch_px"] = args.patch_px
    cfg["jpeg_quality"] = args.jpeg_quality

    metrics, fp = patch_mod.extract(
        out, args.slides, cfg, cap=args.cap,
        blur_min=args.blur_min, ink_max=args.ink_max,
        tissue_min=args.tissue_min, limit=args.limit,
        per_tss=args.per_tss, res_max=args.res_max,
        hema_min=args.hema_min)
    patch_mod.extract_report(metrics)
    print(f"\n  final_patch_manifest.csv : {fp}")
    print(f"  patches/                 : {out / 'patches'}")


# ------------------------------------------------------------------- check

def _load(path):
    p = Path(path).expanduser()
    if not p.exists():
        sys.exit(f"[DUR] dosya yok: {p}")
    return list(csv.DictReader(open(p)))


def cmd_check(args):
    """
    Kohort dogrulamasi. Iddia degil, test.
    Herhangi bir kontrol basarisizsa cikis kodu 1 olur.
    """
    out = Path(args.out).expanduser()
    rows = _load(out / "cohort.csv")
    excl = _load(out / "excluded.csv")

    print(f"{'='*70}\nKOHORT DOGRULAMA: {out}\n{'='*70}")
    problems = []

    def check(name, ok, detail=""):
        print(f"  [{'OK ' if ok else 'HATA'}] {name}{('  ' + detail) if detail else ''}")
        if not ok:
            problems.append(name)

    # 1 - her satirin barkodu ayristirilabiliyor mu
    bad = [r for r in rows if barcode.parse(r["file_name"]) is None]
    check("tum barkodlar ayristirilabilir", not bad,
          f"({len(bad)} bozuk)" if bad else "")

    # 2 - hepsi primary tumor (01)
    non01 = [r for r in rows if r["sample"] != "01"]
    check("hepsi sample_type=01 (primary tumor)", not non01,
          f"({[r['file_name'] for r in non01][:3]})" if non01 else "")

    # 3 - hepsi DX1, DX10 yok
    non_dx1 = [r for r in rows if r["slide"] != "DX1"]
    check("hepsi DX1 (DX10/DX11 yok)", not non_dx1,
          f"({[r['slide'] for r in non_dx1][:5]})" if non_dx1 else "")

    # 4 - hasta basina tek WSI
    dup = [p for p, n in Counter(r["patient"] for r in rows).items() if n > 1]
    check("hasta basina tek WSI", not dup, f"({dup[:5]})" if dup else "")

    # 5 - WSI sayisi = hasta sayisi
    n_wsi, n_pat = len(rows), len({r["patient"] for r in rows})
    check("WSI sayisi = hasta sayisi", n_wsi == n_pat, f"({n_wsi} / {n_pat})")

    # 6 - beklenen sayi
    if args.expect_wsi:
        check(f"beklenen WSI sayisi = {args.expect_wsi}",
              n_wsi == args.expect_wsi, f"(gercek: {n_wsi})")

    # 7 - kohort ile dislanan arasinda kesisim yok
    overlap = {r["file_name"] for r in rows} & {r["file_name"] for r in excl}
    check("kohort/dislanan kesisimi bos", not overlap,
          f"({list(overlap)[:3]})" if overlap else "")

    # 8 - her dislamanin gerekcesi var
    no_reason = [r for r in excl if not r["reason"]]
    check("her dislamanin gerekcesi var", not no_reason,
          f"({len(no_reason)} gerekcesiz)" if no_reason else "")

    # --- bilgi amacli ozet (test degil) ---
    print(f"\n  ozet: {n_wsi} WSI / {n_pat} hasta / "
          f"{len({r['tss'] for r in rows})} merkez / {len(excl)} dislanan")
    print("\n  dislama tablosu (tez icin):")
    for reason, n in Counter(r["reason"] for r in excl).most_common():
        print(f"      {n:>4}  {reason}")

    # --- asama 2: wsi_qc.csv varsa onu da dogrula ---
    qc_path = out / "wsi_qc.csv"
    if qc_path.exists():
        qc = _load(qc_path)
        print(f"\n{'='*70}\nWSI QC DOGRULAMA\n{'='*70}")

        ids_c = {r["file_name"].split(".")[0] for r in rows}
        ids_q = {r["wsi_id"] for r in qc}
        check("wsi_qc kohortla ayni WSI'lari kapsiyor", ids_c == ids_q,
              f"(kohort {len(ids_c)} / qc {len(ids_q)})")

        unread = [r for r in qc if r["readable"] != "True"]
        check("tum WSI openslide ile acilabiliyor", not unread,
              f"({len(unread)} acilamadi)" if unread else "")

        no_mpp = [r for r in qc if not r["base_mpp_x"]]
        check("tum WSI'da MPP okunabiliyor", not no_mpp,
              f"({len(no_mpp)} MPP'siz)" if no_mpp else "")

        bad = [r for r in qc if r["status"] not in ("ok",)]
        check("tum WSI status=ok", not bad,
              f"({Counter(r['status'] for r in bad).most_common()})" if bad else "")

        up = [r for r in qc if r["downscale_factor"]
              and float(r["downscale_factor"]) > 1.05]
        check("hicbir WSI yukari olceklenmiyor", not up,
              f"({len(up)} WSI hedefe ulasamiyor)" if up else "")

        okq = [r for r in qc if r["status"] == "ok"]
        if okq:
            print(f"\n  olculen MPP gruplari:")
            for m, n in Counter(round(float(r["base_mpp_x"]), 2)
                                for r in okq).most_common():
                print(f"      {n:>4}  {m:.2f} um/px")
            print(f"  secilen level: "
                  f"{dict(Counter(r['chosen_level'] for r in okq))}")

    if problems:
        print(f"\n[DUR] {len(problems)} kontrol basarisiz. Sonraki asamaya GECME.")
        sys.exit(1)
    print("\n[GECTI] Tum kontroller gecti. Sonraki asamaya gecebilirsin.")


# ------------------------------------------------------------------- sites

def cmd_sites(args):
    if len(args.cohort) != 2:
        sys.exit("iki cohort.csv gerekli: --cohort A.csv --cohort B.csv")
    a, b = (_load(p) for p in args.cohort)
    la = a[0]["label"] if a else "A"
    lb = b[0]["label"] if b else "B"
    cohort_mod.compare_sites(a, la, b, lb)


# -------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description="Tez pipeline")
    sub = ap.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("cohort", help="kohort kur (asama 1)")
    c.add_argument("--data", required=True, help="ornek: ~/tcga_data/LUAD")
    c.add_argument("--label", required=True, choices=["LUAD", "LUSC"])
    c.add_argument("--primary", default=None,
                   help="kohort adaylarinin bulundugu alt klasor")
    c.add_argument("--out", required=True)
    c.set_defaults(func=cmd_cohort)

    w = sub.add_parser("wsi-qc", help="fiziksel olcek analizi (asama 2)")
    w.add_argument("--slides", required=True,
                   help="ornek: ~/tcga_data/LUAD")
    w.add_argument("--out", required=True, help="cohort.csv'nin bulundugu klasor")
    w.add_argument("--target-mpp", type=float, default=0.50)
    w.add_argument("--patch-px", type=int, default=224)
    w.add_argument("--limit", type=int, default=0, help="ilk N WSI (pilot)")
    w.add_argument("--no-read-test", action="store_true",
                   help="okuma testini atla (daha hizli)")
    w.set_defaults(func=cmd_wsi_qc)

    ps = sub.add_parser("patch-scan", help="tavansiz aday tarama (asama 3a)")
    ps.add_argument("--slides", required=True)
    ps.add_argument("--out", required=True)
    ps.add_argument("--target-mpp", type=float, default=0.50)
    ps.add_argument("--patch-px", type=int, default=224)
    ps.add_argument("--stride-um", type=float, default=None,
                    help="varsayilan: patch_px * target_mpp (overlap yok)")
    ps.add_argument("--limit", type=int, default=0)
    ps.set_defaults(func=cmd_patch_scan)

    pm = sub.add_parser("patch-manifest", help="ortak koordinat kaynagi (asama 3b)")
    pm.add_argument("--out", required=True)
    pm.add_argument("--tissue-min", type=float, default=0.80,
                    help="thumbnail doku orani alt siniri")
    pm.add_argument("--cap", type=int, default=2000, help="WSI basina hedef")
    pm.add_argument("--pen-max", type=float, default=0.01,
                    help="patch'in kalem izi degen max orani")
    pm.add_argument("--oversample", type=float, default=1.5,
                    help="QC elemesi icin fazladan aday orani")
    pm.add_argument("--target-mpp", type=float, default=0.50)
    pm.add_argument("--patch-px", type=int, default=224)
    pm.set_defaults(func=cmd_patch_manifest)

    pe = sub.add_parser("patch-extract", help="JPEG patch yaz (asama 3c)")
    pe.add_argument("--slides", required=True)
    pe.add_argument("--out", required=True)
    pe.add_argument("--cap", type=int, default=2000)
    pe.add_argument("--blur-min", type=float, default=0.0,
                    help="0 = eleme yok, sadece olc")
    pe.add_argument("--ink-max", type=float, default=1.0,
                    help="eski hue-tabanli olcum (onerilmez)")
    pe.add_argument("--res-max", type=float, default=1.0,
                    help="H&E deconvolution kalintisi ust siniri (murekkep)")
    pe.add_argument("--hema-min", type=float, default=0.0,
                    help="hematoksilen/cekirdek orani alt siniri (hucresellik)")
    pe.add_argument("--tissue-min", type=float, default=0.0,
                    help="patch'in KENDI piksellerindeki doku orani")
    pe.add_argument("--jpeg-quality", type=int, default=95)
    pe.add_argument("--target-mpp", type=float, default=0.50)
    pe.add_argument("--patch-px", type=int, default=224)
    pe.add_argument("--limit", type=int, default=0,
                    help="ilk N WSI (alfabetik -- tek merkezden gelir!)")
    pe.add_argument("--per-tss", type=int, default=0,
                    help="her kaynak merkezden N WSI (temsili pilot icin)")
    pe.set_defaults(func=cmd_patch_extract)

    k = sub.add_parser("check", help="asamalari dogrula")
    k.add_argument("--out", required=True)
    k.add_argument("--expect-wsi", type=int, default=None)
    k.set_defaults(func=cmd_check)

    s = sub.add_parser("sites", help="iki kohortun merkez ortusmesi")
    s.add_argument("--cohort", action="append", required=True)
    s.set_defaults(func=cmd_sites)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()