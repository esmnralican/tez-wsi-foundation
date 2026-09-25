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

from tez import barcode, cohort as cohort_mod


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

    if problems:
        print(f"\n[DUR] {len(problems)} kontrol basarisiz. Sonraki asamaya GECME.")
        sys.exit(1)
    print("\n[GECTI] Kohort dogrulandi. Sonraki asamaya gecebilirsin.")


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

    k = sub.add_parser("check", help="kohortu dogrula")
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