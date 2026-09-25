"""
Kohort kurma ve dogrulama.

Bu modul diski tarar, barcode.py'nin kurallarini uygular ve iki dosya uretir:

    cohort.csv      primary analiz setine giren WSI'lar
    excluded.csv    dislanan her WSI + gerekcesi

excluded.csv kritiktir: tezdeki dislama tablosu dogrudan bu dosyadan cikar.
"Neden 468?" sorusunun cevabi tahmine degil kayda dayanir.

Dislama gerekceleri iki kaynaktan gelir:
  1. Barkod kurallari (barcode.cohort_reason)  -- DX10, 02Z, TS/BS ...
  2. Disk uzerindeki _excluded_* klasorleri     -- MPP-unknown gibi,
     onceki asamalarda fiziksel olarak ayrilmis dosyalar
"""
import csv
from collections import Counter, defaultdict
from pathlib import Path

from . import barcode

COHORT_FIELDS = [
    "file_name", "patient", "tss", "sample", "sample_name",
    "slide", "label", "rel_path",
]

EXCLUDED_FIELDS = [
    "file_name", "patient", "tss", "sample", "slide",
    "reason", "source", "rel_path",
]

# _excluded_* klasor adi -> tezde kullanilacak gerekce anahtari
DIR_REASONS = {
    "_excluded_mpp_unknown": "no_reliable_mpp_downsampled_derivative",
    "_excluded_non_primary": "non_primary_sample",
    "_excluded_duplicate_patient": "secondary_diagnostic_slide",
    "_excluded_non_dx1": "secondary_diagnostic_slide",
}


def _dir_reason(path, base):
    """Dosya bir _excluded_* klasoru altindaysa gerekcesini dondurur."""
    for part in Path(path).relative_to(base).parts:
        key = part.lower()
        if key in DIR_REASONS:
            return DIR_REASONS[key]
        if key.startswith("_excluded"):
            return f"excluded_dir({part})"
    return None


def scan(base_dir, label, primary_subdir=None):
    """
    base_dir altindaki tum .svs dosyalarini tarar.

    primary_subdir verilirse yalnizca o klasordeki dosyalar kohort adayi
    sayilir; disindakiler otomatik olarak dislanmis kabul edilir.

    Donus: (cohort_rows, excluded_rows)
    """
    base = Path(base_dir).expanduser().resolve()
    if not base.is_dir():
        raise SystemExit(f"klasor yok: {base}")

    primary = (base / primary_subdir).resolve() if primary_subdir else base

    cohort, excluded = [], []
    for p in sorted(base.rglob("*.svs")):
        rel = str(p.relative_to(base))
        rec = barcode.parse(p.name)
        common = {
            "file_name": p.name,
            "patient": rec["patient"] if rec else "",
            "tss": rec["tss"] if rec else "",
            "sample": rec["sample"] if rec else "",
            "slide": rec["slide"] if rec else "",
            "rel_path": rel,
        }

        # 1) disk uzerinde zaten ayrilmis mi
        dr = _dir_reason(p, base)
        if dr:
            # Barkod kesin bir gerekce veriyorsa ONU kullan. Klasor adlari
            # elle verilir ve yaniltabilir: _excluded_duplicate_patient
            # icinde hem DX10 hem 02Z (nuks tumor) bulunabiliyor.
            ok, brc = barcode.cohort_reason(rec)
            if not ok and brc != "barcode_parse_failed":
                excluded.append({**common, "reason": brc, "source": "barcode"})
            else:
                excluded.append({**common, "reason": dr, "source": "directory"})
            continue

        # 2) primary klasorun disinda mi
        if primary_subdir and primary not in p.parents:
            excluded.append({**common, "reason": "outside_primary_dir",
                             "source": "directory"})
            continue

        # 3) barkod kurallari
        ok, reason = barcode.cohort_reason(rec)
        if not ok:
            excluded.append({**common, "reason": reason, "source": "barcode"})
            continue

        cohort.append({**common, "sample_name": rec["sample_name"],
                       "label": label})

    return cohort, excluded


def dedupe_patients(cohort, excluded):
    """
    Hasta basina tek WSI kurali.

    Barkod kurallari (01 + DX1) gectikten sonra ayni hastadan birden fazla
    dosya kaliyorsa, dosya adina gore deterministik olarak ilki tutulur.
    Kalanlar dislanir ve kaydedilir.

    Not: 468'lik LUAD kohortunda bu durum BEKLENMIYOR. Tetiklenirse
    bilmedigimiz bir sey var demektir -- check bunu uyari olarak basar.
    """
    by_patient = defaultdict(list)
    for r in cohort:
        by_patient[r["patient"]].append(r)

    kept, dropped = [], []
    for pid, rows in by_patient.items():
        rows.sort(key=lambda r: r["file_name"])
        kept.append(rows[0])
        for r in rows[1:]:
            dropped.append({
                **{k: r.get(k, "") for k in
                   ("file_name", "patient", "tss", "sample", "slide", "rel_path")},
                "reason": "duplicate_wsi_same_patient",
                "source": "dedupe",
            })

    kept.sort(key=lambda r: r["file_name"])
    return kept, excluded + dropped


def write(cohort, excluded, out_dir):
    out = Path(out_dir).expanduser()
    out.mkdir(parents=True, exist_ok=True)

    cp = out / "cohort.csv"
    with open(cp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COHORT_FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(cohort)

    ep = out / "excluded.csv"
    with open(ep, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=EXCLUDED_FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(excluded)

    return cp, ep


# ---------------------------------------------------------------- raporlama

def report(cohort, excluded, label):
    print(f"\n{'='*70}\n{label} KOHORT\n{'='*70}")
    patients = {r["patient"] for r in cohort}
    print(f"  kohorta giren WSI : {len(cohort)}")
    print(f"  benzersiz hasta   : {len(patients)}")
    if len(cohort) != len(patients):
        print(f"  [!] WSI != hasta  -- dedupe calismamis olabilir")

    print(f"\n  dislanan          : {len(excluded)}")
    for reason, n in Counter(r["reason"] for r in excluded).most_common():
        print(f"      {n:>4}  {reason}")

    tss = Counter(r["tss"] for r in cohort)
    print(f"\n  kaynak merkez     : {len(tss)} farkli")
    top = tss.most_common(5)
    for t, n in top:
        print(f"      {t}: {n:>4}  (%{n/len(cohort)*100:.1f})")
    if top and top[0][1] / len(cohort) > 0.15:
        print(f"  [!] en buyuk merkez kohortun %{top[0][1]/len(cohort)*100:.0f}'ini tutuyor"
              f" -- site-aware degerlendirme gerekli")


def compare_sites(cohort_a, label_a, cohort_b, label_b):
    """
    Iki kohortun kaynak merkez ortusmesi.
    Ortusmenin SIFIR olmasi, merkez bilgisinin sinifi tam belirledigi
    anlamina gelir -- batch effect acisindan en kotu senaryo.
    """
    a = Counter(r["tss"] for r in cohort_a)
    b = Counter(r["tss"] for r in cohort_b)
    shared = set(a) & set(b)

    print(f"\n{'='*70}\nKAYNAK MERKEZ ORTUSMESI\n{'='*70}")
    print(f"  {label_a}: {len(a)} merkez / {sum(a.values())} WSI")
    print(f"  {label_b}: {len(b)} merkez / {sum(b.values())} WSI")
    print(f"  ortak merkez: {len(shared)}")

    if not shared:
        print("\n  [!!] ORTAK MERKEZ YOK.")
        print("       Kaynak merkez bilgisi sinif etiketini %100 belirliyor.")
        print("       Model morfolojiyi ogrenmeden, boyama/tarayici imzasindan")
        print("       yuksek dogruluk elde edebilir.")
        print("       -> site prediction baseline ve stain ablation zorunlu.")
    else:
        share = sum(min(a[t], b[t]) for t in shared)
        print(f"  ortak merkezlerdeki ortusen WSI: ~{share}")
    return shared