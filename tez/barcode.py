"""
TCGA barkod mantigi.

Bu modul SAF tutulur: dosya okumaz, openslide cagirmaz, disk'e yazmaz.
Sadece dosya ADI alir, yapisal bilgi dondurur.

Neden ayri dosya: kohort tanimindaki iki hatanin ikisi de burada olustu.
  - "DX1" in name  ->  DX10'u da yakaladi
  - sample type kontrol edilmedi  ->  02Z (recurrent) sete sizdi
Saf oldugu icin test edilmesi kolay; degistirirsen once testi kirilir.

TCGA barkod yapisi:
    TCGA-05-4244-01Z-00-DX1
    ---- -- ---- --  --  ---
     |   |   |   |   |    +-- slayt: DX1, DX10, TS1, BS1 ...
     |   |   |   |   +------- portion / analyte
     |   |   |   +----------- sample type (01) + vial (Z)
     |   |   +--------------- hasta numarasi
     |   +------------------- TSS = tissue source site (kaynak merkez)
     +----------------------- proje
"""
import re

# ---------------------------------------------------------------- sabitler

BARCODE_RE = re.compile(
    r"^(?P<patient>TCGA-(?P<tss>[0-9A-Z]{2})-(?P<pid>[0-9A-Z]{4}))"
    r"-(?P<sample>\d{2})(?P<vial>[A-Z])"
    r"-(?P<portion>\d{2})"
    r"-(?P<slide_kind>[A-Z]{2})(?P<slide_no>\d+)"
    r"(?:\.|$)"
)

SAMPLE_TYPES = {
    "01": "Primary Solid Tumor",
    "02": "Recurrent Solid Tumor",
    "03": "Primary Blood Derived Cancer",
    "06": "Metastatic",
    "07": "Additional Metastatic",
    "11": "Solid Tissue Normal",
}

SLIDE_KINDS = {
    "DX": "Diagnostic (FFPE)",
    "TS": "Tissue Slide (frozen, top)",
    "BS": "Tissue Slide (frozen, bottom)",
    "MS": "Tissue Slide (frozen, middle)",
}

# --- kohort kurallari: primary analiz setine girmek icin gereken kosullar ---
COHORT_SAMPLE_TYPE = "01"     # yalnizca primary solid tumor
COHORT_SLIDE_KIND = "DX"      # yalnizca tanisal FFPE slayt
COHORT_SLIDE_NO = 1           # yalnizca birinci tanisal slayt


# ---------------------------------------------------------------- parse

def parse(file_name):
    """
    Dosya adindan barkod alanlarini cikarir.

    Donus: dict veya None (eslesmezse).
      patient     TCGA-05-4244
      tss         05
      sample      01              (string, basindaki sifir korunur)
      sample_name Primary Solid Tumor
      vial        Z
      portion     00
      slide_kind  DX
      slide_no    1               (INT -- DX10 burada 10 olur, 1 degil)
      slide       DX1
    """
    m = BARCODE_RE.match(str(file_name).strip())
    if not m:
        return None
    d = m.groupdict()
    d["slide_no"] = int(d["slide_no"])
    d["slide"] = f"{d['slide_kind']}{d['slide_no']}"
    d["sample_name"] = SAMPLE_TYPES.get(d["sample"], f"UNKNOWN({d['sample']})")
    d["slide_kind_name"] = SLIDE_KINDS.get(d["slide_kind"], d["slide_kind"])
    d["file_name"] = str(file_name)
    return d


# ---------------------------------------------------------------- kurallar

def cohort_reason(rec):
    """
    Bir barkod kaydinin primary kohorta girip girmedigini dondurur.

    Donus: (uygun_mu, gerekce)
      (True,  "ok")
      (False, "recurrent_sample")   gibi

    Gerekce stringleri sabit tutulur; excluded.csv'de raporlanir ve
    tezde dogrudan dislama tablosuna donusur.
    """
    if rec is None:
        return False, "barcode_parse_failed"

    if rec["slide_kind"] != COHORT_SLIDE_KIND:
        return False, f"not_diagnostic_slide({rec['slide_kind']})"

    if rec["slide_no"] != COHORT_SLIDE_NO:
        # DX10, DX11 ... ayni hastanin ek tanisal slaytlari
        return False, f"secondary_diagnostic_slide({rec['slide']})"

    if rec["sample"] != COHORT_SAMPLE_TYPE:
        # 02 = recurrent, 11 = normal, 06 = metastatic
        return False, f"non_primary_sample({rec['sample']}:{rec['sample_name']})"

    return True, "ok"


def is_cohort(file_name):
    """Kisa yol: dosya adi primary kohorta uygun mu?"""
    ok, _ = cohort_reason(parse(file_name))
    return ok


# ---------------------------------------------------------------- kendi testi

SELF_TEST = [
    # (dosya adi, uygun_mu, gerekce_anahtar_kelimesi)
    ("TCGA-05-4244-01Z-00-DX1.d4ff32cd.svs", True, "ok"),
    ("TCGA-49-4488-01Z-00-DX10.930e7dc2.svs", False, "secondary_diagnostic"),
    ("TCGA-50-5066-02Z-00-DX1.A2CA045A.svs", False, "non_primary_sample"),
    ("TCGA-85-A512-01A-01-TSA.C10905AD.svs", False, "parse_failed"),  # TSA: rakam yok
    ("TCGA-44-6146-11A-01-BS1.8db8c567.svs", False, "not_diagnostic"),
    ("bozuk_dosya_adi.svs", False, "parse_failed"),
]


def self_test(verbose=True):
    """
    Modulun kendi testi. run.py check bunu cagirir.
    Yeni bir edge case bulursan SELF_TEST'e ekle.
    """
    fails = []
    for name, want_ok, want_key in SELF_TEST:
        rec = parse(name)
        got_ok, got_reason = cohort_reason(rec)
        ok = (got_ok == want_ok) and (want_key in got_reason or want_key == "ok")
        if not ok:
            fails.append((name, want_ok, want_key, got_ok, got_reason))
        if verbose:
            mark = "OK " if ok else "HATA"
            print(f"  [{mark}] {name[:45]:<45} -> {got_ok} / {got_reason}")

    # DX10 tuzagi ayrica vurgulanir
    assert parse("TCGA-49-4488-01Z-00-DX10.x.svs")["slide_no"] == 10, \
        "DX10 yanlis ayristirildi"
    assert parse("TCGA-05-4244-01Z-00-DX1.x.svs")["slide_no"] == 1

    return fails