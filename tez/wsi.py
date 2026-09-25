"""
WSI teknik QC ve fiziksel olcek analizi.

cohort.csv'deki her WSI icin:
  - openslide ile acilabiliyor mu
  - base MPP (mikron/piksel), vendor, objective power
  - pyramid yapisi: level sayisi, boyutlar, downsample oranlari
  - hedef MPP icin HANGI level'dan okunacagi
  - gercek okuma testi ve suresi

Uretir: wsi_qc.csv

TEMEL ILKE
----------
Patch'ler ayni PIKSEL boyutunda degil, ayni FIZIKSEL alanda olmali.
Iki WSI'in ikisi de "20x" yazsa bile gercek cozunurlukleri farkli olabilir.
Bu yuzden level secimi objective power'a gore degil, olculmus MPP'ye gore
yapilir:

    base_mpp x level_downsample = o level'in gercek MPP'si
                                -> hedefe en yakin uygun level

Level secimi HER ZAMAN asagi yonludur: hedeften daha yuksek cozunurluklu
(daha kucuk MPP) bir level secilir ve yazilimda kucultulur. Yukari
olcekleme yapilmaz -- interpolasyon bilgi kazandirmaz, sahte detay uretir.
"""
import csv
import json
import time
from collections import Counter
from pathlib import Path

try:
    import openslide
except ImportError:
    openslide = None

QC_FIELDS = [
    "wsi_id", "patient", "tss", "label",
    "readable", "vendor", "objective_power",
    "base_mpp_x", "base_mpp_y", "mpp_source", "mpp_consistent",
    "level_count", "level0_w", "level0_h",
    "level_downsamples", "level_mpps",
    "target_mpp", "chosen_level", "chosen_mpp", "downscale_factor",
    "read_size_px", "read_ok", "read_secs",
    "status", "note",
]


# ---------------------------------------------------------------- metadata

def read_mpp(slide):
    """
    Level 0 mikron/piksel.
    Donus: (mpp_x, mpp_y, kaynak) -- bulunamazsa (None, None, "none")
    """
    mx = slide.properties.get("openslide.mpp-x")
    my = slide.properties.get("openslide.mpp-y")
    if mx and my:
        try:
            return float(mx), float(my), "openslide.mpp"
        except ValueError:
            pass

    am = slide.properties.get("aperio.MPP")
    if am:
        try:
            v = float(am)
            return v, v, "aperio.MPP"
        except ValueError:
            pass

    # son care: objective power'dan turet (40x ~ 0.25, 20x ~ 0.50)
    mag = (slide.properties.get("aperio.AppMag")
           or slide.properties.get("openslide.objective-power"))
    if mag:
        try:
            v = 10.0 / float(mag)
            return v, v, "objective_power"
        except ValueError:
            pass

    return None, None, "none"


def choose_level(slide, base_mpp, target_mpp):
    """
    Hedef MPP icin okunacak level.

    Hedefe ESIT veya DAHA YUKSEK cozunurluklu (mpp <= target) en ust level
    secilir. Boylece okunan goruntu her zaman hedefe KUCULTULUR.

    Donus: (level, o level'in mpp'si, downscale_factor)
      downscale_factor = chosen_mpp / target_mpp
        1.0  -> tam hedefte, resize gerekmez
        0.46 -> okunan goruntu ~2.15x kucultulecek
        >1.0 -> hedef yakalanamadi (slayt cok dusuk cozunurluklu)
    """
    best_level, best_mpp = 0, base_mpp
    for lvl, ds in enumerate(slide.level_downsamples):
        mpp = base_mpp * ds
        if mpp <= target_mpp * 1.001:
            best_level, best_mpp = lvl, mpp
    return best_level, best_mpp, best_mpp / target_mpp


def read_test(slide, level, size, n=3):
    """
    Gercek okuma testi ve sure olcumu.
    Slayt ortasindan n bolge okur. Patch cikarmanin ne kadar surecegini
    tahmin etmek icin kullanilir.
    """
    W, H = slide.dimensions
    ds = slide.level_downsamples[level]
    size0 = int(size * ds)
    cx, cy = W // 2, H // 2
    offsets = [(0, 0), (size0, 0), (0, size0)][:n]

    t0 = time.time()
    for dx, dy in offsets:
        x = max(0, min(cx + dx, W - size0 - 1))
        y = max(0, min(cy + dy, H - size0 - 1))
        img = slide.read_region((x, y), level, (size, size))
        img.convert("RGB")
    return (time.time() - t0) / max(1, len(offsets))


# ------------------------------------------------------------------ ana is

def analyze(svs_path, row, target_mpp, patch_px, do_read_test=True):
    """Tek bir WSI'i analiz eder, QC satiri dondurur."""
    rec = {f: "" for f in QC_FIELDS}
    rec.update(wsi_id=row.get("file_name", "").split(".")[0],
               patient=row.get("patient", ""), tss=row.get("tss", ""),
               label=row.get("label", ""), target_mpp=target_mpp,
               readable=False, read_ok=False, status="ok")

    if openslide is None:
        rec.update(status="error", note="openslide kurulu degil")
        return rec

    try:
        slide = openslide.OpenSlide(str(svs_path))
    except Exception as e:
        rec.update(status="error", note=f"acilamadi: {type(e).__name__}: {e}")
        return rec

    rec["readable"] = True
    rec["vendor"] = slide.properties.get("openslide.vendor", "")
    rec["objective_power"] = (slide.properties.get("aperio.AppMag")
                              or slide.properties.get("openslide.objective-power")
                              or "")

    mx, my, src = read_mpp(slide)
    rec["mpp_source"] = src
    if mx is None:
        rec.update(status="excluded", note="MPP hicbir kaynaktan okunamadi")
        slide.close()
        return rec

    rec["base_mpp_x"] = round(mx, 5)
    rec["base_mpp_y"] = round(my, 5)
    # x/y farki %1'den fazlaysa anizotropik piksel -- nadir ama isaretlenir
    rec["mpp_consistent"] = abs(mx - my) / max(mx, my) < 0.01
    if not rec["mpp_consistent"]:
        rec["note"] = f"mpp_x != mpp_y ({mx:.4f} / {my:.4f})"

    base_mpp = (mx + my) / 2
    W, H = slide.dimensions
    rec["level_count"] = slide.level_count
    rec["level0_w"], rec["level0_h"] = W, H
    ds_list = [round(float(d), 3) for d in slide.level_downsamples]
    rec["level_downsamples"] = json.dumps(ds_list)
    rec["level_mpps"] = json.dumps([round(base_mpp * d, 4) for d in ds_list])

    lvl, lvl_mpp, factor = choose_level(slide, base_mpp, target_mpp)
    rec["chosen_level"] = lvl
    rec["chosen_mpp"] = round(lvl_mpp, 4)
    rec["downscale_factor"] = round(factor, 4)

    if factor > 1.05:
        # hedef MPP'ye ulasamiyoruz: slaytin en yuksek cozunurlugu bile yetersiz
        rec.update(status="excluded",
                   note=f"hedef {target_mpp} yakalanamadi (en iyi {lvl_mpp:.3f})")
        slide.close()
        return rec

    # hedef fiziksel FOV'u bu level'da kac piksel eder
    fov_um = patch_px * target_mpp
    read_size = int(round(fov_um / lvl_mpp))
    rec["read_size_px"] = read_size

    if do_read_test:
        try:
            rec["read_secs"] = round(read_test(slide, lvl, read_size), 3)
            rec["read_ok"] = True
        except Exception as e:
            rec.update(status="error", read_ok=False,
                       note=f"okuma testi basarisiz: {type(e).__name__}")

    slide.close()
    return rec


def run(cohort_csv, slides_root, out_dir, target_mpp=0.50, patch_px=224,
        limit=0, read_test_on=True):
    """cohort.csv'deki tum WSI'lari analiz eder, wsi_qc.csv yazar."""
    rows = list(csv.DictReader(open(Path(cohort_csv).expanduser())))
    if limit:
        rows = rows[:limit]
    root = Path(slides_root).expanduser().resolve()

    out = Path(out_dir).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    qc_path = out / "wsi_qc.csv"

    fov_um = patch_px * target_mpp
    print(f"{len(rows)} WSI analiz edilecek")
    print(f"hedef: {target_mpp} um/px  |  patch {patch_px}px "
          f"-> fiziksel FOV {fov_um:.0f} x {fov_um:.0f} um\n")

    results = []
    with open(qc_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=QC_FIELDS, extrasaction="ignore")
        w.writeheader()
        for i, row in enumerate(rows, 1):
            rel = row.get("rel_path") or row["file_name"]
            p = root / rel
            if not p.exists():
                hits = list(root.rglob(row["file_name"]))
                p = hits[0] if hits else p
            rec = analyze(p, row, target_mpp, patch_px, read_test_on)
            w.writerow(rec)
            f.flush()
            results.append(rec)
            if i % 25 == 0 or i == len(rows) or rec["status"] != "ok":
                print(f"[{i}/{len(rows)}] {rec['wsi_id']:32s} "
                      f"{rec['status']:9s} mpp={rec['base_mpp_x'] or '-':<8} "
                      f"lvl={rec['chosen_level']} read={rec['read_size_px']} "
                      f"{rec['note']}")

    return results, qc_path


# ---------------------------------------------------------------- raporlama

def report(results, target_mpp, patch_px):
    ok = [r for r in results if r["status"] == "ok"]
    print(f"\n{'='*70}\nWSI QC OZETI\n{'='*70}")
    print(f"  toplam        : {len(results)}")
    for s, n in Counter(r["status"] for r in results).most_common():
        print(f"      {n:>4}  {s}")

    if not ok:
        return

    print(f"\n  MPP kaynagi:")
    for s, n in Counter(r["mpp_source"] for r in ok).most_common():
        print(f"      {n:>4}  {s}")

    print(f"\n  olculen base MPP dagilimi:")
    for m, n in Counter(round(float(r["base_mpp_x"]), 2) for r in ok).most_common():
        approx = 10.0 / m if m else 0
        print(f"      {n:>4}  {m:.2f} um/px  (~{approx:.0f}x)")

    print(f"\n  secilen level:")
    for lvl, n in Counter(r["chosen_level"] for r in ok).most_common():
        print(f"      {n:>4}  level {lvl}")

    print(f"\n  okuma boyutu (hedef {patch_px}px'e kucultulecek):")
    for rs, n in Counter(r["read_size_px"] for r in ok).most_common(6):
        f = rs / patch_px
        tag = "resize yok" if abs(f - 1) < 0.02 else f"{f:.2f}x kucultme"
        print(f"      {n:>4}  {rs}px  ({tag})")

    print(f"\n  pyramid level sayisi:")
    for lc, n in Counter(r["level_count"] for r in ok).most_common():
        print(f"      {n:>4}  {lc} level")

    secs = [float(r["read_secs"]) for r in ok if r["read_secs"] != ""]
    if secs:
        secs.sort()
        med = secs[len(secs) // 2]
        print(f"\n  patch okuma suresi (ortanca): {med*1000:.0f} ms")
        for npatch in (100, 300, 500, 1000):
            total = med * npatch * len(ok) / 3600
            print(f"      {npatch:>5} patch/WSI -> {len(ok)} WSI icin "
                  f"~{total:.1f} saat")
        print("      (tek islemci; paralel calistirilirsa bolunur)")

    anis = [r for r in ok if r["mpp_consistent"] is False
            or r["mpp_consistent"] == "False"]
    if anis:
        print(f"\n  [!] {len(anis)} WSI'da mpp_x != mpp_y")
        for r in anis[:5]:
            print(f"      {r['wsi_id']}  {r['note']}")