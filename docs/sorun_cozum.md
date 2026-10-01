# SORUN → ÇÖZÜM KAYDI

**TCGA-LUAD / TCGA-LUSC · Foundation Model Karşılaştırması**
Veri hazırlama aşamasında karşılaşılan problemler ve uygulanan çözümler

Depo: `github.com/esmnralican/tez-wsi-foundation`
Tarih: 25 Eylül 2026

---

# BÖLÜM A · ÖZET TABLO

## A.1 · Kohort tanımı

| # | Sorun | Çözüm |
|---|---|---|
| 1 | `"DX1" in dosya_adı` filtresi `DX10`'u da yakalıyordu | Slayt numarası tamsayı olarak ayrıştırılıyor, `slide_no == 1` |
| 2 | Örnek tipi hiç kontrol edilmiyordu; nüks tümör kohorta sızmıştı | `sample_type == 01` kuralı eklendi |
| 3 | Dışlama gerekçesi klasör adından türetiliyordu, iki farklı sebep aynı etiketi alıyordu | Barkod kesin gerekçe verebiliyorsa o kullanılıyor |
| 4 | 10 WSI'ın dışlanma gerekçesi "metadata eksikliği" idi, eksikti | Level-0 boyutları incelendi: tümü tam 10.000 px → küçültülmüş türev görüntü |

## A.2 · Fiziksel ölçek

| # | Sorun | Çözüm |
|---|---|---|
| 5 | Aynı piksel boyutundaki patch'ler aynı doku alanını göstermiyor | Hedef fiziksel ölçek sabitlendi: 0.50 µm/px, 112×112 µm FOV |
| 6 | Metadata'daki 20×/40× bilgisi güvenilmez (20× yazan 43 WSI'ın 19'u gerçekten 0.50) | Level seçimi ölçülmüş MPP'ye göre yapılıyor |
| 7 | Pyramid level sayısı slayttan slayta değişiyor (4, 3, 2) | Sabit level yok; her WSI için `base_mpp × downsample` hesaplanıyor |
| 8 | Yukarı ölçekleme sahte detay üretir | Level seçimi her zaman aşağı yönlü, interpolasyonla büyütme yok |

## A.3 · Patch protokolü

| # | Sorun | Çözüm |
|---|---|---|
| 9 | Doku eşiği karşılaştırması geçersizdi — ölçülen eşik değil tavan etkisiydi | Aday tarama tavansız yeniden çalıştırıldı |
| 10 | `max_patches = 100` "hız için" seçilmiş ama hız ölçülmemişti | Ölçüldü: 15–33 ms/patch. Gerekçe geçersiz kaldı |
| 11 | Blur QC resize öncesi uygulanıyordu, farklı MPP grupları haksız değerlendiriliyordu | Blur QC her zaman 224×224'e küçültme sonrası |
| 12 | Mürekkep maskesinin ton aralığı hematoksileni kapsıyordu | Ton aralığı 125'te kesildi; renk tek başına kullanılmıyor |
| 13 | Uniform turkuaz bölgeler doku maskesini geçiyordu | Blur QC zorunlu hale getirildi (Laplacian varyansı ≈ 0) |
| 14 | Patch seviyesi doku oranı metriği boyama yoğunluğunu ölçüyordu | Metrikle filtreleme yapılmıyor, yalnızca kaydediliyor |
| 15 | Ton aralığını daraltmak mürekkep/hematoksilen ayrımını çözmedi | H&E deconvolution kalıntısına geçildi, eşik 0.55 ölçümle kalibre edildi |
| 16 | İlk kalem eşiği (0.35) 468 slaytın 370'inde yanlış alarm verdi | Kalemli/kalemsiz slayt karşılaştırmasıyla 0.55'e çekildi → 21 WSI |
| 17 | `hash()` ile seed üretimi Python'da süreç başına değişiyor | `hashlib.blake2b` ile sabit hash (LUSC öncesi uygulanacak) |
| 18 | `TCGA-44-7661` sıfır patch üretti (slayt tamamen odak dışı) | Kohorttan çıkarıldı, gerekçe kayda geçti → 467 WSI |

## A.4 · Metodolojik

| # | Sorun | Çözüm |
|---|---|---|
| 19 | Aynı hastanın iki slaytı train/test'e bölünebilir (data leakage) | Hasta başına tek WSI; split `patient` üzerinden |
| 20 | LUAD ve LUSC'un kaynak merkezleri hiç örtüşmüyor (ortak = 0) | Çözülemez; ölçülüyor — `tss` kolonu manifest'te tutuluyor |
| 21 | Patch seviyesinde etiket yok | MIL (bag = hasta, instance = patch), ABMIL |
| 22 | Önceki çalışmada Macenko referansı sınıf-spesifikti (etiket sızıntısı) | Yeni kurguda tek referans, sınıftan bağımsız |

---

# BÖLÜM B · DETAYLI AÇIKLAMA

---

## B.1 · `"DX1" in dosya_adı` filtresi DX10'u da yakalıyordu

### Sorun

TCGA barkodunda son alan slayt numarasını taşır:

```
TCGA-49-4488-01Z-00-DX1     ← birinci tanısal slayt
TCGA-49-4488-01Z-00-DX10    ← AYNI HASTANIN onuncu slaytı
```

Basit bir alt dize kontrolü ikisini de kabul eder:

```python
"DX1" in "TCGA-49-4488-01Z-00-DX10.930e7dc2.svs"   # True
```

Sonuç: aynı hastanın iki kaydı sete girdi. Hasta seviyesinde mükerrerlik,
dolayısıyla split sırasında data leakage riski.

### Çözüm

Slayt numarası tamsayı olarak ayrıştırılıyor ve tam eşitlik aranıyor.

`tez/barcode.py`:

```python
BARCODE_RE = re.compile(
    r"^(?P<patient>TCGA-(?P<tss>[0-9A-Z]{2})-(?P<pid>[0-9A-Z]{4}))"
    r"-(?P<sample>\d{2})(?P<vial>[A-Z])"
    r"-(?P<portion>\d{2})"
    r"-(?P<slide_kind>[A-Z]{2})(?P<slide_no>\d+)"
    r"(?:\.|$)"
)

def parse(file_name):
    m = BARCODE_RE.match(str(file_name).strip())
    if not m:
        return None
    d = m.groupdict()
    d["slide_no"] = int(d["slide_no"])      # "10" -> 10, "1" -> 1
    ...
```

Kural:

```python
COHORT_SLIDE_NO = 1

def cohort_reason(rec):
    ...
    if rec["slide_no"] != COHORT_SLIDE_NO:
        return False, f"secondary_diagnostic_slide({rec['slide']})"
```

### Doğrulama

Modül kendi birim testini içeriyor ve kohort kurulmadan önce otomatik
çalışıyor:

```python
SELF_TEST = [
    ("TCGA-05-4244-01Z-00-DX1.d4ff32cd.svs",  True,  "ok"),
    ("TCGA-49-4488-01Z-00-DX10.930e7dc2.svs", False, "secondary_diagnostic"),
    ...
]

assert parse("TCGA-49-4488-01Z-00-DX10.x.svs")["slide_no"] == 10
assert parse("TCGA-05-4244-01Z-00-DX1.x.svs")["slide_no"] == 1
```

Test başarısızsa `run.py cohort` durur, kohort kurulmaz.

---

## B.2 · Örnek tipi kontrol edilmiyordu

### Sorun

Barkodun beşinci alanı örnek tipini gösterir:

```
TCGA-50-5066-01Z-00-DX1     01 = Primary Solid Tumor
TCGA-50-5066-02Z-00-DX1     02 = Recurrent Solid Tumor  ← AYNI HASTA
              ^^
```

İkisi de `DX1` olduğu için slayt numarası filtresinden geçtiler. Nüks tümör
tedavi sonrası alınmıştır; nekroz, fibrozis ve terapiye bağlı morfolojik
değişiklikler taşıyabilir. Primary kohorta ait değildir.

Başlangıçta bu durum "aynı hastanın ikinci DX1 slaytı, normaldir" şeklinde
yorumlanmıştı. Görsel olarak iki slayt gerçekten farklı doku içeriyordu;
ancak fark görüntüde değil **örnek tipindeydi**.

### Çözüm

`tez/barcode.py`:

```python
SAMPLE_TYPES = {
    "01": "Primary Solid Tumor",
    "02": "Recurrent Solid Tumor",
    "06": "Metastatic",
    "11": "Solid Tissue Normal",
}

COHORT_SAMPLE_TYPE = "01"

def cohort_reason(rec):
    ...
    if rec["sample"] != COHORT_SAMPLE_TYPE:
        return False, f"non_primary_sample({rec['sample']}:{rec['sample_name']})"
    return True, "ok"
```

### Sonuç

```
480  GDC manifest'teki toplam DX kaydı
 -1  DX10 (ikincil tanısal slayt)
 -1  02Z (nüks tümör)
-10  güvenilir MPP bilgisi olmayan WSI
────
468  analiz kohortu · 468 WSI = 468 hasta (birebir)
```

---

## B.3 · Dışlama gerekçesi klasör adından türetiliyordu

### Sorun

Dışlanan dosyalar elle `_excluded_*` klasörlerine taşınmıştı ve gerekçe
klasör adından çıkarılıyordu. Ancak `_excluded_duplicate_patient` klasörü
iki farklı sebeple dışlanmış dosya içeriyordu:

```
_excluded_duplicate_patient/
├── TCGA-49-4488-01Z-00-DX10...svs    ← ikincil tanısal slayt
└── TCGA-50-5066-02Z-00-DX1...svs     ← nüks tümör
```

Her ikisi de `secondary_diagnostic_slide` etiketi alıyordu. Bu dosya tezin
dışlama tablosuna dönüşeceği için etiketlerin doğru olması gerekiyordu.

Benzer şekilde `_excluded_non_DX1` klasör adı da yanıltıcıydı: içindeki
dosya bir DX slaytıydı, sadece onuncu tanısal slayttı.

### Çözüm

Barkod kesin bir gerekçe verebiliyorsa o kullanılıyor; klasör adı yalnızca
barkodun söyleyemediği durumlar için (MPP eksikliği gibi).

`tez/cohort.py`, `scan()` içinde:

```python
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
```

### Sonuç

Öncesi:

```
10  no_reliable_mpp_downsampled_derivative
 2  secondary_diagnostic_slide            ← iki farklı sebep, tek etiket
```

Sonrası:

```
10  no_reliable_mpp_downsampled_derivative
 1  secondary_diagnostic_slide(DX10)
 1  non_primary_sample(02:Recurrent Solid Tumor)
```

`excluded.csv` dosyası doğrudan tezin dışlama tablosuna dönüşüyor.

---

## B.4 · 10 WSI'ın dışlanma gerekçesi eksikti

### Sorun

10 WSI'da `openslide.mpp-x` alanı boştu ve gerekçe "metadata eksikliği"
olarak kaydedilmişti. Ancak bu eksik bir açıklamaydı — dosyaların neden
metadata içermediği bilinmiyordu.

### İnceleme

Vendor-specific alanlar da kontrol edildi:

```bash
python3 -c "
import openslide, glob
for p in sorted(glob.glob('*/*.svs')):
    s = openslide.OpenSlide(p); pr = s.properties
    print(p.split('/')[-1][:30],
          '| mpp-x:', pr.get('openslide.mpp-x'),
          '| AppMag:', pr.get('aperio.AppMag'),
          '| obj-power:', pr.get('openslide.objective-power'),
          '| L0:', s.level_dimensions[0])
    s.close()
"
```

Çıktı:

```
TCGA-05-4384-... | mpp-x: None | AppMag: None | obj-power: None | L0: (10000, 9419)
TCGA-05-4390-... | mpp-x: None | AppMag: None | obj-power: None | L0: (6745, 10000)
TCGA-05-4410-... | mpp-x: None | AppMag: None | obj-power: None | L0: (8369, 10000)
...
```

**Onunda da uzun kenar tam olarak 10.000 piksel.**

Karşılaştırma (aynı merkezden sağlam slaytlar):

| WSI | L0 boyutu | MPP |
|---|---|---|
| TCGA-05-4244 | 46336 × 44288 | 0.2325 |
| TCGA-05-4245 | 37120 × 78080 | 0.2325 |
| TCGA-05-4249 | 68608 × 35072 | 0.2325 |

### Çözüm

MPP eksikliği **semptomdu, sebep değildi.** Bu dosyalar orijinal tarama
değil, uzun kenarı 10.000 piksele sabitlenerek yeniden boyutlandırılmış
türev görüntülerdir. Tahmini fiziksel çözünürlük ~2–2.5 µm/px (≈4–6×),
hedef 0.50 µm/px seviyesinin çok altında.

Hedefe yükseltmek ancak interpolasyonla mümkün olurdu; bu da hücresel
morfoloji içermeyen yapay patch'ler üretirdi.

Kodda çift güvenlik — `tez/wsi.py`:

```python
if not (cfg["min_base_mpp"] <= base_mpp <= cfg["max_base_mpp"]):
    rec.update(status="excluded",
               reason=f"MPP beklenen araligin disinda ({base_mpp:.4f})")
```

### Metodolojik not

Dışlanan 10 WSI'ın **tamamı tek bir kaynak merkezden** (TSS 05) gelmektedir.
Bu sistematik bir kayıptır ve tezde açıkça raporlanacaktır. Aynı merkezden
24 WSI kohortta kalmıştır.

---

## B.5 · Aynı piksel boyutu ≠ aynı doku alanı

### Sorun

Üç foundation model de 224×224 girdi alır. Ancak "hepsi 224 alıyor" demek
"hepsi aynı biyolojik alanı görüyor" anlamına gelmez:

```
224 px @ 0.25 µm/px  →   56 × 56 µm
224 px @ 0.50 µm/px  →  112 × 112 µm
```

Aynı piksel boyutundaki iki patch, slaydın çözünürlüğüne göre dört kat
farklı doku alanı gösterebilir.

### Çözüm

Hedef **fiziksel ölçek** sabitlendi:

```
target_mpp  = 0.50 µm/px
patch_px    = 224
FOV         = 224 × 0.50 = 112 × 112 µm
```

Her WSI kendi ölçülmüş MPP'sine göre farklı piksel boyutunda okunuyor,
ardından 224'e küçültülüyor.

`tez/patches.py`:

```python
fov_um = cfg["patch_px"] * cfg["target_mpp"]      # 112 um
read_size = int(round(fov_um / lvl_mpp))          # her WSI icin farkli
...
if rs != px:
    interp = cv2.INTER_AREA if rs > px else cv2.INTER_LINEAR
    arr = cv2.resize(arr, (px, px), interpolation=interp)
```

### Ölçülen sonuç (468 WSI)

| Base MPP | Okunan | İşlem | WSI |
|---|---|---|---:|
| 0.252 | 444 px | → 224 px (1.98× küçültme) | 308 |
| 0.2527 | 443 px | → 224 px | 37 |
| 0.2325 | 482 px | → 224 px (2.15× küçültme) | 24 |
| 0.2457 | 456 px | → 224 px | 24 |
| 0.500 | 224 px | resize yok | 19 |
| 0.2485 | 451 px | → 224 px | 14 |

Farklı piksel sayısı okunuyor, **aynı 112 × 112 µm fiziksel alan** elde
ediliyor.

### Fairness açısından önemi

```
                  Canonical patch (112 × 112 µm)
                 /              |              \
            DINOv2             UNI          Virchow
```

Üç model aynı WSI'dan, aynı ön işlemeden geçmiş, aynı koordinatlardaki,
aynı fiziksel alanı temsil eden patch'leri görecek. Karşılaştırmada değişen
tek değişken foundation model olacak.

---

## B.6 · Objective power tek başına güvenilmez

### Sorun

Metadata'daki 20×/40× bilgisine güvenilerek sabit level seçmek yanlış
sonuç verir:

```
YANLIŞ:  20× → level 1
         40× → level 2
```

Ölçüm: metadata'da "20×" yazan 43 WSI'ın yalnızca 19'unun gerçek
çözünürlüğü 0.50 µm/px. Kalan 24'ü 0.2325 civarında — yani metadata 20×
dese de fiziksel çözünürlük 40× seviyesinde.

Ayrıca pyramid yapısı slayttan slayta değişiyor: 368 WSI'da 4 level,
99'unda 3 level, 1'inde 2 level. Sabit bir level numarası kullanılamaz.

### Çözüm

Level seçimi ölçülmüş MPP'ye göre, her WSI için ayrı hesaplanıyor.

`tez/wsi.py`:

```python
def choose_level(slide, base_mpp, target_mpp):
    """
    Hedefe ESIT veya DAHA YUKSEK cozunurluklu (mpp <= target) en ust level.
    Boylece okunan goruntu her zaman hedefe KUCULTULUR.
    """
    best_level, best_mpp = 0, base_mpp
    for lvl, ds in enumerate(slide.level_downsamples):
        mpp = base_mpp * ds
        if mpp <= target_mpp * 1.001:
            best_level, best_mpp = lvl, mpp
    return best_level, best_mpp, best_mpp / target_mpp
```

MPP okuma, birden fazla kaynağı sırayla deniyor:

```python
def read_mpp(slide):
    mx = slide.properties.get("openslide.mpp-x")
    my = slide.properties.get("openslide.mpp-y")
    if mx and my:
        return float(mx), float(my), "openslide.mpp"

    am = slide.properties.get("aperio.MPP")
    if am:
        v = float(am)
        return v, v, "aperio.MPP"

    mag = (slide.properties.get("aperio.AppMag")
           or slide.properties.get("openslide.objective-power"))
    if mag:
        v = 10.0 / float(mag)
        return v, v, "objective_power"

    return None, None, "none"
```

### Doğrulama

`run.py check` yukarı ölçekleme olmadığını test ediyor:

```python
up = [r for r in qc if r["downscale_factor"]
      and float(r["downscale_factor"]) > 1.05]
check("hicbir WSI yukari olceklenmiyor", not up, ...)
```

468 WSI'ın tamamı level 0'dan okunuyor — Aperio pyramid'i tipik olarak
`[1, 4, 16, 32]` oranlarında olduğu için level 1 zaten 1.0 µm/px'e düşüyor
(hedefin iki katı kaba). 0.50 arada kalıyor ve level 0'dan küçültülmesi
gerekiyor.

---

## B.7 · Doku eşiği karşılaştırması geçersizdi

### Sorun

Önceki bir denemede 0.70 ve 0.80 doku eşikleri karşılaştırılmış ve şu sonuca
varılmıştı:

> "Aralarında yalnızca 30 patch (%0,064) fark var, 0.80 kullanılabilir."

Ancak o denemede slayt başına 100 patch üst sınırı vardı:

```
469 WSI × 100 max_patches = 46.900   ← tavan
Ölçülen:  0.70 → 46.897   ·   0.80 → 46.867
```

Neredeyse her slayt tavana dayanmıştı. Bir slaytta 5.000 aday varsa, eşiği
sıkılaştırıp 4.000'e düşürsen bile yine 100 seçilir — sayı değişmez.

**Ölçülen şey eşiğin etkisi değil, tavanın etkisiydi.**

### Çözüm

Aday tarama aşaması hiçbir eşik veya tavan uygulamadan yeniden çalıştırıldı.
Her eşik değeri için ayrı ayrı sayım yapılıyor.

`tez/patches.py`, `scan_slide()`:

```python
t_ratio = _window_means(tis, mx, my, wx, wy, mw, mh)
x_ratio = _window_means(tex, mx, my, wx, wy, mw, mh)

for thr in (30, 50, 60, 70, 80, 90):
    rec[f"n_ge_{thr}"] = int((t_ratio >= thr / 100).sum())
```

Grid taraması integral image ile toplu yapılıyor (hızlı):

```python
def _window_means(mask, xs, ys, wx, wy, mw, mh):
    """Integral image ile toplu pencere ortalamasi."""
    ii = cv2.integral(mask.astype(np.uint8))
    ...
    s = (ii[y1, x1] - ii[y0, x1] - ii[y1, x0] + ii[y0, x0])
    return s / area
```

### Gerçek ölçüm (468 WSI)

| Doku eşiği | ortanca/slayt | toplam |
|---|---:|---:|
| ≥ 0.30 | 11.518 | 6,12 M |
| ≥ 0.50 | 10.742 | 5,85 M |
| ≥ 0.70 | 10.126 | 5,54 M |
| ≥ 0.80 | 9.486 | 5,34 M |
| ≥ 0.90 | 8.794 | 5,10 M |

Gerçek fark: 0.70 → 0.80 geçişinde slayt başına ortanca ~640 patch
kayboluyor. Eski audit'in gördüğü "30 patch" değil.

**Bulgu:** 0.30 ile 0.90 arasındaki tüm aralıkta toplam aday sayısı yalnızca
%17 değişiyor. Doku eşiği düşük etkili bir parametre.

---

## B.8 · `max_patches = 100` ölçülmeden gerekçelendirilmişti

### Sorun

Protokolde slayt başına 100 patch belirlenmişti; gerekçe "hız ve disk
kullanımı" idi. Ancak ne hız ne disk kullanımı ölçülmüştü.

Kapsam açısından anlamı:

```
100 patch × (112 µm)² = 1.25 mm²
Tipik WSI doku alanı   = 138 mm² (ortanca)
Kapsanan oran          = %0.9
```

Slaydın binde dokuzu örnekleniyordu.

### Ölçüm

`tez/wsi.py` içine okuma testi eklendi:

```python
def read_test(slide, level, size, n=3):
    """Gercek okuma testi ve sure olcumu."""
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
```

Sonuç: ortanca **15–33 ms/patch**.

| Patch/WSI | 468 WSI için süre | Disk |
|---:|---|---|
| 100 | ~0.4 sa | ~1,4 GB |
| 500 | ~2,1 sa | ~7 GB |
| 1000 | ~4,3 sa | ~14 GB |

### Çözüm

Hız ve disk gerekçesi geçerliliğini yitirdi. Tavan parametresi artık
ölçüme dayanarak seçiliyor ve **geri alınabilir**: tüm aday koordinatlar
`candidates/` altında saklandığı için sayıyı değiştirmek yeniden tarama
gerektirmiyor, yalnızca çıkarma aşaması tekrarlanıyor.

### Literatür bağlamı

Standart yaklaşım tavan koymamaktır — CLAM protokolünde doku konturu
içindeki her koordinat alınır, slayt başına üst sınır yoktur. Bu çalışmada
örnekleme donanım kısıtı nedeniyle zorunludur ve tezde açıkça
belirtilecektir.

---

## B.9 · Blur QC yanlış sırada uygulanıyordu

### Sorun

Laplacian varyansı çözünürlüğe bağlı bir ölçüdür. Farklı base MPP'ye sahip
slaytlar farklı piksel boyutlarında okunur:

```
40× WSI        →  444 px okunur  →  224'e küçültülür
native 20× WSI →  224 px okunur  →  resize yok
```

Blur eşiği resize öncesinde uygulanırsa iki grup farklı ölçekte
değerlendirilir. 19 native-20× slaytın patch'leri sistematik olarak haksız
elenebilir.

### Çözüm

`tez/patches.py`, `extract()` içinde — sıra kesin:

```python
arr = np.array(img.convert("RGB"))

# --- ONCE resize, SONRA QC (kural 1.3) ---
if rs != px:
    interp = cv2.INTER_AREA if rs > px else cv2.INTER_LINEAR
    arr = cv2.resize(arr, (px, px), interpolation=interp)

tis, ink, blur = patch_metrics(arr)
if blur < blur_min:
    qc["n_blur_rej"] += 1
    continue
```

Tüm ölçümler 224×224 üzerinde, yani her slayt için aynı ölçekte yapılıyor.

---

## B.10 · Mürekkep maskesi hematoksileni kapsıyordu

### Sorun

İlk uygulamada ink maskesi HSV ton aralığı `100–140` olarak tanımlanmıştı.
Ancak H&E'de çekirdekleri boyayan hematoksilen **mor-mavidir** ve tam bu
aralıkta yer alır.

Görsel denetimde "mavi oranı ≥ 0.80" kuralıyla ayrılan patch'lerde üç farklı
içerik bulundu:

- gerçek mavi mürekkep / kalem izi
- **çekirdek yoğun gerçek tümör dokusu**
- turkuaz arka plan artefaktı

Ölçülen ton aralıkları (OpenCV 0–179 ölçeği):

| İçerik | Ton | Doygunluk | Laplacian varyansı |
|---|---|---|---|
| Turkuaz artefakt | ~80–95 | orta | ≈ 0 |
| Mavi mürekkep | ~105–125 | çok yüksek | düşük |
| **Hematoksilen** | **~125–150** | orta | **yüksek** |

Uygulansaydı en bilgilendirici patch'ler mürekkep sanılıp silinecekti.

### Çözüm

Ton aralığı 125'te kesildi ve renk tek başına kullanılmıyor.

`tez/patches.py`:

```python
def patch_metrics(rgb):
    """
    NOT: hematoksilen hue ~125-150, mavi murekkep ~105-125.
    Ink maskesi 125'te kesilir; aksi halde cekirdek yogun tumor dokusu
    murekkep sanilip silinir.
    """
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    hue, sat, val = hsv[..., 0], hsv[..., 1], hsv[..., 2]

    ink = ((sat > 90) & (val > 40) &
           (((hue >= 95) & (hue < 125)) |        # mavi murekkep
            ((hue >= 35) & (hue < 95)))).mean()  # yesil/turkuaz marker
```

Asıl ayırt edici özellik renk değil **doku/yapıdır**: mürekkep ve artefakt
düzdür, gerçek doku kenar ve yapı içerir. Bu nedenle Laplacian varyansı
(blur QC) ikinci bir güvence olarak kullanılıyor.

### Literatür doğrulaması

2025 tarihli bir foundation model benchmark çalışmasında UNI, Hibou, Virchow
ve Kaiko gibi modellerin zaman zaman kalem izlerini öne çıkardığı
gözlenmiştir. Aynı çalışmada kalem izlerinin bir kohortun slaytlarının
%90'ında, başkasında %22'sinde, diğerlerinde hiç bulunmadığı ölçülmüştür —
yani mürekkep varlığı da merkeze bağlı sistematik bir değişkendir.

### Ton daraltması yetmedi — nihai çözüm

Ton aralığını 125'te kesmek de işe yaramadı. Görsel denetimde, en yüksek
`ink` skoru alan 15 patch'in **13'ü mürekkep değildi**: TCGA-83-5908'den
gelen çekirdek yoğun gerçek tümör dokusuydu.

Sebep: boyama yoğunluğu merkeze göre kaydığı için sabit bir ton sınırı her
merkezde aynı yere düşmüyor.

**Çözüm: renk tonu yerine H&E stain deconvolution kalıntısı.**

```python
od = -np.log10((rgb + 1.0) / 256.0)
coeffs = od @ HE_PINV.T          # hematoksilen + eozin konsantrasyonu
recon  = coeffs @ HE
residual = ||od - recon|| / ||od||
```

Soru şu: *"Bu rengi hematoksilen ve eozin karışımıyla üretebilir miyim?"*

| | Kalıntı | Neden |
|---|---|---|
| Koyu mor çekirdek | düşük | Rengi zaten hematoksilenden geliyor |
| Kırmızı/yeşil kalem | yüksek | H&E dışı pigment |

Çekirdek ne kadar koyu olursa olsun kalıntı düşük kalır — boyama
yoğunluğundan bağımsız.

Ayrıca tespit **patch seviyesinden slayt seviyesine** taşındı: kalem izleri
belirli bölgelerde yoğunlaştığı için rastgele patch örneklemesinde
görülmeyebiliyor. Artık her slaydın thumbnail'i taranıyor ve işaretli
koordinatlar manifest aşamasında havuzdan çıkarılıyor.

### Eşik kalibrasyonu  `[ölçüm]`

İlk denenen 0.35 değeri 468 slaydın **370'inde** "belirgin kalem" buldu.
ImageScope kontrolünde üstteki 5 slayttan yalnızca 1'inde gerçek kalem vardı.

Üç kalemli ve dört kalemsiz slayt karşılaştırıldı:

| Slayt | Kalem? | >0.35 | >0.45 | **>0.55** | >0.65 |
|---|---|---:|---:|---:|---:|
| TCGA-86-7953 | var | %90.8 | %80.4 | **%63.1** | %28.9 |
| TCGA-86-8056 | var | %55.0 | %44.0 | **%16.5** | %3.1 |
| TCGA-86-7714 | var | %26.5 | %13.6 | **%5.6** | %1.0 |
| TCGA-75-6214 | yok | %57.0 | %3.7 | **%0.1** | %0.0 |
| TCGA-50-6673 | yok | %78.9 | %5.3 | **%0.2** | %0.0 |
| TCGA-50-5045 | yok | %61.7 | %2.0 | **%0.1** | %0.0 |
| TCGA-99-8025 | yok | %54.6 | %0.9 | **%0.0** | %0.0 |

0.55'te iki grup 28 kat ayrışıyor. Eşik burada sabitlendi.

### Sonuç

| | 0.35 | **0.55** |
|---|---:|---:|
| "Belirgin kalem" WSI | 370 | **21** |
| Elenen aday patch | 1.493.976 | **47.057** |

Belirgin olan 21 WSI'ın **tamamı TSS 86'dan**. Merkez ortalamaları: TSS 86
%11.58, ikinci sıradaki TSS 93 %1.65.

> Kalem kullanımı bir hastanenin rutinidir — literatürdeki bulgunun
> doğrudan tekrarı. Merkez sınıfı belirlediği için bu izlerin temizlenmesi
> metodolojik olarak zorunludur.

---

## B.10b · Mürekkep kontrolünün nihai mimarisi — iki aşamalı

Tek filtre yeterli olmadı. Son hal iki farklı ölçüm, iki farklı aşamada.

### Aşama 1 — Slayt seviyesi, manifest'ten önce

`tez/patches.py` → `tissue_and_texture_mask()`:

```python
od     = -np.log10((thumb + 1.0) / 256.0)
coeffs = od @ _HE_PINV.T              # hematoksilen + eozin konsantrasyonu
recon  = coeffs @ _HE
resid  = ||od - recon|| / ||od||      # H&E ile aciklanamayan kisim

pen = (resid > 0.55) & (sat > 60) & (val > 30)
pen = cv2.morphologyEx(pen, MORPH_OPEN, k3)
pen = cv2.dilate(pen, k7)
```

Üç koşul birlikte: H&E ile açıklanamıyor **ve** doygun **ve** yeterince
koyu. Renk listesi yok — yeşil, mavi, kırmızı, siyah fark etmez; H&E'nin
doğal renk uzayı dışındaki her pigment yakalanır.

Manifest aşamasında `--pen-max 0.01` ile uygulanır.

| | Elenen aday |
|---|---:|
| LUAD | 47.057 |
| LUSC | 41.448 |

### Aşama 2 — Patch seviyesi, çıkarmadan sonra

Thumbnail çözünürlüğünde kalemin yayılmış ince gölgesi kaçabiliyor
(2048 px önizlemede 444 px'lik bir patch ~13×13 piksele düşer). Gerçek
patch pikselleri üzerinde ikinci ölçüm:

```python
c   = od @ _HE_PINV.T                    # (N,2) H ve E katsayilari
neg = (c.min(axis=1) < -0.05).mean()     # negatif katsayi orani
```

**Fiziksel mantık:** gerçek H&E dokusunda hem hematoksilen hem eozin
konsantrasyonu pozitif olmalıdır — boya miktarı negatif olamaz. H&E dışı
renkler ancak negatif katsayıyla temsil edilebilir.

| İçerik | negatif% | eozin katsayısı |
|---|---:|---:|
| **Kalem altı turkuaz** | **100,0** | **−0,596** |
| Aynı slayttan ikinci turkuaz | 100,0 | −0,510 |
| Gerçek doku (11 örnek, 7 merkez) | 0,0 – 2,8 | +0,04 … +0,52 |

Ara değer yok; dağılım bimodal. Eşik **0,15** `[ölçüm + görsel]`.

### Neden `he_residual` patch seviyesinde çalışmadı

Aynı kalıntı ölçüsü patch seviyesinde denendiğinde turkuaz bir patch
**0,1346** verdi — kohort ortancasından (0,275) **düşük**.

Sebep: deconvolution matrisi iki vektör içerir (H ve E), renk uzayı üç
boyutludur. Negatif katsayılara izin verildiğinde herhangi bir renk iki
vektörün lineer kombinasyonuyla yaklaşık temsil edilebilir. Kalıntı düşük
çıkar.

Ölçü matematiksel olarak doğru ama **fiziksel olarak anlamsız**. Doğru soru
"temsil edilebilir mi?" değil, **"fiziksel olarak mümkün mü?"**

Ayrıca metrik merkeze göre belirgin değişiyordu (TSS 86 p50 = 0,331;
TSS MP p50 = 0,244) — yani mürekkebi değil boyama yoğunluğunu ölçüyordu.
0,30 eşiği uygulansaydı TSS 86'nın patch'lerinin %74'ü elenecekti.

### Eşik iki kez yanlış seçildi

| Deneme | Sonuç |
|---|---|
| **0,05** | Antrakoz ve yoğun hematoksilenli gerçek dokuyu eliyordu |
| **0,50** | Gevşekti; 0,49 civarında hâlâ mavi/turkuaz kalıntı vardı |
| **0,15** | Görsel denetimle: 0,10–0,15 ağırlıklı temiz, 0,15+ mavi başlıyor |

Antrakoz (akciğerdeki karbon birikimi) gerçek biyolojik yapıdır ve miktarı
hastanın yaşadığı yere, sigara geçmişine, mesleğine bağlıdır — dolayısıyla
merkeze bağlıdır. Elenmesi merkez-bağımlı seçilim yaratır.

### Mimarinin özeti

| | Slayt seviyesi | Patch seviyesi |
|---|---|---|
| Yakaladığı | Kalın kalem konturu | Yayılmış ince gölge |
| Ne zaman | Manifest, çıkarmadan önce | Çıkarmadan sonra |
| Metrik | `he_residual > 0.55` | `neg_frac > 0.15` |
| Sonuç | Aday havuzundan çıkar | Manifest'ten düşer, JPEG diskte kalır |
| Maliyet | Ucuz (thumbnail) | ~25 dk (her patch okunur) |
| Elenen (LUAD) | 47.057 aday | 2.105 patch |

---

## B.10c · Kalem maskesi yalnızca kırmızıyı yakalıyordu

### Nasıl fark edildi

Yoğun kalemli slaytlarda patch sayısı anormal düşüktü. `extract_qc.csv`:

```
TCGA-85-8276:  manifest aday 3.000 -> blur eledi 2.451 -> 549 patch
TCGA-85-8071:  manifest aday 2.928 -> blur eledi 2.355 -> 573 patch
TCGA-86-7701:  manifest aday 3.000 -> blur eledi 2.576 -> 424 patch
```

Mürekkep filtresi **sıfır** eledi, blur filtresi binlerce.

Elenen patch'ler görsel olarak incelendiğinde sebep anlaşıldı: **düz kırmızı
ve düz turkuaz kareler** — doku içermeyen katı mürekkep blokları. Kalem
maskesi onları kaçırmış, blur yakalamıştı.

### Hatalı teşhisler (sırayla)

1. **"Blur eşiği kırmızı patch'leri haksız eliyor"** — kırmızı doku düşük
   kontrastlı diye düşünüldü. Görsel kontrol çürüttü: elenenler doku değil,
   katı mürekkep.
2. **"Bu slaytlar sadece yoğun eozinle boyanmış, dışlamaya gerek yok"** —
   eozin sıralamasına bakıldı, kohort içinde süreklilik görüldü. Ama
   sıralama "neden kırmızı" sorusunu cevaplamıyordu.
3. **"Kalem yayılması var, mesafe filtresi uygulayalım"** — kaleme uzaklık
   ile eozin arasında negatif korelasyon bulundu (−0,42 / −0,28 / −0,42).
   Ancak mesafe kalem *çekirdeğine* göre ölçülüyordu; yayılma bölgesi
   maskede olmadığı için bir patch "uzak" sayılıp tam yayılmanın içinde
   olabiliyordu.
4. **"Kohort geneli eozin eşiği koyalım"** — p95'te `85-8276` 143 patch'e
   düşüyordu. Slayt seviyesi bir sorunu patch seviyesinde çözme denemesi.

Doğru teşhise ancak aday havuzuna bakılınca ulaşıldı: havuzda 8.416 aday
vardı, yani slaytta bol temiz doku vardı. Sorun seçimde değil, **maskede**.

### Ölçüm

`TCGA-85-8276`, doku içinde yakalanan oran:

| Ölçüt | Yakalanan |
|---|---:|
| `resid > 0.55` (eski) | %11,5 |
| `resid > 0.45` | %37,5 |
| `neg < −0,05` | %28,6 |
| **Birleşim** | **%66,0** |

Ayrıca doygunluk koşulu (`sat > 60`) ölçüldü: doku piksellerinin %100'ü
geçiyordu. Hiçbir katkısı yok, kaldırıldı.

### Çözüm

```python
pen = (resid > 0.45) | (neg < -0.05)
```

| Ölçüt | Yakaladığı | Fiziksel sebep |
|---|---|---|
| `resid` | Kırmızı / magenta | Eozin yönüne yakın; temsil edilebiliyor ama kalıntı bırakıyor |
| `neg` | Yeşil / mavi / turkuaz | H&E düzleminin ters yönü; ancak negatif katsayıyla temsil edilebiliyor |

### Sonuç

18 slayt yeniden çıkarıldı (LUAD 15, LUSC 3). Seçim kriteri:
`extract_qc.csv`'de `n_blur_rej > 1000`.

| Slayt | Eski | Yeni |
|---|---:|---:|
| TCGA-86-8055 | 557 | 2.000 |
| TCGA-86-8073 | 620 | 2.000 |
| TCGA-86-6851 | 623 | 2.000 |
| TCGA-86-7954 | 589 | 2.000 |
| TCGA-85-8276 | 549 | 2.000 |
| TCGA-85-8071 | 573 | 2.000 |
| TCGA-85-8052 | 655 | 2.000 |
| TCGA-62-A46P | 452 | 1.989 |

Blur elemesi %80'den %5–30'a düştü.

### Neden tüm kohorta uygulanmadı

Birleşim maskesi kalemsiz slaytlarda normal dokuyu da eliyor
(`TCGA-05-4433`: %19,7 → %8,0). Yalnızca ölçülmüş bir kritere göre seçilen
18 slayda uygulandı.

### Ders

Bir filtrenin çok eleme yapması, o filtrenin hatalı olduğu anlamına gelmez.
Blur doğru çalışıyordu; hatalı olan ondan önceki aşamaydı. Elenen örneklere
bakmadan dört farklı yanlış teşhis kondu.

---

## B.11 · Patch seviyesi doku oranı boyama yoğunluğunu ölçüyordu

### Sorun

Patch'ler thumbnail maskesinde **≥ 0.80** doku oranı ile seçilmişti. Ancak
kendi pikselleri üzerinde ölçüldüğünde ortanca **0.48** çıktı.

Daha kritik olan: en düşük skorlu patch'ler görsel denetimde **en iyi
dokulardı.**

| `tissue_patch` | İçerik |
|---|---|
| 0.08 – 0.12 | Açık boyanmış, yoğun çekirdekli, net morfolojili tümör |
| 0.48 (ortanca) | Koyu boyanmış slaytlar |

Metrik doku varlığını değil **boyama koyuluğunu** ölçüyordu.

Kod tarafında sebep:

```python
tissue = ((sat > 20) & (val > 30) & (val < 235)).mean()
                                     ^^^^^^^^^^^
```

Açık pembe-magenta boyanmış patch'lerin parlaklık değerleri 235'i aşıyor ve
"boş cam" sanılıp eleniyorlardı. Doku orada, sadece açık renkte.

### Neden filtrelemedik

Boyama yoğunluğu **kaynak merkeze bağlı sistematik bir değişkendir.**

```
Açık boyayan merkez  → düşük skor → patch'leri elenir → küçük bag
Koyu boyayan merkez  → yüksek skor → patch'leri kalır → büyük bag
```

Bu metrikle filtreleme **merkez-bağımlı bir seçilim** yaratır ve çalışmanın
başlıca metodolojik riski olan batch effect'i büyütür. Patch seviyesinde
0.80 eşiği uygulansaydı örneklemin **%97,9'u** elenecekti.

### Çözüm

Metrik kaydediliyor ama filtreleme için kullanılmıyor.

```bash
python run.py patch-extract --slides ... --out ... \
  --blur-min 100 --ink-max 0.05
  # --tissue-min verilmiyor (varsayilan 0 = kapali)
```

Doku filtrelemesi thumbnail aşamasında zaten yapılıyor (≥ 0.80). Patch
seviyesinde ikinci bir eşik gerekmiyor.

`tissue_patch` değeri `final_patch_manifest.csv` içinde saklanıyor ve
merkezler arası boyama farkının nicel ölçüsü olarak stain ablation
analizinde kullanılabilir.

### Nihai karar — patch seviyesinde filtre yok

Görsel denetim sonrası iki filtre de devre dışı bırakıldı:

| Filtre | Neden |
|---|---|
| `ink_max` | En yüksek skorlu 15 patch'in 13'ü hematoksilen yoğun gerçek doku |
| `tissue_min` | Boyama yoğunluğunu ölçüyor, merkez-bağımlı seçilim yaratıyor |

Tek aktif patch-seviyesi filtre `blur_min = 100`. Mürekkep manifest
aşamasında slayt seviyesinde zaten elendi.

Her iki metrik `final_patch_manifest.csv` içinde saklanıyor ve merkezler
arası boyama farkının nicel ölçüsü olarak ablation analizlerinde
kullanılabilir.

### Ek not — antrakoz

Görsel denetimde bazı patch'lerde siyah noktacıklar görüldü. Bunlar
**antrakotik pigment** — akciğerde biriken karbon partikülleri. Gerçek
biyolojik yapıdır, artefakt değildir ve filtrelenmemelidir. Doygunlukları
düşük olduğu için ink maskesine takılmıyorlar.

---

## B.12 · Pilot örneklemi temsili değildi

### Sorun

İlk pilot `--limit 10` ile çalıştırıldı. Bu parametre alfabetik sırada ilk
10 WSI'ı alıyor — tamamı `TCGA-05-*`, yani **tek bir kaynak merkezden.**

Ölçülen her şey (mürekkep oranı, blur dağılımı, boyama yoğunluğu) o merkezin
tarayıcısına ve boyama protokolüne ait. Diğer 32 merkezde farklı çıkabilir.

Eşikleri bu ölçüme dayanarak seçmek yanıltıcı olurdu.

### Çözüm

Merkez-katmanlı pilot seçeneği eklendi.

`tez/patches.py`:

```python
if per_tss:
    # MERKEZ-KATMANLI PILOT: her TSS'ten per_tss adet WSI.
    # Alfabetik --limit tek merkezden ornek alir (TCGA-05-* gibi) ve
    # murekkep/boyama dagilimini temsil etmez.
    by_tss = {}
    for w in wids:
        t = by_wsi[w][0].get("tss", "??")
        by_tss.setdefault(t, []).append(w)
    wids = []
    for t in sorted(by_tss):
        wids.extend(by_tss[t][:per_tss])
```

Kullanım:

```bash
python run.py patch-extract --slides ~/tcga_data/LUAD \
  --out ~/tez/out/LUAD --per-tss 1 --cap 200
```

Rapor artık merkez kırılımı da basıyor:

```python
if max(tv) - min(tv) > 0.15:
    print("  [!] tissue_patch merkeze gore belirgin degisiyor "
          "-> bu metrikle filtreleme merkez-bagimli secilim yaratir")
```

---

## B.13 · Hasta sızıntısı (data leakage)

### Sorun

Aynı hastanın iki slaytı train ve test setine ayrı düşerse test doğruluğu
yapay olarak şişer. Daha kötüsü, aynı slayttan gelen patch'ler bölünürse
çok ciddi sızıntı oluşur — patch'ler aynı biyolojik örnekten gelir.

### Çözüm — iki katmanlı

**1. Kohort seviyesinde:** hasta başına tek WSI.

`tez/cohort.py`:

```python
def dedupe_patients(cohort, excluded):
    """
    Hasta basina tek WSI kurali.
    Ayni hastadan birden fazla dosya kaliyorsa dosya adina gore
    deterministik olarak ilki tutulur.
    """
    by_patient = defaultdict(list)
    for r in cohort:
        by_patient[r["patient"]].append(r)

    kept, dropped = [], []
    for pid, rows in by_patient.items():
        rows.sort(key=lambda r: r["file_name"])
        kept.append(rows[0])
        for r in rows[1:]:
            dropped.append({..., "reason": "duplicate_wsi_same_patient"})
    return kept, excluded + dropped
```

**2. Doğrulama seviyesinde:** `run.py check` test ediyor.

```python
dup = [p for p, n in Counter(r["patient"] for r in rows).items() if n > 1]
check("hasta basina tek WSI", not dup, ...)

n_wsi, n_pat = len(rows), len({r["patient"] for r in rows})
check("WSI sayisi = hasta sayisi", n_wsi == n_pat, f"({n_wsi} / {n_pat})")
```

Sonuç: **468 WSI = 468 hasta** (birebir).

**3. Split seviyesinde:** manifest'te `patient` kolonu tutuluyor; split
`GroupKFold` / `StratifiedGroupKFold` ile `groups=patient` üzerinden
yapılacak. Patch seviyesinde split kesinlikle yapılmayacak.

---

## B.14 · Kaynak merkez ile sınıf etiketi arasında %100 ilişki

### Sorun

TCGA barkodunun ikinci alanı dokuyu gönderen hastaneyi gösterir (TSS —
tissue source site).

**İlk ölçüm** GDC manifest'lerinden yapılmıştı. **İkinci ölçüm** her iki
kohort da indirilip QC'den geçtikten sonra, gerçek `cohort.csv` dosyaları
üzerinden tekrarlandı (`run.py sites`):

```
LUAD: 33 merkez / 468 WSI
LUSC: 35 merkez / 478 WSI
ortak merkez: 0
```

| | LUAD | LUSC |
|---|---:|---:|
| Kohort WSI | 468 | 478 |
| Farklı kaynak merkez | 33 | 35 |
| **Ortak merkez** | **0** | **0** |

En büyük merkezler:

| LUAD | WSI | % | | LUSC | WSI | % |
|------|----:|---:|---|------|----:|---:|
| 55 | 91 | 19,4 | | 85 | 50 | 10,5 |
| 44 | 47 | 10,0 | | 77 | 44 | 9,2 |
| 86 | 34 | 7,3 | | 56 | 41 | 8,6 |
| 78 | 31 | 6,6 | | 66 | 38 | 7,9 |
| 49 | 28 | 6,0 | | 22 | 37 | 7,7 |

LUAD'da merkez dağılımı daha yoğunlaşmış (en büyük merkez %19,4), LUSC'ta
daha dengeli (%10,5).

Hiçbir kaynak merkez hem LUAD hem LUSC slaytı göndermemiş. Bu, manifest
seviyesinde tahmin edilen değil, **her iki kohortun tamamı indirilip QC'den
geçtikten sonra ölçülmüş** bir bulgudur. Dolayısıyla
**kaynak merkez bilgisi sınıf etiketini tam olarak belirliyor.**

Merkezler arasında H&E boyama protokolü, tarayıcı modeli, fiksasyon süresi
ve kesit kalınlığı sistematik olarak farklılaşır. Bu koşulda model,
morfolojiyi hiç öğrenmeden yalnızca renk dağılımı ve tarayıcı imzası
üzerinden yüksek doğruluk elde edebilir.

### Neden kritik

1. **Rapor edilen doğruluk yorumlanamaz.** %95'in ne kadarı morfoloji, ne
   kadarı merkez tanıma — bilinmez.
2. **Dış geçerlilik yok.** Başka bir hastanenin slaytında model çöker.
3. **Açıklanabilirlik bozulur.** Attention haritası tümörü değil renk
   özelliğini gösterir.
4. **Model karşılaştırması anlamsızlaşır.** Görev boyama tonundan
   çözülebiliyorsa üç model de kolayca yüksek skor alır ve aralarında
   anlamlı fark çıkmaz — yani tezin araştırma sorusu cevapsız kalır.

### Çözülemez, ama ölçülebilir

Veri böyle; sorun ortadan kaldırılamaz. Kodda alınan önlem, ölçümü mümkün
kılmak: `tss` kolonu kohorttan manifest'e kadar taşınıyor.

`tez/patches.py`, `build_manifest()`:

```python
mw.writerow({
    "wsi_id": wid,
    "patient": row.get("patient", ""),
    "tss": row.get("tss", ""),      # <-- site-aware analizler icin
    "label": row.get("label", ""),
    ...
})
```

Merkez örtüşmesini ölçen komut — `tez/cohort.py`:

```python
def compare_sites(cohort_a, label_a, cohort_b, label_b):
    a = Counter(r["tss"] for r in cohort_a)
    b = Counter(r["tss"] for r in cohort_b)
    shared = set(a) & set(b)
    if not shared:
        print("  [!!] ORTAK MERKEZ YOK.")
        print("       Kaynak merkez bilgisi sinif etiketini %100 belirliyor.")
        print("       -> site prediction baseline ve stain ablation zorunlu.")
```

Kullanım (LUSC hazır olduğunda):

```bash
python run.py sites --cohort out/LUAD/cohort.csv --cohort out/LUSC/cohort.csv
```

### Planlanan analizler

1. **Site prediction baseline** — aynı embedding'lerden TSS kodunu tahmin
   et. Yüksek doğruluk, embedding'lerin merkez bilgisi taşıdığını gösterir.
2. **Leave-two-sites-out** — her sınıftan birer merkez dışarıda bırakılarak
   dengeli değerlendirme. (Merkezler tek sınıflı olduğu için klasik LOSO
   uygulanamaz.)
3. **Stain normalization ablation** — Macenko ile ve olmadan karşılaştırma.
   Doğruluktaki düşüş, renk bağımlılığının doğrudan ölçüsüdür.
4. **Boyut indirgeme görselleştirmesi** — UMAP veya TDA Mapper grafiği bir
   kez sınıfa, bir kez TSS'e göre renklendirilir. Kümelenmenin morfolojiye
   mi merkeze mi göre oluştuğunu doğrudan gösterir.

### Literatür

**Howard FM et al. (2021).** *The impact of site-specific digital histology
signatures on deep learning model accuracy and bias.* Nature Communications
12:4423.

TCGA'da altı kanser tipinde 3.000'den fazla hasta üzerinde, merkezler arası
görüntü farklarının derin öğrenme ile kolayca tespit edilebildiği
gösterilmiştir. Kritik bulgu: **yaygın kullanılan renk normalizasyonu ve
augmentation yöntemlerine rağmen merkez tespiti mümkün olmaya devam
etmektedir.** Yani stain normalization bu sorunu çözmüyor.

Önerdikleri çözüm, modellerin aynı merkezden gelen örneklerle hem eğitilip
hem doğrulanmamasını sağlayan bir kuadratik programlama yöntemidir —
buradaki leave-two-sites-out yaklaşımı bunun uyarlanmış halidir.

---

## B.15 · Patch seviyesinde etiket yok

### Sorun

TCGA hasta düzeyinde tanı veriyor:

```
TCGA-05-4244 → LUAD
```

Patch düzeyinde etiket vermiyor:

```
patch_0001 → tümör mü? stroma mı? normal akciğer mi?   ← bilinmiyor
```

"Hasta LUAD ise tüm patch'leri LUAD" demek yanlış olur. Bir LUAD slaytındaki
2000 patch'in tahminen %20-40'ı tümör, geri kalanı stroma, normal akciğer,
damar, lenfosit ve nekrozdur. Normal akciğer dokusu LUAD hastasında da LUSC
hastasında da aynı görünür:

```
Normal alveol (LUAD hastasından)  → "LUAD"
Normal alveol (LUSC hastasından)  → "LUSC"     ← aynı görüntü, farklı etiket
```

Model bunu çözemez; ezberlemeye çalışırken tutunacağı tek şey boyama tonu
kalır — yani doğrudan B.14'teki merkez imzası problemine düşer.

### Çözüm — MIL

Problem **weakly supervised WSI classification** olarak ele alınıyor:

```
WSI / hasta   = bag (torba)
patch         = instance
etiket        = yalnızca bag seviyesinde
```

Attention mekanizması her patch'e bir ağırlık veriyor ve hangi patch'lerin
kararı belirlediğini model kendi öğreniyor:

```
patch_0001 → 0.003   (normal alveol, önemsiz)
patch_0157 → 0.041   (tümör bölgesi, önemli)
patch_0894 → 0.002   (stroma)
                ↓
   ağırlıklı ortalama → hasta temsili → LUAD
```

### Pipeline

```
patch (224×224)
    ↓  DINOv2 / UNI / Virchow   (frozen backbone, eğitilmiyor)
embedding
    ↓  ABMIL                     (eğitilen tek bileşen)
hasta seviyesi tahmin
```

Foundation model burada yalnızca özellik çıkarıcı. Üç ayrı model
eğitilmiyor; üç ayrı embedding seti çıkarılıp **aynı ABMIL** uygulanıyor.
Değişen tek değişken embedding'in kaynağı.

### Not

Patch seviyesinde çalışmak ancak patch etiketi olan veri setlerinde
mümkündür (LC25000, CAMELYON16 gibi). TCGA LUAD/LUSC ile çalışan literatürün
tamamı MIL kullanır.

---

## B.16 · Önceki Macenko denemesinde sınıf-spesifik referans

### Sorun

Daha önceki bir çalışmada TCGA, LC25000'e Macenko yöntemiyle normalize
edilmişti. Ancak referans görüntü sınıfa göre seçilmişti:

| Sınıf | Referans kaynağı |
|---|---|
| LUAD | LC25000 `lung_aca` klasöründen |
| LUSC | LC25000 `lung_scc` klasöründen |

Bu, etiketi ön işlemeye sızdırır. LUAD görüntüleri LUAD-benzeri bir renk
dağılımına, LUSC görüntüleri LUSC-benzeri bir dağılıma çekilir. Model artık
rengin kendisinden sınıfı çıkarabilir — çünkü renk, etiket kullanılarak
ayarlanmıştır.

Ayrıca elde edilen %93 doğruluk **LC25000 test seti** üzerindedir, TCGA
üzerinde değil. LC25000 sentetik augmentasyonla üretilmiş bir settir ve
TCGA'dan belirgin şekilde kolaydır.

### Çözüm

Yeni kurguda:

- Stain normalization primary deneyde **uygulanmıyor**
- Ablation olarak değerlendirilecek: tek referans görüntü, sınıftan bağımsız
- Amaç batch effect'i çözmek değil (Howard et al. bunun çalışmadığını
  gösteriyor), modelin renge bağımlılığını **ölçmek**

---

## B.17 · `hash()` ile seed üretimi deterministik değil

### Sorun

`build_manifest()` içinde örnekleme seed'i şöyle üretiliyordu:

```python
rng = np.random.default_rng(abs(hash((cfg["seed"], wid))) % (2 ** 32))
```

Python 3'te string hash'i **her süreçte farklı tuzlanır** (PYTHONHASHSEED
varsayılan olarak rastgele). Yani manifest iki kez üretilirse aynı slayttan
farklı patch'ler seçilir.

Üretilen patch'lere zarar vermiyor — çıkanlar geçerli. Ancak "sabit seed
ile deterministik örnekleme" iddiası tutmuyor.

### Çözüm

```python
import hashlib
h = hashlib.blake2b(f"{cfg['seed']}:{wid}".encode(), digest_size=8).digest()
rng = np.random.default_rng(int.from_bytes(h, "big") % (2 ** 32))
```

`blake2b` süreçten bağımsız, her zaman aynı sonucu verir.

LUAD patch'leri yeniden üretilmeyecek; düzeltme LUSC'a geçmeden önce
uygulanacak.

---

## B.18 · Tamamen odak dışı slayt

### Sorun

`TCGA-44-7661` çıkarma aşamasında **sıfır patch** üretti. Manifest'te
koordinatları vardı ancak hiçbiri blur eşiğini geçemedi.

Manuel kontrolde slaydın tamamının odak dışı olduğu doğrulandı — patch
alınabilecek net bölge yok.

### Çözüm

Kohorttan çıkarıldı. Gerekçe `excluded.csv`'ye eklendi:
`no_patch_passed_qc_out_of_focus`.

```
480  GDC manifest DX kaydı
-1   ikincil tanısal slayt (DX10)
-1   nüks tümör (02Z)
-10  küçültülmüş türev görüntü
-1   QC'yi geçen patch yok (odak dışı)
────
467  analiz kohortu · 467 WSI = 467 hasta
```

### Yan bulgu

Bu, blur eşiğinin doğru seçildiğinin kanıtı: kullanılabilir slaytlarda
neredeyse hiç eleme yapmıyor, kullanılamaz slaytta %100 eliyor.

---

# BÖLÜM C · UYGULANAN KONTROL MEKANİZMALARI

Her aşama bir doğrulama adımıyla eşleşiyor. Herhangi bir kontrol
başarısızsa çıkış kodu 1 döner ve sonraki aşamaya geçilmez.

## C.1 · Birim testi

`tez/barcode.py` kendi testini içeriyor, kohort kurulmadan önce otomatik
çalışıyor:

```
barcode.py kendi testi:
  [OK ] TCGA-05-4244-01Z-00-DX1.d4ff32cd.svs   -> True / ok
  [OK ] TCGA-49-4488-01Z-00-DX10.930e7dc2.svs  -> False / secondary_diagnostic_slide(DX10)
  [OK ] TCGA-50-5066-02Z-00-DX1.A2CA045A.svs   -> False / non_primary_sample(02:Recurrent Solid Tumor)
  [OK ] TCGA-85-A512-01A-01-TSA.C10905AD.svs   -> False / barcode_parse_failed
  [OK ] TCGA-44-6146-11A-01-BS1.8db8c567.svs   -> False / not_diagnostic_slide(BS)
  [OK ] bozuk_dosya_adi.svs                    -> False / barcode_parse_failed
  -> gecti
```

## C.2 · Aşama doğrulamaları

`run.py check` 13 kontrol çalıştırıyor:

```
[OK ] tum barkodlar ayristirilabilir
[OK ] hepsi sample_type=01 (primary tumor)
[OK ] hepsi DX1 (DX10/DX11 yok)
[OK ] hasta basina tek WSI
[OK ] WSI sayisi = hasta sayisi  (468 / 468)
[OK ] beklenen WSI sayisi = 468  (gercek: 468)
[OK ] kohort/dislanan kesisimi bos
[OK ] her dislamanin gerekcesi var
[OK ] wsi_qc kohortla ayni WSI'lari kapsiyor
[OK ] tum WSI openslide ile acilabiliyor
[OK ] tum WSI'da MPP okunabiliyor
[OK ] tum WSI status=ok
[OK ] hicbir WSI yukari olceklenmiyor
```

## C.3 · Determinizm

Örnekleme seed ile `wsi_id`'nin birleşimine dayanıyor:

```python
rng = np.random.default_rng(abs(hash((cfg["seed"], wid))) % (2 ** 32))
```

Aynı WSI her çalıştırmada aynı patch'leri veriyor — işlem sırası değişse,
pipeline yarıda kesilse bile.

## C.4 · Parametre kaynağı etiketleme

Her eşik ve parametre kaynağıyla birlikte kaydediliyor:

| Etiket | Anlamı |
|---|---|
| `[ölçüm]` | Bu veride hesaplandı |
| `[literatür]` | Kaynak var |
| `[tahmin]` | Doğrulanması gerekiyor |

`[tahmin]` etiketli hiçbir parametre kalibre edilmeden kullanılmıyor.

## C.5 · Geri alınabilirlik

Tüm aday koordinatlar `candidates/` altında saklanıyor. Patch sayısı veya
doku eşiği değişirse yeniden tarama gerekmiyor, yalnızca manifest ve çıkarma
aşamaları tekrarlanıyor.

---

# BÖLÜM D · KOD HARİTASI

```
~/tez/
├── run.py                    ← tek giriş noktası, alt komutlar
├── tez/
│   ├── barcode.py            ← TCGA barkod mantığı (saf, I/O yok, self-test)
│   ├── cohort.py             ← kohort kurma + dışlama kaydı + site karşılaştırma
│   ├── wsi.py                ← MPP ölçümü, level seçimi, okuma testi
│   └── patches.py            ← aday tarama, manifest, çıkarma, QC
├── docs/
│   ├── sorun_cozum.md        ← bu dosya
│   └── danisman_raporu.md    ← teknik durum raporu
└── out/LUAD/                 ← üretilen dosyalar (git'e girmez)
    ├── cohort.csv
    ├── excluded.csv          ← tezin dışlama tablosu
    ├── wsi_qc.csv
    ├── scan_summary.csv
    ├── candidates/*.csv      ← tüm aday koordinatlar
    ├── canonical_patch_manifest.csv   ← üç modelin ortak kaynağı
    ├── final_patch_manifest.csv       ← QC sonrası yazılanlar
    └── patches/<hasta>/<wsi>/*.jpg
```

## D.1 · Sorunların kod karşılığı

| Sorun | Modül | Fonksiyon / alan |
|---|---|---|
| B.1 DX10 filtresi | `barcode.py` | `parse()`, `BARCODE_RE`, `COHORT_SLIDE_NO` |
| B.2 Örnek tipi | `barcode.py` | `cohort_reason()`, `COHORT_SAMPLE_TYPE` |
| B.3 Dışlama gerekçesi | `cohort.py` | `scan()` |
| B.4 Türev görüntüler | `wsi.py` | `read_mpp()`, `min_base_mpp` / `max_base_mpp` |
| B.5 Fiziksel ölçek | `patches.py` | `scan_slide()`, `extract()` |
| B.6 Level seçimi | `wsi.py` | `choose_level()` |
| B.7 Tavansız tarama | `patches.py` | `scan_slide()`, `_window_means()` |
| B.8 Okuma süresi | `wsi.py` | `read_test()` |
| B.9 Blur sırası | `patches.py` | `extract()` |
| B.10 Ink / hematoksilen | `patches.py` | `patch_metrics()` |
| B.11 Doku metriği | `patches.py` | `patch_metrics()`, `extract_report()` |
| B.12 Temsili pilot | `patches.py` | `extract(per_tss=...)` |
| B.13 Hasta sızıntısı | `cohort.py` | `dedupe_patients()` |
| B.14 Merkez etkisi | `cohort.py` / `patches.py` | `compare_sites()`, manifest `tss` kolonu |

## D.2 · Çalıştırma sırası

```bash
python run.py cohort         --data ~/tcga_data/LUAD --label LUAD \
                             --primary <slayt_klasoru> --out out/LUAD
python run.py check          --out out/LUAD --expect-wsi 468

python run.py wsi-qc         --slides ~/tcga_data/LUAD --out out/LUAD
python run.py check          --out out/LUAD --expect-wsi 468

python run.py patch-scan     --slides ~/tcga_data/LUAD --out out/LUAD
python run.py patch-manifest --out out/LUAD --tissue-min 0.80 --cap 2000
python run.py patch-extract  --slides ~/tcga_data/LUAD --out out/LUAD \
                             --per-tss 1 --cap 200          # pilot
python run.py patch-extract  --slides ~/tcga_data/LUAD --out out/LUAD \
                             --cap 2000 --blur-min 100 --ink-max 0.05
```

---

# BÖLÜM E · İKİ KOHORT — NİHAİ DURUM

| | WSI | Patch |
|---|---:|---:|
| **LUAD** | **466** | **900.828** |
| **LUSC** | **478** | **934.875** |
| Toplam | 944 | 1.835.703 |

## Dışlamalar

**LUAD** (480 GDC kaydı → 466)

| n | Gerekçe |
|---:|---|
| 10 | Küçültülmüş türev görüntü (MPP yok, L0 uzun kenar = 10.000 px) |
| 1 | İkincil tanısal slayt (DX10) |
| 1 | Nüks tümör (02Z) |
| 1 | Odak dışı — QC'yi geçen patch üretilemedi |
| 1 | Trikrom boyama (H&E değil) |

**LUSC** (512 tanısal slayt → 478)

34 ikincil tanısal slayt (DX2–DX9) indirme aşamasında filtrelendi, dışlama
gerekmedi.

## Kesinleşen protokol

| Parametre | Değer | Kaynak |
|---|---|---|
| target_mpp | 0,50 µm/px | `[literatür]` |
| patch_size | 224 px (112 × 112 µm) | `[literatür]` |
| stride | 112 µm, overlap yok | `[literatür]` CLAM |
| tissue_min (thumbnail) | 0,80 | `[ölçüm]` |
| kalem maskesi (standart) | `resid > 0.55` | `[ölçüm]` |
| kalem maskesi (18 slayt) | `resid > 0.45 \| neg < −0.05` | `[ölçüm]` |
| pen_max | 0,01 | `[ölçüm]` |
| max_patches | 2.000 | `[ölçüm]` |
| blur_min | 100 (resize sonrası) | `[ölçüm + görsel]` |
| neg_frac (patch) | 0,15 | `[ölçüm + görsel]` |
| stain normalization | yok (ablation planlı) | `[literatür]` |

## Düşük bag'ler

| Slayt | Patch | Sebep |
|---|---:|---|
| TCGA-55-8207 | 78 | 2 mm² doku |
| TCGA-67-4679 | 289 | Doku %17, blur 5.776 eledi |
| TCGA-98-8022 | 290 | Tarama kalitesi (blur p50 = 869 vs kohort 1.089) |
| TCGA-55-8301 | 326 | Az doku |
| TCGA-56-8082 | 381 | Parçalı doku |

---

# BÖLÜM F · AÇIK KONULAR

| # | Konu | Durum |
|---|---|---|
| 1 | Deterministik hash düzeltmesi (B.17) | LUSC öncesi uygulanacak |
| 2 | Bag boyutu dengesizliği: 20 WSI 2000'in altında (78-890) | Karar embedding aşamasında |
| 3 | 112 µm FOV yeterli mi | Sonuçlar görüldükten sonra |
| 4 | Nekroz nasıl ele alınacak | Açık |
| 5 | Macenko ablation: `patches_macenko/` ikinci set | Embedding sonrası |
| 6 | Açıklanabilirlik yöntemi (attention / UMAP / TDA Mapper) | Sonuçlar görüldükten sonra |
| 7 | Hasta başına tek slayt mı, tüm DX slaytları mı | Danışman leakage olmadığını belirtti; yeniden değerlendirilecek |
| 8 | Site prediction baseline, leave-two-sites-out | Embedding sonrası |
| 9 | Native 20× vs 40× grup farkı (LUSC'ta %12, LUAD'da %4) | Değerlendirmede kontrol edilecek |

---

*Son güncelleme: 27 Eylül 2026 · LUAD patch çıkarma tamamlandı, embedding aşamasına geçiliyor*
