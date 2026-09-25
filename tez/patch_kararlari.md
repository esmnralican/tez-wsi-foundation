# PATCH ÇIKARMA — KARAR VE HATA KAYDI

> Bu dosya her patch kararında güncellenir. Amaç: "neden böyle yaptık?"
> sorusunun cevabını hatırlamaya değil kayda dayandırmak.
>
> Üç bölüm: **neye dikkat ediyoruz** · **ne değiştirdik** · **nerede hata yapıldı**

---

# BÖLÜM 1 — NEYE DİKKAT EDİYORUZ

Patch çıkarırken sürekli kontrol edilmesi gereken 10 nokta.

## 1.1 · Fiziksel ölçek, piksel boyutu değil

Aynı piksel boyutundaki iki patch aynı doku alanını göstermeyebilir.

```
224 px @ 0.25 µm/px  →   56 × 56 µm
224 px @ 0.50 µm/px  →  112 × 112 µm     ← bizim hedefimiz
```

**Kural:** Her WSI aynı pyramid level'ından değil, aynı **fiziksel ölçeğe**
(0.50 µm/px) normalize edilerek örneklenir.

Ölçülen durum (468 LUAD WSI):

| Base MPP | WSI | Okunan | Sonuç |
|---|---:|---|---|
| 0.252 | 414 | 444 px | → 224 px |
| 0.233 | 35 | 482 px | → 224 px |
| 0.500 | 19 | 224 px | resize yok |

Farklı piksel sayısı okunuyor, aynı 112 µm alan elde ediliyor.

## 1.2 · Level seçimi her zaman aşağı yönlü

Hedeften **daha yüksek çözünürlüklü** (daha küçük MPP) level seçilir,
yazılımda küçültülür. Yukarı ölçekleme yapılmaz.

Sebep: interpolasyon bilgi kazandırmaz, sahte detay üretir. Modelin
gerçekte var olmayan hücresel yapı görmesine yol açar.

Ölçülen: 468 WSI'ın **tamamı level 0**'dan okunuyor. Aperio pyramid'i
tipik olarak `[1, 4, 16, 32]` — level 1 zaten 1.0 µm/px, hedefin iki katı
kaba. Yani 0.50 arada kalıyor, level 0'dan okuyup küçültmek zorunlu.

## 1.3 · Blur QC her zaman resize SONRASI

Laplacian varyansı çözünürlüğe bağlı bir ölçüdür.

```
40x WSI        →  444 px okunur  →  224'e küçültülür
native 20x WSI →  224 px okunur  →  resize yok
```

Blur eşiği resize öncesi uygulanırsa iki grup farklı ölçekte değerlendirilir.
19 native-20x slaytın patch'leri sistematik olarak haksız elenebilir.

**Kural:** önce resize, sonra blur.

## 1.4 · Mürekkep ≠ hematoksilen

H&E'de çekirdekleri boyayan hematoksilen **mor-mavidir**; RGB'de mavi kanalı
baskındır. "Mavi oranı yüksek" kuralı, çekirdek yoğun tümör dokusunu
mürekkeple aynı kefeye koyar — yani en değerli patch'leri siler.

Ölçülen hue aralıkları (OpenCV 0-179 ölçeği):

| | Hue | Doygunluk | Laplacian var. |
|---|---|---|---|
| Turkuaz artefakt | ~80-95 | orta | **≈ 0** |
| Mavi mürekkep | ~105-125 | çok yüksek | düşük |
| **Hematoksilen** | **~125-150** | orta | **yüksek** |

**Asıl ayırt edici renk değil, dokudur.** Mürekkep ve artefakt düzdür;
gerçek doku kenar ve yapı içerir.

Daha sağlam ölçü: **H&E stain deconvolution kalıntısı.** Patch'in rengi H&E
boya vektörleriyle açıklanabiliyor mu? Doku → düşük kalıntı. Mürekkep →
yüksek kalıntı.

## 1.5 · Doku oranı nerede ölçülüyor

Thumbnail üzerinde ölçmek kabadır:

```
Slayt 68.000 px → thumbnail 2048 px → ölçek ~33x
444 px'lik patch → thumbnail'de 13 × 13 piksel
```

13×13 pencerede %80 doku şart koşmak gürültülüdür ve **doku kenarlarını
sistematik eler**. Tümör-stroma sınırı LUAD/LUSC ayrımının en bilgilendirici
bölgesi olabilir.

**Çözüm:** iki aşamalı ölçüm.
- Manifest aşaması: thumbnail oranı, gevşek eşik → ucuz aday eleme
- Extract aşaması: patch'in kendi pikselleri → asıl karar

İkisi de manifest'te saklanır (`tissue_ratio_thumb`, `tissue_ratio_patch`),
aralarındaki fark raporlanabilir.

## 1.6 · Örnekleme deterministik olmalı

Seed, `wsi_id` ile birleştirilir. Aynı WSI her çalıştırmada aynı patch'leri
verir — işlem sırası değişse, pipeline yarıda kesilse bile.

Yeniden üretilebilirlik iddiası için gerekli.

## 1.7 · Patch sayısı ile doku kapsamı

```
100 patch × (112 µm)² ≈ 1.25 mm²
Tipik WSI doku alanı   ≈ 200–600 mm²
Kapsanan oran          ≈ %0.2–0.6
```

MIL'in mantığı "çok sayıda instance arasından modelin önemlileri bulması".
100 instance ile attention mekanizmasının seçecek pek bir şeyi olmaz.

Karşılaştırma: ProtoMIL (TCGA-NSCLC) ~3.961 patch/WSI kullanıyor.

**Ölçülen süre** (468 WSI, tek işlemci, 15 ms/patch okuma):

| Patch/WSI | Süre | Disk |
|---:|---|---|
| 100 | 0.2 sa | ~1.4 GB |
| 500 | 1.0 sa | ~7 GB |
| 1000 | 1.9 sa | ~14 GB |

> Bu ölçüm iyimser: test slayt ortasından üç bitişik bölge okuyor, disk
> önbelleğinden faydalanıyor. Gerçek çıkarma dağınık koordinatlardan okur,
> üstüne JPEG kodlama ve maske işi biner. **3-5 katı bekle.**

## 1.8 · Üç model, tek patch seti

```
YANLIŞ                          DOĞRU
WSI                             WSI
 ├── DINOv2 → patch set A        ↓
 ├── UNI    → patch set B       MASTER MANIFEST
 └── Virchow→ patch set C        ├── DINOv2
                                 ├── UNI
                                 └── Virchow
```

Patch seçimi deneysel bir confounder haline gelmemeli. Değişen tek şey
foundation model olmalı.

Bu yüzden pipeline iki aşamalı: manifest (koordinatlar) → extract (görüntü).
Model değişirse manifest **aynen kalır**.

## 1.9 · Hasta sızıntısı

Aynı hastanın patch'leri train ve test'e bölünürse sonuç şişer. Patch
seviyesinde split **asla** yapılmaz.

Split her zaman `patient` kolonu üzerinden, `GroupKFold` /
`StratifiedGroupKFold` ile.

Kohort aşamasında hasta başına tek WSI kuralı uygulandı (468 WSI = 468 hasta),
bu riski baştan azaltıyor.

## 1.10 · Kaynak merkez (TSS) kaydı

LUAD 33 merkez, LUSC 35 merkez, **ortak sıfır**. Merkez bilgisi sınıf
etiketini %100 belirliyor.

Manifest'te `tss` kolonu **mutlaka** olmalı. Yoksa şu analizlerin hiçbiri
sonradan yapılamaz:

- site prediction baseline
- leave-two-sites-out değerlendirme
- stain normalization ablation
- TDA Mapper'da sınıf vs merkez renklendirmesi

---

# BÖLÜM 2 — NE DEĞİŞTİRDİK

Tarih sırasıyla. Her satır: ne → neye, neden.

## 2.1 · Kohort: 479 → 468

| Adım | Sayı | Gerekçe |
|---|---:|---|
| GDC manifest DX kayıtları | 480 | |
| − `TCGA-49-4488-...-DX10` | 479 | Aynı hastanın 2. tanısal slaytı |
| − 10 MPP'siz WSI | 469 | L0 = tam 10.000 px, küçültülmüş türev |
| − `TCGA-50-5066-02Z` | **468** | Nüks tümör (recurrent), primary değil |

468 WSI = 468 hasta, birebir.

## 2.2 · Dışlama gerekçesi: klasör adı → barkod

**Önce:** gerekçe `_excluded_*` klasör adından alınıyordu.
**Sonra:** barkod kesin bir gerekçe veriyorsa o kullanılıyor.

Sebep: `_excluded_duplicate_patient` klasöründe hem DX10 hem 02Z vardı;
ikisi de "secondary_diagnostic_slide" etiketi alıyordu. Farklı şeyler.

Klasör adları elle verilir ve yanılır. Barkod yanılmaz.

## 2.3 · 10 MPP'siz WSI'ın gerekçesi güçlendi

**Önce:** "metadata eksikliği"
**Sonra:** "düşük çözünürlüklü türev görüntü"

Bulgu: onunun da L0 uzun kenarı **tam 10.000 piksel**. Sağlam slaytlar
37.000–78.000 aralığında. Yani orijinal tarama değil, küçültülmüş kopya.
Tahmini çözünürlük ~2–2.5 µm/px (4–6x), hedefin çok altında.

MPP eksikliği semptomdu, sebep değildi.

## 2.4 · max_patches: 100 → (karar bekliyor, 500 öneriliyor)

**Önce:** 100, gerekçe "hız ve disk".
**Sonra:** Hız ölçüldü — 15 ms/patch, 500 patch için ~1 saat. Disk 7 GB.

Gerekçe geçersiz kaldı. Nihai karar manifest dağılımı görüldükten sonra.

## 2.5 · tissue_threshold: 0.80 → (karar bekliyor)

Audit geçersizdi (bkz. 3.3). Manifest aşamasında tavan koymadan tüm adaylar
ve doku oranları kaydedilecek, karar dağılıma bakarak verilecek.

## 2.6 · Blur QC sırası: resize öncesi → resize sonrası

Bkz. 1.3. Farklı MPP'li slaytların aynı ölçekte değerlendirilmesi için.

## 2.7 · Doku oranı ölçüm yeri: sadece thumbnail → thumbnail + patch

Bkz. 1.5. Manifest aşamasında gevşek thumbnail eşiği, extract aşamasında
patch'in kendi piksellerinden asıl karar.

---

# BÖLÜM 3 — NEREDE HATA YAPILDI

Dürüst kayıt. Tekrarlanmasın diye.

## 3.1 · `"DX1" in filename` → DX10'u da yakaladı

**Etki:** Aynı hastanın iki slaytı sete girdi, hasta sızıntısı riski.
**Nerede:** Erken kohort filtresi.
**Düzeltme:** `barcode.py` — slayt numarası int olarak ayrıştırılıyor,
`slide_no == 1` kontrolü yapılıyor. `self_test()` bunu test ediyor.

## 3.2 · Sample type hiç kontrol edilmedi → 02Z (nüks tümör) sızdı

**Etki:** Nüks tümör primary kohortta sayıldı. Farklı biyoloji, tedavi
etkisi taşıyabilir. Ayrıca aynı hastanın ikinci kaydı.
**Nerede:** Kohort tanımı. Notlarda "bir hastada iki DX1 olması normaldir"
diye yanlış yorumlanmıştı.
**Düzeltme:** `cohort_reason()` sample type'ı kontrol ediyor.

## 3.3 · Tissue threshold audit'i matematiksel olarak geçersizdi

```
469 WSI × 100 max_patches = 46.900
Ölçülen: 0.70 → 46.897  ·  0.80 → 46.867
```

**Etki:** "0.70 ve 0.80 arasında %0.064 fark var, 0.80 kullanalım" sonucu
hatalı. Neredeyse her slayt tavana dayanmış; **bağlayıcı kısıt eşik değil,
max_patches tavanıydı.** Eşik ne olursa olsun sayı değişmezdi.
**Düzeltme:** Manifest aşaması tavansız çalışacak, gerçek aday dağılımı
ölçülecek.

## 3.4 · [Claude] Ink hue aralığı hematoksileni kapsıyordu

`build_patches.py` içinde ink maskesi `hue 100-140` olarak yazılmıştı.
Hematoksilen **125-150** aralığında.

**Etki:** Uygulansaydı çekirdek yoğun tümör dokusunu mürekkep sanıp silerdi —
yani en bilgilendirici patch'leri.
**Nerede:** Benim yazdığım ilk patch script'i.
**Düzeltme:** Yeni pipeline'da renk tek başına kullanılmayacak; doku
(Laplacian) ve H&E deconvolution kalıntısı birlikte kullanılacak.

## 3.5 · [Claude] 0.60 doku eşiği uydurulmuştu

"0.80 yerine 0.60 kullan" önerisi ne literatürden ne ölçümden geliyordu.
Kullanıcı sordu, doğrulanamadı.

**Düzeltme kuralı:** Bundan sonra her eşik ve parametre kaynağıyla birlikte
verilecek:
- `[ölçüm]` — bu veride hesaplandı
- `[literatür]` — kaynak var
- `[tahmin]` — doğrulanması gerekiyor

## 3.6 · max_patches = 100 ölçülmeden gerekçelendirildi

"Hızlı olsun" diye seçilmişti ama süre hiç ölçülmemişti. Ölçüldüğünde
15 ms/patch çıktı — 500 patch için bile sadece ~1 saat.

**Ders:** Performans gerekçesiyle verilen kararlar ölçülmeden verilmemeli.

## 3.7 · Uniform turkuaz patch'ler tissue_ratio ≥ 0.80'i geçti

`green_ge_080` klasöründeki 20 patch'in hepsi düz turkuaz — hiç doku yok,
ama doku maskesini geçmişler.

**Etki:** Renk tabanlı thumbnail maskesi uniform renkli background'ı doku
sanıyor.
**Not:** Blur QC devrede olsaydı hepsi elenirdi (Laplacian var ≈ 0). O audit
sırasında blur filtresi yoktu.
**Düzeltme:** Blur QC extract aşamasının zorunlu parçası.

## 3.8 · Klasör adı `_excluded_non_DX1` yanıltıcıydı

İçindeki dosya bir **DX slaytıydı** — sadece 10. tanısal slayt. "non-DX"
değildi.

**Ders:** Klasör adları açıklayıcı olmalı, ama gerekçe kaynağı kod olmalı.

---

# BÖLÜM 4 — AÇIK KONULAR

Henüz karar verilmemiş, manifest dağılımı görülünce karara bağlanacak.

- [ ] `max_patches` — 500 öneriliyor, dağılıma bakılacak
- [ ] `tissue_threshold` — thumbnail için gevşek, patch için asıl eşik
- [ ] `blur_threshold` — 100.0 `[tahmin]`, pilot patch'lerde kalibre edilecek
- [ ] İnk tespiti — deconvolution kalıntısı eşiği `[ölçülecek]`
- [ ] Nekroz nasıl ele alınacak (artefakt değil, gerçek biyoloji)
- [ ] Patch sayısını doku alanına orantılamak mantıklı mı
- [ ] FOV 112 µm yeterli mi — LUAD/LUSC ayrımı doku mimarisine dayanıyor,
      112 µm ≈ 10-15 hücre çapı. Zayıf sonuç çıkarsa ikinci deneme burası.

---

*Son güncelleme: 25 Eylül 2026 — aşama 2 (WSI QC) tamamlandı*
