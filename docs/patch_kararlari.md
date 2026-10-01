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

## 1.4 · Mürekkep ≠ hematoksilen  [ÇÖZÜLDÜ]

H&E'de çekirdekleri boyayan hematoksilen **mor-mavidir**; RGB'de mavi kanalı
baskındır. "Mavi oranı yüksek" kuralı, çekirdek yoğun tümör dokusunu
mürekkeple aynı kefeye koyar — yani en değerli patch'leri siler.

Ölçülen hue aralıkları (OpenCV 0-179 ölçeği):

| | Hue | Doygunluk | Laplacian var. |
|---|---|---|---|
| Turkuaz artefakt | ~80-95 | orta | **≈ 0** |
| Mavi mürekkep | ~105-125 | çok yüksek | düşük |
| **Hematoksilen** | **~125-150** | orta | **yüksek** |

Ton aralığını daraltmak **yetmedi**: boyama yoğunluğu merkeze göre kaydığı
için sabit bir ton sınırı her merkezde aynı yere düşmüyor.

**Nihai çözüm: H&E stain deconvolution kalıntısı.**

Renge değil boya kimyasına bakılıyor:

```
OD = -log10(I/256)
     ↓ hematoksilen + eozin vektörlerinin gerdiği düzleme izdüşüm
kalıntı = açıklanamayan kısım
```

| | Kalıntı | Neden |
|---|---|---|
| Koyu mor çekirdek | düşük | Rengi zaten hematoksilenden geliyor |
| Kırmızı/yeşil kalem | yüksek | H&E dışı pigment |

Boyama yoğunluğundan bağımsız çalışıyor — ton eşiğinin çözemediği şey buydu.

**Eşik: 0.55**  `[ölçüm]` — bkz. 2.8

## 1.4b · Mürekkep kontrolü nasıl çalışıyor (iki aşamalı)

Tek bir filtre yok; iki farklı ölçüm, iki farklı aşamada.

### Aşama 1 — Slayt seviyesi (manifest'ten önce)

Slaydın küçültülmüş önizlemesinde kalem maskesi çıkarılır.

`tez/patches.py` → `tissue_and_texture_mask()`:

```python
od     = -np.log10((thumb + 1.0) / 256.0)
coeffs = od @ _HE_PINV.T              # H ve E konsantrasyonu
recon  = coeffs @ _HE
resid  = ||od - recon|| / ||od||      # H&E ile aciklanamayan kisim

pen = (resid > 0.55) & (sat > 60) & (val > 30)
pen = morfolojik temizlik + dilatasyon
```

Sorulan soru: *"Bu rengi hematoksilen ve eozin karışımıyla üretebilir
miyim?"* Kalem H&E dışı bir pigment olduğu için yüksek kalıntı verir.

Üç koşul birlikte: H&E ile açıklanamıyor **ve** doygun **ve** yeterince koyu.
Renk listesi yok — yeşil, mavi, kırmızı, siyah fark etmez.

Manifest aşamasında `--pen-max 0.01`: patch'in %1'inden fazlası kalem
değiyorsa aday havuzundan çıkar.

| | Elenen aday |
|---|---:|
| LUAD | 47.057 |
| LUSC | 41.448 |

**Eşik 0.55 `[ölçüm]`** — bkz. 2.8. Üç kalemli / dört kalemsiz slayt
karşılaştırmasıyla seçildi; 0.55'te ayrım 28 kat.

### Aşama 2 — Patch seviyesi (çıkarmadan sonra)

Thumbnail çözünürlüğünde kalemin **yayılmış ince gölgesi** kaçabiliyor.
Gerçek patch piksellerinde ikinci ölçüm yapılır: **negatif boya
konsantrasyonu**.

```python
c   = od @ _HE_PINV.T                    # (N,2) H ve E katsayilari
neg = (c.min(axis=1) < -0.05).mean()     # negatif katsayi orani
```

Fiziksel mantık: gerçek H&E dokusunda **hem hematoksilen hem eozin
konsantrasyonu pozitif olmalıdır.** Boya miktarı negatif olamaz. Turkuaz
gibi H&E dışı renkler ancak negatif katsayıyla temsil edilebilir.

Ölçülen ayrım:

| İçerik | negatif% | eozin katsayısı |
|---|---:|---:|
| **Kalem altı turkuaz** | **100,0** | **−0,596** |
| Aynı slayttan ikinci turkuaz | 100,0 | −0,510 |
| Gerçek doku (11 örnek) | 0,0 – 2,8 | +0,04 … +0,52 |

Ara değer yok, bimodal dağılım.

**Eşik 0.15 `[ölçüm + görsel]`** — bkz. 3.12.

### Neden iki aşama

| | Slayt seviyesi | Patch seviyesi |
|---|---|---|
| Yakaladığı | Kalın kalem konturu | Yayılmış ince gölge |
| Ne zaman | Manifest, çıkarmadan önce | Çıkarmadan sonra |
| Metrik | `he_residual > 0.55` | `neg_frac > 0.15` |
| Sonuç | Aday havuzundan çıkar | Manifest'ten düşer, JPEG diskte kalır |
| Maliyet | Ucuz (thumbnail) | Pahalı (her patch okunur, ~25 dk) |

### `he_residual` neden patch seviyesinde çalışmadı

İlk denemede aynı kalıntı ölçüsü patch seviyesinde de kullanılmak istendi.
Ölçüm gösterdi ki turkuaz bir patch **0,1346** kalıntı veriyor — tüm veri
setinin ortancasından (0,275) **düşük**.

Sebep: deconvolution matrisi iki vektör içeriyor (H ve E), üçüncü boyut
serbest. Herhangi bir renk iki vektörün lineer kombinasyonuyla — negatif
katsayılarla bile — yaklaşık temsil edilebiliyor. Yani kalıntı düşük
çıkıyor.

Ölçü matematiksel olarak doğru ama **fiziksel olarak anlamsız**. Doğru soru
"temsil edilebilir mi?" değil, **"fiziksel olarak mümkün mü?"** — negatif
konsantrasyon metriği bunu sorar.

Ayrıca `he_residual` merkeze göre belirgin değişiyordu (TSS 86 p50=0,331 vs
TSS MP p50=0,244); yani boyama yoğunluğunu ölçüyordu, mürekkebi değil.

## 1.4c · Kalem maskesi tek ölçütle çalışmıyor — birleşim gerekti

Aşama 1'deki maske (`resid > 0.55`) yalnızca **kırmızı/magenta** kalemi
yakalıyordu. Yeşil, mavi ve turkuaz kalem ile mürekkep yayılma bölgeleri
aday havuzuna sızıyordu.

### Nasıl fark edildi

Yoğun kalemli slaytlarda patch sayısı anormal düşük kalmıştı. Çıkarma
kayıtlarına bakıldığında sebep göründü:

```
TCGA-85-8276:  8.416 aday -> 3.000 secildi -> 2.451'i BLUR'a takildi
                                           -> 549 patch
```

Elenen patch'ler görsel olarak incelendiğinde **düz kırmızı ve düz turkuaz
kareler** çıktı — doku içermeyen katı mürekkep blokları. Yani blur filtresi,
kalem maskesinin kaçırdığı mürekkebi eliyordu.

### Ölçüm (TCGA-85-8276, doku içinde yakalanan oran)

| Ölçüt | Yakalanan |
|---|---:|
| `resid > 0.55` (eski) | %11,5 |
| `resid > 0.45` | %37,5 |
| `neg < −0,05` | %28,6 |
| **Birleşim** | **%66,0** |

Doygunluk koşulu (`sat > 60`) ölçüldü: doku piksellerinin **%100'ü** zaten
geçiyordu — hiçbir katkısı yok, kaldırıldı.

### Neden iki ölçüt

| Ölçüt | Yakaladığı renk | Fiziksel sebep |
|---|---|---|
| `resid` | Kırmızı / magenta | Eozin yönüne yakın; temsil edilebiliyor ama kalıntı bırakıyor |
| `neg` | Yeşil / mavi / turkuaz | H&E düzleminin ters yönü; ancak negatif katsayıyla temsil edilebiliyor |

Hiçbiri tek başına yetmiyor.

**Nihai maske:**

```python
pen = (resid > 0.45) | (neg < -0.05)
```

### Sonuç

18 slayt birleşim maskesiyle yeniden çıkarıldı (LUAD 15, LUSC 3):

| Slayt | Eski | Yeni |
|---|---:|---:|
| TCGA-86-8055 | 557 | 2.000 |
| TCGA-86-8073 | 620 | 2.000 |
| TCGA-86-6851 | 623 | 2.000 |
| TCGA-85-8276 | 549 | 2.000 |
| TCGA-85-8071 | 573 | 2.000 |
| TCGA-85-8052 | 655 | 2.000 |
| TCGA-62-A46P | 452 | 1.989 |

Blur elemesi 2.451'den 586/107/636'ya düştü — mürekkep artık aday havuzuna
girmiyor.

> **Not:** Bu maske yalnızca etkilenen slaytlara uygulandı, tüm kohorta
> değil. Sebebi 2.11'de.

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

**Ölçüm doğrulandı** (27 Eylül, `run.py sites` — her iki kohort da indirilip
QC'den geçtikten sonra):

```
LUAD: 33 merkez / 468 WSI
LUSC: 35 merkez / 478 WSI
ortak merkez: 0
```

En büyük merkez LUAD'da %19,4 (TSS 55), LUSC'ta %10,5 (TSS 85).

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

## 2.11 · Birleşim maskesi neden tüm kohorta uygulanmadı

Birleşim maskesi (`resid>0.45 | neg<−0.05`) mevcut maskeden belirgin daha
geniş. Kalemsiz slaytlarda test edildiğinde normal dokuyu da elemeye
başlıyor (ör. `TCGA-05-4433`: OD maskesi %19,7 → %8,0).

Bu nedenle yalnızca **blur elemesi 1.000'in üzerinde olan** 18 slayda
uygulandı. Seçim kriteri `extract_qc.csv` üzerinden ölçüldü, keyfi değil.

Kalan slaytlarda mevcut maske (`resid > 0.55`) geçerliliğini koruyor;
onlarda mürekkep zaten yok veya ihmal edilebilir.

## 2.8 · Kalem tespiti: renk tonu → H&E deconvolution, eşik 0.35 → 0.55

**Önce:** HSV ton aralığı (100-140), hematoksileni kapsıyordu.
**Sonra:** H&E deconvolution kalıntısı, eşik 0.55.

Eşik üç kalemli / dört kalemsiz slaytın karşılaştırmasıyla seçildi. İlk
denenen 0.35 değerinde iki grup üst üste biniyordu:

| Slayt | Kalem? | >0.35 | >0.45 | **>0.55** | >0.65 |
|---|---|---:|---:|---:|---:|
| TCGA-86-7953 | var | %90.8 | %80.4 | **%63.1** | %28.9 |
| TCGA-86-8056 | var | %55.0 | %44.0 | **%16.5** | %3.1 |
| TCGA-86-7714 | var | %26.5 | %13.6 | **%5.6** | %1.0 |
| TCGA-75-6214 | yok | %57.0 | %3.7 | **%0.1** | %0.0 |
| TCGA-50-6673 | yok | %78.9 | %5.3 | **%0.2** | %0.0 |
| TCGA-50-5045 | yok | %61.7 | %2.0 | **%0.1** | %0.0 |
| TCGA-99-8025 | yok | %54.6 | %0.9 | **%0.0** | %0.0 |

0.55'te kalemsiz slaytlar ≤ %0.2, kalemliler %5.6-63.1 — 28 kat ayrım.

Etki: "belirgin kalem" sayısı 370 WSI → **21 WSI**. Elenen aday
1.49 M → **47.057**.

## 2.9 · Patch seviyesi ink ve doku eşikleri kaldırıldı

Görsel denetim sonrası iki filtre de devre dışı:

| Filtre | Neden kaldırıldı |
|---|---|
| `ink_max` | İncelenen 15 patch'in 13'ü mürekkep değil, hematoksilen yoğun tümör dokusuydu |
| `tissue_min` | Metrik boyama yoğunluğunu ölçüyor, merkez-bağımlı seçilim yaratıyor |

Tek aktif patch-seviyesi filtre: `blur_min = 100`.

## 2.10 · Kohort 468 → 467

`TCGA-44-7661` çıkarma aşamasında **sıfır patch** üretti. Manuel kontrolde
slaydın tamamının odak dışı olduğu doğrulandı — QC eşiğini geçen patch yok.

Gerekçe `excluded.csv`'ye eklendi: `no_patch_passed_qc_out_of_focus`.

Bu aynı zamanda blur eşiğinin doğru seçildiğinin kanıtı: kullanılabilir
slaytlarda neredeyse hiç eleme yapmıyor, kullanılamaz slaytta hepsini eliyor.

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

## 3.8 · [Claude] Kalem eşiği 0.35 tahminle kondu, 370 slaytta yanlış alarm verdi

İlk tarama 468 slaydın **370'inde** "belirgin kalem" buldu; ortanca %7.4,
max %92. Gerçek olsa TCGA kullanılamaz olurdu.

Kullanıcı ImageScope'ta kontrol etti: metriğin üstte gösterdiği 5 slayttan
yalnızca 1'inde gerçek kalem vardı. Diğerleri yoğun/atipik boyanmış normal
dokuydu.

**Ders:** `[tahmin]` etiketli bir eşik, doğrulanmadan çalıştırılmamalı.
Düzeltme: 2.8'deki kalibrasyon.

## 3.9 · [Claude] "Mürekkep ihmal edilebilir" sonucu bozuk metrikle verildi

Merkez-katmanlı pilotta `ink_patch` p95 değeri 0.000-0.081 çıkmıştı ve
"mürekkep az" denmişti. Ancak o metrik hematoksileni mürekkep sanıyordu —
yani ölçüm aleti bozuktu.

Ayrıca pilot slayt başına 200 patch alıyordu (%2 örnekleme); kalem izleri
belirli bölgelerde yoğunlaştığı için rastgele örneklemede denk gelinmemişti.

**Ders:** Bozuk aletin düşük okuma vermesi "sorun yok" demek değil.
Düzeltme: tespit patch seviyesinden slayt seviyesine taşındı.

## 3.10 · [Claude] Deterministik örnekleme iddiası yanlıştı

`build_manifest()` içinde:

```python
rng = np.random.default_rng(abs(hash((cfg["seed"], wid))) % (2 ** 32))
```

Python'da string hash'i her süreçte farklı tuzlanır (PYTHONHASHSEED).
Yani manifest iki kez üretilirse farklı patch'ler seçilir.

Üretilen patch'lere zarar vermiyor ama "sabit seed ile deterministik"
iddiası tutmuyor.

**Düzeltme:** `hashlib.blake2b` ile sabit hash. LUSC'a geçmeden önce
uygulanmalı.

## 3.12 · [Claude] `he_residual` patch seviyesinde mürekkebi ölçmüyordu

Slayt seviyesinde çalışan kalıntı ölçüsü patch seviyesinde de kullanılmak
istendi. Ancak turkuaz (kalem altı) bir patch **0,1346** kalıntı verdi —
kohort ortancasından (0,275) düşük.

Sebep: H ve E iki vektör, renk uzayı üç boyutlu. Negatif katsayılara izin
verildiğinde her renk yaklaşık temsil edilebiliyor.

Ayrıca metrik merkeze göre değişiyordu (TSS 86 %74 elenecekti) — yani
boyama yoğunluğunu ölçüyordu.

**Düzeltme:** negatif konsantrasyon oranı (`neg_frac`). Gerçek H&E dokusunda
her iki boya konsantrasyonu pozitif olmalı; turkuaz patch'te eozin −0,596
çıktı. Ayrım 100% vs 0–2,8%.

## 3.13 · [Claude] `neg_frac` eşiği iki kez yanlış seçildi

**İlk deneme 0.05** — antrakoz (akciğerdeki karbon birikimi) ve yoğun
hematoksilenli gerçek dokuyu eliyordu. Siyah karbon partikülleri H&E
vektörlerine dik durduğu için hafif negatif katsayı üretiyor.

**İkinci deneme 0.50** — bu sefer gevşekti; 0,49 civarındaki patch'lerde
hâlâ mavi/turkuaz kalıntı görülüyordu.

**Nihai 0.15** — aralık aralık görsel denetimle bulundu: 0,10–0,15 ağırlıklı
temiz (antrakoz, yoğun çekirdek), 0,15 üstünde mavi/turkuaz başlıyor.

**Ders:** Antrakoz gerçek biyolojik yapıdır ve miktarı hastanın yaşadığı
yere/sigara geçmişine bağlıdır — yani merkeze bağlı. Elenmesi merkez-bağımlı
seçilim yaratır.

## 3.11 · Klasör adı `_excluded_non_DX1` yanıltıcıydı

İçindeki dosya bir **DX slaytıydı** — sadece 10. tanısal slayt. "non-DX"
değildi.

**Ders:** Klasör adları açıklayıcı olmalı, ama gerekçe kaynağı kod olmalı.

---

# BÖLÜM 4 — AÇIK KONULAR

- [x] ~~`max_patches`~~ → **2000** `[ölçüm]`. Ortanca slaytın %19'u.
      Geri alınabilir: tüm adaylar `candidates/` altında saklı
- [x] ~~`tissue_threshold`~~ → thumbnail'de **0.80**, patch seviyesinde
      filtre **yok** `[ölçüm]`
- [x] ~~`blur_threshold`~~ → **100** `[ölçüm + görsel]`. 884.987 patch'te
      eleme oranı düşük; tamamen odak dışı slaytta %100 eliyor
- [x] ~~İnk tespiti~~ → H&E deconvolution kalıntısı, eşik **0.55** `[ölçüm]`
- [ ] Deterministik hash düzeltmesi (3.10) — LUSC öncesi
- [ ] Bag boyutu dengesizliği: 20 WSI 2000'in altında (78-890 arası).
      Olduğu gibi bırak / dışla / MIL'de normalize et
- [ ] Nekroz nasıl ele alınacak (artefakt değil, gerçek biyoloji)
- [ ] FOV 112 µm yeterli mi — LUAD/LUSC ayrımı doku mimarisine dayanıyor
- [ ] Macenko ablation: `patches/` → `patches_macenko/`, ikinci set

---

# BÖLÜM 5 — LUAD PATCH ÇIKARMA SONUÇLARI

| | |
|---|---|
| İşlenen WSI | 468 |
| Toplam patch | **884.987** |
| Tam (2000 patch) | 447 WSI |
| Eksik (78-890) | 20 WSI |
| Sıfır patch | 1 WSI → kohorttan çıkarıldı |
| **Analiz kohortu** | **467 WSI / 467 hasta** |
| Disk | ~19 GB |
| Süre | ~9 saat (tek işlemci) |

## Uygulanan protokol

| Parametre | Değer | Kaynak |
|---|---|---|
| target_mpp | 0.50 µm/px | `[literatür]` UNI/Virchow eğitim ölçeği |
| patch_size | 224 px | `[literatür]` üç modelin ortak girdisi |
| FOV | 112 × 112 µm | türetilmiş |
| stride | 112 µm (overlap yok) | `[literatür]` CLAM |
| tissue_min (thumbnail) | 0.80 | `[ölçüm]` eşik düşük etkili (%17 aralık) |
| pen_max | 0.01 | `[ölçüm]` deconvolution kalıntısı > 0.55 |
| max_patches | 2000 | `[ölçüm]` süre ve disk ölçüldü |
| blur_min | 100 | `[ölçüm + görsel]` |
| ink / tissue (patch) | filtre yok | `[görsel]` merkez-bağımlı seçilim riski |

## Kalem izi dağılımı

| | WSI |
|---|---:|
| Hiç yok | 413 |
| Eser (<%1) | 34 |
| Belirgin (>%1) | 21 |

Belirgin olan 21 WSI'ın **tamamı TSS 86'dan**. Merkez ortalamaları:
TSS 86 %11.58, ikinci sıradaki TSS 93 %1.65 — yirmi kat fark.

> **Tez için bulgu:** Kalem kullanımı bir hastanenin rutini. Merkez sınıfı
> belirlediği için (LUAD/LUSC merkezleri hiç örtüşmüyor), kalem izi
> temizlenmezse model onu sinyal olarak kullanabilirdi. Literatürde de
> foundation modellerin kalem izlerine attention verdiği gösterilmiştir.

---

*Son güncelleme: 27 Eylül 2026 — LUAD patch çıkarma tamamlandı*

---

# BÖLÜM 6 — LUSC KOHORTU

| | |
|---|---|
| GDC tanısal slayt | 512 |
| İndirilen (01 + DX1) | **478 WSI / 478 hasta** |
| Atlanan | 34 (DX2–DX9) |
| Dışlanan | **0** |

LUAD'dan farklı olarak hiç dışlama gerekmedi: nüks tümör yok, MPP'si
okunamayan WSI yok, mükerrer hasta yok.

## MPP dağılımı — kohortlar arası fark  `[ölçüm]`

| | LUAD | LUSC |
|---|---:|---:|
| ~0,25 (40×) | 414 | 414 |
| ~0,23 (43×) | 35 | 5 |
| **0,50 (native 20×)** | **19 (%4)** | **59 (%12)** |

LUSC'ta native 20× oranı üç kat fazla. Bu slaytlar küçültülmüyor, 40×
olanlar `INTER_AREA` ile ~2× küçültülüyor (hafif yumuşatma).

Fark sınıflar arasında dengesiz dağılmış — teorik olarak sinyal olarak
kullanılabilir. `base_mpp` manifest'te saklı, değerlendirmede kontrol
edilecek.

## Okuma süresi

Ortanca 12 ms/patch (LUAD'da 15–33 ms). 478 × 2000 patch ≈ 3 saat.

---

*Son güncelleme: 28 Eylül 2026 — LUSC kohortu ve QC tamamlandı*

---

# BÖLÜM 7 — İKİ KOHORT TAMAMLANDI

| | WSI | Patch |
|---|---:|---:|
| **LUAD** | **466** | **900.828** |
| **LUSC** | **478** | **934.875** |
| Toplam | 944 | 1.835.703 |

## Dışlamalar

**LUAD** (480 → 466)

| n | Gerekçe |
|---:|---|
| 10 | Küçültülmüş türev görüntü (MPP yok, L0 = 10.000 px) |
| 1 | İkincil tanısal slayt (DX10) |
| 1 | Nüks tümör (02Z) |
| 1 | Odak dışı — QC'yi geçen patch yok |
| 1 | Trikrom boyama (H&E değil) |

**LUSC** (512 tanısal → 478, dışlama yok)

34 ikincil tanısal slayt (DX2–DX9) indirme aşamasında filtrelendi.

## Düşük bag'ler

| Slayt | Patch | Sebep |
|---|---:|---|
| TCGA-55-8207 | 78 | 2 mm² doku, küçük biyopsi |
| TCGA-67-4679 | 289 | Doku %17, blur 5.776 eledi |
| TCGA-98-8022 | 290 | Genel tarama kalitesi düşük (blur p50 = 869 vs kohort 1.089) |
| TCGA-55-8301 | 326 | Az doku |
| TCGA-56-8082 | 381 | Parçalı doku |

Bag boyutu kararı (olduğu gibi bırak / dışla / MIL'de normalize et) embedding
aşamasında verilecek.

---

*Son güncelleme: 30 Eylül 2026 — her iki kohortun patch çıkarması tamamlandı*
