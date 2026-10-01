# TEZ TEKNİK DURUM RAPORU

**TCGA-LUAD / TCGA-LUSC · Foundation Model Karşılaştırması**
DINOv2 · UNI · Virchow · ABMIL

Tarih: 30 Eylül 2026
Aşama: Her iki kohortun patch çıkarması tamamlandı · embedding öncesi

---

## 0 · YAŞANAN PROBLEMLER VE ÇÖZÜMLERİ

Bu bölüm, veri hazırlama sürecinde karşılaşılan somut problemleri, nasıl
fark edildiklerini ve nasıl çözüldüklerini kaydeder. Amaç şeffaflık:
"neden böyle yapıldı?" sorusunun cevabı hatırlamaya değil kayda dayanır.

### 0.1 · Kohort tanımındaki iki sessiz hata

| Problem | Nasıl fark edildi | Çözüm |
|---|---|---|
| `"DX1" in dosya_adı` filtresi `DX10`'u da yakalıyordu | Dosya listesi barkod alanlarına ayrıştırılınca | Slayt numarası tamsayı olarak ayrıştırılıyor, `slide_no == 1` kontrolü |
| Sample type hiç kontrol edilmiyordu; `TCGA-50-5066-02Z` (nüks tümör) kohorta sızmıştı | Barkodun 4. alanı (`02Z`) incelenince | `sample_type == 01` kuralı eklendi |

Her iki hata da aynı hastanın iki kaydının sete girmesine yol açıyordu —
yani hasta seviyesinde mükerrerlik ve potansiyel data leakage.

> Not: `02Z` dosyası başlangıçta "aynı hastanın ikinci DX1 slaytı, bu
> normaldir" şeklinde yorumlanmıştı. Görsel olarak iki slayt gerçekten
> farklı doku içeriyordu; ancak fark görüntüde değil **örnek tipindeydi**
> (primary vs recurrent).

**Alınan önlem:** Barkod mantığı ayrı bir modüle (`tez/barcode.py`) taşındı
ve kendi birim testini içeriyor (`self_test()`). Kohort kurulmadan önce
otomatik çalışıyor; test başarısızsa pipeline duruyor.

### 0.2 · Dışlama gerekçesi klasör adından alınıyordu

Dışlanan dosyalar `_excluded_*` klasörlerine elle taşınmıştı ve gerekçe
klasör adından türetiliyordu. Ancak `_excluded_duplicate_patient` klasörü
hem `DX10` hem `02Z` içeriyordu — ikisi farklı sebeple dışlanmıştı ama
aynı etiketi alıyorlardı.

**Çözüm:** Barkod kesin bir gerekçe verebiliyorsa o kullanılıyor; klasör
adı yalnızca barkodun söyleyemediği durumlar için (MPP eksikliği gibi).
Klasör adları elle verilir ve yanılabilir; barkod yanılmaz.

Benzer şekilde `_excluded_non_DX1` klasör adı da yanıltıcıydı: içindeki
dosya bir DX slaytıydı, sadece 10. tanısal slayttı.

### 0.3 · Doku eşiği karşılaştırması matematiksel olarak geçersizdi

Önceki bir denemede 0.70 ve 0.80 doku eşikleri karşılaştırılmış ve
"aralarında yalnızca 30 patch (%0,064) fark var" sonucuna varılmıştı.

```
469 WSI × 100 max_patches = 46.900   ← tavan
Ölçülen:  0.70 → 46.897   ·   0.80 → 46.867
```

Neredeyse her slayt tavana dayanmıştı. Ölçülen şey eşiğin etkisi değil,
**tavanın etkisiydi.** Eşik ne seçilirse seçilsin sonuç değişmezdi.

**Çözüm:** Aday tarama aşaması hiçbir eşik veya tavan uygulamadan yeniden
çalıştırıldı (Bölüm 4). Gerçek fark şimdi ölçülebiliyor: 0.70 → 0.80
geçişinde slayt başına ortanca ~700 patch kayboluyor, %0,064 değil.

### 0.4 · 10 WSI'ın dışlanma gerekçesi eksikti

Başlangıçta gerekçe "metadata eksikliği" olarak kaydedilmişti. Vendor
metadata'sı (`aperio.MPP`, `aperio.AppMag`) kontrol edildiğinde bunların
da boş olduğu doğrulandı — ancak asıl bulgu level-0 boyutlarındaydı:
onunun da uzun kenarı **tam 10.000 piksel**.

MPP eksikliği semptomdu; sebep, dosyaların yeniden boyutlandırılmış türev
görüntüler olmasıydı (Bölüm 2.4).

### 0.5 · Uniform renkli bölgeler doku maskesini geçiyordu

Görsel denetimde, tamamen düz turkuaz görünen patch'lerin doku oranı
eşiğini (≥ 0.80) geçtiği görüldü. Bu bölgeler biyolojik doku içermiyor;
tarayıcı artefaktı veya renk kayması.

Renk tabanlı doygunluk maskesi bu bölgeleri doku sanıyor. Laplacian
varyansı ≈ 0 olduğu için blur QC bunları eleyecektir; ancak o denetim
sırasında blur filtresi devrede değildi.

**Durum:** Çıkarma aşamasında blur QC zorunlu hale getirildi.

### 0.6 · Mürekkep tespiti gerçek dokuyu siliyordu

"Mavi oranı ≥ 0.80" kuralıyla ayrılan patch'ler görsel olarak
incelendiğinde üç farklı içerik bulundu:

- gerçek mavi mürekkep / kalem izi
- **çekirdek yoğun gerçek tümör dokusu** (yoğun hematoksilen boyanması)
- turkuaz arka plan artefaktı

H&E'de çekirdekleri boyayan hematoksilen mor-mavidir ve RGB'de mavi kanalı
baskındır. Renk kanalı oranına dayalı bir kural, en bilgilendirici
patch'leri mürekkep sanıp siler.

Ölçülen HSV ton aralıkları (OpenCV 0–179 ölçeği):

| İçerik | Ton | Doygunluk | Laplacian varyansı |
|---|---|---|---|
| Turkuaz artefakt | ~80–95 | orta | ≈ 0 |
| Mavi mürekkep | ~105–125 | çok yüksek | düşük |
| **Hematoksilen** | **~125–150** | orta | **yüksek** |

**Çözüm yönü:** Ayırt edici özellik renk değil **doku/yapıdır**. Mürekkep
ve artefakt düzdür; gerçek doku kenar ve yapı içerir. Daha sağlam bir ölçü
olarak H&E stain deconvolution kalıntısı (patch'in rengi H&E boya
vektörleriyle açıklanabiliyor mu) değerlendirilmektedir.

Bu endişe literatürde de doğrulanmıştır (Bölüm 5, Soru 6 notu).

**Uygulanan çözüm ve kalibrasyon.** Ton aralığını daraltmak yetmedi; boyama
yoğunluğu merkeze göre kaydığı için sabit bir ton sınırı her merkezde aynı
yere düşmüyor. Renk tonu yerine **H&E stain deconvolution kalıntısı**
kullanıldı: *"Bu rengi hematoksilen ve eozin karışımıyla üretebilir miyim?"*
Çekirdek ne kadar koyu olursa olsun kalıntı düşük kalır (rengi zaten
hematoksilenden gelir); kalem H&E dışı bir pigment olduğu için yüksek
kalıntı verir.

Tespit ayrıca patch seviyesinden **slayt seviyesine** taşındı — kalem izleri
belirli bölgelerde yoğunlaştığı için rastgele patch örneklemesinde
görülmeyebiliyordu.

Eşik, üç kalemli ve dört kalemsiz slaydın karşılaştırılmasıyla belirlendi:

| Slayt | Kalem | >0.35 | >0.45 | **>0.55** | >0.65 |
|---|---|---:|---:|---:|---:|
| TCGA-86-7953 | var | %90,8 | %80,4 | **%63,1** | %28,9 |
| TCGA-86-8056 | var | %55,0 | %44,0 | **%16,5** | %3,1 |
| TCGA-86-7714 | var | %26,5 | %13,6 | **%5,6** | %1,0 |
| TCGA-75-6214 | yok | %57,0 | %3,7 | **%0,1** | %0,0 |
| TCGA-50-6673 | yok | %78,9 | %5,3 | **%0,2** | %0,0 |
| TCGA-50-5045 | yok | %61,7 | %2,0 | **%0,1** | %0,0 |
| TCGA-99-8025 | yok | %54,6 | %0,9 | **%0,0** | %0,0 |

İlk denenen 0,35 eşiğinde iki grup üst üste biniyordu ve 468 WSI'ın 370'inde
yanlış alarm veriyordu. 0,55'te ayrım 28 kata çıkmaktadır.

| | eşik 0,35 | **eşik 0,55** |
|---|---:|---:|
| "Belirgin kalem" WSI | 370 | **21** |
| Elenen aday patch | 1.493.976 | **47.057** |

Belirgin kalem izi bulunan 21 WSI'ın **tamamı tek bir kaynak merkezden**
(TSS 86) gelmektedir. Merkez ortalamaları: TSS 86 %11,58; ikinci sıradaki
TSS 93 %1,65 — yirmi kat fark. Kalem kullanımı bir hastanenin rutinidir ve
literatürdeki bulgunun doğrudan tekrarıdır (bkz. Bölüm 5, Soru 6).

Kalem izi tespit edilen koordinatlar manifest aşamasında aday havuzundan
çıkarılmıştır; slayt dışlanmamıştır.

### 0.7 · Blur QC'nin yanlış sırada uygulanması

Laplacian varyansı çözünürlüğe bağlı bir ölçüdür. Farklı base MPP'ye sahip
slaytlar farklı piksel boyutlarında okunmaktadır:

```
40× WSI        →  444 px okunur  →  224'e küçültülür
native 20× WSI →  224 px okunur  →  resize yok
```

Blur eşiği resize öncesinde uygulanırsa iki grup farklı ölçekte
değerlendirilir ve 19 native-20× slaytın patch'leri sistematik olarak
haksız elenebilir.

**Kural:** Blur QC her zaman 224×224'e küçültme **sonrasında** uygulanır.

### 0.8 · Ölçülmeden verilen kararlar

İki parametre gerekçelendirilirken ölçüm yapılmamıştı:

| Parametre | Verilen gerekçe | Ölçüm sonrası |
|---|---|---|
| `max_patches = 100` | "hız ve disk kullanımı" | 15 ms/patch → 500 patch için ~1 sa, 7 GB. Gerekçe geçersiz. |
| Texture eşiği | (tahmin) | Aday sayısını yalnızca %0,6 değiştiriyor — kalibre edilmeden kullanılmıyor |

**Alınan önlem:** Bundan sonra her eşik ve parametre kaynağıyla
kaydedilmektedir: `[ölçüm]` (bu veride hesaplandı), `[literatür]`
(kaynak var), `[tahmin]` (doğrulanması gerekiyor).

### 0.9 · Pilot örneklem temsili değildi

İlk pilot `--limit 10` ile çalıştırıldı; bu parametre alfabetik sırada ilk
10 WSI'ı alır ve tamamı tek bir kaynak merkezdendi (TSS 05). Ölçülen
mürekkep, blur ve boyama dağılımları o merkezin tarayıcı ve boyama
protokolüne aitti.

**Çözüm:** Merkez-katmanlı pilot seçeneği eklendi (her merkezden birer WSI).
33 merkezden alınan örneklemde `tissue_patch` ortancası **0,58 – 0,96**
arasında değişti — yani ilk pilotla eşik seçmek yanlış olurdu.

### 0.10 · Patch seviyesindeki doku ve mürekkep eşikleri kaldırıldı

Görsel denetim sonrası iki filtre de devre dışı bırakıldı:

| Filtre | Gerekçe |
|---|---|
| Mürekkep oranı | En yüksek skorlu 15 patch'in 13'ü mürekkep değil, çekirdek yoğun gerçek tümör dokusuydu |
| Doku oranı | Metrik boyama yoğunluğunu ölçüyor; merkeze göre 0,58–0,96 değişiyor, filtreleme merkez-bağımlı seçilim yaratır |

Tek aktif patch seviyesi filtre **blur ≥ 100**'dür. Her iki metrik yine de
kaydedilmekte, merkezler arası boyama farkının nicel ölçüsü olarak ablation
analizlerinde kullanılabilir durumdadır.

### 0.11 · Örnekleme deterministik değildi

Manifest üretiminde seed, `hash()` fonksiyonuyla üretiliyordu. Python 3'te
string hash'i her süreçte farklı tuzlanır (PYTHONHASHSEED); dolayısıyla
manifest iki kez üretilirse farklı patch'ler seçilirdi.

Üretilen patch'lere zarar vermemektedir; ancak "sabit seed ile deterministik
örnekleme" iddiası geçersizdi. `hashlib.blake2b` ile düzeltilmiştir; LUSC
kohortu bu düzeltmeyle işlenecektir.

---

## 1 · ÖZET

| Aşama | Durum | Çıktı |
|---|---|---|
| Kohort tanımı | ✅ Tamamlandı | **467 WSI / 467 hasta** |
| WSI teknik QC | ✅ Tamamlandı | Tümü okunabilir, MPP ölçüldü |
| Fiziksel ölçek standardizasyonu | ✅ Tamamlandı | 0.50 µm/px, 112×112 µm FOV |
| Patch aday taraması | ✅ Tamamlandı | ~5,3 M aday koordinat |
| Kalem izi tespiti | ✅ Tamamlandı | 21 WSI'da belirgin, tümü TSS 86 |
| **Patch çıkarma (LUAD)** | ✅ **Tamamlandı** | **900.828 patch** |
| LUSC indirme | ✅ Tamamlandı | **478 WSI / 478 hasta · 376 GB** |
| LUSC kohort + QC | ✅ Tamamlandı | 13 kontrol geçti, dışlama yok |
| **Merkez örtüşmesi ölçümü** | ✅ **Tamamlandı** | **Ortak merkez: 0** |
| LUSC patch çıkarma | ✅ Tamamlandı | **934.875 patch** |
| Kalem maskesi düzeltmesi | ✅ Tamamlandı | 18 slayt yeniden çıkarıldı |
| Embedding çıkarma | ⬜ Başlanmadı | DINOv2 / UNI / Virchow |
| MIL eğitimi | ⬜ Başlanmadı | ABMIL |

Tüm kod ve çıktılar versiyonlanmıştır:
`github.com/esmnralican/tez-wsi-foundation` (private)

**Danışman görüşü beklenen kararlar Bölüm 5'tedir.**

---

## 2 · KOHORT TANIMI

### 2.1 · Ne yapıldı

TCGA-LUAD projesinden tanısal (diagnostic, FFPE) WSI'lar indirildi. Kohort,
üç kurala göre tanımlandı:

1. `sample_type = 01` → yalnızca **primary solid tumor**
2. `slide = DX1` → yalnızca **birinci tanısal slayt**
3. Hasta başına tek WSI

### 2.2 · Sonuç

```
480   GDC manifest'teki toplam DX kaydı
 −1   ikincil tanısal slayt (DX10)
 −1   nüks tümör örneği (sample_type = 02)
−10   güvenilir fiziksel çözünürlük bilgisi olmayan WSI
 −1   QC'yi geçen patch üretilemeyen WSI (tamamen odak dışı)
────
467   analiz kohortu  ·  467 WSI = 467 hasta (birebir)
```

### 2.3 · Dışlama gerekçeleri

| n | Gerekçe | Açıklama |
|---:|---|---|
| 10 | Güvenilir MPP/büyütme metadata'sı yok | Aşağıda ayrıntı |
| 1 | İkincil tanısal slayt (`TCGA-49-4488-...-DX10`) | Aynı hastanın 10. slaytı |
| 1 | Nüks tümör (`TCGA-50-5066-02Z`) | Primary değil, recurrent |
| 1 | Odak dışı (`TCGA-44-7661`) | QC eşiğini geçen patch üretilemedi |

Son satır çıkarma aşamasında ortaya çıkmıştır: slayt teknik olarak
okunabilir ancak tamamı odak dışıdır; 2000 aday koordinatın hiçbiri blur
eşiğini (Laplacian varyansı ≥ 100) geçememiştir. Manuel görsel kontrolle
doğrulanmıştır. Bu aynı zamanda blur eşiğinin doğru kalibre edildiğini
gösterir: kullanılabilir slaytlarda eleme oranı ihmal edilebilir düzeydeyken
bu slaytta %100'dür.

### 2.4 · 10 WSI'ın dışlanma gerekçesi (ayrıntı)

Bu WSI'lar teknik olarak **okunabilir** durumdadır; bozuk değillerdir. Ancak:

- `openslide.mpp-x`, `aperio.MPP`, `aperio.AppMag`, `objective-power`
  alanlarının **tamamı boş**
- Level-0 boyutlarının uzun kenarı **onunda da tam 10.000 piksel**

Karşılaştırma (aynı kaynak merkezden, TSS 05):

| | Level-0 boyutu | MPP |
|---|---|---|
| Sağlam slaytlar | 37.120 – 78.080 px | 0.2325 |
| Dışlanan 10 slayt | **10.000 px (sabit)** | yok |

Bu, orijinal tarama değil **yeniden boyutlandırılmış türev görüntü** olduğunu
gösterir. Tahmini fiziksel çözünürlük ~2–2.5 µm/px (≈4–6× büyütme), hedef
0.50 µm/px seviyesinin çok altında. Hedefe yükseltmek ancak interpolasyonla
mümkün olurdu; bu da hücresel morfoloji içermeyen yapay patch'ler üretirdi.

> **Metodolojik not:** Dışlanan 10 WSI'ın tamamı **tek bir kaynak merkezden**
> (TSS 05) gelmektedir. Bu sistematik bir kayıptır ve tezde açıkça
> raporlanacaktır. Aynı merkezden 24 WSI kohortta kalmıştır.

### 2.5 · Literatürdeki standart kohort ile karşılaştırma

TCGA-NSCLC alt tip sınıflandırması için yaygın olarak kullanılan kurulum,
478 vakadan 512 LUSC ve 478 vakadan 541 LUAD olmak üzere toplam **1.053
tanısal WSI** içermektedir.

Görüldüğü üzere literatür **hasta başına birden fazla tanısal slayt**
kullanmaktadır (541 LUAD slaytı / 478 hasta).

Bu çalışmada hasta başına tek WSI (DX1) tercih edilmiştir. Bu tercih:

- **Avantaj:** Hasta seviyesinde mükerrerlik ve data leakage riski sıfıra
  iner; bag'ler arası bağımsızlık garanti altına alınır.
- **Dezavantaj:** Örneklem büyüklüğü literatürden düşük olur (468 vs 541).
  Ayrıca aynı hastanın diğer slaytlarındaki doku bilgisi kullanılmaz.

**Alternatif:** Aynı hastanın tüm DX slaytlarının patch'leri **tek bir
bag'de** birleştirilebilir. Bu, hasta seviyesinde etiket kullanıldığı için
leakage oluşturmaz ve daha fazla doku görülmesini sağlar.

**Görüş isteniyor:** Hasta başına tek slayt mı, yoksa tüm DX slaytları tek
bag'de birleştirilerek mi ilerlenmeli?

### 2.6 · Tez metni için önerilen ifade

> Initial DX1 cohort consisted of 479 WSIs. Ten WSIs were excluded from the
> primary analysis: although technically readable, these files contained no
> resolution metadata and their level-0 dimensions were capped at exactly
> 10,000 pixels on the long axis, indicating downsampled derivative images
> rather than original scans (comparable slides from the same source site
> ranged from 37,120 to 78,080 pixels). One additional WSI was excluded as a
> recurrent-tumour sample, and one further WSI was removed after patch
> extraction because no patch passed the focus quality threshold. The final
> cohort comprised **467 primary-tumour WSIs from 467 unique patients**.

---

## 3 · FİZİKSEL ÖLÇEK STANDARDİZASYONU

### 3.1 · Problem

Aynı piksel boyutundaki iki patch, aynı fiziksel doku alanını
göstermeyebilir:

```
224 px @ 0.25 µm/px  →   56 × 56 µm
224 px @ 0.50 µm/px  →  112 × 112 µm
```

Üç foundation model de 224×224 girdi alır. Ancak "hepsi 224 alıyor" demek,
"hepsi aynı biyolojik alanı görüyor" anlamına gelmez.

Ayrıca metadata'daki objective power (20×/40×) tek başına güvenilir değildir:
metadata'da 20× yazan 43 WSI'ın yalnızca 19'unun gerçek çözünürlüğü
0.50 µm/px'tir; kalan 24'ü 0.2325 civarındadır.

### 3.2 · Benimsenen yaklaşım

Level seçimi objective power'a göre değil, **ölçülmüş MPP'ye** göre yapılır:

```
base_mpp × level_downsample = o level'in gerçek MPP'si
                            → hedefe en yakın uygun level
```

Level seçimi **her zaman aşağı yönlüdür**: hedeften daha yüksek çözünürlüklü
bir level seçilir ve yazılımda küçültülür. Yukarı ölçekleme yapılmaz —
interpolasyon bilgi kazandırmaz, sahte detay üretir.

### 3.3 · Ölçüm sonuçları (468 WSI)

**Base MPP dağılımı:**

| Ölçülen MPP | WSI | Yaklaşık büyütme |
|---|---:|---|
| 0.25 | 414 | ~40× |
| 0.23 | 35 | ~43× |
| 0.50 | 19 | 20× |

**Pyramid yapısı:** 368 WSI'da 4 level, 99'unda 3 level, 1'inde 2 level.
Bu nedenle sabit bir level numarası (`level = 2` gibi) kullanılamaz.

**Seçilen level:** 468 WSI'ın tamamında **level 0**. Aperio pyramid'i
tipik olarak `[1, 4, 16, 32]` oranlarındadır; level 1 zaten 1.0 µm/px'e
düşmektedir (hedefin iki katı kaba). Dolayısıyla 0.50 µm/px arada kalmakta
ve level 0'dan okunup küçültülmesi gerekmektedir.

**Okuma boyutları:**

| Base MPP | Okunan | İşlem | WSI |
|---|---|---|---:|
| 0.252 | 444 px | → 224 px (1.98× küçültme) | 308 |
| 0.2527 | 443 px | → 224 px | 37 |
| 0.2325 | 482 px | → 224 px (2.15× küçültme) | 24 |
| 0.2457 | 456 px | → 224 px | 24 |
| 0.500 | 224 px | resize yok | 19 |
| 0.2485 | 451 px | → 224 px | 14 |

Farklı piksel sayısı okunmakta, **aynı 112 × 112 µm fiziksel alan** elde
edilmektedir.

### 3.4 · Fairness açısından önemi

```
                  Canonical patch (112 × 112 µm)
                 /              |              \
            DINOv2             UNI          Virchow
```

Üç model aynı WSI'dan, aynı ön işlemeden geçmiş, aynı koordinatlardaki,
aynı fiziksel alanı temsil eden patch'leri görecektir. Karşılaştırmada
değişen tek değişken **foundation model** olacaktır.

---

## 4 · PATCH ADAY TARAMASI

### 4.1 · Ne yapıldı

468 WSI üzerinde, **hiçbir eşik veya üst sınır uygulanmadan** tüm aday
patch koordinatları ve doku oranları ölçüldü.

> **Neden tavansız:** Önceki bir denemede doku eşiği karşılaştırması
> (0.70 vs 0.80) yapılmış ve "aralarında yalnızca %0.064 fark var" sonucuna
> varılmıştı. Bu sonuç geçersizdir: o denemede slayt başına 100 patch üst
> sınırı vardı ve 469 × 100 = 46.900 tavanına neredeyse her slayt dayanmıştı.
> Ölçülen şey eşiğin etkisi değil, tavanın etkisiydi. Bu nedenle bu kez
> önce ölçüm yapılmış, karar sonraya bırakılmıştır.

### 4.2 · Doku alanı dağılımı (468 WSI)

| | mm² |
|---|---:|
| min | 2 |
| p25 | 59 |
| **ortanca** | **138** |
| p75 | 232 |
| max | 668 |

### 4.3 · Eşiğe göre aday patch sayısı

| Doku eşiği | ortanca/slayt | toplam |
|---|---:|---:|
| ≥ 0.30 | 11.842 | 6,07 M |
| ≥ 0.50 | 11.255 | 5,81 M |
| ≥ 0.70 | 10.414 | 5,50 M |
| **≥ 0.80** | **9.692** | **5,31 M** |
| ≥ 0.90 | 9.130 | 5,07 M |

**Bulgu:** 0.30 ile 0.90 arasındaki tüm aralıkta toplam aday sayısı yalnızca
%17 değişmektedir. Doku eşiği düşük etkili bir parametredir.

### 4.4 · Doku kapsamı

Bir patch 112 × 112 µm = 0.0125 mm² alanı temsil eder. Ortanca slayt
(138 mm² doku) için:

| Patch/WSI | Kapsanan doku |
|---:|---|
| 100 | %0.9 |
| 300 | %2.7 |
| 500 | %4.5 |
| 1000 | %9.1 |

### 4.5 · Yetersiz aday içeren WSI'lar

Doku eşiği 0.80'de:

| Hedef | Ulaşamayan WSI |
|---|---:|
| < 100 aday | 1 |
| < 500 aday | 4 |
| < 1000 aday | 8 |

Minimum doku alanı 2 mm² (küçük biyopsi materyali).

### 4.6 · Zaman ve depolama ölçümleri

- Aday tarama: ortanca 0.8 s/WSI, toplam ~9 dakika
- Patch okuma: ortanca **15 ms/patch** (tek işlemci)

| Patch/WSI | Çıkarma süresi | Disk |
|---:|---|---|
| 100 | ~0.2 sa | ~1,4 GB |
| 500 | ~1,0 sa | ~7 GB |
| 1000 | ~1,9 sa | ~14 GB |

> Bu ölçüm iyimserdir; gerçek çıkarma dağınık koordinatlardan okuyacak ve
> JPEG kodlama eklenecektir. Gerçekte 3–5 katı beklenmektedir.

---

## 5 · DANIŞMAN GÖRÜŞÜ BEKLENEN KARARLAR

### Soru 1 — Slayt başına patch sayısı  ⚠️ KARAR VERİLDİ, ONAY BEKLİYOR

Önceki protokolde **100 patch/WSI** belirlenmişti; gerekçe "hız ve disk
kullanımı" idi. Ölçüm sonrası bu gerekçe geçerliliğini yitirdi
(15–33 ms/patch).

#### Literatür ne yapıyor

**Standart yaklaşım tavan koymamaktır — tüm doku patch'leri kullanılır.**

CLAM (Lu et al., *Nat Biomed Eng* 2021) protokolü: 20× seviyesinde
256×256 piksel, step 256 (overlap yok). Slayt başına üst sınır yoktur; doku
konturu içindeki her koordinat alınır. Foundation model karşılaştırmalarında
da aynı yaklaşım izlenir.

#### Donanım kısıtı

| Yaklaşım | Patch/WSI | LUAD toplam | 3 model × LUAD+LUSC |
|---|---:|---:|---:|
| Literatür (tümü) | ~9.486 | 4,4 M | **~27 M forward pass** |
| **2000 (seçilen)** | **2.000** | **885 K** | **5,3 M** |
| 1000 | 1.000 | 467 K | 2,8 M |
| 100 (önceki) | 100 | 47 K | 281 K |

Mevcut donanım (RTX 3050 Ti, 4 GB VRAM) ile ViT-L/14 çıkarımı kabaca
40–60 patch/sn hızındadır. Literatür standardı ~130–190 saat sürer.

#### Alınan karar: 2000 patch/WSI

Gerekçe:

1. **Geri alınamazlık.** `.svs` dosyaları disk kısıtı nedeniyle
   kaldırılacaktır (LUAD + LUSC ≈ 706 GB, kullanılabilir alan ~530 GB).
   Patch sayısını sonradan **azaltmak** mümkün, **artırmak** değildir.
2. **Maliyet düşük.** 2000 patch = ~19 GB disk, ~9 saat çıkarma.
3. **Kapsam.** Ortanca slaytın **%19'u** (100 patch %0,9 idi).

Embedding aşamasında alt küme kullanılabilir; manifest'ten seçim tek
satırlık işlemdir.

**Gerçekleşen:** 884.987 patch / 468 WSI. 447 WSI tam 2000'e ulaştı,
20 WSI 78–890 arasında kaldı (küçük biyopsi materyali), 1 WSI sıfır patch
üretti ve kohorttan çıkarıldı.

**Görüş isteniyor:** Donanım kısıtı nedeniyle örnekleme zorunludur ve tezde
açıkça belirtilecektir. 2000 patch/WSI (ortanca slaytın %19'u) metodolojik
olarak savunulabilir midir?

> Tüm aday koordinatlar saklandığı için bu karar **geri alınabilir**;
> sayıyı azaltmak yeniden tarama gerektirmez.

---

### Soru 2 — Doku oranı eşiği  ✅ KARAR VERİLDİ

Ölçüm, eşiğin düşük etkili olduğunu gösterdi (0,30–0,90 aralığında toplam
aday sayısı yalnızca %17 değişiyor).

**Karar: thumbnail seviyesinde 0,80, patch seviyesinde eşik yok.**

Patch seviyesinde ikinci bir doku eşiği uygulanmamasının gerekçesi Bölüm
0.10'dadır: metrik doku varlığını değil boyama yoğunluğunu ölçmekte ve
merkeze göre 0,58–0,96 arasında değişmektedir. Bu metrikle filtreleme
merkez-bağımlı bir seçilim yaratır ve Soru 5'teki batch effect riskini
büyütür.

#### Literatür ne yapıyor

CLAM'da patch seçimi **oran tabanlı değil, kontur tabanlıdır**: görüntü HSV
uzayına çevrilir, kernel 7 ile median blur uygulanır, sabit doygunluk eşiği
(20) ile doku bölgeleri belirlenir; ardından patch merkezinin doku konturu
içinde olup olmadığına bakılır. Patch içindeki doku oranı hesaplanmaz.

Bazı çalışmalar ek bir texture filtresi kullanır: Otsu eşikleme ve
morfolojik işlemlerle arka plan ve küçük doku parçaları (alan < 20)
çıkarılır, 40× slaytlar 20×'e indirilir, son olarak gradyan tabanlı bir
entropi eşiği (= 4) ile düşük dokulu patch'ler elenir.

**Bu çalışmadaki fark:** Oran tabanlı bir kriter (≥ 0.80) kullanılmaktadır.
Ölçüm, eşiğin düşük etkili olduğunu gösterdiği için bu tercih maliyetsizdir;
ancak yöntem farkı tezde belirtilmelidir. Gradyan entropisi tabanlı texture
filtresi, mevcut kalibre edilmemiş yerel-std yaklaşımının yerine
değerlendirilebilir.

**Not — metodolojik sınırlılık:** Doku oranı şu anda 2048 piksellik
thumbnail maskesi üzerinde hesaplanmaktadır. 68.000 piksellik bir slaytta
ölçek ~33× olduğundan, 444 piksellik bir patch thumbnail'de yalnızca
~13 × 13 piksele karşılık gelmektedir. Bu çözünürlükte %80 doku şart koşmak
gürültülüdür ve doku kenarlarını sistematik olarak eleyebilir. Tümör–stroma
sınırı LUAD/LUSC ayrımı açısından bilgilendirici olabileceğinden, çıkarma
aşamasında patch'in kendi pikselleri üzerinde ikinci bir ölçüm yapılması
planlanmaktadır.

---

### Soru 3 — Fiziksel görüş alanı yeterli mi?

Mevcut protokol: **224 px @ 0.50 µm/px = 112 × 112 µm**
(yaklaşık 10–15 hücre çapı)

LUAD/LUSC ayrımı büyük ölçüde **doku mimarisine** dayanır. Glandüler
organizasyonu değerlendirmek için 112 µm dar kalabilir.

Alternatifler:

| Seçenek | FOV | Model girdisi |
|---|---|---|
| Mevcut | 112 µm | 224 px @ 0.50 |
| Geniş bağlam | 224 µm | 224 px @ 1.00 |
| Çok ölçekli | 112 + 224 µm | iki ayrı embedding |

**Görüş isteniyor:** İlk denemede 112 µm ile devam edip, sonuç zayıf çıkarsa
mı değiştirilmeli; yoksa baştan çok ölçekli mi kurulmalı?

---

### Soru 4 — Küçük materyalli WSI'lar  ⚠️ AÇIK

Gerçekleşen çıkarmada **20 WSI** 2000 patch hedefine ulaşamamıştır:

| Patch sayısı | WSI |
|---|---:|
| 78 – 500 | 6 |
| 501 – 1000 | 8 |
| 1001 – 2000 | 6 |

En küçüğü 78 patch (`TCGA-55-8207`, 2 mm² doku — küçük biyopsi materyali).

MIL değişken bag boyutuyla çalışabilir; ancak 78 instance'lık bir bag ile
2000'lik bir bag'in aynı ağırlıkta değerlendirilmesi sorunlu olabilir.

Seçenekler:

1. Olduğu gibi bırak, bag boyutunu manifest'te kaydet (mevcut durum)
2. Minimum patch sayısının altındaki WSI'ları dışla
3. MIL'de bag boyutuna göre normalizasyon uygula

**Görüş isteniyor:** Hangi yaklaşım tercih edilmeli? Dışlama seçilirse
minimum patch sayısı ne olmalı?

---

### Soru 5 — Kaynak merkez (TSS) kaynaklı batch effect

**Bu, çalışmanın en ciddi metodolojik riski olarak görülmektedir.**

Bu bulgu artık **her iki kohortun tamamı indirilip QC'den geçtikten sonra
ölçülmüştür** (`run.py sites`), manifest üzerinden tahmin değildir:

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

Hiçbir kaynak merkez hem LUAD hem LUSC slaytı göndermemiştir. Dolayısıyla
**kaynak merkez bilgisi sınıf etiketini tam olarak belirlemektedir.**

Merkezler arasında H&E boyama protokolü, tarayıcı modeli, fiksasyon süresi
ve kesit kalınlığı sistematik olarak farklılaşır. Bu koşulda bir model,
morfolojiyi hiç öğrenmeden yalnızca renk dağılımı ve tarayıcı imzası
üzerinden yüksek doğruluk elde edebilir.

Bu durum TCGA LUAD/LUSC ile çalışan tüm literatür için geçerlidir; ancak
çalışmaların çoğu bunu ölçmemektedir.

#### Literatür ne diyor — bu problemin adı var

**Howard FM et al. (2021). "The impact of site-specific digital histology
signatures on deep learning model accuracy and bias." *Nature
Communications* 12:4423. DOI: 10.1038/s41467-021-24698-1**

Bu çalışma TCGA'daki 6 kanser tipinde 3.000'den fazla hasta üzerinde,
tahmin edilen özelliklerin doku gönderen merkezler arasında ciddi biçimde
farklılaştığını göstermekte; ayrıca merkezler arası görüntü farklarının
derin öğrenme ile kolayca tespit edilebildiğini ortaya koymaktadır.

**Kritik bulgu:** Yaygın olarak kullanılan renk normalizasyonu ve
augmentation yöntemlerine rağmen merkez tespiti mümkün olmaya devam
etmektedir; yazarlar bu merkez imzasını oluşturan görüntü özelliklerini
nicelleştirmiş ve bunların sağkalım, genomik mutasyon ve tümör evresi
tahmininde yanlı doğruluğa yol açtığını göstermişlerdir.

Önerdikleri çözüm, modellerin aynı merkezden gelen örneklerle hem eğitilip
hem doğrulanmamasını sağlayan bir kuadratik programlama yöntemidir.

#### Bu çalışma için anlamı

1. **Stain normalization bu sorunu çözmüyor.** Soru 6'daki karar bu bulgu
   ışığında değerlendirilmelidir: normalization uygulansa bile merkez
   imzası kalmaktadır.
2. **Önerilen leave-two-sites-out yaklaşımı literatürle uyumludur** —
   Howard et al.'ın split kısıtı fikrinin, merkezlerin tek sınıflı olduğu
   duruma uyarlanmış halidir.
3. **Site prediction baseline doğrudan bu makalenin yöntemidir.** Aynı
   embedding'lerden TSS tahmin etmek, "merkez imzası embedding'lerde var
   mı?" sorusunun standart testidir.
4. Bu bir kişisel endişe değil, **belgelenmiş ve yaygın atıf alan bir
   metodolojik problemdir**. Tezde bu referansla ele alınması, sonuçların
   savunulabilirliğini artırır.

**Önerilen ek analizler** (çalışmanın özgün katkısı olabilir):

1. **Site prediction baseline** — aynı embedding'lerden TSS kodunu tahmin
   et. Yüksek doğruluk, embedding'lerin merkez bilgisi taşıdığını gösterir.
2. **Leave-two-sites-out** — her sınıftan birer merkez dışarıda bırakılarak
   dengeli değerlendirme. (Merkezler tek sınıflı olduğu için klasik LOSO
   uygulanamaz.)
3. **Stain normalization ablation** — Macenko/Vahadane ile ve olmadan
   karşılaştırma. Doğruluktaki düşüş, renk bağımlılığının doğrudan ölçüsüdür.
4. **TDA Mapper renklendirmesi** — Mapper grafiği bir kez sınıfa, bir kez
   TSS'e göre renklendirilir. Kümelenmenin morfolojiye mi merkeze mi göre
   oluştuğunu doğrudan gösterir.

**Görüş isteniyor:** Bu analizler ana deneyin parçası mı olmalı, yoksa ayrı
bir bölüm/ablation olarak mı ele alınmalı?

---

### Soru 6 — Stain normalization

Mevcut karar: **primary deneyde uygulanmıyor**, ayrı bir ablation olarak
değerlendirilecek.

Gerekçe: foundation model karşılaştırmasında ön işlemenin gereksiz
değiştirilmesi modellerin karşılaştırılabilirliğini etkileyebilir. Nitekim
2025 tarihli bir UNI/Virchow2/CONCH karşılaştırması minimal preprocessing
(stain normalization yok, color augmentation yok) yaklaşımını benimsemiştir.

Ancak Soru 5'teki bulgu ışığında bu kararın yeniden değerlendirilmesi
gerekebilir. Howard et al. (2021) renk normalizasyonunun merkez imzasını
ortadan kaldırmadığını göstermektedir; dolayısıyla normalization'ı batch
effect çözümü olarak konumlandırmak yanıltıcı olur. Ablation olarak ele
alınması, "renge ne kadar bağımlıyız?" sorusunu ölçmek açısından yine de
değerlidir.

#### Mürekkep/kalem izi — literatür bu endişeyi doğruluyor

2025 tarihli bir foundation model benchmark çalışmasında, UNI, Hibou,
Virchow ve Kaiko gibi modellerin zaman zaman kalem izlerini öne çıkardığı
gözlenmiştir; yazarlar bunu, tahminlerin altta yatan biyolojiyi anlamak
yerine bir tür örüntü ilişkilendirmesiyle yapıldığını düşündüren istenmeyen
bir davranış olarak nitelendirmektedir.

Aynı çalışmada kalem izlerinin bir kohortun slaytlarının %90'ında, bir
başkasında %22'sinde bulunduğu, diğerlerinde ise hiç bulunmadığı
ölçülmüştür — yani mürekkep varlığı da merkeze bağlı sistematik bir
değişkendir.

Bu, Bölüm 0.6'da tarif edilen ink filtreleme çalışmasının gerekliliğini
desteklemektedir. Ancak filtrenin hematoksilen yoğun gerçek dokuyu
elememesi kritik önemdedir.

**Bu çalışmadaki ölçüm literatürü doğrulamaktadır.** H&E deconvolution
tabanlı tarama sonucunda 468 WSI'ın 21'inde belirgin kalem izi bulunmuş ve
**bu 21 WSI'ın tamamı tek bir kaynak merkezden (TSS 86) gelmiştir.** Merkez
ortalamaları: TSS 86 %11,58; ikinci sıradaki TSS 93 %1,65.

Yani mürekkep varlığı merkez imzasının bir bileşenidir. Merkez sınıf
etiketini belirlediği için (Soru 5), kalem izleri temizlenmeseydi model
bunları sınıf sinyali olarak kullanabilirdi. Kalem izi tespit edilen
koordinatlar aday havuzundan çıkarılmıştır (47.057 aday).

---

## 6 · KOD HARİTASI

Depo: `github.com/esmnralican/tez-wsi-foundation` (private)

```
~/tez/
├── run.py                    ← tek giriş noktası, alt komutlar
├── tez/
│   ├── barcode.py            ← TCGA barkod mantığı (saf, I/O yok)
│   ├── cohort.py             ← kohort kurma + dışlama kaydı
│   ├── wsi.py                ← MPP ölçümü, level seçimi, okuma testi
│   └── patches.py            ← aday tarama, manifest, çıkarma
├── docs/
│   └── patch_kararlari.md    ← karar ve hata kaydı
└── out/LUAD/                 ← üretilen dosyalar (git'e girmez)
    ├── cohort.csv
    ├── excluded.csv
    ├── wsi_qc.csv
    ├── scan_summary.csv
    └── candidates/*.csv
```

### 6.1 · Rapordaki bölümlerin kod karşılığı

| Rapor bölümü | Modül | Fonksiyon / alan |
|---|---|---|
| 2.1 Kohort kuralları | `tez/barcode.py` | `cohort_reason()`, `COHORT_SAMPLE_TYPE`, `COHORT_SLIDE_NO` |
| 2.2 Barkod ayrıştırma | `tez/barcode.py` | `parse()`, `BARCODE_RE` |
| 2.3 Dışlama tablosu | `tez/cohort.py` | `scan()` → `excluded.csv` |
| 3.2 Level seçimi | `tez/wsi.py` | `choose_level()` |
| 3.3 MPP ölçümü | `tez/wsi.py` | `read_mpp()` |
| 4.6 Okuma süresi | `tez/wsi.py` | `read_test()` |
| 4.1 Aday tarama | `tez/patches.py` | `scan_slide()` |
| 4.2–4.5 Dağılımlar | `tez/patches.py` | `scan_report()` |
| Doku maskesi | `tez/patches.py` | `tissue_and_texture_mask()` |

### 6.2 · Çalıştırma sırası

```bash
python run.py cohort     --data ~/tcga_data/LUAD --label LUAD \
                         --primary <slayt_klasoru> --out out/LUAD
python run.py check      --out out/LUAD --expect-wsi 468
python run.py wsi-qc     --slides ~/tcga_data/LUAD --out out/LUAD
python run.py check      --out out/LUAD --expect-wsi 468
python run.py patch-scan --slides ~/tcga_data/LUAD --out out/LUAD
```

### 6.3 · Doğrulama yaklaşımı

Her aşama bir doğrulama adımıyla eşleşir. `check` komutu 13 kontrol
çalıştırır (barkod ayrıştırma, sample type, slayt numarası, hasta
mükerrerliği, beklenen sayı, kohort/dışlanan kesişimi, MPP okunabilirliği,
yukarı ölçekleme kontrolü, vb.). Herhangi biri başarısızsa çıkış kodu 1
döner ve sonraki aşamaya geçilmez.

`barcode.py` ayrıca kendi birim testini içerir (`self_test()`), kohort
kurulmadan önce otomatik çalışır.

---

## 7 · TESPİT EDİLEN METODOLOJİK RİSKLER

| # | Risk | Durum |
|---|---|---|
| 1 | Kaynak merkez ile sınıf etiketi arasında %100 ilişki | **Açık** — Soru 5. Ölçülecek, çözülemez |
| 2 | Düşük doku kapsamı morfolojiyi kaçırabilir | Azaltıldı — 100 → 2000 patch (%0,9 → %19 kapsam) |
| 3 | 112 µm FOV doku mimarisi için dar olabilir | **Açık** — Soru 3 |
| 4 | Doku oranı thumbnail'de ölçülüyor (kaba) | Çözüldü — patch seviyesinde ikinci ölçüm eklendi; metriğin boyamaya bağımlı olduğu görülünce filtre olarak kullanılmadı |
| 5 | Hasta sızıntısı | Çözüldü — 467 WSI = 467 hasta, split `patient` üzerinden |
| 6 | Mürekkep tespitinin hematoksileni elemesi | Çözüldü — H&E deconvolution kalıntısı, eşik 0,55 `[ölçüm]` |
| 7 | Uniform renkli background'ın doku sanılması | Çözüldü — blur QC zorunlu (Laplacian ≥ 100) |
| 8 | Nekroz artefakt sanılıp silinmesi | Nekroz gerçek biyolojik yapıdır; otomatik silinmiyor |
| 9 | Dışlanan 10 WSI'ın tek merkezden olması | Tezde açıkça raporlanacak |
| 10 | Kalem izi varlığının merkeze bağlı olması | Çözüldü — 21 WSI'da tespit, koordinatlar havuzdan çıkarıldı. **Bulgu tezde raporlanacak** |
| 11 | Bag boyutu dengesizliği (78–2000 patch) | **Açık** — Soru 4 |
| 12 | Örnekleme determinizmi (`hash()` süreç başına değişiyor) | Çözüldü — `blake2b`. LUSC bu düzeltmeyle işlenecek |
| 13 | Boyama yoğunluğunun merkeze göre değişmesi (0,58–0,96) | Ölçüldü. Filtreleme kararlarında dikkate alındı; ablation'da nicel ölçü olarak kullanılabilir |
| 15 | Kalem maskesi tek ölçütle yeşil/mavi kalemi kaçırıyordu | Çözüldü — birleşim maskesi, 18 slayt yeniden çıkarıldı |
| 14 | Native 20× oranı sınıflar arası dengesiz (LUAD %4, LUSC %12) | **Açık** — `base_mpp` manifest'te saklı, değerlendirmede kontrol edilecek |

---

## 8 · SONRAKİ ADIMLAR

**Tamamlananlar**

1. ✅ Kohort tanımı ve doğrulama — 467 WSI / 467 hasta
2. ✅ WSI teknik QC — MPP ölçümü, level seçimi
3. ✅ Patch aday taraması (tavansız) — ~5,3 M aday
4. ✅ Kalem izi tespiti ve kalibrasyonu
5. ✅ `canonical_patch_manifest.csv` — üç modelin ortak koordinat kaynağı
6. ✅ Patch çıkarma + görsel QC — 884.987 patch

7. ✅ LUSC indirme, kohort ve WSI QC — 478 WSI / 478 hasta
8. ✅ Merkez örtüşmesi ölçümü — ortak merkez 0

**Sıradakiler**

9. ✅ LUSC patch çıkarma — 934.875 patch
10. ✅ Kalem maskesi düzeltmesi — 18 slayt yeniden çıkarıldı

**Sıradakiler**

11. ⬜ Hasta seviyesinde split (`GroupKFold`, TSS stratifikasyonu dahil)
12. ⬜ DINOv2 / UNI / Virchow embedding çıkarma (frozen backbone)
13. ⬜ ABMIL eğitimi ve hasta seviyesinde değerlendirme
14. ⬜ Site-aware analizler: site prediction baseline, leave-two-sites-out
15. ⬜ Stain normalization ablation (`patches/` → `patches_macenko/`)
16. ⬜ Açıklanabilirlik katmanı (yöntem henüz seçilmedi)

---

## 8.1 · LUAD PATCH ÇIKARMA SONUÇLARI

| | |
|---|---|
| **Analiz kohortu** | **466 WSI / 466 hasta** |
| Toplam patch | **900.828** |
| Disk | ~20 GB |

### Kesinleşen patch protokolü

| Parametre | Değer | Kaynak |
|---|---|---|
| Hedef çözünürlük | 0,50 µm/px | `[literatür]` UNI/Virchow eğitim ölçeği |
| Patch boyutu | 224 × 224 px | `[literatür]` üç modelin ortak girdisi |
| Fiziksel görüş alanı | 112 × 112 µm | türetilmiş |
| Stride | 112 µm (overlap yok) | `[literatür]` CLAM |
| Doku oranı (thumbnail) | ≥ 0,80 | `[ölçüm]` eşik düşük etkili |
| Kalem izi | kalıntı > 0,55 → dışla | `[ölçüm]` kalibre edildi |
| Slayt başına tavan | 2000 | `[ölçüm]` süre ve disk ölçüldü |
| Blur (resize sonrası) | Laplacian ≥ 100 | `[ölçüm + görsel]` |
| Patch seviyesi ink/doku | filtre yok | `[görsel]` merkez-bağımlı seçilim riski |
| Stain normalization | yok (ablation planlı) | `[literatür]` |

### Tez metni için önerilen ifade

> Patches were extracted at a target resolution of 0.50 µm/px
> (224 × 224 px, 112 × 112 µm field of view) without overlap. Pyramid level
> selection was based on the measured MPP of each slide rather than the
> reported objective power, and reading was always performed from a level at
> or above the target resolution followed by software downsampling; no
> upsampling was applied. Candidate coordinates required a tissue fraction
> ≥ 0.80 on a downsampled tissue mask; regions containing pen markings
> (H&E stain deconvolution residual > 0.55) were excluded. Up to 2,000
> patches per slide were sampled with a fixed seed. Patches with a Laplacian
> variance < 100 after resizing to 224 × 224 were discarded as out of focus.
> This yielded 887,665 patches from 466 WSIs.

---

## 8.1b · KALEM MASKESİ DÜZELTMESİ

Patch çıkarma sonrası, yoğun kalemli slaytlarda patch sayısının anormal
düşük kaldığı görülmüştür. Çıkarma kayıtları sebebi göstermiştir:

```
TCGA-85-8276:  3.000 aday -> blur eledi 2.451 -> 549 patch
TCGA-86-7701:  3.000 aday -> blur eledi 2.576 -> 424 patch
```

Mürekkep filtresi hiçbir patch elememiş, blur filtresi binlercesini
elemiştir. Elenen patch'lerin görsel incelemesinde bunların **düz kırmızı ve
düz turkuaz kareler** olduğu, yani doku içermeyen katı mürekkep blokları
olduğu görülmüştür.

Yani kalem maskesi bu bölgeleri kaçırmış, blur filtresi onun yerine
temizlemiştir.

### Ölçüm

`TCGA-85-8276` üzerinde, doku içinde yakalanan oran:

| Ölçüt | Yakalanan |
|---|---:|
| `residual > 0,55` (mevcut) | %11,5 |
| `residual > 0,45` | %37,5 |
| `negatif konsantrasyon < −0,05` | %28,6 |
| **Birleşim** | **%66,0** |

Ayrıca doygunluk koşulunun (`saturation > 60`) doku piksellerinin %100'ünü
geçirdiği, dolayısıyla maskeye hiçbir katkısı olmadığı ölçülmüş ve
kaldırılmıştır.

### Neden iki ölçüt gerekiyor

| Ölçüt | Yakaladığı renk | Fiziksel sebep |
|---|---|---|
| Kalıntı | Kırmızı / magenta | Eozin yönüne yakın; H&E ile temsil edilebiliyor ancak kalıntı bırakıyor |
| Negatif konsantrasyon | Yeşil / mavi / turkuaz | H&E düzleminin ters yönünde; ancak negatif katsayıyla temsil edilebiliyor |

Hiçbiri tek başına yeterli değildir.

### Uygulama

```python
pen = (residual > 0.45) | (negatif_konsantrasyon < -0.05)
```

Birleşim maskesi kalemsiz slaytlarda normal dokuyu da elediğinden
(`TCGA-05-4433`: %19,7 → %8,0) tüm kohorta uygulanmamıştır. Yalnızca
`extract_qc.csv`'de blur elemesi 1.000'in üzerinde olan **18 slayt** (LUAD
15, LUSC 3) yeniden çıkarılmıştır.

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

Blur eleme oranı %80'den %5–30 aralığına inmiştir.

---

## 8.2 · LUSC KOHORTU

| | |
|---|---|
| GDC tanısal slayt kaydı | 512 |
| İndirilen (01 + DX1) | **478** |
| Atlanan | 34 (tümü ikincil tanısal slayt, DX2–DX9) |
| **Analiz kohortu** | **478 WSI / 478 hasta** |
| Dışlanan | **0** |
| İndirme boyutu | 376 GB |

LUAD'dan farklı olarak hiç dışlama gerekmemiştir: nüks tümör örneği yok,
MPP'si okunamayan WSI yok, hasta mükerrerliği yok. Barkod filtresi indirme
aşamasında uygulandığından ikincil tanısal slaytlar (DX2–DX9) hiç
indirilmemiş, yalnızca CSV'ye kaydedilmiştir.

`run.py check` 13 kontrolün tamamını geçmiştir. Tüm WSI okunabilir, tümünde
MPP mevcut, tümü level 0'dan okunacak, hiçbiri yukarı ölçeklenmiyor.

### MPP dağılımı — kohortlar arası fark

| Ölçülen MPP | LUAD | LUSC |
|---|---:|---:|
| ~0,25 (40×) | 414 (%88) | 414 (%87) |
| ~0,23 (43×) | 35 (%7) | 5 (%1) |
| **0,50 (native 20×)** | **19 (%4)** | **59 (%12)** |

**Metodolojik not.** LUSC'ta native 20× tarama oranı LUAD'ın üç katıdır.
Native 20× slaytlar hedef çözünürlükte olduğundan küçültülmez; 40× slaytlar
~2× küçültülür ve `INTER_AREA` interpolasyonu hafif bir yumuşatma etkisi
yaratır.

Bu, iki grup arasında küçük bir doku karakteri farkı doğurmakta ve bu fark
**sınıflar arasında dengesiz dağılmaktadır**. Tarayıcı seçimi de merkeze
bağlı olduğundan bu, Soru 5'teki merkez imzası probleminin bir bileşenidir.

`base_mpp` alanı manifest'te saklandığından, sonuçların 20× ve 40×
gruplarında farklılaşıp farklılaşmadığı değerlendirme aşamasında kontrol
edilebilir.

### Patch çıkarma sonuçları

| | |
|---|---|
| **Analiz kohortu** | **478 WSI / 478 hasta** |
| Toplam patch | **934.875** |
| Disk | ~30 GB |
| Okuma süresi | ortanca 12 ms/patch |

### Düşük bag'ler (her iki kohort)

| Slayt | Patch | Sebep |
|---|---:|---|
| TCGA-55-8207 (LUAD) | 78 | 2 mm² doku |
| TCGA-67-4679 (LUAD) | 289 | Doku %17, blur 5.776 eledi |
| TCGA-98-8022 (LUSC) | 290 | Tarama kalitesi (blur p50 = 869 vs kohort 1.089) |
| TCGA-55-8301 (LUAD) | 326 | Az doku |
| TCGA-56-8082 (LUSC) | 381 | Parçalı doku |

Soru 4'teki bag boyutu kararı bu slaytlar için geçerlidir.

---

## 9 · KAYNAKLAR

1. **Howard FM, Dolezal J, Kochanny S, et al.** The impact of site-specific
   digital histology signatures on deep learning model accuracy and bias.
   *Nature Communications* 2021;12:4423.
   https://doi.org/10.1038/s41467-021-24698-1
   → Bölüm 5, Soru 5 ve Soru 6'nın temel referansı.

2. **Lu MY, Williamson DFK, Chen TY, et al.** Data-efficient and weakly
   supervised computational pathology on whole-slide images (CLAM).
   *Nature Biomedical Engineering* 2021;5:555–570.
   → Patch çıkarma ve doku segmentasyonu standardı.

3. **Chen RJ, Ding T, Lu MY, et al.** Towards a general-purpose foundation
   model for computational pathology (UNI). *Nature Medicine* 2024;30:850–862.
   → ABMIL konfigürasyonu ve embedding protokolü.

4. **Vorontsov E, Bozkurt A, Casson A, et al.** A foundation model for
   clinical-grade computational pathology and rare cancers detection
   (Virchow). *Nature Medicine* 2024;30:2924–2935.

5. **Oquab M, Darcet T, Moutakanni T, et al.** DINOv2: Learning robust
   visual features without supervision. *TMLR* 2024.

6. **Ilse M, Tomczak JM, Welling M.** Attention-based deep multiple
   instance learning (ABMIL). *ICML* 2018.

7. Foundation model benchmark (2025) — kalem izi/attention bulgusu için;
   Bölüm 5, Soru 6 notu.

> Kaynak 1 mutlaka okunmalıdır: bu çalışmanın en ciddi metodolojik riskinin
> (Soru 5) literatürdeki tanımı ve önerilen çözümü bu makaledir.

---

## 10 · İKİ KOHORT — ÖZET

| | WSI | Hasta | Patch | Merkez |
|---|---:|---:|---:|---:|
| **LUAD** | 466 | 466 | 900.828 | 33 |
| **LUSC** | 478 | 478 | 934.875 | 35 |
| **Toplam** | **944** | **944** | **1.835.703** | **68 (ortak: 0)** |

Her iki kohort da aynı protokolle işlenmiştir: 0,50 µm/px hedef çözünürlük,
224 × 224 px (112 × 112 µm), örtüşmesiz, slayt başına en fazla 2.000 patch.

Üç foundation model (DINOv2, UNI, Virchow) aynı koordinatlardan, aynı
fiziksel alanı temsil eden, aynı ön işlemeden geçmiş patch'leri görecektir.
Karşılaştırmada değişen tek değişken model olacaktır.
