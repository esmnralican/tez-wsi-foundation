# DİKKAT AĞIRLIKLARININ ANALİZİ

Tarih: 8 Ekim 2026
Kapsam: DINOv2 + ABMIL, merkez-korumalı bölme, fold 0 (190 test slaydı)

---
~/tez/docs/dikkat_analizi.md olarak kaydet.

Görseller zaten out/dikkat/dinov2_site_k5_fold0/ altında. Tezde şekil olarak kullanacağın üç tanesi:

Dosya	Tezdeki rolü
TCGA-85-A50M_dogru.jpg	Şekil A — LUSC, morfolojiye odaklanma
TCGA-55-A57B_dogru.jpg	Şekil B — LUAD, morfolojiye odaklanma
TCGA-37-3783_yanlis.jpg	Şekil C — başarısızlık modu

Üçünü yan yana koyarsan hem doğru davranışı hem hata mekanizmasını tek sayfada gösterirsin.


## 1 · GEREKÇE

Bölüm 10'da tespit edilen açık soru, modelin kurumlar arası aktardığı
sinyalin morfolojiye mi yoksa birçok kurumda ortak olan bir edinim
karakteristiğine mi karşılık geldiğidir. Nicel ölçümler bu iki olasılığı
ayırt edememektedir.

Dikkat tabanlı çoklu örnek öğrenmesi, karara katkı veren örnekleri
açıkça ağırlıklandırdığından, bu ağırlıkların incelenmesi doğrudan nitel
kanıt sağlar: modelin en yüksek ağırlık verdiği bölgeler morfolojik yapı
mı yoksa artefakt mı içermektedir?

## 2 · YÖNTEM

Merkez-korumalı bölmenin ilk fold'u üzerinde model eğitilmiş, aynı
fold'un 190 test slaydında dikkat ağırlıkları çıkarılmıştır. Her slayt
için en yüksek ve en düşük ağırlıklı on patch yan yana yerleştirilerek
kontak sayfaları üretilmiştir.

Dikkat ağırlıkları softmax çıktısıdır ve toplamları bire eşittir.
Görsellerdeki etiketler düzgün dağılıma göre normalize edilmiştir:
`3,0x` ifadesi, o patch'in düzgün dağılımda alacağı ağırlığın üç katını
aldığını gösterir.

İki yoğunlaşma ölçüsü hesaplanmıştır:

**Üst %1 dikkat payı.** En yüksek ağırlıklı %1'lik patch diliminin
(2.000 patch için yaklaşık 20 patch) toplam dikkat içindeki payı.
Düzgün dağılımda bu değer 0,01'dir.

**Etkin patch sayısı.** Dikkat dağılımının perplexity değeri,
exp(−Σ aᵢ log aᵢ). Düzgün dağılımda patch sayısına, tek bir patch'e
yoğunlaşmada bire eşittir. Karara fiilen katılan patch sayısının ölçüsüdür.

## 3 · BULGULAR
### 3.1 · Dikkatin yoğunlaşma derecesi

| Ölçü | Ortanca | Düzgün dağılımda |
|---|---|---|
| Üst %1 dikkat payı | **0,596** | 0,010 |
| Etkin patch sayısı | **66** | 2.000 |

En yüksek ağırlıklı yirmi patch, toplam dikkatin yaklaşık %60'ını
almaktadır; bu, düzgün dağılıma göre altmış kat yoğunlaşmaya karşılık
gelir. İki bin patch'ten yalnızca yaklaşık altmış altısı karara fiilen
katılmaktadır.

Bu bulgu, Bölüm 9.2'deki örnekleme yoğunluğu ablasyonuyla tutarlıdır:
hasta-katmanlı değerlendirmenin 1.000 patch civarında doyuma ulaşması,
modelin zaten sınırlı sayıda patch kullanmasıyla açıklanabilir.

### 3.2 · Dikkatin yöneldiği morfoloji

Üç temsili örnek incelenmiştir.

**TCGA-85-A50M (LUSC, doğru sınıflandırıldı, skor 0,944).**
En yüksek ağırlıklı patch'ler yoğun tümör hücre adaları, büyük ve
hiperkromatik çekirdekler, yüksek çekirdek/sitoplazma oranı ve kısmen
eozinofilik sitoplazma göstermektedir — skuamöz hücreli karsinom ile
uyumlu morfoloji. En düşük ağırlıklı patch'ler antrakotik pigment
birikimi, fibröz stroma ve hücreden yoksun alanlardan oluşmaktadır.

**TCGA-55-A57B (LUAD, doğru sınıflandırıldı, skor 0,162).**
En yüksek ağırlıklı patch'ler iri veziküler çekirdekli, belirgin
nükleollü ve berrak sitoplazmalı hücre grupları içermektedir —
adenokarsinom ile uyumlu morfoloji. En düşük ağırlıklı patch'ler
kollajen demetleri ve hücresiz bağ dokusudur.

**TCGA-37-3783 (LUSC, yanlış sınıflandırıldı, skor 0,001).**
En yüksek ağırlıklı patch 1.688 kat ağırlık almış, etkin patch sayısı
2,8'e düşmüştür. İlgili bölgeler alveolar yapı, lenfositik infiltrasyon
ve ince septalar göstermektedir; skuamöz tümör morfolojisi
içermemektedir. Model temsil edici olmayan birkaç bölgeye yoğunlaşmış ve
yanılmıştır.

**Değerlendirme.** İncelenen slaytlarda yüksek ağırlık hücresel tümör
bölgelerine yönelmekte, stroma ve antrakotik pigment düşük ağırlık
almaktadır. Üst sıralarda mürekkep kenarı, arka plan veya kesit
artefaktı gözlenmemiştir.

Bu, modelin kurumlar arası aktardığı sinyalin morfolojik olduğu yorumunu
destekleyen nitel kanıttır. Nicel bir ayrıştırma sağlamaz; gri tonlama
ablasyonu bu amaçla planlanmaktadır.

### 3.3 · Yoğunlaşma ve doğruluk ilişkisi

Dikkat yoğunlaşmasının sınıflandırma doğruluğunu öngörüp öngörmediği,
her iki bölme şemasının beş fold'unda ayrı ayrı sınanmıştır (toplam on
model, 1.888 test slaydı).

Havuzlanmış analizde doğru sınıflandırılan slaytlarda dikkat belirgin
biçimde daha yoğundur (üst %1 payı 0,514'e karşı 0,448; etkin patch
109'a karşı 149; Mann-Whitney p = 0,0003 ve p < 0,0001).

**Ancak bu sonuç fold'lar arası farklardan kaynaklanmaktadır.** Fold'lar
hem doğruluk hem yoğunlaşma bakımından sistematik olarak
farklılaşmaktadır: dikkatin dağınık olduğu fold'lar aynı zamanda düşük
doğruluğa sahiptir (örneğin merkez-korumalı fold 2'de ortanca etkin
patch 388, doğruluk 0,762; hasta-katmanlı fold 4'te sırasıyla 137 ve
0,899). Slaytların havuzlanması bu fold seviyesi ilişkiyi slayt seviyesi
bir ilişki gibi göstermektedir.

Fold içi analiz — her fold içinde yoğunlaşmanın doğru sınıflandırmayı
öngörme gücü, ROC eğrisi altındaki alan olarak — bu yorumu
desteklememektedir:

| Ölçü | Ortalama AUC | Std | Şans üstü fold |
|---|---|---|---|
| Üst %1 dikkat payı | 0,568 | 0,094 | 6/10 |
| Etkin patch (ters) | 0,581 | 0,081 | 6/10 |

Fold değerleri 0,458 ile 0,725 arasında değişmektedir; dört fold'da
ilişki şans seviyesinin altındadır. Ortalamanın yaklaşık %95 güven
aralığı 0,51–0,63'tür ve sıfır etkiyi ancak sınırda dışlamaktadır.

**Sonuç:** Dikkat yoğunlaşması, tek tek slaytlar düzeyinde güvenilir bir
belirsizlik göstergesi sağlamamaktadır. Bu çalışmada bir bulgu olarak
raporlanmamaktadır.

#### Keşifsel gözlem

Bölme şemasına göre ayrıldığında ilişki merkez-korumalı şemada daha
güçlüdür (ortalama AUC 0,615'e karşı 0,521). Merkez kestirmesinin
bulunmadığı koşulda dikkat yoğunlaşmasının modelin gerçek kanıt bulup
bulmadığını daha iyi yansıttığı düşünülebilir. Her grupta beş fold
bulunduğundan ve fold'lar arası sapma yüksek olduğundan bu gözlem
hipotez üretme amaçlıdır; doğrulanması bağımsız bir çalışma gerektirir.

#### Dilimlere göre doğruluk

Havuzlanmış veride yoğunlaşma dilimlerine göre doğruluk:

| Etkin patch | n | Doğruluk |
|---|---|---|
| < 20 | 150 | 0,853 |
| 20–60 | 420 | **0,905** |
| 60–150 | 556 | 0,840 |
| 150–300 | 412 | 0,828 |
| > 300 | 350 | 0,774 |

İlişki monoton değildir: en yüksek doğruluk orta düzeyde yoğunlaşmada
(20–60 etkin patch) gözlenmekte, aşırı yoğunlaşma (< 20) ve aşırı
dağınıklık (> 300) birlikte daha düşük doğruluk vermektedir. Aşırı
yoğunlaşmanın hataya yol açtığı bir örnek Bölüm 3.2'de verilmiştir
(TCGA-37-3783, etkin patch 2,8).

Fold'lar arası karışım nedeniyle bu tablo da nedensel olarak
yorumlanmamalıdır.
### 3.3 · Yoğunlaşma ve doğruluk ilişkisi   ESKİ HALİİ

| | n | Üst %1 dikkat payı | Etkin patch |
|---|---|---|---|
| Doğru sınıflandırılan | 145 | **0,619** | **60** |
| Yanlış sınıflandırılan | 45 | 0,534 | 95 |

| Dikkat dağılımı | Doğruluk |
|---|---|
| Yoğun (etkin patch < 20) | **16/19 = %84** |
| Dağınık (etkin patch > 150) | 28/42 = **%67** |

Doğru sınıflandırılan slaytlarda dikkat belirgin biçimde daha
yoğunlaşmıştır. Dikkatin iki binden fazla patch'e yayıldığı slaytlarda
doğruluk düşmektedir.

**Yorum.** Model ayırt edici morfolojik kanıt bulduğunda sınırlı sayıda
bölgeye yoğunlaşmakta ve doğru karar vermektedir. Böyle bir kanıt
bulunmadığında dikkat dağılmakta ve karar güvenilmez hale gelmektedir.

TCGA-37-3783 örneği bu eğilimin istisnasıdır: aşırı yoğunlaşma, seçilen
bölgelerin temsil edici olmaması durumunda hataya yol açmaktadır.

**Uygulama değeri.** Dikkat yoğunlaşması, model çıktısının güvenilirliği
için ölçülebilir bir gösterge sunmaktadır. Klinik karar destek
bağlamında, dikkati dağınık olan slaytların insan incelemesine
yönlendirilmesi için kullanılabilir. Bu, tahmin olasılığından bağımsız
ve doğrudan modelin iç davranışından türetilen bir ölçüdür.

### 3.4 · Tutarlılık kontrolü

Bu analizde elde edilen fold doğruluğu 145/190 = 0,763 olup Tablo 3'te
raporlanan merkez-korumalı fold 0 doğruluğuyla (0,7632) birebir
örtüşmektedir.

## 4 · SINIRLILIKLAR

Nitel değerlendirme patolog onayı almamıştır; morfolojik tanımlamalar
genel histopatolojik örüntülere dayanmaktadır.

Analiz tek bir fold üzerinde yürütülmüştür. Yoğunlaşma–doğruluk
ilişkisinin diğer fold'larda ve hasta-katmanlı bölmede tekrarlanması
gerekmektedir.

Dikkat ağırlıkları modelin hangi bölgeleri kullandığını gösterir, o
bölgelerin **hangi özelliğinin** kullanıldığını göstermez. Yüksek
ağırlıklı bir patch'in tümör morfolojisi içermesi, modelin o morfolojiyi
kullandığını kanıtlamaz; aynı patch kuruma özgü bir boyama
karakteristiği de taşıyor olabilir. Gri tonlama ablasyonu bu ayrımı
hedeflemektedir.

## 5 · PLANLANAN DOĞRULAMA

| İş | Ne cevaplar |
|---|---|
| Hasta-katmanlı fold 4 (AUC 0,962) üzerinde aynı analiz | İyi çalışan bir fold'da yoğunlaşma daha mı yüksek |
| Tüm fold'larda yoğunlaşma–doğruluk ilişkisi | Bulgu tekrarlanabilir mi |
| Patolog değerlendirmesi | Morfolojik tanımlamaların doğrulanması |
| Gri tonlama ablasyonu | Yüksek ağırlıklı bölgelerde kullanılan özellik renk mi |