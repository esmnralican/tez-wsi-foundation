# MERKEZ (TSS) ETKİSİ — ÖLÇÜM VE LİTERATÜR

Tarih: 7 Ekim 2026
Kapsam: Kaynak merkez imzasının sınıflandırmaya etkisi, ölçüm sonuçları,
literatürdeki durum ve bu çalışmada izlenecek yol.

---

## 1 · PROBLEMİN YAPISI

| | LUAD | LUSC |
|---|---:|---:|
| WSI | 466 | 478 |
| Farklı kaynak merkez | 33 | 35 |
| **Ortak merkez** | **0** | **0** |

Hiçbir kaynak merkez hem LUAD hem LUSC slaytı göndermemiştir.
**Kaynak merkez bilgisi sınıf etiketini tam olarak belirlemektedir.**

Merkezler arasında H&E boyama protokolü, tarayıcı modeli, fiksasyon
süresi ve kesit kalınlığı sistematik olarak farklılaşır. Bu koşulda bir
model morfolojiyi hiç öğrenmeden, yalnızca merkez imzası üzerinden
yüksek doğruluk elde edebilir.

### 1.1 · Neden bölme stratejisi tek başına çözmüyor

Site-preserved cross-validation (hiçbir merkez hem eğitimde hem testte
olmayacak şekilde bölme) **ölçümü dürüstleştirir, temsili düzeltmez.**
Merkez bilgisi embedding'lerde durmaya devam eder.

Bu çalışmada daha temel bir kilit vardır: merkezler tek sınıflı olduğu
için model hiçbir zaman *"aynı merkez içinde LUAD ile LUSC'u ne ayırır"*
sorusunu öğrenemez. Her merkezde tek sınıf bulunduğundan
"bu merkezin görünümü = bu sınıf" kuralı eğitim verisinde daima
doğrudur.

Merkezler karışık sınıflı olsaydı merkez-içi karşılaştırma mümkün
olurdu. Burada yapısal olarak mümkün değildir.

**Sonuç:** Hiçbir bölme stratejisi bu konfounding'i ortadan kaldıramaz.
Gerçek çözümler yalnızca (a) harici kohort ve (b) ölçüp raporlamaktır.

---

## 2 · ÖLÇÜM (DINOv2, 100 patch/WSI pilot)

5-fold çapraz doğrulama, ABMIL, iki farklı bölme ile:

| | patient bölme | site bölme | fark |
|---|---:|---:|---:|
| AUC | 0,9168 | 0,8676 | **−0,049** |
| Doğruluk | 0,8337 | 0,7565 | **−0,077** |
| Dengeli doğruluk | 0,8329 | 0,7553 | −0,078 |
| std (AUC) | 0,021 | 0,025 | |

Düşüş beş fold'un hepsinde tutarlıdır.

### 2.1 · Val–test uçurumu

Site bölmesinde doğrulama skorları test skorlarının belirgin üzerindedir:

| fold | val AUC | test AUC | fark |
|---|---:|---:|---:|
| 0 | 0,963 | 0,832 | −0,131 |
| 1 | 0,925 | 0,854 | −0,071 |
| 2 | 0,929 | 0,873 | −0,057 |
| 3 | 0,932 | 0,906 | −0,026 |
| 4 | 0,939 | 0,874 | −0,065 |

Patient bölmesinde böyle sistematik bir açık yoktur; test skoru zaman
zaman val'in üzerindedir (fold 0: val 0,893 → test 0,949).

**Sebep:** İç doğrulama kümesi eğitim fold'larından ayrılmaktadır, yani
aynı merkezlerden gelmektedir. Test ise site bölmesinde hiç görülmemiş
merkezlerdendir. Model tanıdığı merkezlerde iyi, yabancı merkezlerde
belirgin biçimde kötüdür.

> **Düzeltilecek:** Erken durdurma, merkez bilgisine yaslanan bir modeli
> "iyi" sanarak seçmektedir. İç doğrulama kümesi de merkez seviyesinde
> ayrılmalıdır.

### 2.2 · Yorumlama sınırı

Mutlak değerler tezde kullanılmayacaktır. 100 patch, ortanca slaydın
%0,9'udur; bu kadar az instance ile attention mekanizması anlamlı bir
şey öğrenemez. Bu koşunun kanıtladığı tek şey boru hattının uçtan uca
çalıştığıdır.

Anlamlı olan **iki bölme arasındaki farktır** ve bu farkın 2000 patch'te
nasıl değiştiği ölçülecektir.

---

## 3 · LİTERATÜR

### 3.1 · Howard et al. 2021 — problemin tanımı

**Howard FM, Dolezal J, Kochanny S, et al.** The impact of site-specific
digital histology signatures on deep learning model accuracy and bias.
*Nature Communications* 2021;12:4423.
https://doi.org/10.1038/s41467-021-24698-1

TCGA'daki 6 kanser tipinde 3.000'den fazla hasta üzerinde, tahmin edilen
özelliklerin doku gönderen merkezler arasında ciddi biçimde
farklılaştığını; merkezler arası görüntü farklarının derin öğrenme ile
kolayca tespit edilebildiğini göstermiştir.

Yaygın renk normalizasyonu ve augmentation yöntemlerine rağmen merkez
tespiti mümkün olmaya devam etmektedir. Önerdikleri çözüm, modellerin
aynı merkezden gelen örneklerle hem eğitilip hem doğrulanmamasını
sağlayan bir kuadratik programlama yöntemidir.

### 3.2 · Foundation model'ler batch effect'i ortadan kaldırıyor mu

**Do Histopathological Foundation Models Eliminate Batch Effects?
A Comparative Study.** arXiv:2411.05489
https://arxiv.org/pdf/2411.05489

Bu çalışma doğrudan bu tezin konusuyla örtüşmektedir: TCGA-LUSC
(5 hastane) ve CAMELYON16 üzerinde UNI, Virchow, Prov-GigaPath ve
Phikon incelenmiştir.

Bulgular:

- Kaynak merkez, foundation model embedding'lerinden **%80–99
  doğrulukla** tahmin edilebilmektedir
- **Daha iyi performans gösteren modeller daha güçlü merkez imzası
  kodlamaktadır**
- Merkez bilgisi öznitelik uzayındaki mesafelere hâkim olmakta, mesafe
  tabanlı yöntemleri bozmaktadır
- Batch effect ana temel bileşenlerde yoğunlaşmakta, bu da alt görev
  modelleri için kolayca erişilebilir hale getirmektedir

**Stain normalization sonucu:** Reinhard ve Macenko yalnızca mütevazı
bir azalma sağlamıştır (TCGA-LUSC-5'te ~%90'dan ~%80'e) ve dengesiz
veri kümelerinde alt görev performansını anlamlı biçimde
iyileştirmemiştir.

### 3.3 · Diğer kaynaklar

- **Dawood et al. (ICLR-W 2023)** — Do Tissue Source Sites Leave
  Identifiable Signatures in Whole Slide Images Beyond Staining?
  https://mlanthology.org/iclrw/2023/dawood2023iclrw-tissue
  → İmzanın yalnızca boyamadan kaynaklanmadığını gösterir.

- **Physical Color Calibration of Digital Pathology Scanners.**
  arXiv:2307.05519 → Donanım seviyesinde kalibrasyon; mevcut veriye
  geriye dönük uygulanamaz.

- **Beyond Diagnostic Performance: Revealing and Quantifying Ethical
  Risks in Pathology Foundation Models.** arXiv:2502.16889

---

## 4 · NEDEN MACENKO / VAHADANE ÇÖZMÜYOR

Stain deconvolution rengi ayrıştırır; Macenko (SVD) ve Vahadane
(sparse NMF) bunun üstüne normalizasyon ekleyerek boya miktarlarını
referans slayta göre yeniden ölçekler.

Ancak merkez imzası yalnızca renk değildir:

| Bileşen | Normalizasyon etkiler mi |
|---|---|
| Boyama yoğunluğu, renk tonu | Evet |
| Kesit kalınlığı | Hayır |
| Fiksasyon süresi | Hayır |
| Tarayıcı optiği, odaklama | Hayır |
| JPEG sıkıştırma artefaktı | Hayır |
| Kalem/mürekkep kullanım rutini | Hayır |

Ölçüm de bunu doğruluyor: normalizasyon sonrası merkez tahmini
%90'dan %80'e iniyor — azalıyor ama kaybolmuyor.

**Karar:** Macenko bu çalışmada **ablation olarak** ele alınacaktır,
batch effect çözümü olarak değil. Ölçtüğü şey "modelin renge bağımlılığı",
"merkez etkisinin giderilmesi" değildir.

---

## 5 · BU ÇALIŞMA İÇİN YOL HARİTASI

| # | Adım | Ne ölçer | Durum |
|---|---|---|---|
| 1 | Site-preserved CV | Yabancı merkezde gerçek performans | ✅ Pilot yapıldı |
| 2 | Site prediction baseline | İmza embedding'de var mı, ne kadar | ⬜ Sıradaki |
| 3 | İç doğrulamayı da site-aware yapmak | Erken durdurmanın yanlılığı | ⬜ Düzeltilecek |
| 4 | Macenko ablation | Renge bağımlılık | ⬜ Embedding sonrası |
| 5 | Harici kohort (CPTAC) | Gerçek genelleme | ⬜ Takvim uygunsa |
| 6 | TDA Mapper / UMAP — TSS'e göre renklendirme | Kümelenme morfolojiye mi merkeze mi | ⬜ |

### 5.1 · Site prediction baseline — neden önemli

Aynı embedding'lerden sınıf yerine TSS tahmin edilir. Yüksek doğruluk,
merkez bilgisinin embedding'lerde kodlandığını doğrudan gösterir.
Howard et al.'ın ve arXiv:2411.05489'un kullandığı standart testtir.

**Bu çalışma için özgün katkı potansiyeli:** Üç modelde (DINOv2, UNI,
Virchow) ayrı ayrı çalıştırıldığında, literatürdeki *"daha iyi model =
daha güçlü merkez imzası"* iddiası kendi verimizde test edilmiş olur.

Bu, tezin ana sorusuyla doğrudan bağlantılıdır: **UNI/Virchow DINOv2'yi
geçerse, bunun ne kadarı morfoloji, ne kadarı daha iyi merkez ayrımı?**

### 5.2 · Raporlama ilkesi

Tezde **her iki bölmenin sonucu da** raporlanacaktır. Yalnızca patient
bölmesini raporlamak, literatürün çoğunun yaptığı gibi, yanıltıcı olur.

---

## 6 · AÇIKLANABİLİRLİK KATMANINA ETKİSİ

arXiv:2411.05489'un bulgusu — merkez bilgisinin öznitelik uzayındaki
mesafelere hâkim olması — planlanan TDA Mapper / UMAP analizlerini
doğrudan etkiler. Kümelenme morfolojiye değil merkeze göre oluşabilir.

**Önlem:** Mapper grafiği bir kez sınıfa, bir kez TSS'e göre
renklendirilecek. İki renklendirme aynı yapıyı gösteriyorsa kümelenme
morfolojik değil, merkez kaynaklıdır.cd ~/tez
python3 -c "
import csv
from collections import defaultdict
d = defaultdict(lambda: [0,0])
for r in csv.DictReader(open('out/mil/dinov2_n100_site_k5/predictions.csv')):
    d[r['tss']][0] += 1
    d[r['tss']][1] += (r['tahmin'] != r['gercek'])
print(f\"{'tss':>5} {'n':>5} {'hata':>5} {'oran':>7}\")
for t,(n,h) in sorted(d.items(), key=lambda x: -x[1][1]/x[1][0])[:12]:
    print(f'{t:>5} {n:>5} {h:>5} {h/n:>6.0%}')
"

---

## 7 · MERKEZ BAZINDA HATA — ⚠ TEZE KESİN EKLENECEK

**Bu tablo, ortalama metriklerin gizlediği etkiyi doğrudan gösterir ve
tez metnine mutlaka girecektir.**

DINOv2, 100 patch/WSI, 5-fold. Aynı merkezlerin iki bölmedeki hata
oranları:

| TSS | Sınıf | n | patient bölme | site bölme | artış |
|---|---|---:|---:|---:|---:|
| **66** | LUSC | 38 | **%5** | **%58** | **11×** |
| 55 | LUAD | 91 | %33 | %60 | 1,8× |
| 05 | LUAD | 24 | %25 | %50 | 2,0× |
| MP | LUAD | 19 | %16 | %47 | 2,9× |
| 38 | LUAD | 11 | %45 | %45 | — |

### 7.1 · TSS 66 — ezberin kanıtı

Merkez eğitim kümesindeyken model 38 slaydın 36'sını doğru
sınıflandırmaktadır (%5 hata). Aynı merkez eğitimden çıkarıldığında
38 slaydın 22'sini yanlış sınıflandırmaktadır (%58 hata).

İkili sınıflandırmada şans seviyesi %50'dir. **%58 hata, şanstan
kötüdür** — model o merkezin tamamını sistematik olarak ters sınıfa
atamaktadır. Bu gürültü değil, yön hatasıdır.

Yorum: Model bu merkezin morfolojisini öğrenmemiş, *"bu merkezin
görünümü = LUSC"* kuralını ezberlemiştir. Kural eğitim verisinde %100
doğru olduğu için hiçbir ceza almamıştır. Yabancı merkez olarak
karşısına çıktığında elinde kullanacak morfolojik bilgi kalmamıştır.

### 7.2 · TSS 55 — kohortun %10'u

91 slayt (kohortun yaklaşık %10'u), site bölmesinde %60 hata.
Tek başına ortalama doğruluğu aşağı çeken en büyük kalem.

### 7.3 · TSS 38 — karşı örnek

İki bölmede de %45. Ezberlenmemiş; muhtemelen gerçekten zor slaytlar.
n=11 olduğu için yorum temkinli yapılmalıdır. Bu örnek, etkinin her
merkezde aynı olmadığını gösterir.

### 7.4 · Neden ortalama metrik yetersiz

Ortalama AUC farkı −0,049 olarak görünür ve küçük sanılabilir. Aynı
veri merkez bazında incelendiğinde tablo değişir: tek bir merkezde hata
oranı %5'ten %58'e çıkmaktadır.

**Tezde yalnızca ortalama metrik raporlamak, bu tabloyu gizlemek
anlamına gelir.**

> **Yapılacak:** Bu analiz tüm merkezler için üretilecek, dağılım
> (ör. hata oranı histogramı) tezde verilecektir. Tam embedding
> (2000 patch) ile tekrarlanacaktır.

---

## 8 · "MODELİN MORFOLOJİ ÖĞRENMESİ" NEDEN DOĞRUDAN SAĞLANAMIYOR

Modele "merkez bilgisini kullanma, morfolojiyi kullan" demek bu veri
setinde mümkün değildir.

Ortak merkez sayısı sıfır olduğundan, eğitim verisinde
*"bu slayt TSS 66'dan gelmektedir"* ile *"bu slayt LUSC'tir"*
ifadeleri **birebir aynı bilgidir.** Merkez bilgisini çıkaran herhangi
bir yöntem, sınıf bilgisini de çıkarır.

### 8.1 · Domain-adversarial eğitim neden çalışmaz

Gradyan tersine çevirme ile merkez bilgisinin bastırılması, literatürde
merkezlerin **karışık sınıflı** olduğu veri setlerinde işe yarar: model
"merkezi unut, sınıfı tut" diyebilir.

Bu veri setinde böyle bir ayrım yoktur. Adversarial kayıp, sınıflandırma
kaybıyla doğrudan çelişir ve sınıflandırıcıyı çökertir.

### 8.2 · Uygulanabilir yollar

| Yöntem | Ne yapar | Bu çalışmada |
|---|---|---|
| Domain-adversarial | Merkez bilgisini bastırır | ❌ Sınıfı da siler |
| Stain normalization | Rengi hizalar | ◐ Kısmi (%90→%80) |
| Ağır augmentation | Kestirme yolu güvenilmez kılar | ⬜ Denenecek |
| Gri tonlama ablasyonu | Rengi tamamen çıkarır | ⬜ Denenecek |
| Harici kohort (CPTAC) | Konfoundu baypas eder | ⬜ Takvime bağlı |
| Ölçüp raporlamak | Çözmez, dürüstleştirir | ✅ Zorunlu |

### 8.3 · Benimsenen duruş

Bu çalışmanın iddiası **batch effect'i çözmek değildir.** İddia, üç
foundation model'in merkez etkisi kontrol altında karşılaştırılmasıdır.

Literatürün büyük bölümü TCGA LUAD/LUSC üzerinde çalışırken bu kontrolü
yapmamaktadır. Her iki bölmenin de raporlanması ve merkez bazında hata
tablosunun verilmesi, bu çalışmanın metodolojik katkısıdır.

### 8.4 · Planlanan ek deneyler

**Augmentation ablasyonu** — eğitim sırasında renk sarsması, hafif
bulanıklık, JPEG gürültüsü. Site bölmesindeki skoru iyileştirirse bu
başlı başına bir bulgudur.

**Gri tonlama ablasyonu** — patch'ler griye çevrilip aynı deney
koşulur. Performans korunursa morfoloji taşıyor, çökerse renk
yapıyordu demektir. Düşük maliyetli, yüksek bilgi değerli.