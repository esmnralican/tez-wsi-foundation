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
morfolojik değil, merkez kaynaklıdır.
