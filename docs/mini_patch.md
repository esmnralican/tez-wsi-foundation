# MİNİ PİLOT — DINOv2 / 100 patch per WSI

Tarih: 7 Ekim 2026
Amaç: Embedding boru hattını küçük veriyle kurmak ve doğrulamak.

---

## 1 · NEDEN PİLOT

Danışman önerisi: tüm veriyi bir anda eğitime sokma. Önce az veriyle
çalış, küçük batch ile test et, sonra büyüt.

Bu aşamanın amacı **sonuç almak değil**, boru hattının çalıştığını
görmektir:

- Embedding'ler doğru boyutta mı çıkıyor
- Patch → vektör eşlemesi kayıt altına alınıyor mu
- Donanım sınırları ne (VRAM, hız, disk)
- Manifest sırası örneklemeyi yanlı hale getiriyor mu

100 instance'lık bag'de MIL attention anlamlı bir şey öğrenemez;
bu aşamadan çıkacak doğruluk değeri yorumlanmayacaktır.

---

## 2 · YAPILAN

| | WSI | Patch | Süre |
|---|---:|---:|---:|
| LUAD | 466 | 46.578 | ~14 dk |
| LUSC | 478 | 47.800 | 18,3 dk |
| **Toplam** | **944** | **94.378** | **~33 dk** |

Model: **DINOv2 ViT-L/14**, donuk (frozen) backbone, CLS token.
Çıktı: slayt başına bir `.npy`, `(100, 1024)`, fp16.

Donanım: RTX 3050 Ti Laptop, 4 GB VRAM. Ölçülen hız **~57 patch/sn**
(xFormers etkin).

Değer aralığı doğrulandı: LUAD `[-6,53 , 5,05]`, LUSC `[-7,27 , 5,23]`.

---

## 3 · PATCH SEÇİMİ

`embed.py` manifest'teki **ilk N satırı** alır, ayrıca rastgeleleştirme
yapmaz. Bu, manifest sırasının rastgele olmasına dayanır.

**Doğrulandı** — `TCGA-05-4244` manifest sırası:

```
x28920_y18798
x27474_y13496
x14942_y11086
x29884_y17352
...
```

Koordinatlar dağınık. Manifest, çıkarma aşamasındaki seed=42 (blake2b)
seçim sırasını koruyor. Dolayısıyla ilk 100 satır slaydın her yerinden
gelen deterministik bir rastgele örnektir.

**Sonucu:** `--limit 500` ile çalıştırıldığında ilk 100'ü aynı
patch'ler olur. Alt küme ablasyonunda tek değişken patch sayısı kalır,
örnekleme rastgeleliği karışmaz.

> Eğer manifest koordinat sırasında olsaydı ilk 100 patch slaydın
> sol üst köşesinden gelecek ve örneklem temsil edici olmayacaktı.
> Bu kontrol yapılmadan devam edilmedi.

---

## 4 · ÇIKTI YAPISI

```
out/<LABEL>/embeddings/dinov2_n100/
├── TCGA-05-4244-01Z-00-DX1.npy      (100, 1024) fp16
├── TCGA-05-4244-01Z-00-DX1.json     patch sırası + hasta kimliği
├── ...
└── _log.csv                          slayt, süre, hız
```

`.json` dosyası her satırın hangi patch'e karşılık geldiğini saklar.
Bu sayede alt küme ablasyonu **embedding yeniden çıkarmadan** yapılabilir:
`.npy`'den ilk N satır alınır, karşılığı json'dan okunur.

Klasör adında `_n100` eki var; tam çıkarma `dinov2/` klasörüne gidecek,
karışmayacak.

Disk: 100 × 1024 × 2 bayt ≈ 200 KB/slayt → 944 slayt ≈ **190 MB**.

---

## 5 · KURULUM SORUNLARI VE ÇÖZÜMLERİ

### 5.1 · NumPy sürüm çakışması

PyTorch NumPy 1.x'e göre derlenmişti, ortamda NumPy 2.2.6 vardı.

```
RuntimeError: Numpy is not available
UserWarning: Failed to initialize NumPy: _ARRAY_API not found
```

Model yükleniyor, veri yükleyici ilk batch'te çöküyordu.

**Çözüm:** `pip install "numpy<2"` → 1.26.4

> Yan etki: `opencv-python-headless 5.0` NumPy ≥ 2 istiyor.
> Patch çıkarma ve ileride Macenko için `opencv-python-headless<5`
> gerekebilir.

### 5.2 · WSL'de libcuda bulunamıyor

```
Could not load library libcudnn_cnn_infer.so.8.
Error: libcuda.so: cannot open shared object file
Aborted (core dumped)
```

WSL'de CUDA sürücü kütüphaneleri `/usr/lib/wsl/lib` altındadır ve
varsayılan arama yoluna dahil değildir.

**Çözüm:**

```bash
echo 'export LD_LIBRARY_PATH=/usr/lib/wsl/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
```

### 5.3 · tmux oturumu açılmıyordu

`tmux new -s embed` komutu `duplicate session` döndürüyor, komut
tmux dışında çalışıyordu. Terminal kapanırsa iş ölürdü.

**Çözüm:** `tmux kill-server` ile temizlik, sonra yeni oturum.

---

## 6 · ÇALIŞTIRMA

```bash
conda activate esma_dinov2
cd ~/tez

tmux new -s embed
python embed.py --model dinov2 --label LUAD --limit 100 --batch 32
python embed.py --model dinov2 --label LUSC --limit 100 --batch 32
# Ctrl+B sonra D ile ayril
```

Doğrulama:

```bash
python3 -c "
import numpy as np, glob
for lab in ('LUAD','LUSC'):
    f = sorted(glob.glob(f'out/{lab}/embeddings/dinov2_n100/*.npy'))
    a = np.load(f[0])
    print(f'{lab}: {len(f)} slayt | {a.shape} {a.dtype}')
"
```

---

## 7 · ÖLÇÜLEN HIZ VE TAM ÇIKARMA TAHMİNİ

| Yapılandırma | Patch | Tahmini süre |
|---|---:|---:|
| Pilot (100/WSI) | 94.378 | ~33 dk `[ölçüm]` |
| Tam (2000/WSI) | 1.835.703 | ~9 sa `[tahmin]` |

Tahmin 57 patch/sn ölçümünden doğrusal ölçeklenmiştir; RTX 5000
(17 GB) üzerinde daha büyük batch ile daha kısa sürmesi beklenir.

UNI benzer boyutta (ViT-L, 1024), benzer sürede biter.
Virchow daha ağır (ViT-H, 632M parametre) ve 4 GB VRAM'de küçük batch
gerektirir.

---

## 8 · SONRAKİ ADIMLAR

1. `split.py` — hasta seviyesinde bölme, **TSS katmanlı**
   (merkez sınıf etiketini %100 belirlediği için rastgele bölme
   yanıltıcı sonuç verir)
2. `train_mil.py` — ABMIL, k-fold çapraz doğrulama, hasta seviyesinde
   değerlendirme
3. Tam embedding çıkarma (2000 patch/WSI)
4. Alt küme ablasyonu: 10 / 50 / 100 / 500 / 1000 / 2000 —
   embedding yeniden çıkarılmadan
5. UNI ve Virchow (HuggingFace erişim onayı bekleniyor)

---

## 9 · AÇIK KONULAR

| # | Konu | Durum |
|---|---|---|
| 1 | Bag boyutu dengesizliği (78–2000 patch) | Karar eğitim aşamasında |
| 2 | Native 20× / 40× grup farkı (LUAD %4, LUSC %12) | Değerlendirmede kontrol |
| 3 | Patch verisinin sunucuya aktarımı (~50 GB) | Manuel aktarım sürüyor |
| 4 | UNI / Virchow HuggingFace erişimi | Başvuru bekleniyor |
| 5 | Macenko ablation | Embedding sonrası |