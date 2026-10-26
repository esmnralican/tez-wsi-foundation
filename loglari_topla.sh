#!/usr/bin/env bash
# Bildirideki sayilarin kaynagi olan script ciktilarini tek klasorde
# toplar. Iki makinede ayri ayri calistirilir; her makine kendi
# bulabildigini uretir, bulamadigini ATLANDI diye isaretler.
#
# Depo kokunde (out/ klasorunun yaninda) calistirilir:
#     bash loglari_topla.sh
#
# Cikti: loglar/<tarih>_<makine>/ ve icinde ozet.txt

set -u

KOK="$(pwd)"
[ -d out ] || { echo "HATA: out/ klasoru yok. Depo kokunde calistir."; exit 1; }

MAKINE="$(hostname | tr -cd '[:alnum:]._-')"
DIZIN="loglar/$(date +%Y%m%d_%H%M)_${MAKINE}"
mkdir -p "$DIZIN"
OZET="$DIZIN/ozet.txt"

# ----------------------------------------------------------- ortam
{
  echo "tarih   : $(date -Iseconds)"
  echo "makine  : $(hostname)"
  echo "dizin   : $KOK"
  echo "python  : $(python -V 2>&1)"
  echo -n "git     : "
  if git rev-parse --git-dir > /dev/null 2>&1; then
      echo -n "$(git rev-parse --short HEAD)"
      # Yalnizca IZLENEN dosyalardaki degisiklige bakilir: kaynak kaydi
      # icin onemli olan kodun commit'ten farkli olup olmamasidir.
      # Izlenmeyen dosyalar (bu scriptin kendi loglar/ klasoru dahil)
      # sayilmaz.
      if [ -n "$(git status --porcelain --untracked-files=no 2>/dev/null)" ]; then
          echo "  (KIRLI -- izlenen dosyalarda kaydedilmemis degisiklik var)"
          echo "degisen dosyalar:"
          git status --porcelain --untracked-files=no | sed 's/^/  /'
      else
          echo "  (temiz -- kod commit ile ayni)"
      fi
  else
      echo "depo degil"
  fi
  echo -n "paketler: "
  python - <<'PY' 2>/dev/null || echo "okunamadi"
import numpy, scipy
print(f"numpy {numpy.__version__}  scipy {scipy.__version__}")
PY
  echo
  echo "veri sayimi:"
  echo "  fold_metrics.csv : $(find out/mil -name fold_metrics.csv 2>/dev/null | wc -l)"
  echo "  splits_*.csv     : $(ls out/splits_*.csv 2>/dev/null | wc -l)"
  echo "  dikkat kompakt   : $(ls out/dikkat_kompakt_*.csv 2>/dev/null | wc -l)"
  echo
} | tee "$OZET"

# ----------------------------------------------------------- kosucu
SAYAC_OK=0
SAYAC_HATA=0
SAYAC_ATLA=0

kos() {                       # kos <cikti_adi> <script> [arg...]
    local ad="$1"; shift
    local script="$1"
    local hedef="$DIZIN/${ad}.txt"

    if [ ! -f "$script" ]; then
        printf '  %-26s ATLANDI  (%s yok)\n' "$ad" "$script" | tee -a "$OZET"
        SAYAC_ATLA=$((SAYAC_ATLA + 1))
        return
    fi

    { echo "# komut : python $*"; echo "# tarih : $(date -Iseconds)"; echo; } > "$hedef"
    if python "$@" >> "$hedef" 2>&1; then
        printf '  %-26s TAMAM    (%s satir)\n' "$ad" "$(wc -l < "$hedef")" | tee -a "$OZET"
        SAYAC_OK=$((SAYAC_OK + 1))
    else
        local kod=$?
        printf '  %-26s HATA %-3s (bkz. %s.txt)\n' "$ad" "$kod" "$ad" | tee -a "$OZET"
        echo "# cikis kodu: $kod" >> "$hedef"
        SAYAC_HATA=$((SAYAC_HATA + 1))
    fi
}

echo "calistirilan kontroller:" | tee -a "$OZET"

# Istatistik dogrulama -- Tablo 1, Tablo 2, NB, Holm, fold boyutlari.
# Dokuz kosulun tamami bu makinede degilse [0] butunluk hatasi verip
# durur; o durumda yalnizca fold boyutlari ([6]) alinir. Hangi
# makinede oldugumuzu sormuyoruz, deniyoruz.
for m in site patient; do
    if python dogrula.py --mod "$m" > /dev/null 2>&1; then
        kos "dogrula_$m" dogrula.py --mod "$m"
    else
        printf '  %-26s tam kontrol gecersiz -> yalniz fold boyutlari\n' \
            "dogrula_$m" | tee -a "$OZET"
        kos "foldboyut_$m" dogrula.py --mod "$m" --yalniz-fold
    fi
done

# Eslestirilmis karsilastirmalar ve guven araliklari (§3.1, §3.3)
kos karsilastir      karsilastir.py
kos karsilastir_abl  karsilastir.py --abl dinov2

# Dikkat yogunlasmasi (§3.4) -- 2. asama, kompakt dosyalar gerekli
if [ "$(ls out/dikkat_kompakt_*.csv 2>/dev/null | wc -l)" -gt 0 ]; then
    kos dikkat dikkat_test.py
else
    printf '  %-26s ATLANDI  (kompakt dosya yok: once --etiket)\n' dikkat \
        | tee -a "$OZET"
    SAYAC_ATLA=$((SAYAC_ATLA + 1))
fi

# Merkez basina hata orani (§3.5)
for m in dinov2 uni virchow; do
    kos "merkez_${m}" merkez_hata.py --model "$m"
done

# ----------------------------------------------------------- kapanis
{
  echo
  echo "tamam: $SAYAC_OK   hata: $SAYAC_HATA   atlanan: $SAYAC_ATLA"
  echo
  echo "Atlananlar bu makinede o verinin olmadigini gosterir; digerinde"
  echo "ayni scripti calistirip iki klasoru birlestir."
} | tee -a "$OZET"

echo
echo "cikti: $DIZIN"
echo "arsiv icin:  tar czf ${DIZIN//\//_}.tar.gz $DIZIN"