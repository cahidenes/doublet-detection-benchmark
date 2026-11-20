#!/usr/bin/env bash
# Download Xi & Li (2021) benchmark data from Zenodo 4062232.
#   real_datasets.zip        ~750 MB  (16 real datasets with experimental doublet labels)
#   simulation_datasets.zip  ~2.7 GB  (synthetic; only needed for stretch goals)
# Optional:
#   --with-sim          Xi & Li synthetic datasets (not used by the report)
#   --with-singletcode  singletCode data (Zhang et al. 2024, figshare 25478680, 5.5 GB zip):
#                       Watermelon T47D matrices + singlet lists -> data/raw/singletcode/wm/,
#                       benchmark plot data -> data/raw/singletcode/plotData/.
#                       Then: python scripts/prepare_singletcode.py, and simulate_doublets.py per
#                       sample (see scripts/run_final.sh). The 5.5 GB zip can be deleted afterwards.
# Always: Xi & Li per-dataset results (Zenodo 4444303, result.xlsx) -> data/raw/xili/, used by
#         analysis/xili_comparison.py via results/reference/xili2021_real.tsv.
# Usage: scripts/download_data.sh [--with-sim] [--with-singletcode]
set -euo pipefail
cd "$(dirname "$0")/.."
RAW=data/raw
BASE=https://zenodo.org/api/records/4062232/files
mkdir -p "$RAW"

fetch() {
  local f=$1
  if [[ -f "$RAW/$f" ]]; then echo "skip $f (exists)"; else
    curl -L --fail --retry 3 -C - -o "$RAW/$f.part" "$BASE/$f/content" && mv "$RAW/$f.part" "$RAW/$f"
  fi
  unzip -n -q "$RAW/$f" -d "$RAW/${f%.zip}"
}

fetch real_datasets.zip
for a in "$@"; do [[ $a == --with-sim ]] && fetch simulation_datasets.zip; done

mkdir -p "$RAW/xili"
[[ -f $RAW/xili/result.xlsx ]] || curl -L --fail --retry 3 -o "$RAW/xili/result.xlsx" \
    https://zenodo.org/api/records/4444303/files/result.xlsx/content

if [[ " $* " == *" --with-singletcode "* ]]; then
  SC=$RAW/singletcode; mkdir -p "$SC"
  if [[ ! -f $SC/singletCodeData.zip ]]; then
    curl -L --fail --retry 3 -C - -o "$SC/singletCodeData.zip.part" https://ndownloader.figshare.com/files/47877253
    echo "1e8251a4f8a7d8b6029d01e0e2a884a0  $SC/singletCodeData.zip.part" | md5sum -c -
    mv "$SC/singletCodeData.zip.part" "$SC/singletCodeData.zip"
  fi
  unzip -n -q "$SC/singletCodeData.zip" watermelon.zip 'plotData/benchmarking/*' -x '__MACOSX/*' -d "$SC"
  unzip -n -q "$SC/watermelon.zip" -x '__MACOSX/*' -d "$SC/wm" && rm "$SC/watermelon.zip"
fi
echo "Done. Contents:"; find "$RAW" -maxdepth 3 -name '*.rds' | head -40
