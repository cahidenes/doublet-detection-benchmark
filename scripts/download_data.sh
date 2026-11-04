#!/usr/bin/env bash
# Download Xi & Li (2021) benchmark data from Zenodo 4062232.
#   real_datasets.zip        ~750 MB  (16 real datasets with experimental doublet labels)
#   simulation_datasets.zip  ~2.7 GB  (synthetic; only needed for stretch goals)
# Usage: scripts/download_data.sh [--with-sim]
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
[[ "${1:-}" == "--with-sim" ]] && fetch simulation_datasets.zip
echo "Done. Contents:"; find "$RAW" -maxdepth 3 -name '*.rds' | head -40
