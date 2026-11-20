"""Xi & Li (2021) per-dataset AUPRC/AUROC on the 16 real datasets, for analysis/xili_comparison.py.

Source: Zenodo 4444303, result.xlsx, sheet real_data (downloaded by scripts/download_data.sh to
data/raw/xili/). Needs openpyxl (not in the pinned envs: pip install --target <dir> openpyxl and
put <dir> on PYTHONPATH). The output is committed, so this only has to be rerun to check it.

Usage: python scripts/prepare_xili_reference.py
Writes: results/reference/xili2021_real.tsv
"""
from pathlib import Path

import pandas as pd

x = pd.read_excel("data/raw/xili/result.xlsx", sheet_name="real_data")
x = x[["dataset", "method", "auprc", "auroc"]].dropna(subset=["dataset", "method"])
Path("results/reference").mkdir(parents=True, exist_ok=True)
x.to_csv("results/reference/xili2021_real.tsv", sep="\t", index=False)
print(len(x), "rows:", x.method.value_counts().to_dict())
