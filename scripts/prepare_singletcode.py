"""singletCode (Zhang et al. 2024): independent validation data.

The Watermelon-barcoded T47D datasets from Zhang et al. (2024) come with the singlets that
singletCode identified from lineage barcodes (singlets_all.txt). These cells were never used to
develop or tune any of the benchmarked tools. For each of the 6 samples this script writes a
singlets-only base dataset, data/processed/wm-<sample>/ (all is_doublet = 0); doublets are then
made with scripts/simulate_doublets.py (random pairing), as in Zhang et al.

Source: figshare 25478680 (singletCodeData.zip -> watermelon.zip), unpacked by
scripts/download_data.sh --with-singletcode into data/raw/singletcode/wm/watermelon/.

Usage: python scripts/prepare_singletcode.py
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io

SRC = Path("data/raw/singletcode/wm/watermelon")
for d in sorted(SRC.glob("T47D-*")):
    X = scipy.io.mmread(d / "matrix.mtx.gz").tocsc()                # genes x droplets
    bc = [b.split("-")[0] for b in gzip.open(d / "barcodes.tsv.gz", "rt").read().split()]
    genes = [l.split("\t")[1] for l in gzip.open(d / "features.tsv.gz", "rt").read().strip().split("\n")]
    genes = pd.Index(genes)
    if genes.has_duplicates:                                        # make symbols unique
        dup = genes.duplicated(keep=False)
        genes = pd.Index([f"{g}.{i}" if dup[i] else g for i, g in enumerate(genes)])
    singlets = set((d / "singlets_all.txt").read_text().split())
    keep = np.array([b in singlets for b in bc])
    X = X[:, np.flatnonzero(keep)]
    gk = np.asarray((X > 0).sum(1)).ravel() >= 3                    # genes in >= 3 cells
    X = X[gk]
    out = Path("data/processed") / f"wm-{d.name}"
    out.mkdir(parents=True, exist_ok=True)
    with gzip.open(out / "matrix.mtx.gz", "wb") as fh:
        scipy.io.mmwrite(fh, X.tocoo())
    (out / "genes.tsv").write_text("\n".join(genes[gk]) + "\n")
    names = [b for b, k in zip(bc, keep) if k]
    (out / "barcodes.tsv").write_text("\n".join(names) + "\n")
    pd.DataFrame({"barcode": names, "is_doublet": 0}).to_csv(out / "labels.tsv", sep="\t", index=False)
    print(f"{out.name}: {len(bc)} droplets, {keep.sum()} singletCode singlets, {X.shape[0]} genes")
