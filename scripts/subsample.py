"""Nested, stratified droplet subsamples of one processed dataset (scaling experiment, metric 9).

Each size keeps the dataset's doublet rate, and smaller subsets are contained in larger ones
(one random permutation per class, then take prefixes), so the size is the only thing that varies.
Writes data/processed/scale-<ds>-n<size>/ in the same format as scripts/convert_rds.R.

Usage: python scripts/subsample.py <ds> <size> [<size> ...] [--seed 0]
"""
import argparse
import gzip
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io

p = argparse.ArgumentParser()
p.add_argument("ds"); p.add_argument("sizes", nargs="+", type=int); p.add_argument("--seed", type=int, default=0)
a = p.parse_args()

src = Path("data/processed") / a.ds
m = scipy.io.mmread(src / "matrix.mtx.gz").tocsc()          # genes x droplets
lab = pd.read_csv(src / "labels.tsv", sep="\t")
rng = np.random.default_rng(a.seed)
perm = {c: rng.permutation(np.flatnonzero(lab.is_doublet.values == c)) for c in (0, 1)}
rate = lab.is_doublet.mean()

for n in a.sizes:
    n = min(n, len(lab))
    k1 = round(n * rate)
    idx = np.sort(np.concatenate([perm[1][:k1], perm[0][:n - k1]]))
    out = Path("data/processed") / f"scale-{a.ds}-n{n}"
    out.mkdir(parents=True, exist_ok=True)
    with gzip.open(out / "matrix.mtx.gz", "wb") as fh:
        scipy.io.mmwrite(fh, m[:, idx])
    shutil.copy(src / "genes.tsv", out / "genes.tsv")
    lab.iloc[idx].to_csv(out / "labels.tsv", sep="\t", index=False)
    (out / "barcodes.tsv").write_text("\n".join(lab.barcode.iloc[idx]) + "\n")
    print(f"{out.name}: {n} droplets, {k1} doublets")
