"""Shared I/O for Python tool wrappers.

Every wrapper (R or Python) writes the SAME output schema, one row per droplet:

    barcode  score  call  tool  dataset  seed  expected_rate

- score: higher = more doublet-like (flip the sign if a tool reports the opposite)
- call:  1/0 from the tool's default thresholding; NA if the tool has none
Runtime and peak memory come from Snakemake's `benchmark:` file, not from here.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def parse_args(extra=None):
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="processed dataset dir")
    p.add_argument("--output", required=True, help="output .tsv")
    p.add_argument("--dataset", required=True)
    p.add_argument("--tool", required=True)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--expected-rate", type=float, default=None,
                   help="expected doublet rate; default = labeled rate (the '1x' setting)")
    for a, kw in (extra or {}).items():
        p.add_argument(a, **kw)
    return p.parse_args()


def load(input_dir):
    """Return (X droplets x genes CSR int matrix, barcodes, genes, labels Series)."""
    d = Path(input_dir)
    X = sp.csr_matrix(scipy.io.mmread(d / "matrix.mtx.gz").T)
    barcodes = (d / "barcodes.tsv").read_text().split()
    genes = (d / "genes.tsv").read_text().split()
    labels = pd.read_csv(d / "labels.tsv", sep="\t").set_index("barcode")["is_doublet"]
    return X, barcodes, genes, labels


def resolve_rate(args, labels):
    return args.expected_rate if args.expected_rate is not None else float(labels.mean())


def write(args, barcodes, score, call, rate):
    out = pd.DataFrame({
        "barcode": barcodes,
        "score": np.asarray(score, dtype=float),
        "call": pd.array(call if call is not None else [pd.NA] * len(barcodes), dtype="Int64"),
        "tool": args.tool, "dataset": args.dataset,
        "seed": args.seed, "expected_rate": rate,
    })
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)
