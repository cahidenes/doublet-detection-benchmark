"""Trivial baselines (Xi & Li 2021): library size and number of expressed genes.
No default call; evaluation uses scores (AUPRC/AUROC) and top-k identification rates."""
import numpy as np
from common import load, parse_args, resolve_rate, write

args = parse_args()
X, barcodes, _, labels = load(args.input)
if args.tool == "libsize":
    score = np.asarray(X.sum(axis=1)).ravel()
elif args.tool == "ngenes":
    score = np.diff(X.indptr)  # nonzeros per row
else:
    raise SystemExit(f"unknown baseline {args.tool}")
write(args, barcodes, score, None, resolve_rate(args, labels))
