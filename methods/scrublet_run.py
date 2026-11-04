"""Scrublet with default parameters except the expected doublet rate."""
import numpy as np
import scrublet as scr
from common import load, parse_args, resolve_rate, write

args = parse_args()
X, barcodes, _, labels = load(args.input)
rate = resolve_rate(args, labels)
np.random.seed(args.seed)
s = scr.Scrublet(X, expected_doublet_rate=rate, random_state=args.seed)
score, call = s.scrub_doublets(verbose=False)
# Automatic threshold can fail (no bimodality) -> predicted_doublets is None
write(args, barcodes, score, None if call is None else call.astype(int), rate)
