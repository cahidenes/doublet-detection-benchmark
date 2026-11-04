"""DoubletDetection 4.3 with package defaults (BoostClassifier(), predict()).
score = doublet_score(): mean -log p-value over iterations (higher = more doublet-like).
Does not take an expected doublet rate (excluded from the rate-robustness experiment)."""
import numpy as np
import doubletdetection
from common import load, parse_args, resolve_rate, write

args = parse_args()
X, barcodes, _, labels = load(args.input)
clf = doubletdetection.BoostClassifier(random_state=args.seed)
call = clf.fit(X).predict()  # 1 doublet, 0 singlet, NaN if undetermined
score = clf.doublet_score()
call = np.where(np.isnan(call), 0, call).astype(int)
write(args, barcodes, score, call, resolve_rate(args, labels))
