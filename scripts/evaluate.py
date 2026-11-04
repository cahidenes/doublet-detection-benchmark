"""Compute accuracy metrics for one tool output against labels.
Usage: python scripts/evaluate.py <pred.tsv> <labels.tsv> <out.tsv>
Metrics: AUPRC, AUROC; precision/recall/TNR at 10/20/40% identification rates;
precision/recall/F1 and estimated-rate error at the tool's default call."""
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

pred = pd.read_csv(sys.argv[1], sep="\t")
if "status" in pred.columns:  # tool crashed: pass the failure through to metrics.tsv
    pred.to_csv(sys.argv[3], sep="\t", index=False)
    sys.exit(0)
lab = pd.read_csv(sys.argv[2], sep="\t")
df = pred.merge(lab, on="barcode", validate="one_to_one")
y, s = df["is_doublet"].to_numpy(), df["score"].to_numpy()


def prt(y, c):
    tp = int(((c == 1) & (y == 1)).sum()); fp = int(((c == 1) & (y == 0)).sum())
    fn = int(((c == 0) & (y == 1)).sum()); tn = int(((c == 0) & (y == 0)).sum())
    p = tp / (tp + fp) if tp + fp else np.nan
    r = tp / (tp + fn) if tp + fn else np.nan
    f1 = 2 * p * r / (p + r) if p + r else np.nan
    return p, r, tn / (tn + fp) if tn + fp else np.nan, f1


row = {k: df[k].iloc[0] for k in ("tool", "dataset", "seed", "expected_rate")}
row["status"] = "ok"
row.update(n=len(df), true_rate=y.mean(), auprc=average_precision_score(y, s), auroc=roc_auc_score(y, s))
order = np.argsort(-s, kind="stable")
for ir in (0.1, 0.2, 0.4):
    c = np.zeros_like(y); c[order[: int(round(ir * len(y)))]] = 1
    p, r, tnr, _ = prt(y, c)
    row.update({f"precision@{ir:.0%}": p, f"recall@{ir:.0%}": r, f"tnr@{ir:.0%}": tnr})
if df["call"].notna().any():
    c = df["call"].fillna(0).astype(int).to_numpy()
    p, r, tnr, f1 = prt(y, c)
    row.update(precision_default=p, recall_default=r, f1_default=f1,
               n_removed=int(c.sum()), est_rate=c.mean(), rate_error=c.mean() - y.mean())
pd.DataFrame([row]).to_csv(sys.argv[3], sep="\t", index=False)
