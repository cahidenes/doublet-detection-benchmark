"""Simple ensembles: rank-average the scores of 2-3 tools (seed 0, labeled rate).

Each tool's score is turned into its within-dataset rank (0-1), and the ranks are averaged.
- Fixed ensembles chosen a priori from the top of the ranking and from different families.
- Leave-one-dataset-out (LODO) selection: for each held-out dataset, pick the 2- or 3-tool
  combination with the highest mean AUPRC on the other 15 datasets, and score it on the held-out
  one. This is the honest estimate; picking the best combination on all 16 is optimistic.

Usage: python analysis/ensemble.py
Writes: results/stats/ensemble.tsv, results/stats/ensemble_lodo.tsv
"""
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import rankdata, wilcoxon
from sklearn.metrics import average_precision_score

from style import BASELINES, NAMES, ROOT, STATS, main_metrics

m = main_metrics()
TOOLS = sorted(t for t in m.tool.unique() if t not in BASELINES)
DATASETS = sorted(m.dataset.unique())

ranks, y = {}, {}
for ds in DATASETS:
    bc = (ROOT / f"data/processed/{ds}/barcodes.tsv").read_text().split()
    lab = pd.read_csv(ROOT / f"data/processed/{ds}/labels.tsv", sep="\t").set_index("barcode")
    y[ds] = lab.loc[bc, "is_doublet"].to_numpy()
    for t in TOOLS:
        s = pd.read_csv(ROOT / f"results/pred/{t}/{ds}/s0_r1.0.tsv", sep="\t").set_index("barcode").loc[bc, "score"]
        ranks[t, ds] = rankdata(s.to_numpy()) / len(s)

combos = [c for k in (1, 2, 3) for c in combinations(TOOLS, k)]
auprc = pd.DataFrame({ds: {"+".join(c): average_precision_score(y[ds], np.mean([ranks[t, ds] for t in c], 0))
                           for c in combos} for ds in DATASETS})       # combos x datasets

best_single = auprc.loc[TOOLS].mean(1).idxmax()
fixed = ["scdblfinder+vaeda", "scdblfinder+doubletdetection", "doubletdetection+scdblfinder+vaeda",
         "cxds+scdblfinder+vaeda", "bcds+cxds+doubletfinder"]  # last = Chord's components
fixed = ["+".join(sorted(f.split("+"))) for f in fixed]
rows = []
for c in [best_single] + fixed:
    d = auprc.loc[c] - auprc.loc[best_single]
    rows.append({"ensemble": c, "mean_auprc": auprc.loc[c].mean(), "median_auprc": auprc.loc[c].median(),
                 "wins_vs_best_single": int((d > 0).sum()), "losses": int((d < 0).sum()),
                 "p_vs_best_single": wilcoxon(d).pvalue if (d != 0).any() else 1.0})
rows.append({"ensemble": "chord (GBM on its components)",
             "mean_auprc": m[m.tool == "chord"].set_index("dataset").auprc.loc[DATASETS].mean(),
             "median_auprc": m[m.tool == "chord"].auprc.median()})
top_all = auprc.mean(1).sort_values(ascending=False)
rows.append({"ensemble": f"best on all 16 (optimistic): {top_all.index[0]}",
             "mean_auprc": top_all.iloc[0], "median_auprc": auprc.loc[top_all.index[0]].median()})

lodo = []
for k in (1, 2, 3):
    cand = [c for c in auprc.index if c.count("+") == k - 1]
    for ds in DATASETS:
        pick = auprc.loc[cand, [x for x in DATASETS if x != ds]].mean(1).idxmax()
        lodo.append({"k": k, "held_out": ds, "picked": pick, "auprc": auprc.loc[pick, ds],
                     "best_single_auprc": auprc.loc[best_single, ds]})
lodo = pd.DataFrame(lodo)
for k, g in lodo.groupby("k"):
    d = g.auprc - g.best_single_auprc
    rows.append({"ensemble": f"LODO-selected, {k} tool(s)", "mean_auprc": g.auprc.mean(),
                 "median_auprc": g.auprc.median(), "wins_vs_best_single": int((d > 0).sum()),
                 "losses": int((d < 0).sum()),
                 "p_vs_best_single": wilcoxon(d).pvalue if (d != 0).any() else 1.0,
                 "picked": "; ".join(f"{c} ({n})" for c, n in g.picked.value_counts().items())})

out = pd.DataFrame(rows)
STATS.mkdir(parents=True, exist_ok=True)
out.to_csv(STATS / "ensemble.tsv", sep="\t", index=False, float_format="%.4f")
lodo.to_csv(STATS / "ensemble_lodo.tsv", sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 90)
print(out.round(4).to_string(index=False))
print("\nTop 8 combinations on all 16 datasets:\n", top_all.head(8).round(4).to_string())
