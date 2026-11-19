"""Realistic expected-rate setting, and rate robustness on simulated doublets.

1. The main run gives rate-taking tools (Scrublet, scDblFinder, DoubletFinder, Chord) the
   labeled rate, which no user knows. Here they get the 10x loading guideline instead,
   0.8% per 1,000 recovered droplets (results/heuristic_metrics.tsv). The default-call
   comparison (F1, droplets removed) is redone with these runs; tools without a rate input
   are unchanged.
2. Rate robustness (0.5x / 1x / 2x) on the 16 seed-0 simulated datasets, where the labels are
   complete (results/sim_robustness_metrics.tsv).

Usage: python analysis/rate_settings.py
Writes: results/stats/heuristic_rate.tsv, results/stats/f1_realistic.tsv,
        results/stats/sim_robustness.tsv
"""
import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare

from style import NAMES, ROOT, STATS, main_metrics

RATE_TOOLS = ["scrublet", "scdblfinder", "doubletfinder", "chord"]
m = main_metrics()
h = pd.read_csv(ROOT / "results/heuristic_metrics.tsv", sep="\t")
print("heuristic runs:", h.status.value_counts().to_dict())
if "error" in h:
    print(h[h.status != "ok"][["tool", "dataset", "error"]].to_string(index=False))
h = h[h.status == "ok"].copy()
# F1 is undefined (0/0) when precision and recall are both 0: e.g. Chord calls 2 droplets on
# J293t-dm (heuristic rate 0.004), both singlets. That is a complete miss, so F1 = 0.
h.loc[(h.precision_default == 0) & (h.recall_default == 0), "f1_default"] = 0.0
for x in (m, h):
    x["removed_over_true"] = x.n_removed / (x.true_rate * x.n)

# 1a. per rate tool: labeled vs heuristic
rows = []
for t in RATE_TOOLS:
    for name, x in (("labeled", m), ("heuristic", h)):
        g = x[x.tool == t]
        rows.append({"tool": t, "setting": name, "n_ok": len(g), "median_auprc": g.auprc.median(),
                     "median_precision": g.precision_default.median(),
                     "median_recall": g.recall_default.median(), "median_f1": g.f1_default.median(),
                     "median_removed_over_true": g.removed_over_true.median()})
ht = pd.DataFrame(rows)
rate_cmp = (h.drop_duplicates("dataset").set_index("dataset")[["expected_rate", "true_rate"]]
            .sort_values("true_rate"))
print("\nheuristic vs labeled rate per dataset:\n", rate_cmp.round(3).to_string())
ht.to_csv(STATS / "heuristic_rate.tsv", sep="\t", index=False, float_format="%.4f")
print("\n", ht.round(3).to_string(index=False))

# 1b. F1 ranking of all callers under the realistic setting
callers = [t for t in m.tool.unique() if t not in ("libsize", "ngenes", "dbldensity")]
real = pd.concat([m[m.tool.isin(callers) & ~m.tool.isin(RATE_TOOLS)], h])
f1 = real.pivot(index="dataset", columns="tool", values="f1_default")
f1 = f1.dropna()                                            # paired: datasets where every tool ran
stat, p = friedmanchisquare(*[f1[c] for c in f1.columns])
out = pd.DataFrame({"median_f1_realistic": f1.median(),
                    "mean_rank_f1_realistic": f1.rank(axis=1, ascending=False).mean(),
                    "median_f1_labeled": m[m.tool.isin(callers)].groupby("tool").f1_default.median()})
out = out.sort_values("mean_rank_f1_realistic")
out.attrs = {}
out.index.name = "tool"
out.to_csv(STATS / "f1_realistic.tsv", sep="\t", float_format="%.4f")
print(f"\nF1 with realistic rates ({len(f1)} datasets): Friedman chi2 = {stat:.1f}, p = {p:.2g}")
print(out.rename(index=NAMES).round(3).to_string())

# 2. rate robustness on simulated data
sr = ROOT / "results/sim_robustness_metrics.tsv"
if sr.exists():
    s = pd.read_csv(sr, sep="\t")
    s = s[s.status == "ok"].copy()
    s["removed_over_true"] = s.n_removed / (s.true_rate * s.n)
    g = s.groupby(["tool", "rate_mult"])
    rob = g[["precision_default", "recall_default", "f1_default", "removed_over_true", "auprc"]].median()
    rob["n"] = g.size()
    rob.to_csv(STATS / "sim_robustness.tsv", sep="\t", float_format="%.4f")
    print("\nRate robustness on simulated doublets (medians over 16 datasets):\n", rob.round(3).to_string())
