"""Stretch goal: stability across random seeds (real data, labeled rate).

Seed 0 is the main run (results/metrics.tsv); seeds 1 and 2 come from the stability run
(results/stability_metrics.tsv) for every tool with randomness (cxds and the baselines are
deterministic). Per tool and dataset:
  auprc_sd       standard deviation of AUPRC over the 3 seeds
  call_jaccard   mean pairwise Jaccard of the called droplets between seeds (tools with calls)
  rank changes   does the per-seed ranking (mean AUPRC) of the tools change?
DoubletFinder comes out identical across seeds: Seurat's RunPCA calls set.seed(42) internally,
which overrides the wrapper's seed before DoubletFinder draws its artificial doublets.

Usage: python analysis/stability.py   (call_jaccard needs the local results/pred/ files)
Writes: results/stats/stability.tsv, results/stats/stability_ranks.tsv,
        results/figures/stability.png
"""
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from style import CAT, FIG, INK2, NAMES, ROOT, STATS, SURF, clean, main_metrics, setup

m0 = main_metrics()
st = pd.read_csv(ROOT / "results/stability_metrics.tsv", sep="\t")
st = st[(st.status == "ok") & (st.rate_mult == 1.0)]
tools = sorted(st.tool.unique())
allm = pd.concat([m0[m0.tool.isin(tools)], st])
allm["seed"] = allm.seed.astype(int)
# keep only (tool, dataset) pairs present for all 3 seeds, so per-seed means are comparable
allm = allm.groupby(["tool", "dataset"]).filter(lambda g: g.seed.nunique() == 3)

rows = []
for (t, ds), g in allm.groupby(["tool", "dataset"]):
    r = {"tool": t, "dataset": ds, "n_seeds": len(g), "auprc_mean": g.auprc.mean(),
         "auprc_sd": g.auprc.std(ddof=1), "auprc_range": g.auprc.max() - g.auprc.min(),
         "f1_sd": g.f1_default.std(ddof=1)}
    calls = {}
    for s in sorted(g.seed):
        p = pd.read_csv(ROOT / f"results/pred/{t}/{ds}/s{s}_r1.0.tsv", sep="\t")
        if p["call"].notna().any():
            calls[s] = set(p.barcode[p["call"] == 1])
    if len(calls) >= 2:
        js = [len(calls[a] & calls[b]) / max(1, len(calls[a] | calls[b]))
              for a, b in combinations(calls, 2)]
        r["call_jaccard"] = float(np.mean(js))
    rows.append(r)
r = pd.DataFrame(rows)
STATS.mkdir(parents=True, exist_ok=True)
r.to_csv(STATS / "stability.tsv", sep="\t", index=False, float_format="%.4f")
summ = r.groupby("tool")[["auprc_sd", "auprc_range", "f1_sd", "call_jaccard"]].median() \
        .sort_values("auprc_sd")
print(summ.round(4).to_string())

# ranking per seed (deterministic tools carry their seed-0 value)
det = m0[~m0.tool.isin(tools)]
per_seed = {}
for s in (0, 1, 2):
    x = pd.concat([allm[allm.seed == s], det])
    per_seed[f"seed{s}"] = x.groupby("tool").auprc.mean()
ps = pd.DataFrame(per_seed)
ranks = ps.rank(ascending=False)
ps.join(ranks.add_suffix("_rank")).sort_values("seed0_rank") \
  .to_csv(STATS / "stability_ranks.tsv", sep="\t", float_format="%.4f")
print(ps.join(ranks.add_suffix("_rank")).sort_values("seed0_rank").round(3).to_string())

# ---------------------------------------------------------------- figure
setup()
order = list(summ.index)[::-1]
fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)
for ax, (col, title) in zip(axes, [("auprc_sd", "AUPRC standard deviation over 3 seeds"),
                                    ("call_jaccard", "Overlap of called doublets between seeds (Jaccard)")]):
    for i, t in enumerate(order):
        v = r[r.tool == t][col].dropna()
        if v.empty:
            ax.text(0.02, i, "no default call", transform=ax.get_yaxis_transform(), va="center",
                    fontsize=8, color=INK2)
            continue
        ax.scatter(v, np.full(len(v), i) + np.linspace(-.2, .2, len(v)), s=10, color="#b9b8b2", zorder=2)
        ax.scatter([v.median()], [i], s=60, color=CAT[0], edgecolor=SURF, lw=1.2, zorder=3)
    ax.set_title(title, loc="left", fontsize=10, pad=8)
    clean(ax, grid="x")
axes[0].set_yticks(range(len(order)), [NAMES[t] for t in order])
axes[1].set_xlim(0, 1.02)
fig.suptitle("Stability across random seeds (16 datasets; blue = median, grey = datasets)",
             x=0.01, ha="left", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=3)
fig.savefig(FIG / "stability.png", dpi=200)
