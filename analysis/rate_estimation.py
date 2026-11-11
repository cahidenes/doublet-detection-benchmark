"""Metric 4: doublet-rate estimation from each tool's default calls (main run, 1x).

Estimated rate = called doublets / droplets. Two groups are reported separately:
- rate supplied: Scrublet, scDblFinder, DoubletFinder, Chord got the labeled rate as input
  (DoubletFinder and Chord call exactly round(rate * n), so their error is 0 by construction);
- self-estimated: DoubletDetection, Solo, vaeda, cxds/bcds/hybrid (estNdbl) choose the number.
Tools without a call (baselines, DoubletDensity) are not included.

Usage: python analysis/rate_estimation.py
Writes: results/stats/rate_estimation.tsv, results/figures/rate_estimation.png
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from style import CAT, FIG, GRID, INK2, NAMES, STATS, SURF, clean, main_metrics, setup

SUPPLIED = ["scrublet", "scdblfinder", "doubletfinder", "chord"]
d = main_metrics()
d = d[d.est_rate.notna()].copy()
d["ratio"] = d.est_rate / d.true_rate

rows = []
for t, g in d.groupby("tool"):
    rows.append({"tool": t, "rate_input": "supplied" if t in SUPPLIED else "self-estimated",
                 "median_error": g.rate_error.median(),
                 "mean_abs_error": g.rate_error.abs().mean(),
                 "median_ratio_est_true": g.ratio.median(),
                 "spearman_est_true": spearmanr(g.est_rate, g.true_rate).statistic,
                 "within_25pct": float((g.ratio.between(0.75, 1.25)).mean())})
s = pd.DataFrame(rows).sort_values(["rate_input", "mean_abs_error"], ascending=[False, True])
STATS.mkdir(parents=True, exist_ok=True)
s.to_csv(STATS / "rate_estimation.tsv", sep="\t", index=False, float_format="%.4f")
print(s.round(3).to_string(index=False))

# figure: small multiples, estimated vs labeled rate, one panel per tool
setup()
order = list(s.tool)
fig, axes = plt.subplots(2, 5, figsize=(14, 6.2), sharex=True, sharey=True)
lim = 0.45  # all labeled rates are < 0.4; points above are drawn on the edge with their value
for ax, t in zip(axes.flat, order):
    g = d[d.tool == t]
    col = CAT[1] if t in SUPPLIED else CAT[0]
    ax.plot([0, lim], [0, lim], color=INK2, lw=1, ls=(0, (3, 3)), zorder=1)
    inside = g.est_rate <= lim
    ax.scatter(g.true_rate[inside], g.est_rate[inside], s=36, color=col, edgecolor=SURF, lw=1,
               zorder=3)
    for x, y in zip(g.true_rate[~inside], g.est_rate[~inside]):   # off-scale: edge marker
        ax.scatter([x], [lim * 0.985], marker="^", s=50, color=col, edgecolor=SURF, zorder=3,
                   clip_on=False)
        ax.annotate(f"{y:.2f}", (x, lim * 0.985), xytext=(6, -3), textcoords="offset points",
                    fontsize=7.5, color=INK2)
    r = s.set_index("tool").loc[t]
    ax.set_title(f"{NAMES[t]}", loc="left", fontsize=10, pad=4)
    ax.text(0.97, 0.03, f"MAE {r.mean_abs_error:.3f}\nmedian est/true {r.median_ratio_est_true:.2f}",
            transform=ax.transAxes, va="bottom", ha="right", fontsize=7.5, color=INK2)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
    clean(ax, grid="both")
for ax in axes[1]:
    ax.set_xlabel("Labeled doublet rate")
for ax in axes[:, 0]:
    ax.set_ylabel("Called / droplets")
fig.suptitle("Estimated vs labeled doublet rate, 16 datasets (dashed = perfect; ▲ = off scale). "
             "Orange: labeled rate given as input; blue: tool chooses the number itself",
             x=0.01, ha="left", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.94), h_pad=2.5)
fig.savefig(FIG / "rate_estimation.png", dpi=200)
