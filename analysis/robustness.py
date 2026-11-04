"""Metric 5: robustness to a misspecified expected doublet rate.

Tools that take an expected rate (rate: true in config) are run with 0.5x, 1x and 2x the labeled
rate. Their default calls are compared on precision, recall, F1 and the number of droplets
removed relative to the true number of doublets. AUPRC/AUROC are reported too: they should
barely move, because the rate mostly sets the threshold (for DoubletFinder pANN does not
depend on nExp at all).

Usage: python analysis/robustness.py   (reads results/metrics.tsv)
Writes: results/robustness_summary.tsv, results/figures/robustness.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
cfg = yaml.safe_load(open(ROOT / "config/config.yaml"))
RATE_TOOLS = [t for t, s in cfg["tools"].items() if s.get("rate")]

d = pd.read_csv(ROOT / "results/metrics.tsv", sep="\t")
d = d[(d.status == "ok") & d.tool.isin(RATE_TOOLS)].copy()
d["removed_over_true"] = d.n_removed / (d.true_rate * d.n)
# paired comparison: keep a (tool, dataset) only if it succeeded at every multiplier
# (Chord fails at 2x on HMEC-rep-MULTI and mkidney-ch: 2x their rate exceeds 0.6)
n_mult = d.rate_mult.nunique()
d = d[d.groupby(["tool", "dataset"]).rate_mult.transform("nunique") == n_mult]

METRICS = {"precision_default": "Precision", "recall_default": "Recall", "f1_default": "F1",
           "removed_over_true": "Removed ÷ true doublets", "auprc": "AUPRC"}

# --- summary table: median and IQR across datasets, per tool x multiplier ---
g = d.groupby(["tool", "rate_mult"])
summ = pd.concat({f"{m}_{q}": g[m].quantile(p) for m in METRICS
                  for q, p in (("q25", .25), ("median", .5), ("q75", .75))}, axis=1)
summ["n_datasets"] = g.size()
summ = summ.reset_index()
summ.to_csv(ROOT / "results/robustness_summary.tsv", sep="\t", index=False, float_format="%.4f")

# --- figure: small multiples, one panel per metric, median across datasets ---
# Categorical slots 1-4 of the reference palette (validated: scripts/validate_palette.js).
# Yellow/aqua are < 3:1 on the surface -> every line is direct-labeled, distinct markers.
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
MARKERS = ["o", "s", "^", "D"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
NAMES = {"scrublet": "Scrublet", "scdblfinder": "scDblFinder",
         "doubletfinder": "DoubletFinder", "chord": "Chord"}
tools = [t for t in RATE_TOOLS if t in d.tool.unique()]
panels = ["precision_default", "recall_default", "f1_default", "removed_over_true"]

plt.rcParams.update({"font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
fig, axes = plt.subplots(1, 4, figsize=(13, 3.6), facecolor=SURF)
for ax, m in zip(axes, panels):
    ax.set_facecolor(SURF)
    med = summ.pivot(index="rate_mult", columns="tool", values=f"{m}_median")
    for i, t in enumerate(tools):
        ax.plot(med.index, med[t], color=COLORS[i], lw=2, marker=MARKERS[i], ms=7,
                mec=SURF, mew=1.5, zorder=3)
    # direct labels at the right end; nudge apart if they collide
    ends = sorted(((med[t].iloc[-1], t, i) for i, t in enumerate(tools)))
    span = ax.get_ylim()[1] - ax.get_ylim()[0]
    gap = 0.075 * span
    ys = [e[0] for e in ends]
    for k in range(1, len(ys)):                      # push apart upward...
        ys[k] = max(ys[k], ys[k - 1] + gap)
    shift = (sum(ys) - sum(e[0] for e in ends)) / len(ys)  # ...then re-center the group
    for (y0, t, i), y in zip(ends, ys):
        ax.annotate(NAMES[t], (2.0, y0), xytext=(2.08, y - shift), fontsize=8,
                    color=INK, va="center", annotation_clip=False)
    if m == "removed_over_true":
        ax.axhline(1, color=INK2, lw=1, ls=(0, (3, 3)), zorder=1)
        ax.annotate("calls = true doublets", (1.4, 1), xytext=(0, -11), textcoords="offset points",
                    fontsize=7.5, color=INK2)
    ax.set_xscale("log", base=2)
    ax.set_xticks([0.5, 1, 2], ["0.5×", "1×", "2×"])
    ax.set_xlim(0.45, 2.15)
    ax.set_title(METRICS[m], loc="left", fontsize=10, color=INK, pad=18)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.tick_params(length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
nd = summ.groupby("tool").n_datasets.first()
fig.suptitle("Default-call metrics under a misspecified expected doublet rate (median over "
             f"{int(nd.max())} datasets; Chord {int(nd['chord'])}, failed at 2× on 2)",
             x=0.01, ha="left", fontsize=11, color=INK)
fig.supxlabel("Expected rate given to tool (× true rate)", fontsize=9, color=INK2)
fig.tight_layout(rect=(0, 0, 0.97, 0.93), w_pad=6)
out = ROOT / "results/figures"; out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "robustness.png", dpi=200, facecolor=SURF)

cols = ["tool", "rate_mult"] + [f"{m}_median" for m in METRICS] + ["n_datasets"]
print(summ[cols].round(3).to_string(index=False))
