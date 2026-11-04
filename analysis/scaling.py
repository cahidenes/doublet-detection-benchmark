"""Metric 9: runtime and peak memory vs. number of droplets.

Nested stratified subsamples of HMEC-orig-MULTI (scripts/subsample.py), every tool at the
labeled rate, one job at a time (--cores 4, 4 threads/job) so timings are not affected by
concurrent jobs. Empirical scaling exponent = slope of log(time) on log(n) over n >= 5k
(small n is dominated by fixed start-up costs: package loading, model setup).

Usage: python analysis/scaling.py   (reads results/scaling_metrics.tsv)
Writes: results/scaling_summary.tsv, results/figures/scaling.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "results/scaling_metrics.tsv", sep="\t")
d["n_cells"] = d.dataset.str.extract(r"-n(\d+)$")[0].astype(int)
d["gb"] = d.max_rss_mb / 1024

ok = d[d.status == "ok"]
rows = []
for t, g in ok.groupby("tool"):
    big = g[g.n_cells >= 5000]
    slope = np.polyfit(np.log(big.n_cells), np.log(big.wall_s), 1)[0] if len(big) >= 2 else np.nan
    top = g.loc[g.n_cells.idxmax()]
    rows.append({"tool": t, "time_exponent": slope, "n_max_ok": int(top.n_cells),
                 "wall_s_at_max": top.wall_s, "peak_gb_at_max": top.gb,
                 "auprc_at_max": top.auprc,
                 "n_failed": int((d[d.tool == t].status != "ok").sum())})
summ = pd.DataFrame(rows).sort_values("wall_s_at_max")
summ.to_csv(ROOT / "results/scaling_summary.tsv", sep="\t", index=False, float_format="%.3f")
print(summ.round(2).to_string(index=False))
print("\nwall time (s) by n:")
print(ok.pivot(index="tool", columns="n_cells", values="wall_s").round(0).to_string())

# --- figure: 13 tools are too many hues for one panel -> small multiples by family,
# every line direct-labeled; grey context lines show all other tools in each panel.
NAMES = {"libsize": "Library size", "ngenes": "# genes", "scrublet": "Scrublet",
         "doubletdetection": "DoubletDetection", "solo": "Solo", "vaeda": "vaeda",
         "scdblfinder": "scDblFinder", "dbldensity": "DoubletDensity", "cxds": "cxds",
         "bcds": "bcds", "hybrid": "hybrid", "doubletfinder": "DoubletFinder", "chord": "Chord"}
GROUPS = [("Python", ["scrublet", "doubletdetection", "libsize"]),
          ("Neural (GPU)", ["solo", "vaeda"]),
          ("scDblFinder / scds", ["scdblfinder", "dbldensity", "cxds", "bcds"]),
          ("Seurat-based / ensembles", ["doubletfinder", "chord", "hybrid", "ngenes"])]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
MARKERS = ["o", "s", "^", "D"]
INK, INK2, GRID, SURF, CTX = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb", "#d6d5d0"

plt.rcParams.update({"font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
fig, axes = plt.subplots(2, 4, figsize=(14, 6.6), facecolor=SURF, sharex=True)
for row, (col, ylab) in enumerate((("wall_s", "Wall time (s)"), ("gb", "Peak memory (GB)"))):
    piv = ok.pivot(index="n_cells", columns="tool", values=col)
    lo, hi = piv.min().min(), piv.max().max()
    for ax, (title, tools) in zip(axes[row], GROUPS):
        ax.set_facecolor(SURF)
        for t in piv.columns:                           # context: all tools, grey
            if t not in tools:
                ax.plot(piv.index, piv[t], color=CTX, lw=1, zorder=1)
        ends = []
        for i, t in enumerate(tools):
            s = piv[t].dropna()
            ax.plot(s.index, s.values, color=COLORS[i], lw=2, marker=MARKERS[i], ms=6,
                    mec=SURF, mew=1.2, zorder=3)
            ends.append((np.log10(s.values[-1]), t, s.index[-1], s.values[-1]))
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_ylim(lo / 1.5, hi * 1.5)
        # direct labels: spread in log space, then re-center
        ends.sort(); gap = 0.075 * np.log10(hi * 1.5 / (lo / 1.5))
        ys = [e[0] for e in ends]
        for k in range(1, len(ys)):
            ys[k] = max(ys[k], ys[k - 1] + gap)
        shift = (sum(ys) - sum(e[0] for e in ends)) / len(ys)
        for (y0, t, x, v), y in zip(ends, ys):
            ax.annotate(NAMES[t], (x, v), xytext=(x * 1.12, 10 ** (y - shift)), fontsize=8,
                        va="center", annotation_clip=False)
        ax.set_xlim(800, 6e4)
        ax.set_xticks([1e3, 2e3, 5e3, 1e4, 26426], ["1k", "2k", "5k", "10k", "26k"])
        ax.minorticks_off()
        if row == 0:
            ax.set_title(title, loc="left", fontsize=10, pad=8)
        if ax is axes[row][0]:
            ax.set_ylabel(ylab)
        ax.grid(color=GRID, lw=0.8); ax.tick_params(length=0)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
fig.suptitle("Runtime and peak memory vs. droplets (nested subsamples of HMEC-orig-MULTI, "
             "one job at a time, 4 threads; grey = other tools)", x=0.01, ha="left", fontsize=11)
fig.supxlabel("Droplets (log scale)", fontsize=9, color=INK2)
fig.tight_layout(rect=(0, 0, 1, 0.95), w_pad=4)
out = ROOT / "results/figures"; out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "scaling.png", dpi=200, facecolor=SURF)
