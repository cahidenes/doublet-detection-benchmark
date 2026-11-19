"""Per-dataset comparison with Xi & Li (2021) on the same 16 datasets.

Reference values: their Zenodo record 4444303, result.xlsx, sheet real_data
(saved as results/reference/xili2021_real.tsv). doubletCells (scran) is the predecessor of
computeDoubletDensity (DoubletDensity here).

For every tool and dataset: our AUPRC - their AUPRC. A large gap is either a change in the tool
or its dependencies, a difference in our wrapper, or run-to-run noise (cf. the seed spread).
For DoubletFinder, the pK picked by BCmvn in our run is listed next to the gap; its fixed-pK
sweep (results/sensitivity/df_pk.tsv) shows what the tool reaches with a small pK.

Usage: python analysis/xili_comparison.py
Writes: results/stats/xili_comparison.tsv, results/stats/xili_summary.tsv, results/figures/xili.png
"""
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from style import CAT, FIG, INK2, NAMES, ROOT, STATS, clean, main_metrics, setup

MAP = {"Scrublet": "scrublet", "doubletCells": "dbldensity", "cxds": "cxds", "bcds": "bcds",
       "hybrid": "hybrid", "DoubletDetection": "doubletdetection", "DoubletFinder": "doubletfinder",
       "solo": "solo", "scDblFinder": "scdblfinder"}
ref = pd.read_csv(ROOT / "results/reference/xili2021_real.tsv", sep="\t")
ref["tool"] = ref.method.map(MAP)
ours = main_metrics()[["tool", "dataset", "auprc", "auroc"]]
c = ref.merge(ours, on=["tool", "dataset"], suffixes=("_xili", "_ours"))
assert len(c) == 9 * 16, len(c)   # Xi & Li have no DoubletDetection value for pdx-MULTI (NaN)
c["d_auprc"] = c.auprc_ours - c.auprc_xili
c["d_auroc"] = c.auroc_ours - c.auroc_xili

# DoubletFinder: chosen pK in our run, and the best fixed pK from the sweep
pk = {}
for ds in c.dataset.unique():
    log = (ROOT / f"logs/doubletfinder/{ds}/s0_r1.0.log").read_text()
    pk[ds] = float(re.search(r"DoubletFinder: pK = ([\d.e-]+), nExp", log).group(1))
c["df_pk_bcmvn"] = np.where(c.tool == "doubletfinder", c.dataset.map(pk), np.nan)
sweep = ROOT / "results/sensitivity/df_pk.tsv"
if sweep.exists():
    sw = pd.read_csv(sweep, sep="\t")
    small = sw[sw.pK == 0.005].set_index("dataset").auprc
    c["df_auprc_pk0.005"] = np.where(c.tool == "doubletfinder", c.dataset.map(small), np.nan)
STATS.mkdir(parents=True, exist_ok=True)
c.drop(columns="method").to_csv(STATS / "xili_comparison.tsv", sep="\t", index=False, float_format="%.4f")

summ = c.groupby("tool").agg(mean_auprc_xili=("auprc_xili", "mean"), mean_auprc_ours=("auprc_ours", "mean"),
                             median_d_auprc=("d_auprc", "median"),
                             max_abs_d_auprc=("d_auprc", lambda x: x.abs().max()),
                             n_abs_d_gt_0_1=("d_auprc", lambda x: int((x.abs() > 0.1).sum())))
summ["spearman_per_dataset"] = {t: spearmanr(g.auprc_xili, g.auprc_ours, nan_policy="omit").statistic
                                for t, g in c.groupby("tool")}
summ["rank_xili"] = summ.mean_auprc_xili.rank(ascending=False)
summ["rank_ours"] = summ.mean_auprc_ours.rank(ascending=False)
summ = summ.sort_values("mean_auprc_xili", ascending=False)
summ.to_csv(STATS / "xili_summary.tsv", sep="\t", float_format="%.4f")
pd.set_option("display.width", 200)
print(summ.round(3).to_string())
print("\nLargest gaps (|d AUPRC| > 0.1):")
print(c[c.d_auprc.abs() > 0.1].sort_values("d_auprc")[
    ["tool", "dataset", "auprc_xili", "auprc_ours", "d_auprc", "df_pk_bcmvn"]].round(3).to_string(index=False))

# figure: ours vs theirs, one panel per tool
setup()
tools = list(summ.index)
fig, axes = plt.subplots(2, 5, figsize=(14, 6), sharex=True, sharey=True)
for ax, t in zip(axes.flat, tools):
    g = c[c.tool == t]
    ax.plot([0, 1], [0, 1], color=INK2, lw=0.8, ls=(0, (3, 3)))
    ax.scatter(g.auprc_xili, g.auprc_ours, s=16, color=CAT[0], lw=0)
    ax.set_title(f"{NAMES[t]}  (median Δ {g.d_auprc.median():+.2f})", loc="left", fontsize=9)
    clean(ax, grid="both")
axes.flat[-1].axis("off")
for ax in axes[1]:
    ax.set_xlabel("AUPRC, Xi & Li (2021)")
for ax in axes[:, 0]:
    ax.set_ylabel("AUPRC, this benchmark")
fig.suptitle("Per-dataset AUPRC: this benchmark vs Xi & Li (2021), same 16 datasets "
             "(DoubletDensity vs its predecessor doubletCells)", x=0.01, ha="left", fontsize=10.5)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(FIG / "xili.png", dpi=200)

# DoubletFinder: how much of the gap does pK explain? (fixed-pK sweep on all 16 datasets)
if sweep.exists():
    sw = pd.read_csv(sweep, sep="\t").pivot(index="dataset", columns="pK", values="auprc")
    dfr = c[c.tool == "doubletfinder"].set_index("dataset")
    rows = [{"setting": f"fixed pK = {k:g}", "mean_auprc": v} for k, v in sw.mean().items()]
    rows += [{"setting": "best fixed pK per dataset (oracle)", "mean_auprc": sw.max(1).mean()},
             {"setting": "automatic BCmvn (main run)", "mean_auprc": dfr.auprc_ours.mean()},
             {"setting": "Xi & Li (2021)", "mean_auprc": dfr.auprc_xili.mean()}]
    pk_s = pd.DataFrame(rows)
    pk_s.to_csv(STATS / "df_pk_summary.tsv", sep="\t", index=False, float_format="%.4f")
    print("\n", pk_s.round(3).to_string(index=False))
