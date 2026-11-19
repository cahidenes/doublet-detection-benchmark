"""Independent validation on singletCode data (Zhang et al. 2024).

Base cells: singlets identified by lineage barcodes (singletCode) in 6 Watermelon-barcoded T47D
samples (one breast-cancer cell line; naive, lag and late time points, 2 replicates each).
None of the tools was developed or tuned on these cells. Doublets were simulated from them with
scripts/simulate_doublets.py: random, homotypic and heteroweighted pairing, f = 1 and 0.5,
doublet rate 0.1 (= chance-level AUPRC).

The question is whether the real-data ranking, in particular scDblFinder's lead (its parameters
were tuned on the 16 real datasets), holds on data no tool has seen.

Usage: python analysis/singletcode.py
Writes: results/stats/singletcode_auprc.tsv, results/stats/singletcode_stats.tsv,
        results/figures/singletcode.png
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests

from style import FIG, INK2, NAMES, ROOT, SEQ, STATS, SURF, main_metrics, setup, text_on

s = pd.read_csv(ROOT / "results/singletcode_metrics.tsv", sep="\t")
print("runs:", s.status.value_counts().to_dict())
s = s[s.status == "ok"].copy()
parts = s.dataset.str.extract(r"^sim-wm-(?P<sample>T47D-[a-z]+-\d)-(?P<pairing>random|homotypic|heteroweighted)"
                              r"-f(?P<f>[\d.]+)$")
s = pd.concat([s, parts], axis=1)
s["scenario"] = s.pairing + ", f=" + s.f

order = list(pd.read_csv(STATS / "ranks.tsv", sep="\t", index_col=0).index)
SCEN = [f"{p}, f={f}" for f in ("1", "0.5") for p in ("heteroweighted", "random", "homotypic")]
SCEN = [c for c in SCEN if c in set(s.scenario)]
mean = s.pivot_table(index="tool", columns="scenario", values="auprc").loc[order, SCEN]
mean["overall"] = mean.mean(axis=1)
mean.to_csv(STATS / "singletcode_auprc.tsv", sep="\t", float_format="%.4f")
print(mean.round(3).to_string())

# statistics over the 6 samples, per scenario of interest (random f=1 = natural mix)
rows = []
for sc_ in SCEN + ["all"]:
    g = s if sc_ == "all" else s[s.scenario == sc_]
    m = g.pivot_table(index="dataset", columns="tool", values="auprc").dropna(axis=1)
    stat, p = friedmanchisquare(*[m[t] for t in m.columns])
    r = m.rank(axis=1, ascending=False).mean()
    best = r.idxmin()
    pw = []
    for t in m.columns:
        if t != best and (m[best] - m[t]).abs().sum() > 0:
            pw.append((t, wilcoxon(m[best], m[t]).pvalue))
    padj = dict(zip([t for t, _ in pw], multipletests([p_ for _, p_ in pw], method="holm")[1]))
    rows.append({"scenario": sc_, "n_datasets": len(m), "friedman_p": p, "best_mean_rank": best,
                 "scdblfinder_rank": r.rank().get("scdblfinder"),
                 "top3": ", ".join(r.sort_values().index[:3]),
                 "not_sig_worse_than_best": ", ".join(sorted(t for t, q in padj.items() if q >= 0.05))})
st = pd.DataFrame(rows)
st.to_csv(STATS / "singletcode_stats.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80)
print(st.to_string(index=False))

real = main_metrics().groupby("tool").auprc.mean()
rho = spearmanr(real.loc[order], mean.loc[order, "overall"]).statistic
print(f"\nSpearman(real mean AUPRC, singletCode overall) = {rho:.2f}")
pd.Series({"spearman_real_vs_singletcode": rho}).to_csv(STATS / "singletcode_vs_real.tsv", sep="\t",
                                                         header=False, float_format="%.4f")

setup()
mm = mean.loc[order, SCEN]
fig, ax = plt.subplots(figsize=(9.5, 6.4))
im = ax.imshow(mm.values, cmap=SEQ, vmin=0, vmax=1, aspect="auto")
for i in range(mm.shape[0]):
    for j in range(mm.shape[1]):
        v = mm.values[i, j]
        ax.text(j, i, f"{v:.2f}".lstrip("0"), ha="center", va="center", fontsize=8, color=text_on(v, 0, 1))
ax.set_yticks(range(len(order)), [f"{NAMES[t]}  (real {real[t]:.2f})" for t in order])
ax.set_xticks(range(len(SCEN)), [c.split(",")[0].replace("heteroweighted", "heterotypic") for c in SCEN])
ax.tick_params(length=0)
for sp_ in ax.spines.values():
    sp_.set_visible(False)
nf1 = sum(c.endswith("f=1") for c in SCEN)
ax.axvline(nf1 - .5, color=SURF, lw=4)
for x, t in (((nf1 - 1) / 2, "summed library (f = 1)"), (nf1 + (len(SCEN) - nf1 - 1) / 2, "singlet-like library (f = 0.5)")):
    ax.text(x, -0.9, t, ha="center", va="bottom", fontsize=9, color=INK2)
cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.outline.set_visible(False)
cb.set_label("AUPRC (mean of 6 T47D samples)", color=INK2)
fig.suptitle("Independent validation: singletCode singlets (Watermelon T47D) + simulated doublets, rate 0.10\n"
             f"Rows in real-data order; rank correlation with real data = {rho:.2f}",
             x=0.01, ha="left", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(FIG / "singletcode.png", dpi=200)
