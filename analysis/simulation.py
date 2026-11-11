"""Stretch goal: accuracy on synthetic doublets with known composition (scripts/simulate_doublets.py).

12 datasets = 2 bases (pbmc-2ctrl-dm, pdx-MULTI singlets) x 3 pairings (random, heterotypic,
homotypic) x 2 depths (f = 1: summed library; f = 0.5: singlet-like library). Every tool at its
defaults with the true rate (0.1).

Usage: python analysis/simulation.py
Writes: results/stats/simulation_auprc.tsv, results/figures/simulation.png
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from style import FIG, INK2, NAMES, ROOT, SEQ, STATS, SURF, main_metrics, setup, text_on

s = pd.read_csv(ROOT / "results/sim_metrics.tsv", sep="\t")
s = s[s.status == "ok"].copy()
parts = s.dataset.str.extract(r"^sim-(?P<base>.+)-(?P<pairing>random|heterotypic|homotypic)-f(?P<f>[\d.]+)$")
s = pd.concat([s, parts], axis=1)
s["scenario"] = s.pairing + ", f=" + s.f

order = list(pd.read_csv(STATS / "ranks.tsv", sep="\t", index_col=0).index)
SCEN = [f"{p}, f={f}" for f in ("1", "0.5") for p in ("heterotypic", "random", "homotypic")]

tab = s.pivot_table(index="tool", columns=["scenario", "base"], values="auprc")
mean = s.pivot_table(index="tool", columns="scenario", values="auprc").loc[order, SCEN]
mean["overall"] = mean.mean(axis=1)
STATS.mkdir(parents=True, exist_ok=True)
mean.to_csv(STATS / "simulation_auprc.tsv", sep="\t", float_format="%.4f")
tab.loc[order].to_csv(STATS / "simulation_auprc_by_base.tsv", sep="\t", float_format="%.4f")
print(mean.round(3).to_string())

real = main_metrics().groupby("tool").auprc.mean()
rho = spearmanr(real.loc[order], mean.loc[order, "overall"]).statistic
rho_r = spearmanr(real.loc[order], mean.loc[order, "random, f=1"]).statistic
print(f"\nSpearman(real mean AUPRC, sim overall) = {rho:.2f}; vs random f=1 only = {rho_r:.2f}")
pd.Series({"spearman_real_vs_sim_overall": rho, "spearman_real_vs_sim_random_f1": rho_r}) \
  .to_csv(STATS / "simulation_vs_real.tsv", sep="\t", header=False, float_format="%.4f")

# ---------------------------------------------------------------- figure
setup()
m = mean.loc[order, SCEN]
fig, ax = plt.subplots(figsize=(9.5, 6.4))
im = ax.imshow(m.values, cmap=SEQ, vmin=0, vmax=1, aspect="auto")
for i in range(m.shape[0]):
    for j in range(m.shape[1]):
        v = m.values[i, j]
        ax.text(j, i, f"{v:.2f}".lstrip("0"), ha="center", va="center", fontsize=8,
                color=text_on(v, 0, 1))
ax.set_yticks(range(len(order)), [f"{NAMES[t]}  (real {real[t]:.2f})" for t in order])
ax.set_xticks(range(len(SCEN)), [c.split(",")[0] for c in SCEN])
ax.tick_params(length=0)
for sp_ in ax.spines.values():
    sp_.set_visible(False)
ax.axvline(2.5, color=SURF, lw=4)
for x, t in ((1, "summed library (f = 1, ≈2× a singlet)"), (4, "singlet-like library (f = 0.5)")):
    ax.text(x, -0.9, t, ha="center", va="bottom", fontsize=9, color=INK2)
cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.outline.set_visible(False)
cb.set_label("AUPRC (mean of 2 base datasets)", color=INK2)
fig.suptitle("AUPRC on synthetic doublets by pairing and depth (doublet rate 0.10 = chance level)\n"
             "Rows in real-data order; bracket = real-data mean AUPRC; "
             f"rank correlation real vs simulated = {rho:.2f}",
             x=0.01, ha="left", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(FIG / "simulation.png", dpi=200)
