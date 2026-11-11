"""Metric 6: agreement between tools (main run, labeled rate, seed 0).

Per dataset, then averaged over the 16 datasets:
- Spearman correlation of the per-droplet scores (all 13 methods);
- Jaccard overlap of the top-k droplets, k = number of labeled doublets. Every tool has a
  score, so this compares all tools at the same budget, independent of default thresholds.
Plus "hard doublets": labeled doublets that no non-baseline tool puts in its top k, and
"oracle union": share of doublets found by at least one tool (upper bound for an ensemble
that could pick the right tool per droplet).

Usage: python analysis/agreement.py   (needs the local results/pred/ files)
Writes: results/stats/agreement_{spearman,topk_jaccard}.tsv, results/stats/hard_doublets.tsv,
        results/figures/agreement.png
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

from style import (BASELINES, FIG, INK2, NAMES, ROOT, SEQ, STATS, main_metrics, setup,
                   text_on)

d = main_metrics()
tools = sorted(d.tool.unique())
datasets = sorted(d.dataset.unique())
STATS.mkdir(parents=True, exist_ok=True)

rho_all, jac_all, hard = [], [], []
for ds in datasets:
    lab = pd.read_csv(ROOT / f"data/processed/{ds}/labels.tsv", sep="\t").set_index("barcode")
    S = pd.DataFrame({t: pd.read_csv(ROOT / f"results/pred/{t}/{ds}/s0_r1.0.tsv", sep="\t")
                      .set_index("barcode").score for t in tools}).loc[lab.index]
    rho_all.append(S.rank().corr().values)          # Spearman = Pearson on ranks
    k = int(lab.is_doublet.sum())
    top = {t: set(S[t].rank(ascending=False, method="first").nsmallest(k).index) for t in tools}
    J = np.array([[len(top[a] & top[b]) / len(top[a] | top[b]) for b in tools] for a in tools])
    jac_all.append(J)
    dbl = set(lab.index[lab.is_doublet == 1])
    real = [t for t in tools if t not in BASELINES]
    found_any = set().union(*(top[t] & dbl for t in real))
    best = max(len(top[t] & dbl) for t in real)
    hard.append({"dataset": ds, "n_doublets": k,
                 "best_single_recall": best / k, "oracle_union_recall": len(found_any) / k,
                 "hard_fraction": 1 - len(found_any) / k})

rho = pd.DataFrame(np.mean(rho_all, axis=0), index=tools, columns=tools)
jac = pd.DataFrame(np.mean(jac_all, axis=0), index=tools, columns=tools)
rho.to_csv(STATS / "agreement_spearman.tsv", sep="\t", float_format="%.3f")
jac.to_csv(STATS / "agreement_topk_jaccard.tsv", sep="\t", float_format="%.3f")
hard = pd.DataFrame(hard)
hard.to_csv(STATS / "hard_doublets.tsv", sep="\t", index=False, float_format="%.4f")

# order by average-linkage clustering on 1 - Spearman
order = [tools[i] for i in leaves_list(linkage(squareform(1 - rho.values, checks=False), "average"))]
off = ~np.eye(len(tools), dtype=bool)
print("mean pairwise Spearman (non-baseline):",
      rho.drop(index=BASELINES, columns=BASELINES).values[off[:-2, :-2]].mean().round(3))
print(rho.loc[order, order].round(2).to_string())
print(jac.loc[order, order].round(2).to_string())
print(hard.round(3).to_string(index=False))
print("median best single / oracle union recall:",
      hard.best_single_recall.median().round(3), hard.oracle_union_recall.median().round(3))

# ------------------------------------------------------------------ figure
setup()
fig, axes = plt.subplots(1, 2, figsize=(15, 6.8))
for ax, (M, title, vmin) in zip(axes, [(rho, "Spearman correlation of scores", 0),
                                       (jac, "Jaccard overlap of top-k droplets "
                                             "(k = labeled doublets)", 0)]):
    M = M.loc[order, order]
    im = ax.imshow(M.values, cmap=SEQ, vmin=vmin, vmax=1)
    for i in range(len(order)):
        for j in range(len(order)):
            if i != j:
                v = M.values[i, j]
                ax.text(j, i, f"{v:.2f}".replace("0.", "."), ha="center", va="center",
                        fontsize=7, color=text_on(v, vmin, 1))
    labels = [NAMES[t] for t in order]
    ax.set_xticks(range(len(order)), labels, rotation=45, ha="right")
    ax.set_yticks(range(len(order)), labels)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title, loc="left", fontsize=10, pad=8)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.outline.set_visible(False)
fig.suptitle("Agreement between tools, mean over 16 datasets (order: hierarchical clustering "
             "on score correlation)", x=0.01, ha="left", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.95), w_pad=3)
fig.savefig(FIG / "agreement.png", dpi=200)
