"""Metrics 1-3 + statistics: AUPRC/AUROC across the 16 datasets.

- Friedman test over tools (blocks = datasets), average ranks (1 = best).
- All pairwise two-sided Wilcoxon signed-rank tests (paired by dataset), Holm-corrected
  within each metric; effect sizes: Cliff's delta over the 16 dataset values (unpaired, as
  committed in the proposal) and two paired ones: the matched-pairs rank-biserial r and the
  median paired difference with a 95% bootstrap CI (resampling datasets).
- Default-call F1: same tests over the tools that make calls (no baselines, no DoubletDensity).
- Sensitivity: datasets from the same study are not independent (5 PBMC demuxlet lanes, 3 HMEC
  MULTI-seq, 2 species mixes), so the Friedman test and ranks are repeated after averaging
  within study first (9 blocks).
- Median precision / recall / TNR at the 10/20/40% identification rates.

Usage: python analysis/accuracy_stats.py
Writes: results/stats/{friedman,ranks,pairwise_auprc,pairwise_auroc,pairwise_f1_default,
        identification,ranks_by_study}.tsv,
        results/figures/accuracy_heatmap.png, results/figures/pairwise_auprc.png
"""
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, wilcoxon
from statsmodels.stats.multitest import multipletests

from style import (DIV, FIG, INK, INK2, NAMES, SEQ, STATS, SURF, TECH, main_metrics, setup,
                   text_on)


def cliffs_delta(a, b):
    """P(a > b) - P(a < b) over all value pairs; +1 = a always higher."""
    a, b = np.asarray(a)[:, None], np.asarray(b)[None, :]
    return float((a > b).mean() - (a < b).mean())


def rank_biserial(diff):
    """Matched-pairs rank-biserial correlation: (sum of positive - negative signed ranks) / total."""
    diff = np.asarray(diff); diff = diff[diff != 0]
    if not len(diff):
        return 0.0
    r = pd.Series(np.abs(diff)).rank().to_numpy()
    return float((r[diff > 0].sum() - r[diff < 0].sum()) / r.sum())


def boot_ci(diff, n=5000, seed=0):
    rng = np.random.default_rng(seed); diff = np.asarray(diff)
    meds = np.median(diff[rng.integers(0, len(diff), (n, len(diff)))], axis=1)
    return np.quantile(meds, [0.025, 0.975])


# studies for the dependence sensitivity analysis (Xi & Li 2021, Table 1)
STUDY = {**{d: "PBMC demuxlet (Kang 2018)" for d in
            ["pbmc-1A-dm", "pbmc-1B-dm", "pbmc-1C-dm", "pbmc-2ctrl-dm", "pbmc-2stim-dm"]},
         **{d: "HMEC MULTI-seq" for d in ["HMEC-orig-MULTI", "HMEC-rep-MULTI", "HEK-HMEC-MULTI"]},
         "hm-12k": "Species mix (10x)", "hm-6k": "Species mix (10x)"}


def magnitude(d):  # Romano et al. (2006) thresholds
    d = abs(d)
    return "negligible" if d < .147 else "small" if d < .33 else "medium" if d < .474 else "large"


d = main_metrics()
STATS.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)

CALLERS = [t for t in d.tool.unique() if t not in ("libsize", "ngenes", "dbldensity")]
fried, ranks, mats = [], {}, {}
for metric in ("auprc", "auroc", "f1_default"):
    m = d.pivot(index="dataset", columns="tool", values=metric)       # datasets x tools
    if metric == "f1_default":
        m = m[sorted(CALLERS)]
    assert m.notna().all().all(), f"missing {metric} values"
    mats[metric] = m
    stat, p = friedmanchisquare(*[m[t] for t in m.columns])
    fried.append({"metric": metric, "chi2": stat, "df": m.shape[1] - 1, "p": p,
                  "n_datasets": m.shape[0], "n_tools": m.shape[1]})
    r = m.rank(axis=1, ascending=False)                                  # 1 = best
    ranks[metric] = r.mean().rename(f"mean_rank_{metric}")

    rows = []
    for a, b in combinations(m.columns, 2):
        diff = m[a] - m[b]
        p = wilcoxon(m[a], m[b]).pvalue if (diff != 0).any() else 1.0
        cd = cliffs_delta(m[a], m[b])
        lo, hi = boot_ci(diff)
        rows.append({"tool_a": a, "tool_b": b, "median_diff": diff.median(),
                     "median_diff_lo": lo, "median_diff_hi": hi,
                     "rank_biserial": rank_biserial(diff),
                     "wins_a": int((diff > 0).sum()), "wins_b": int((diff < 0).sum()),
                     "p": p, "cliffs_delta": cd, "magnitude": magnitude(cd)})
    pw = pd.DataFrame(rows)
    pw["p_holm"] = multipletests(pw.p, method="holm")[1]
    pw["significant"] = pw.p_holm < 0.05
    pw.sort_values("p").to_csv(STATS / f"pairwise_{metric}.tsv", sep="\t", index=False,
                               float_format="%.4g")

# sensitivity: average within study first, then Friedman + mean ranks over 9 blocks
d["study"] = d.dataset.map(lambda x: STUDY.get(x, x))
by_study = []
for metric in ("auprc", "auroc"):
    m = d.groupby(["study", "tool"])[metric].mean().unstack()
    stat, p = friedmanchisquare(*[m[t] for t in m.columns])
    fried.append({"metric": f"{metric} (study-averaged)", "chi2": stat, "df": m.shape[1] - 1,
                  "p": p, "n_datasets": m.shape[0], "n_tools": m.shape[1]})
    by_study.append(m.rank(axis=1, ascending=False).mean().rename(f"mean_rank_{metric}_study"))
    by_study.append(m.mean().rename(f"mean_{metric}_study"))

pd.DataFrame(fried).to_csv(STATS / "friedman.tsv", sep="\t", index=False, float_format="%.4g")
summ = pd.concat([mats["auprc"].mean().rename("mean_auprc"),
                  mats["auprc"].median().rename("median_auprc"), ranks["auprc"],
                  mats["auroc"].mean().rename("mean_auroc"),
                  mats["auroc"].median().rename("median_auroc"), ranks["auroc"]], axis=1)
summ = summ.sort_values(["mean_rank_auprc", "mean_auprc"], ascending=[True, False])  # ties -> higher mean
summ.index.name = "tool"
summ["median_f1_default"] = mats["f1_default"].median()
summ["mean_rank_f1_default"] = ranks["f1_default"]
summ.to_csv(STATS / "ranks.tsv", sep="\t", float_format="%.4f")
bs = pd.concat(by_study, axis=1).loc[summ.index]
bs.index.name = "tool"
bs.to_csv(STATS / "ranks_by_study.tsv", sep="\t", float_format="%.4f")

ident = d.groupby("tool")[[c for c in d.columns if "@" in c]].median().loc[summ.index]
ident.to_csv(STATS / "identification.tsv", sep="\t", float_format="%.4f")

print(pd.DataFrame(fried).to_string(index=False))
print(summ.round(3).to_string())
pw = pd.read_csv(STATS / "pairwise_auprc.tsv", sep="\t")
print(f"\nAUPRC: {pw.significant.sum()} of {len(pw)} pairs significant after Holm")
top = summ.index[0]
print(pw[(pw.tool_a == top) | (pw.tool_b == top)].round(4).to_string(index=False))
pf = pd.read_csv(STATS / "pairwise_f1_default.tsv", sep="\t")
print(f"\nF1: {pf.significant.sum()} of {len(pf)} pairs significant after Holm")
print(pf[pf.significant].round(4).to_string(index=False))
print("\nStudy-averaged ranks\n", bs.round(3).to_string())

# ------------------------------------------------------------------ figures
setup()
order = list(summ.index)
cols = sorted(mats["auprc"].index, key=lambda s: (list(dict.fromkeys(TECH.values())).index(TECH[s]), s))

# Fig: AUPRC heatmap, tools (rows, by mean rank) x datasets (columns, by technology)
m = mats["auprc"].T.loc[order, cols]
fig, ax = plt.subplots(figsize=(12.5, 6.4))
im = ax.imshow(m.values, cmap=SEQ, vmin=0, vmax=1, aspect="auto")
for i in range(m.shape[0]):
    for j in range(m.shape[1]):
        v = m.values[i, j]
        ax.text(j, i, f"{v:.2f}".lstrip("0"), ha="center", va="center", fontsize=7.5,
                color=text_on(v, 0, 1))
ax.set_yticks(range(len(order)),
              [f"{NAMES[t]}  ({summ.loc[t, 'mean_rank_auprc']:.1f})" for t in order])
ax.set_xticks(range(len(cols)), cols, rotation=40, ha="right")
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
# separators + labels between technologies (light lines on the surface)
techs = [TECH[c] for c in cols]
starts = [0] + [j for j in range(1, len(cols)) if techs[j] != techs[j - 1]] + [len(cols)]
for a, b in zip(starts[:-1], starts[1:]):
    if a:
        ax.axvline(a - 0.5, color=SURF, lw=3)
    ax.text((a + b - 1) / 2, -0.9, techs[a], ha="center", va="bottom", fontsize=8.5, color=INK2)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
cb.outline.set_visible(False); cb.set_label("AUPRC", color=INK2)
fig.suptitle("AUPRC per tool and dataset (rows ordered by mean rank, shown in brackets; "
             "1 = best)", x=0.01, ha="left", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(FIG / "accuracy_heatmap.png", dpi=200)

# Fig: pairwise Cliff's delta (row tool vs column tool), Holm-significant pairs marked
pw = pd.read_csv(STATS / "pairwise_auprc.tsv", sep="\t")
n = len(order); cd = np.full((n, n), np.nan); sig = np.zeros((n, n), bool)
idx = {t: i for i, t in enumerate(order)}
for r in pw.itertuples():
    i, j = idx[r.tool_a], idx[r.tool_b]
    cd[i, j], cd[j, i] = r.cliffs_delta, -r.cliffs_delta
    sig[i, j] = sig[j, i] = r.significant
fig, ax = plt.subplots(figsize=(8.6, 7.4))
im = ax.imshow(cd, cmap=DIV, vmin=-1, vmax=1)
for i in range(n):
    for j in range(n):
        if i != j:
            ax.text(j, i, f"{cd[i, j]:+.2f}" + ("*" if sig[i, j] else ""), ha="center",
                    va="center", fontsize=6.8, color="white" if abs(cd[i, j]) > .6 else INK,
                    fontweight="bold" if sig[i, j] else "normal")
labels = [NAMES[t] for t in order]
ax.set_xticks(range(n), labels, rotation=45, ha="right"); ax.set_yticks(range(n), labels)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
cb.outline.set_visible(False); cb.set_label("Cliff's δ (row vs column)", color=INK2)
fig.suptitle("Pairwise AUPRC comparison: Cliff's δ over 16 datasets\n"
             "red = row tool higher; * bold = Wilcoxon signed-rank p < 0.05 after Holm",
             x=0.01, ha="left", fontsize=10.5)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(FIG / "pairwise_auprc.png", dpi=200)
