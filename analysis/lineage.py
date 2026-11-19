"""Dimension 15: lineage co-expression in the PBMC datasets (label-free check of predicted doublets).

A droplet co-expresses two incompatible immune lineages if it passes the marker threshold for
two of: T cell (CD3E + CD3D >= 2 UMIs), B cell (CD79A + MS4A1 >= 2), myeloid (CD14 + LYZ >= 5;
higher because LYZ is abundant in ambient RNA). Such droplets are heterotypic doublets with
high confidence; the check cannot confirm homotypic doublets (e.g. T + T), which are the
majority in PBMC.

Per tool and dataset, on the top-k droplets (k = labeled doublets) and on the default calls:
- coexp_frac: fraction of the tool's doublets that co-express (vs. the labeled doublets and
  the labeled singlets as references);
- coexp_recall: fraction of all co-expressing droplets the tool flags;
- coexp_unlabeled_found: co-expressing droplets labeled singlet (likely missed by the labels,
  e.g. same-donor doublets) that the tool flags.

Usage: python analysis/lineage.py
Writes: results/stats/lineage.tsv, results/stats/lineage_summary.tsv, results/figures/lineage.png
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io

from style import CAT, FIG, INK2, NAMES, ROOT, STATS, clean, main_metrics, setup

PBMC = ["pbmc-ch", "pbmc-1A-dm", "pbmc-1B-dm", "pbmc-1C-dm", "pbmc-2ctrl-dm", "pbmc-2stim-dm"]
MARKERS = {  # lineage: (symbols, Ensembl IDs, min UMIs)
    "T": (["CD3E", "CD3D"], ["ENSG00000198851", "ENSG00000167286"], 2),
    "B": (["CD79A", "MS4A1"], ["ENSG00000105369", "ENSG00000156738"], 2),
    "Myeloid": (["CD14", "LYZ"], ["ENSG00000170458", "ENSG00000090382"], 5),
}


def coexpression(ds):
    d = ROOT / "data/processed" / ds
    X = scipy.io.mmread(d / "matrix.mtx.gz").tocsr()            # genes x droplets
    genes = np.array((d / "genes.tsv").read_text().split())
    pos = []
    for sym, ens, th in MARKERS.values():
        rows = np.isin(genes, sym + ens)
        assert rows.sum() == 2, (ds, sym)
        pos.append(np.asarray(X[rows].sum(0)).ravel() >= th)
    pos = np.array(pos)
    return pos.sum(0) >= 2


rows = []
tools = main_metrics().tool.unique()
for ds in PBMC:
    d = ROOT / "data/processed" / ds
    bc = (d / "barcodes.tsv").read_text().split()
    lab = pd.read_csv(d / "labels.tsv", sep="\t").set_index("barcode").loc[bc, "is_doublet"].to_numpy()
    co = coexpression(ds)
    k = int(lab.sum())
    ref = {"labeled doublets": lab == 1, "labeled singlets": lab == 0}
    for name, sel in ref.items():
        rows.append({"dataset": ds, "tool": name, "set": "reference", "n": int(sel.sum()),
                     "coexp_frac": co[sel].mean(), "coexp_recall": (co & sel).sum() / co.sum(),
                     "coexp_unlabeled_found": int((co & sel & (lab == 0)).sum())})
    for tool in tools:
        p = pd.read_csv(ROOT / f"results/pred/{tool}/{ds}/s0_r1.0.tsv", sep="\t").set_index("barcode").loc[bc]
        top = np.zeros(len(bc), bool); top[np.argsort(-p.score.to_numpy(), kind="stable")[:k]] = True
        sets = {"top-k": top}
        if p["call"].notna().any():
            sets["default calls"] = p["call"].fillna(0).to_numpy().astype(bool)
        for name, sel in sets.items():
            rows.append({"dataset": ds, "tool": tool, "set": name, "n": int(sel.sum()),
                         "coexp_frac": co[sel].mean() if sel.any() else np.nan,
                         "coexp_recall": (co & sel).sum() / co.sum(),
                         "coexp_unlabeled_found": int((co & sel & (lab == 0)).sum())})
    print(f"{ds}: {co.sum()} co-expressing droplets, {int((co & (lab == 0)).sum())} labeled singlet")

out = pd.DataFrame(rows)
STATS.mkdir(parents=True, exist_ok=True)
out.to_csv(STATS / "lineage.tsv", sep="\t", index=False, float_format="%.4f")

summ = (out.groupby(["tool", "set"])
        .agg(median_coexp_frac=("coexp_frac", "median"), median_coexp_recall=("coexp_recall", "median"),
             unlabeled_found=("coexp_unlabeled_found", "sum"))
        .reset_index())
# total co-expressing droplets labeled singlet, over the 6 datasets
n_unl = int(out[out.tool == "labeled singlets"].coexp_unlabeled_found.sum())
summ["unlabeled_total"] = n_unl
summ = summ.sort_values(["set", "median_coexp_recall"], ascending=[True, False])
summ.to_csv(STATS / "lineage_summary.tsv", sep="\t", index=False, float_format="%.4f")
print(summ.round(3).to_string(index=False))

# figure: recall of co-expressing droplets at top-k, per tool (median + per-dataset points)
setup()
t = out[out.set == "top-k"]
order = t.groupby("tool").coexp_recall.median().sort_values().index
fig, ax = plt.subplots(figsize=(7.2, 4.6))
for i, tool in enumerate(order):
    v = t[t.tool == tool].coexp_recall
    ax.scatter(v, np.full(len(v), i) + np.linspace(-.18, .18, len(v)), s=14, color=CAT[0], alpha=.55, lw=0)
    ax.plot([v.median()] * 2, [i - .3, i + .3], color=CAT[1], lw=2.2)
ref = out[out.tool == "labeled doublets"].coexp_recall.median()
ax.axvline(ref, color=INK2, lw=1, ls="--")
ax.text(ref, len(order) - .3, " labeled doublets", color=INK2, fontsize=8, va="bottom")
ax.set_yticks(range(len(order)), [NAMES[x] for x in order])
ax.set_xlabel("Fraction of lineage-co-expressing droplets ranked in the tool's top k")
ax.set_xlim(0, 1)
clean(ax, grid="x")
fig.suptitle("Lineage co-expression (T / B / myeloid markers), 6 PBMC datasets\n"
             "orange = median; dots = datasets; k = number of labeled doublets",
             x=0.01, ha="left", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(FIG / "lineage.png", dpi=200)
