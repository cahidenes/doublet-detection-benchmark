"""Stretch goal: downstream impact of each tool's doublet removal (real data, main run).

For every dataset, droplets are removed according to each tool, then a standard scanpy
analysis is run on what is left:
  normalize_total, log1p, 2000 HVG (seurat), PCA 30, kNN 15, Leiden res 1.0, t-test markers.
Removal: the tool's default calls; for score-only methods (baselines, DoubletDensity), the
top round(labeled rate x n) droplets by score. References:
  oracle   remove exactly the labeled doublets (the target)
  none     remove nothing (what happens without doublet detection)
  reseed   oracle again with another clustering seed: the noise floor of the comparison

Metrics vs the oracle analysis:
  hvg_jaccard       overlap of the 2000 HVGs
  ari               adjusted Rand index of Leiden clusters, on the singlets kept by both
  marker_jaccard    each oracle cluster is matched to the tool cluster holding most of its
                    cells; mean Jaccard of their top-50 markers
  doublet_clusters  clusters (>= 20 droplets) made of >= 50% labeled doublets: spurious
                    "doublet cell types" left in the data
Trajectory inference is not evaluated: none of these datasets has a ground-truth trajectory.

Usage: python analysis/downstream.py [--jobs 6]
Writes: results/downstream_metrics.tsv, results/stats/downstream_summary.tsv,
        results/figures/downstream.png
"""
import argparse
import os
import warnings

os.environ.setdefault("OMP_NUM_THREADS", "2")
warnings.filterwarnings("ignore")

from multiprocessing import Pool

import numpy as np
import pandas as pd

from style import FIG, NAMES, ROOT, STATS, main_metrics

NMARK = 50


def _patch_scanpy():
    """scanpy 1.11 calls .getnnz() on the kNN distance matrix, which newer scipy/sklearn
    return as a sparse *array* (no getnnz) on larger datasets: convert it first."""
    import scipy.sparse as sp
    from scanpy.neighbors import _common
    if getattr(_common, "_patched", False):
        return
    orig = _common._ind_dist_shortcut
    _common._ind_dist_shortcut = lambda d: orig(sp.csr_matrix(d))
    _common._patched = True


def pipeline(X, keep, seed=0):
    import anndata as ad
    import scanpy as sc
    _patch_scanpy()
    A = ad.AnnData(X[keep].astype(np.float32))
    sc.pp.filter_genes(A, min_cells=3)
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    sc.pp.highly_variable_genes(A, n_top_genes=2000)
    hvg = set(A.var_names[A.var.highly_variable])
    sc.pp.pca(A, n_comps=30, mask_var="highly_variable", random_state=seed)
    sc.pp.neighbors(A, n_neighbors=15, random_state=seed)
    sc.tl.leiden(A, resolution=1.0, random_state=seed, flavor="igraph", n_iterations=2)
    cl = A.obs.leiden.astype(int).values
    markers = {}
    if len(np.unique(cl)) > 1:
        sc.tl.rank_genes_groups(A, "leiden", method="t-test", n_genes=NMARK)
        names = A.uns["rank_genes_groups"]["names"]
        markers = {int(g): set(names[g]) for g in names.dtype.names}
    return {"keep": np.flatnonzero(keep), "cl": cl, "hvg": hvg, "markers": markers}


def compare(ref, res, is_dbl):
    from sklearn.metrics import adjusted_rand_score
    out = {"hvg_jaccard": len(ref["hvg"] & res["hvg"]) / len(ref["hvg"] | res["hvg"])}
    r = pd.Series(ref["cl"], index=ref["keep"]); t = pd.Series(res["cl"], index=res["keep"])
    common = r.index.intersection(t.index)
    out["ari"] = adjusted_rand_score(r[common], t[common])
    js = []
    for c, mk in ref["markers"].items():
        cells = r.index[r.values == c].intersection(t.index)
        if len(cells) == 0:
            js.append(0.0); continue
        m = t[cells].value_counts().idxmax()
        tm = res["markers"].get(m, set())
        js.append(len(mk & tm) / len(mk | tm) if (mk | tm) else 1.0)
    out["marker_jaccard"] = float(np.mean(js)) if js else np.nan
    dfrac = pd.Series(is_dbl[res["keep"]]).groupby(res["cl"]).agg(["mean", "size"])
    out["doublet_clusters"] = int(((dfrac["mean"] >= 0.5) & (dfrac["size"] >= 20)).sum())
    out["n_clusters"] = int(len(dfrac))
    return out


def run_dataset(ds):
    import scipy.io
    d = main_metrics(); d = d[d.dataset == ds]
    src = ROOT / "data/processed" / ds
    import scipy.sparse as sp
    X = sp.csr_matrix(scipy.io.mmread(src / "matrix.mtx.gz").T)  # scanpy needs spmatrix, not sparray
    lab = pd.read_csv(src / "labels.tsv", sep="\t")
    is_dbl = lab.is_doublet.values.astype(bool); n = len(lab); k = int(is_dbl.sum())
    ref = pipeline(X, ~is_dbl)
    rows = []

    def record(name, keep, seed=0):
        res = pipeline(X, keep, seed)
        r = compare(ref, res, is_dbl)
        r.update({"dataset": ds, "tool": name, "n_removed": int((~keep).sum()),
                  "residual_doublets": int((keep & is_dbl).sum()),
                  "singlets_removed": int((~keep & ~is_dbl).sum())})
        rows.append(r)

    record("reseed", ~is_dbl, seed=1)
    record("none", np.ones(n, bool))
    for t in sorted(d.tool):
        p = pd.read_csv(ROOT / f"results/pred/{t}/{ds}/s0_r1.0.tsv", sep="\t")
        p = p.set_index("barcode").loc[lab.barcode]
        if p["call"].notna().any():
            remove = p["call"].fillna(0).values.astype(bool)
        else:
            remove = np.zeros(n, bool)
            remove[np.argsort(-p.score.values, kind="stable")[:k]] = True
        record(t, ~remove)
    print(f"{ds} done", flush=True)
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--plot-only", action="store_true", help="reuse results/downstream_metrics.tsv")
    args = ap.parse_args()
    if args.plot_only:
        m = pd.read_csv(ROOT / "results/downstream_metrics.tsv", sep="\t")
    else:
        datasets = sorted(main_metrics().dataset.unique())
        with Pool(args.jobs) as pool:
            rows = [r for rs in pool.map(run_dataset, datasets, chunksize=1) for r in rs]
        m = pd.DataFrame(rows)
        m.to_csv(ROOT / "results/downstream_metrics.tsv", sep="\t", index=False, float_format="%.4f")

    cols = ["ari", "hvg_jaccard", "marker_jaccard", "doublet_clusters", "residual_doublets"]
    summ = m.groupby("tool")[cols].median()
    summ["doublet_clusters_total"] = m.groupby("tool").doublet_clusters.sum()
    summ = summ.sort_values("ari", ascending=False)
    STATS.mkdir(parents=True, exist_ok=True)
    summ.to_csv(STATS / "downstream_summary.tsv", sep="\t", float_format="%.4f")
    print(summ.round(3).to_string())

    # ------------------------------------------------------------ figure
    import matplotlib.pyplot as plt
    from style import CAT, INK2, SURF, clean, setup
    setup()
    REF = {"reseed": "Oracle, other seed", "none": "No removal"}
    rank_order = list(pd.read_csv(STATS / "ranks.tsv", sep="\t", index_col=0).index)
    order = rank_order[::-1]          # same order as the accuracy figures (best on top)
    panels = [("ari", "Clustering ARI vs oracle"), ("hvg_jaccard", "HVG Jaccard vs oracle"),
              ("marker_jaccard", "Marker Jaccard vs oracle"),
              ("doublet_clusters", "Doublet clusters left (total)")]
    fig, axes = plt.subplots(1, 4, figsize=(15, 5.6), sharey=True)
    for ax, (col, title) in zip(axes, panels):
        for i, t in enumerate(order):
            v = m[m.tool == t][col]
            if col == "doublet_clusters":
                ax.barh(i, v.sum(), height=0.6, color=CAT[0])
                ax.text(v.sum() + 0.3, i, f"{v.sum():d}", va="center", fontsize=8, color=INK2)
            else:
                ax.scatter(v, np.full(len(v), i) + np.linspace(-.18, .18, len(v)), s=10,
                           color="#b9b8b2", zorder=2)
                ax.scatter([v.median()], [i], s=60, color=CAT[0], edgecolor=SURF, lw=1.2, zorder=3)
        for ref, ls in (("none", (0, (3, 3))), ("reseed", (0, (1, 2)))):
            v = m[m.tool == ref][col]
            x = v.sum() if col == "doublet_clusters" else v.median()
            ax.axvline(x, color=INK2, lw=1.2, ls=ls, zorder=1)
            if col == "doublet_clusters" and ref == "none":
                ax.annotate(f"no removal: {x:d}", (x, len(order) - 0.5), xytext=(-4, 0),
                            textcoords="offset points", ha="right", fontsize=8, color=INK2)
        ax.set_title(title, loc="left", fontsize=10, pad=8)
        clean(ax, grid="x")
    axes[0].set_yticks(range(len(order)), [NAMES[t] for t in order])
    fig.suptitle("Downstream impact of each tool's doublet removal, 16 real datasets. Blue = median "
                 "(grey = datasets); dashed = no removal; dotted = oracle re-clustered with another seed",
                 x=0.01, ha="left", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=2)
    fig.savefig(FIG / "downstream.png", dpi=200)
