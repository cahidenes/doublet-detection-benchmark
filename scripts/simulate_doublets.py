"""Stretch goal: synthetic doublets built from labeled singlets, with known composition.

From one processed real dataset:
  1. keep only labeled singlets; draw N of them;
  2. cluster them (scanpy: normalize, log1p, 2000 HVG, PCA 30, kNN 15, Leiden res 0.5)
     to get cell types;
  3. draw D = rate*N/(1+rate) pairs of parent cells without replacement and remove the
     parents, so the final doublet rate is exactly `rate` and no parent is also a singlet;
  4. doublet counts = Binomial(parent1 + parent2, f): f = 1 keeps the summed library
     (about 2x a singlet); f = 0.5 gives singlet-like depth, so library size is uninformative.

Pairing scenarios:
  random        any two cells (the natural mix)
  heterotypic   parents from different clusters
  homotypic     parents from the same cluster (hard: the profile looks like one cell type)

This is our own procedure, not any tool's, but like most tools it is additive (sum of two
cells), so it is "independent" only in implementation, not in its generative assumption.

Usage: python scripts/simulate_doublets.py <base_ds> [--n 6000] [--rate 0.1] [--seed 0]
Writes data/processed/sim-<base>-<pairing>-f<f>/ for all 3 pairings x f in {1, 0.5}.
"""
import argparse
import gzip
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp

p = argparse.ArgumentParser()
p.add_argument("base"); p.add_argument("--n", type=int, default=6000)
p.add_argument("--rate", type=float, default=0.1); p.add_argument("--seed", type=int, default=0)
a = p.parse_args()
rng = np.random.default_rng(a.seed)

src = Path("data/processed") / a.base
X = scipy.io.mmread(src / "matrix.mtx.gz").tocsc().T.tocsr()      # droplets x genes
lab = pd.read_csv(src / "labels.tsv", sep="\t")
sing = np.flatnonzero(lab.is_doublet.values == 0)
pick = np.sort(rng.choice(sing, size=min(a.n, len(sing)), replace=False))
X = X[pick].astype(np.int64); bc = lab.barcode.values[pick]; N = X.shape[0]

# cell types of the singlets
A = ad.AnnData(X.astype(np.float32))
sc.pp.filter_genes(A, min_cells=3)
sc.pp.normalize_total(A); sc.pp.log1p(A)
sc.pp.highly_variable_genes(A, n_top_genes=2000)
sc.pp.pca(A, n_comps=30, mask_var="highly_variable", random_state=a.seed)
sc.pp.neighbors(A, n_neighbors=15, random_state=a.seed)
sc.tl.leiden(A, resolution=0.5, random_state=a.seed, flavor="igraph", n_iterations=2)
cl = A.obs.leiden.astype(int).values
print(f"{a.base}: {N} singlets, {cl.max() + 1} clusters, sizes {np.bincount(cl).tolist()}")

D = int(round(a.rate * N / (1 + a.rate)))


def draw_pairs(mode):
    """D disjoint parent pairs under the pairing rule."""
    avail = np.ones(N, bool); pairs = []
    order = rng.permutation(N)
    by_cl = {c: list(rng.permutation(np.flatnonzero(cl == c))) for c in np.unique(cl)}
    if mode == "random":
        idx = order[:2 * D]; return idx[:D], idx[D:]
    for i in order:
        if len(pairs) == D:
            break
        if not avail[i]:
            continue
        if mode == "homotypic":
            cand = [j for j in by_cl[cl[i]] if avail[j] and j != i]
        else:
            # heterotypic: take the next available cell of a random other cluster
            others = [c for c in by_cl if c != cl[i] and any(avail[j] for j in by_cl[c])]
            if not others:
                continue
            c = others[rng.integers(len(others))]
            cand = [j for j in by_cl[c] if avail[j]]
        if not cand:
            continue
        j = cand[0]; avail[i] = avail[j] = False; pairs.append((i, j))
    assert len(pairs) == D, f"{mode}: only {len(pairs)} of {D} pairs"
    p1, p2 = map(np.array, zip(*pairs)); return p1, p2


for mode in ("random", "heterotypic", "homotypic"):
    p1, p2 = draw_pairs(mode)
    parents = np.concatenate([p1, p2])
    keep = np.setdiff1d(np.arange(N), parents)
    summed = (X[p1] + X[p2]).tocsr()
    for f in (1.0, 0.5):
        dbl = summed.copy()
        if f < 1:
            dbl.data = rng.binomial(dbl.data, f)
            dbl.eliminate_zeros()
        M = sp.vstack([X[keep], dbl]).tocsr()
        is_d = np.r_[np.zeros(len(keep), int), np.ones(D, int)]
        names = np.r_[bc[keep], [f"dbl{k:05d}_{bc[i]}+{bc[j]}" for k, (i, j) in enumerate(zip(p1, p2))]]
        perm = rng.permutation(M.shape[0])                          # shuffle droplet order
        M, is_d, names = M[perm], is_d[perm], names[perm]
        out = Path("data/processed") / f"sim-{a.base}-{mode}-f{f:g}"
        out.mkdir(parents=True, exist_ok=True)
        with gzip.open(out / "matrix.mtx.gz", "wb") as fh:
            scipy.io.mmwrite(fh, M.T.tocoo())                         # genes x droplets
        shutil.copy(src / "genes.tsv", out / "genes.tsv")
        (out / "barcodes.tsv").write_text("\n".join(names) + "\n")
        pd.DataFrame({"barcode": names, "is_doublet": is_d}).to_csv(out / "labels.tsv", sep="\t", index=False)
        same = float((cl[p1] == cl[p2]).mean())
        print(f"  {out.name}: {M.shape[0]} droplets, {D} doublets, same-cluster pairs {same:.2f}")
