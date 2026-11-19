"""Label noise on the species mixtures: are a tool's "false positives" unlabeled doublets?

In hm-12k / hm-6k only human-mouse doublets are labeled; human-human and mouse-mouse doublets
are labeled singlets. If a tool's top-scored labeled singlets are such within-species doublets,
they should carry about twice the RNA of a singlet of their own species.

For each tool: take its top-k droplets (k = labeled doublets), keep the labeled singlets among
them (false positives), assign each to a species, and divide
its library size by the median library size of labeled singlets of that species.
A ratio near 2 = likely unlabeled homotypic doublet; near 1 = genuine false positive.

The processed matrices of Xi & Li map mouse genes to human orthologs, so species cannot be read
from gene IDs. Both datasets mix two cell lines (HEK293T, NIH3T3), so species = 2-means
cluster on the expression PCA of labeled singlets (each droplet assigned to the nearest centre).

Usage: python analysis/label_noise.py
Writes: results/stats/label_noise.tsv
"""
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
from sklearn.cluster import KMeans

from style import NAMES, ROOT, STATS, main_metrics

warnings.filterwarnings("ignore", category=RuntimeWarning)
rows = []
for ds in ("hm-12k", "hm-6k"):
    d = ROOT / "data/processed" / ds
    X = scipy.io.mmread(d / "matrix.mtx.gz").tocsc()            # genes x droplets
    genes = np.array((d / "genes.tsv").read_text().split())
    bc = (d / "barcodes.tsv").read_text().split()
    lib = np.asarray(X.sum(0)).ravel()
    lab = pd.read_csv(d / "labels.tsv", sep="\t").set_index("barcode").loc[bc, "is_doublet"].to_numpy()
    sing = lab == 0
    ad = sc.AnnData(X.T.tocsr().astype(np.float32))
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, n_top_genes=2000)
    sc.pp.pca(ad, n_comps=20, mask_var="highly_variable", random_state=0)
    km = KMeans(2, n_init=10, random_state=0).fit(ad.obsm["X_pca"][sing])
    sp = np.array(["A", "B"])[km.predict(ad.obsm["X_pca"])]
    med = {s: np.median(lib[sing & (sp == s)]) for s in ("A", "B")}
    rel = lib / np.vectorize(med.get)(sp)                       # library size / own-species median
    k = int(lab.sum())
    # expected unlabeled within-species doublets: with human fraction q among cells, a random
    # doublet is cross-species with prob 2q(1-q) and within-species with prob q^2 + (1-q)^2
    q = (sp[sing] == "A").mean()
    exp_within = k * (q ** 2 + (1 - q) ** 2) / (2 * q * (1 - q))
    print(f"{ds}: cluster-A fraction {q:.2f}, labeled (cross-species) doublets {k}, "
          f"expected unlabeled within-species doublets ~{exp_within:.0f}")
    for tool in main_metrics().tool.unique():
        p = pd.read_csv(ROOT / f"results/pred/{tool}/{ds}/s0_r1.0.tsv", sep="\t").set_index("barcode")
        s = p.loc[bc, "score"].to_numpy()
        top = np.argsort(-s, kind="stable")[:k]
        fp = top[lab[top] == 0]
        rows.append({"dataset": ds, "tool": tool, "k": k, "expected_unlabeled": round(exp_within),
                     "n_fp": len(fp),
                     "precision_topk": 1 - len(fp) / k,
                     "fp_median_rel_lib": np.median(rel[fp]) if len(fp) else np.nan,
                     "fp_frac_rel_lib_gt_1.6": (rel[fp] > 1.6).mean() if len(fp) else np.nan,
                     "singlet_frac_rel_lib_gt_1.6": (rel[sing] > 1.6).mean()})

out = pd.DataFrame(rows).sort_values(["dataset", "precision_topk"], ascending=[True, False])
STATS.mkdir(parents=True, exist_ok=True)
out.to_csv(STATS / "label_noise.tsv", sep="\t", index=False, float_format="%.4f")
out["tool"] = out.tool.map(NAMES)
print(out.round(3).to_string(index=False))
