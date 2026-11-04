"""vaeda 0.0.30 (GitHub) with package defaults; raw counts in adata.X.
score = adata.obs['vaeda_scores'], call = adata.obs['vaeda_calls'].
Runs in the dbl-vaeda env (TensorFlow). No expected-rate input."""
import anndata as ad
import numpy as np
import vaeda
from common import load, parse_args, resolve_rate, write

args = parse_args()
X, barcodes, genes, labels = load(args.input)
adata = ad.AnnData(X.astype(np.float32))
adata.obs_names, adata.var_names = barcodes, genes
adata.var_names_make_unique()

adata = vaeda.vaeda(adata, seed=args.seed)
obs = adata.obs.loc[barcodes]
write(args, barcodes, obs["vaeda_scores"].to_numpy(),
      (obs["vaeda_calls"] == "doublet").astype(int).to_numpy(), resolve_rate(args, labels))
