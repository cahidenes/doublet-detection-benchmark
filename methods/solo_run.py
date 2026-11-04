"""Solo via scvi-tools (the original calico/Solo repo is deprecated in favor of scvi.external.SOLO).
Follows the scvi-tools doublet tutorial: 2000 HVGs (seurat_v3 on counts) -> SCVI -> SOLO.
score = P(doublet) (soft prediction); call = SOLO's hard prediction. No expected-rate input."""
import numpy as np
import anndata as ad
import scanpy as sc
import scvi
from scvi.external import SOLO
from common import load, parse_args, resolve_rate, write

args = parse_args()
X, barcodes, genes, labels = load(args.input)
scvi.settings.seed = args.seed

adata = ad.AnnData(X.astype(np.float32))
adata.obs_names, adata.var_names = barcodes, genes
adata.var_names_make_unique()
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", subset=True)

scvi.model.SCVI.setup_anndata(adata)
vae = scvi.model.SCVI(adata)
vae.train()
solo = SOLO.from_scvi_model(vae)
solo.train()

soft = solo.predict(soft=True)
hard = solo.predict(soft=False)
soft, hard = soft.loc[barcodes], hard.loc[barcodes]
write(args, barcodes, soft["doublet"].to_numpy(), (hard == "doublet").astype(int).to_numpy(),
      resolve_rate(args, labels))
