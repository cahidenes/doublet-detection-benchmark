"""Sensitivity analysis: the original Solo (solo-sc 0.6 on scvi 0.6.8), as Xi & Li (2021) ran it.

Why: on the species mixes hm-12k / hm-6k, Xi & Li report Solo AUPRC 0.995 / 0.972; the
scvi-tools reimplementation (scvi.external.SOLO, main run) reaches 0.56 / 0.61, and changing
its HVG count does not close the gap (results/sensitivity/solo_hvg.tsv). If the original
package reproduces their values on our matrices, the gap is a change in the tool, not in the
data or the labels.

Setup: separate env dbl-solo0g: python 3.9, scvi 0.6.8 (the last version before scvi-tools),
solo-sc 0.6 (copied into site-packages: its setup.py no longer builds), torch 2.7.1+cu128 (the
torch 1.x of 2020 has no kernels for current GPUs, and on CPU one dataset takes days),
anndata 0.7.8, scanpy 1.8.2, numpy 1.23.5, pandas 1.3.5, numba 0.56.4, matplotlib 3.5.3,
get_version 2.1, umap-learn 0.5.3. Solo's code is unchanged; only torch is newer than in 2020.
Parameters: the repository's solo_params_example.json at tag 0.6
(config/solo0_params.json). All genes, CLI defaults otherwise (doublet ratio 2, multinomial
doublets). Score = softmax_scores.npy, the file Xi & Li evaluated. The CLI has no random-seed
option, so runs are not seeded.

Usage: ~/miniforge3/envs/dbl-solo0g/bin/python analysis/solo_original.py hm-12k hm-6k ...
Writes: results/sensitivity/solo0/<dataset>.tsv (standard prediction format)
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/sensitivity/solo0"
OUT.mkdir(parents=True, exist_ok=True)
PARAMS = ROOT / "config/solo0_params.json"

for ds in sys.argv[1:]:
    if (OUT / f"{ds}.tsv").exists():                # resumable: skip finished datasets
        print(ds, "exists, skipped", flush=True)
        continue
    d = ROOT / "data/processed" / ds
    X = sp.csr_matrix(scipy.io.mmread(str(d / "matrix.mtx.gz")).T).astype(np.float32)
    barcodes = (d / "barcodes.tsv").read_text().split()
    genes = (d / "genes.tsv").read_text().split()
    labels = pd.read_csv(d / "labels.tsv", sep="\t").set_index("barcode")["is_doublet"]
    with tempfile.TemporaryDirectory() as tmp:
        a = anndata.AnnData(X, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=genes))
        a.var_names_make_unique()
        a.write(f"{tmp}/data.h5ad")
        cmd = ("import sys; from solo.solo import main; "
               f"sys.argv = ['solo', '-l', '-o', '{tmp}/out', '{PARAMS}', '{tmp}/data.h5ad']; main()")
        with open(ROOT / f"logs/sensitivity/solo0_{ds}.log", "w") as log:
            subprocess.run([sys.executable, "-c", cmd], check=True, stdout=log, stderr=subprocess.STDOUT,
                           cwd=tmp)
        score = np.load(f"{tmp}/out/softmax_scores.npy")
    assert len(score) == len(barcodes)
    pd.DataFrame({"barcode": barcodes, "score": score, "call": pd.NA, "tool": "solo0", "dataset": ds,
                  "seed": 0, "expected_rate": float(labels.mean())}).to_csv(OUT / f"{ds}.tsv", sep="\t",
                                                                            index=False)
    print(ds, "done", flush=True)
