# Usability log (metric 10)

Fill this in while installing and running each tool; don't rely on memory later.
Score each 1–5 at the end: installation, documentation, ease of use, failure rate.

| Tool | Version | Language | Install time | Install problems | Docs quality | Failed runs (dataset: error) | Notes |
|------|---------|----------|--------------|------------------|--------------|------------------------------|-------|
| Scrublet | 0.2.3 | Python | | | | | |
| DoubletFinder | 2.0.6 @1b244d8 | R | | | | | |
| scDblFinder | 1.23.4 (bioconda) | R | | | | | |
| scds (cxds/bcds/hybrid) | 1.22.0 (bioconda) | R | | bcds/hybrid incompatible with current xgboost: 3.x breaks `xgboost()` call; 2.x moved `xgb.cv` fields `best_iteration`, `pred`. Pinned xgboost 2.1.4 + shim (methods/common.R). | | | cxds unaffected (no xgboost). |
| Solo (scvi-tools) | 1.4.2 | Python | | | | | Original repo deprecated → scvi.external.SOLO. On J293t-dm (500 droplets) the classifier stopped at epoch 31 with no improvement after epoch 1 → called 467/500 as doublets. |
| DoubletDetection | 4.3.0.post1 | Python | | | | | |
| vaeda | 0.0.30 @7b4979c (GitHub) | Python | | PyPI 0.0.1 unusable: bogus deps (random, math, pathlib, sklearn); absolute intra-package imports; plot_results.py missing from wheel. GitHub version fixes imports but omits tensorflow_probability from install_requires and imports seaborn undeclared. | | | No expected-rate parameter → excluded from metric 5. Needs python-igraph + leidenalg (undeclared). |
| Chord | 2.0.1 @00ac9a1 | R | | DESCRIPTION has no Imports field → install_github fails (adabag, gbm missing); deps must be installed by hand. | | | Unmaintained since 2022: incompatible with DoubletFinder ≥2.0.4 (`_v3` names, `reuse.pANN=FALSE`), Seurat v5 assays, and (via bcds) xgboost ≥2. Runs only with 4 shims (methods/chord_run.R). Writes PDFs/.robj/CSVs into the working dir. |

## Dated notes
- 2026-11-03: project scaffold created.
- 2026-11-03: conda envs dbl-r, dbl-py created without errors (Seurat 5.5.1, snakemake 8.30.0).
- 2026-11-03: vaeda from PyPI (only release: 0.0.1) fails at `pip install` (declares stdlib modules as deps);
  with --no-deps it fails at import (`from vae import ...`, missing `plot_results`). Needs GitHub install.
- 2026-11-03: Smoke test on J293t-dm passed for libsize, ngenes, Scrublet, scDblFinder, computeDoubletDensity.
  Scrublet's automatic threshold called only 5/500 doublets (true 42) despite expected_doublet_rate=0.084.
- 2026-11-04: DoubletFinder 2.0.6 installed from GitHub without problems.
- 2026-11-04: Chord install failed: `there is no package called 'adabag'`. Cause: DESCRIPTION declares no
  Imports although NAMESPACE imports 9 packages. Fixed by installing r-adabag, r-gbm from conda-forge.
  Code also references xgboost/lightgbm (optional boosting methods; default is gbm).
- 2026-11-04: vaeda 0.0.30 import failed: no tensorflow_probability (not in setup.py), seaborn imported
  undeclared. Fixed with tensorflow-probability==0.23.0 + seaborn. TF 2.15.1 sees 1 GPU.
- 2026-11-04: Wrappers written for all tools. Smoke test on J293t-dm surfaced:
  - vaeda: `ImportError: Please install the igraph package` (sc.tl.leiden) → added python-igraph, leidenalg.
  - bcds/hybrid/Chord with xgboost 3.4.2: `argument "y" is missing` (xgboost 3 API) → pinned 2.1.4.
  - with 2.1.4: `vector size cannot be infinite` / `'response' must have two levels` — xgb.cv result fields
    moved in xgboost 2 → shim restoring `best_iteration` and `pred` (methods/common.R).
  - Chord: `paramSweep_v3`/`doubletFinder_v3` missing (DoubletFinder 2.0.4 renamed them), `@counts`
    slot missing on Seurat v5 assays, `reuse.pANN=FALSE` misread by DoubletFinder 2.0.6 → shims.
  - scDblFinder result on J293t-dm changed with xgboost 3.4.2 → 2.1.4 (AUPRC 0.222 → 0.211, calls 38 → 20):
    report the xgboost version used.

### 2026-11-04 — DoubletFinder pK sensitivity (to revisit in the report)
`analysis/df_pk_sensitivity.sh` → `results/sensitivity/df_pk.tsv`. AUPRC falls monotonically
with pK on pbmc-ch, pbmc-2stim-dm, pbmc-1C-dm (e.g. pbmc-ch 0.50 at pK=0.005 vs 0.29 at 0.3).
The automatic BCmvn choice picked pK=0.3 (pbmc-ch) and 0.13 (pbmc-2stim-dm), roughly halving
AUPRC there. DoubletFinder's low rank is mainly its pK-selection step. Main ranking keeps the
default procedure (picking pK with labels would be tuning on the test set).

### 2026-11-04 — Rate robustness (metric 5)
Chord fails at 2× on HMEC-rep-MULTI and mkidney-ch (2× rate = 0.62 / 0.75):
`overkillDB2 -> sample(): cannot take a sample larger than the population`. Chord cannot
accept an expected rate above ~0.5. Recorded as failed; excluded from paired medians.
