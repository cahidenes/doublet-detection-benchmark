# Usability scores (metric 10)

**Provisional:** these scores are a proposal drawn from the facts in
[`usability_log.md`](usability_log.md). The *Documentation* column is the most
subjective and should be checked by a second reader.

Each criterion is scored from 1 (poor) to 5 (good). Every score cites the evidence behind it.

| Criterion | 5 means | 1 means |
|---|---|---|
| **Installation** | One conda/pip/Bioconductor command | Broken package; manual dependency hunting; a separate environment |
| **Compatibility** | Runs on current dependency versions as released | Needs several wrapper shims to run at all |
| **Ease of use** | One call with sensible defaults | Multi-step workflow with decisions the user must make (parameters, preprocessing) |
| **Reliability** | No failures; sensible default calls | Crashes in normal settings or default calls far off |
| **Documentation** | Vignette or docs site with a worked example | README only, or outdated |

| Tool | Install | Compat. | Ease | Reliab. | Docs | **Mean** | Evidence |
|---|---|---|---|---|---|---|---|
| scDblFinder | 5 | 4 | 5 | 5 | 5 | **4.8** | Bioconductor; one call. Scores change with the xgboost version (AUPRC 0.222 → 0.211 on J293t-dm). |
| DoubletDetection | 5 | 5 | 5 | 4 | 4 | **4.6** | pip; one call. Peak memory up to 19 GB. |
| Scrublet | 5 | 5 | 5 | 3 | 4 | **4.4** | pip; one call. The automatic threshold calls about 0.23× the true number of doublets, even when given the true rate. |
| DoubletDensity | 5 | 5 | 5 | 3 | 4 | **4.4** | Part of scDblFinder. Scores only, no default call. |
| cxds | 5 | 5 | 4 | 4 | 4 | **4.4** | Bioconductor; no xgboost dependency. Over-calls slightly (median est/true 1.25). |
| bcds / hybrid | 5 | 2 | 4 | 4 | 4 | **3.8** | Broken with xgboost ≥ 2: needs xgboost pinned to 2.1.4 plus the field shim (`common.R`). |
| Solo | 4 | 5 | 3 | 3 | 4 | **3.8** | Now part of scvi-tools (the original repo is deprecated). Needs an SCVI model and HVG choices first; GPU advisable. Over-calls (median est/true 1.9; 467/500 on J293t-dm). |
| vaeda | 1 | 3 | 4 | 4 | 2 | **2.8** | The PyPI release cannot be imported. The GitHub version leaves out 4 dependencies, and TensorFlow forces a separate env. No expected-rate input. |
| DoubletFinder | 3 | 3 | 2 | 3 | 3 | **2.8** | GitHub only. Functions renamed in 2.0.4. Multi-step workflow (Seurat preprocessing → pK sweep → nExp → optional homotypic adjustment). The automatic pK choice can halve AUPRC. |
| Chord | 1 | 1 | 3 | 2 | 2 | **1.8** | GitHub only, with no declared dependencies; unmaintained since 2022. Needs 4 shims (DoubletFinder renames, Seurat v5, xgboost). Crashes when the expected rate is above about 0.5, and writes files into the working directory. |

Baselines (library size, # genes) are not scored.
