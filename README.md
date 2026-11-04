# Benchmarking computational doublet-detection methods for scRNA-seq

A neutral, reproducible benchmark of current doublet detectors for droplet-based single-cell
RNA-seq. It extends Xi & Li (2021) with newer tools (Chord, vaeda), evaluates on the 16 real
datasets with experimentally annotated doublets, and measures accuracy, robustness, speed,
memory and usability. The full plan is in [`docs/proposal.pdf`](docs/proposal.pdf).

## Methods

| Tool | Language | Env | Version | Expected rate input | GPU |
|---|---|---|---|---|---|
| Scrublet | Python | `dbl-py` | 0.2.3 | yes | |
| DoubletDetection | Python | `dbl-py` | 4.3.0.post1 | no | |
| Solo (`scvi.external.SOLO`) | Python | `dbl-py` | scvi-tools 1.4.2 | no | yes |
| vaeda | Python | `dbl-vaeda` | 0.0.30 (GitHub `7b4979c`) | no | yes |
| scDblFinder | R | `dbl-r` | 1.23.4 | yes | |
| computeDoubletDensity (scDblFinder pkg) | R | `dbl-r` | 1.23.4 | no | |
| cxds, bcds, hybrid (scds) | R | `dbl-r` | 1.22.0 | no | |
| DoubletFinder | R | `dbl-r` | 2.0.6 (GitHub `1b244d8`) | yes | |
| Chord | R | `dbl-r` | 2.0.1 (GitHub `00ac9a1`) | yes | |
| Baselines: library size, # expressed genes | Python | `dbl-py` | – | no | |

Every tool runs with its **default parameters**. The expected doublet rate is set to the
labeled rate (1×) in the main run. Only tools that accept a rate are rerun at 0.5×/2× for the
robustness experiment. Wrapper-specific choices are documented at the top of each file in
`methods/`.

## Data

The 16 real datasets assembled by Xi & Li (2021), from Zenodo
[4062232](https://zenodo.org/records/4062232). Doublets are labeled by cell hashing, demuxlet,
MULTI-seq or species mixing.

| Technology | Datasets |
|---|---|
| Cell hashing | `pbmc-ch`, `cline-ch`, `mkidney-ch` |
| Species mixture | `hm-12k`, `hm-6k` |
| demuxlet | `pbmc-1A-dm`, `pbmc-1B-dm`, `pbmc-1C-dm`, `pbmc-2ctrl-dm`, `pbmc-2stim-dm`, `J293t-dm` |
| MULTI-seq | `pdx-MULTI`, `HMEC-orig-MULTI`, `HMEC-rep-MULTI`, `HEK-HMEC-MULTI`, `nuc-MULTI` |

These range from 500 to 26,426 droplets, with doublet rates of 2.5–37.3%. Droplet counts and
rates were checked against the proposal after download. Simulated data (2.7 GB) is optional:
`scripts/download_data.sh --with-sim`.

## Setup

Tested on Linux x86-64 with an NVIDIA GPU. You need about 10 GB for the conda envs and 2 GB
for the real data.

```bash
# 1. Miniforge (conda + mamba)
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3
M=~/miniforge3/bin/mamba

# 2. Environments
$M env create -f envs/py-env.yaml       # dbl-py: Python tools + snakemake
$M env create -f envs/r-env.yaml        # dbl-r: R tools
$M env create -f envs/vaeda-env.yaml    # dbl-vaeda: TensorFlow (conflicts with PyTorch in dbl-py)

# 3. GitHub-only packages, pinned to exact commits
$M run -n dbl-r Rscript envs/install_r_github.R    # DoubletFinder, Chord
$M run -n dbl-vaeda pip install --no-deps git+https://github.com/kostkalab/vaeda@7b4979c12d

# 4. Data
scripts/download_data.sh
```

`install_r_github.R` prints OK/FAIL for each package. Do not install vaeda from PyPI: the only
release there (0.0.1) cannot be imported (see `docs/usability_log.md`).

## Running

The pipeline is a Snakemake workflow:

```
.rds → convert → data/processed/<ds>/   (matrix.mtx.gz, genes, barcodes, labels)
                  → run_tool → results/pred/<tool>/<ds>/s<seed>_r<mult>.tsv
                             + results/bench/...   (wall time, peak RSS)
                  → evaluate → results/eval/...
                  → collect  → results/metrics.tsv
```

```bash
SMK="$HOME/miniforge3/envs/dbl-py/bin/snakemake -s workflow/Snakefile --resources gpu=1 heavy=1"

$SMK --cores 16 --config datasets='[J293t-dm]'      # smoke test (~3 min)
$SMK --cores 8                                       # all datasets (see memory note)
$SMK --cores 8 --config only_tools='[chord,vaeda]'   # subset of tools
```

- Each tool runs with its own env's interpreter (`interpreters:` in `config/config.yaml`), so
  you don't need to activate any env.
- `--resources gpu=1 heavy=1` runs at most one GPU tool and one memory-heavy tool
  (`heavy: 1` in the config) at a time.
- If a tool crashes, its output records `status=failed` with the last lines of the error, and
  the workflow continues. Failures appear in `metrics.tsv`, because they count toward the
  usability metric. Full logs are in `logs/<tool>/<ds>/`.
- **Memory:** some tools need a lot of RAM. DoubletDetection peaked at 18.7 GB and Chord at
  11 GB on only 6.8k droplets. Keep `--cores` low enough that two of these jobs fit in memory
  together (each job uses 4 threads).
- Seeds, the rate multipliers for robustness (`[0.5, 1.0, 2.0]`) and the dataset list are set
  in `config/config.yaml`.

### Output schema

Every wrapper, in R or Python, writes one row per droplet:

| column | meaning |
|---|---|
| `barcode` | droplet barcode |
| `score` | doublet score, **higher = more doublet-like** |
| `call` | 1/0 from the tool's default threshold; `NA` if the tool has none |
| `tool`, `dataset`, `seed`, `expected_rate` | run identifiers |

`scripts/evaluate.py` turns these into AUPRC and AUROC; precision, recall and TNR at 10/20/40%
identification rates; precision, recall and F1 at the default call; the number of droplets
removed; and the error in the estimated doublet rate. `collect` adds the wall time and peak
memory.

## Compatibility shims

Several tools no longer work with current versions of their dependencies. Tool source code is
**never edited**. Each fix lives in a wrapper and is reported as a usability finding:

| Where | Problem | Fix |
|---|---|---|
| `methods/common.R::patch_scds_xgboost` | scds 1.22 (bcds, hybrid, Chord) reads `xgb.cv()$best_iteration` and `$pred`, which xgboost ≥ 2 moved. xgboost 3 also breaks its `xgboost()` call. | Pin xgboost 2.1.4 and copy both fields back to where scds expects them |
| `methods/chord_run.R` | Chord 2.0.1 calls `paramSweep_v3` and `doubletFinder_v3` (renamed in DoubletFinder 2.0.4), reads `@counts` (missing from Seurat v5 assays), and passes `reuse.pANN = FALSE` (misread by DoubletFinder 2.0.6) | Aliases, v3 assays, argument translation |
| `envs/r-env.yaml` | Chord's DESCRIPTION declares no Imports | Its dependencies are installed explicitly |
| `envs/vaeda-env.yaml` | vaeda leaves out tensorflow-probability, seaborn, igraph and leidenalg from its dependencies | Installed explicitly |

Note: scDblFinder also uses xgboost, and its scores changed slightly between xgboost 3.4.2 and
2.1.4. Always report the versions you used.

## Repository layout

```
config/config.yaml   tools, interpreters, datasets, seeds, rate multipliers
envs/                conda env files + pinned GitHub installs
methods/             one wrapper per tool; common.py / common.R hold the shared I/O
scripts/             download_data.sh, convert_rds.R, evaluate.py
workflow/Snakefile   the pipeline
docs/                proposal, usability log (metric 10)
data/, results/, logs/   generated, not tracked in git
```

## Adding a tool

1. Write `methods/<tool>_run.py`, using `methods/common.py`
   (`parse_args`, `load`, `resolve_rate`, `write`), or `methods/<tool>_run.R`, using
   `methods/common.R` (`parse_cli`, `load_counts`, `write_out`).
2. Register it in `config/config.yaml` with `env`, `script`, `rate` and optionally `gpu`.
3. Run the smoke test on `J293t-dm`, then a dataset with a clear signal (`hm-6k`, where every
   tool should reach AUROC above 0.98).
4. Record install and run problems in `docs/usability_log.md`.

## Status (2026-11-04)

- [x] Environments, data download and conversion (all 16 datasets)
- [x] Wrappers for all 13 methods; smoke test passes on `J293t-dm`, `hm-6k`, `pbmc-1A-dm`
- [ ] Full run on all 16 datasets
- [ ] Rate robustness (0.5×/2×), scaling with droplet number, tool agreement
- [ ] Independent doublet simulation from labeled singlets
- [ ] Statistics (Friedman + Wilcoxon/Holm, Cliff's δ), figures, report

## References

- Xi, N. M. & Li, J. J. Benchmarking computational doublet-detection methods for single-cell
  RNA sequencing data. *Cell Systems* 12, 176–194 (2021).
- Zhang, Z. et al. *Cell Genomics* 4, 100592 (2024).
- Tool references are listed in the proposal.
