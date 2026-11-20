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
rates were checked against the proposal after download. Xi & Li's own simulated data (2.7 GB)
is not used: `scripts/download_data.sh --with-sim`.

**Simulated doublets** (`scripts/simulate_doublets.py`): synthetic doublets made from the
labeled singlets of `pbmc-2ctrl-dm` and `pdx-MULTI` (random, homotypic, heterotypic and
abundance-weighted heterotypic pairing; summed or singlet-depth library; 3 draws).

**Independent validation: singletCode** (Zhang et al. 2024, figshare
[25478680](https://figshare.com/articles/dataset/25478680)). Six Watermelon-barcoded T47D
samples whose true singlets were identified from lineage barcodes. No tool was tuned on them.
`scripts/download_data.sh --with-singletcode` fetches them (5.5 GB zip; about 0.7 GB is kept,
and the zip can be deleted afterwards). `scripts/prepare_singletcode.py` writes the
singlets-only bases `data/processed/wm-T47D-*`, and `simulate_doublets.py` adds the doublets
(`sim-wm-*`, 36 datasets).

**Reference results:** Xi & Li's per-dataset AUPRC/AUROC (Zenodo 4444303, `result.xlsx`,
downloaded by default) → `results/reference/xili2021_real.tsv`
(`scripts/prepare_xili_reference.py`).

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
# (to get the exact versions used for the results: envs/lock/<env>.yaml instead)

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
  in `config/config.yaml`. The multiplier `h` gives the rate-taking tools the 10x loading
  guideline (0.8% per 1,000 droplets) instead of the labeled rate.

Beyond the main run, the experiments write separate metrics files (`results/metrics.tsv` is
never touched):

| Script | Runs | Output |
|---|---|---|
| `scripts/run_stretch.sh` | seeds 1–2 for every tool with randomness | `results/stability_metrics.tsv` |
| `scripts/run_final.sh` | realistic rate (`h`); simulated doublets; singletCode; 0.5×/1×/2× on simulations | `results/{heuristic,sim,singletcode,sim_robustness}_metrics.tsv` |
| `analysis/df_pk_sensitivity.sh` | DoubletFinder with fixed pK (options `DATASETS`, `PKS`, `TAG`) | `results/sensitivity/df_pk*.tsv` |
| `analysis/solo_original.py` | the original Solo (solo-sc 0.6), in its own env `dbl-solo0g` (setup in the docstring) | `results/sensitivity/solo0/` |

`run_final.sh` passes `--rerun-incomplete`, so it can simply be restarted after an
interruption (run `snakemake --unlock` first if a crash left a lock).

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

**Evidence that the shims preserve each tool's behaviour.** On the same 16 datasets, bcds,
hybrid and cxds reproduce Xi & Li's (2021) per-dataset AUPRC to within 0.005 on average
(per-dataset Spearman 0.98–1.00), and so does scDblFinder (report §3.12). Chord has no earlier
reference; its components (cxds, bcds, DoubletFinder) are checked individually, and its
inverted scores on some datasets come from its own training design, not from a shim
(report §3.10).

## Repository layout

```
config/config.yaml   tools, interpreters, datasets, seeds, rate multipliers
envs/                conda env files + pinned GitHub installs; envs/lock/ = exact versions used
analysis/            statistics and figures (run_all.sh rebuilds everything)
report/              final report (report.md) and every number behind it (data.md)
methods/             one wrapper per tool; common.py / common.R hold the shared I/O
scripts/             download_data.sh, convert_rds.R, evaluate.py, subsample.py,
                     simulate_doublets.py, prepare_singletcode.py, prepare_xili_reference.py,
                     run_stretch.sh (seeds), run_final.sh (final-report runs)
workflow/Snakefile   the pipeline
docs/                proposal, usability log (metric 10)
results/             *_metrics.tsv, stats/, figures/, bench/, sensitivity/, reference/ are
                     tracked; pred/ and eval/ are not
data/, logs/         generated, not tracked in git
```

## Adding a tool

1. Write `methods/<tool>_run.py`, using `methods/common.py`
   (`parse_args`, `load`, `resolve_rate`, `write`), or `methods/<tool>_run.R`, using
   `methods/common.R` (`parse_cli`, `load_counts`, `write_out`).
2. Register it in `config/config.yaml` with `env`, `script`, `rate` and optionally `gpu`.
3. Run the smoke test on `J293t-dm`, then a dataset with a clear signal (`hm-6k`, where every
   tool should reach AUROC above 0.98).
4. Record install and run problems in `docs/usability_log.md`.

## Status (2026-11-20)

- [x] Environments, data download and conversion (all 16 datasets)
- [x] Wrappers for all 13 methods; full run on all 16 datasets (208/208 ok)
- [x] Rate robustness (0.5×/2×) and a realistic-rate setting (10x guideline, 64/64 ok)
- [x] Scaling with droplet number (1k–26k): `scripts/subsample.py`, `analysis/scaling.py`
- [x] Tool agreement, rate estimation, usability scores
- [x] Statistics: Friedman + Wilcoxon/Holm, Cliff's δ, paired effect sizes, F1 tests,
      study-averaged sensitivity (`analysis/accuracy_stats.py`)
- [x] Simulated doublets, 4 pairings × 2 depths × 3 draws (624/624 ok): `analysis/simulation.py`
- [x] Stability across seeds 0–2 (320/320 ok): `analysis/stability.py`
- [x] Downstream impact with noise-floor references: `analysis/downstream.py`
- [x] Independent validation on singletCode (468/468 ok): `analysis/singletcode.py`
- [x] Comparison with Xi & Li (2021) per dataset; DoubletFinder pK sweep on all 16 datasets;
      original Solo: `analysis/xili_comparison.py`, `analysis/solo_original.py`
- [x] Label noise, lineage co-expression, ensembles: `analysis/label_noise.py`, `lineage.py`,
      `ensemble.py`
- [x] Final report: [`report/report.md`](report/report.md), all numbers in [`report/data.md`](report/data.md)
- [ ] Zenodo deposit of `results/pred/` (needed to rerun the analyses that use per-droplet
      predictions without rerunning the tools)
- [ ] Scaling to 100k–200k droplets (stated as a limitation)

## Review issues (2026-11-17) and how they were resolved

### Fix first

- [x] **Independent validation of scDblFinder.** On singletCode data (6 T47D samples, 36
      datasets with simulated doublets; no tool was tuned on them), scDblFinder still ranks
      first overall (report §3.10). Its lead depends on doublet depth: with singlet-sized
      doublets, Scrublet, vaeda and DoubletDensity are better. This also reconciles our
      result with Zhang et al. (2024), who used averaged doublets and found scDblFinder last.
      The abstract states both.
- [x] **Label noise is not tool-neutral.** Confirmed on `hm-12k`/`hm-6k`: about as many
      same-species doublets are expected as labeled ones, and Solo's and bcds' top "false
      positives" have a median 1.9–2.0× the library of their own species (85–90% above 1.6×,
      vs 8.5% of singlets). Lineage co-expression shows that the PBMC labels miss about a third
      of heterotypic doublets (report §3.12–3.13).
- [x] **Per-dataset comparison with Xi & Li (2021).** Six of eight comparable tools reproduce
      (within 0.005 mean AUPRC for scds, Scrublet, scDblFinder). DoubletFinder's gap (0.42 vs
      0.54) is about half pK: BCmvn's choice (logged for all 16 datasets) costs 0.06, and the
      best fixed pK per dataset still stays 0.05 below. The pK claim in the report was softened
      accordingly. Solo's gap on the species mixtures is not the HVG choice; see the original
      Solo run (report §3.12).
- [x] **Realistic expected-rate setting.** With the 10x guideline, scDblFinder's F1 is
      unchanged (0.49 → 0.50); DoubletFinder (0.37 → 0.32), Chord and Scrublet lose. The
      labeled rate flattered them, not scDblFinder (report §3.3).
- [x] **Paired effect sizes and study dependence.** Matched-pairs rank-biserial and bootstrap
      CIs of the median paired difference added (scDblFinder vs vaeda: r = 0.40, 0.015
      [−0.011, 0.026]). Averaging within study first gives the same ranking (p = 3.5·10⁻⁸).

### Smaller fixes

- [x] Downstream noise floor: "remove k random droplets" and "oracle minus k random
      singlets" references. HVG differences between most tools are within that noise.
- [x] Simulation replicates (3 draws; SD mostly 0.01–0.04) and the real doublet/singlet
      library ratio (median 1.66, range 0.89–1.91). Abundance-weighted heterotypic pairing
      added.
- [x] DoubletDensity dropped from the F1 Friedman test (10 callers; p = 3.1·10⁻⁵, 4 of 45
      pairs significant).
- [ ] Milestone table in `report/milestone2_design.docx`: singletCode ticks and "sel" (to be
      edited by hand).
- [x] References: Chord is *Communications Biology* 5, 510 (2022). Zhang et al. (2024) cited
      and discussed. All references checked against the proposal.
- [ ] Reproducibility: deposit `results/pred/` on Zenodo.
- [x] Compatibility shims: evidence that they preserve behaviour added above.

### Worth adding

1. [x] **Simple ensemble** with leave-one-dataset-out selection: no gain over scDblFinder
   (0.561 vs 0.563). Chord's GBM is worse than a rank average of its own components.
2. [ ] **Scaling to 100k–200k droplets**: not done; stated as a limitation.
3. [x] **Lineage co-expression (dimension 15)**: report §3.13.
4. [x] **Label-free pK rule for DoubletFinder**, tested on singletCode data (report §3.10).

## References

- Xi, N. M. & Li, J. J. Benchmarking computational doublet-detection methods for single-cell
  RNA sequencing data. *Cell Systems* 12, 176–194 (2021).
- Zhang, Z., Melzer, M. E. et al. Synthetic DNA barcodes identify singlets in scRNA-seq
  datasets and evaluate doublet algorithms. *Cell Genomics* 4, 100592 (2024).
- Tool references are listed in the proposal.
