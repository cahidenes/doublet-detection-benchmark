#!/usr/bin/env bash
# Final-report benchmark runs (separate metrics files; results/metrics.tsv is untouched).
#   1. realistic expected rate: rate-taking tools with the 10x guideline (0.8% per 1k droplets)
#                                                         -> results/heuristic_metrics.tsv
#   2. simulated doublets incl. abundance-weighted heterotypic pairing and replicate draws r1, r2
#      (scripts/simulate_doublets.py --seed 1/2)          -> results/sim_metrics.tsv
#   3. rate robustness (0.5x/2x) on the seed-0 simulations -> results/sim_robustness_metrics.tsv
#   4. singletCode / Watermelon T47D singlets + simulated doublets (independent validation;
#      scripts/prepare_singletcode.py, then simulate_doublets.py per sample)
#                                                         -> results/singletcode_metrics.tsv
set -uo pipefail
cd "$(dirname "$0")/.."
SMK=(~/miniforge3/envs/dbl-py/bin/snakemake -s workflow/Snakefile --cores 8
     --resources gpu=1 heavy=1 --keep-going --rerun-incomplete)
mkdir -p logs/runs
lst() { ls -d data/processed/$1 | xargs -n1 basename | paste -sd,; }

"${SMK[@]}" --config 'rate_multipliers=["h"]' 'only_tools=[scrublet,scdblfinder,doubletfinder,chord]' \
    metrics_out=results/heuristic_metrics.tsv > logs/runs/heuristic.log 2>&1
echo "heuristic exit $?"
"${SMK[@]}" --config "datasets=[$(lst 'sim-pbmc-2ctrl-dm-*'),$(lst 'sim-pdx-MULTI-*')]" 'rate_multipliers=[1.0]' \
    metrics_out=results/sim_metrics.tsv > logs/runs/sim.log 2>&1
echo "sim exit $?"
"${SMK[@]}" --config "datasets=[$(lst 'sim-wm-*')]" 'rate_multipliers=[1.0]' \
    metrics_out=results/singletcode_metrics.tsv > logs/runs/singletcode.log 2>&1
echo "singletcode exit $?"
S0=$(ls -d data/processed/sim-{pbmc-2ctrl-dm,pdx-MULTI}-*-f{1,0.5} | xargs -n1 basename | paste -sd,)
"${SMK[@]}" --config "datasets=[$S0]" 'only_tools=[scrublet,scdblfinder,doubletfinder,chord]' \
    metrics_out=results/sim_robustness_metrics.tsv > logs/runs/sim_robustness.log 2>&1
echo "sim robustness exit $?"
