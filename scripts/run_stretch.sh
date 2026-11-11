#!/usr/bin/env bash
# Stretch-goal benchmark runs (separate metrics files; results/metrics.tsv is untouched).
#   1. simulated doublets (scripts/simulate_doublets.py) -> results/sim_metrics.tsv
#   2. seeds 1 and 2 for every tool with randomness   -> results/stability_metrics.tsv
set -uo pipefail
cd "$(dirname "$0")/.."
SMK=(~/miniforge3/envs/dbl-py/bin/snakemake -s workflow/Snakefile --cores 8
     --resources gpu=1 heavy=1 --keep-going)
SIM=$(ls -d data/processed/sim-* | xargs -n1 basename | paste -sd,)
"${SMK[@]}" --config "datasets=[$SIM]" 'rate_multipliers=[1.0]' metrics_out=results/sim_metrics.tsv \
    > logs/sim_run.log 2>&1
echo "sim exit $?"; tail -2 logs/sim_run.log
"${SMK[@]}" --config 'seeds=[1,2]' 'rate_multipliers=[1.0]' \
    'only_tools=[scrublet,doubletdetection,solo,vaeda,scdblfinder,dbldensity,bcds,hybrid,doubletfinder,chord]' \
    metrics_out=results/stability_metrics.tsv > logs/stability_run.log 2>&1
echo "stability exit $?"; tail -2 logs/stability_run.log
