#!/usr/bin/env bash
# Rebuild every table and figure in the report from the pipeline outputs.
#   needs: results/metrics.tsv, results/scaling_metrics.tsv, results/sensitivity/df_pk.tsv,
#          results/sim_metrics.tsv, results/stability_metrics.tsv (scripts/run_stretch.sh),
#          results/{heuristic,sim_robustness,singletcode}_metrics.tsv (scripts/run_final.sh),
#          results/reference/xili2021_real.tsv, and
#          (agreement, stability, downstream) the per-droplet predictions in results/pred/ and
#          data/processed/. downstream.py takes ~10 min.
# Usage: analysis/run_all.sh
set -euo pipefail
cd "$(dirname "$0")"
PY=~/miniforge3/envs/dbl-py/bin/python
for s in accuracy_stats.py agreement.py rate_estimation.py robustness.py scaling.py \
         simulation.py stability.py downstream.py label_noise.py lineage.py ensemble.py \
         xili_comparison.py rate_settings.py singletcode.py; do
  echo "== $s"; $PY "$s" > /dev/null
done
$PY make_report_data.py
