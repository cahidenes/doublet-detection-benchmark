#!/usr/bin/env bash
# Rebuild every table and figure in the report from the pipeline outputs.
#   needs: results/metrics.tsv, results/scaling_metrics.tsv, results/sensitivity/df_pk.tsv,
#          and (for agreement.py only) the per-droplet predictions in results/pred/.
# Usage: analysis/run_all.sh
set -euo pipefail
cd "$(dirname "$0")"
PY=~/miniforge3/envs/dbl-py/bin/python
for s in accuracy_stats.py agreement.py rate_estimation.py robustness.py scaling.py; do
  echo "== $s"; $PY "$s" > /dev/null
done
$PY make_report_data.py
