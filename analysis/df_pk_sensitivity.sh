#!/usr/bin/env bash
# DoubletFinder pK sensitivity: rerun with a fixed pK grid on all 16 real datasets.
# Not part of the main ranking, which keeps the tool's default (BCmvn) procedure.
# The pK chosen by BCmvn in the main run is read from logs/doubletfinder/<ds>/s0_r1.0.log.
# Output: results/sensitivity/df_pk${TAG}.tsv
# Options (env): DATASETS, PKS (pK grid), TAG (output suffix), e.g. the label-free pK test on
# singletCode data: DATASETS="$(ls data/processed | grep ^sim-wm-)" PKS="0.005 0.01" TAG=_singletcode
set -euo pipefail
cd "$(dirname "$0")/.."
R=~/miniforge3/envs/dbl-r/bin/Rscript
PY=~/miniforge3/envs/dbl-py/bin/python
OUT=results/sensitivity; mkdir -p $OUT/pred $OUT/eval logs/sensitivity
DATASETS=${DATASETS:-$(ls data/processed | grep -v -e "^scale-" -e "^sim-" -e "^sc-" -e "^wm-")}
for ds in $DATASETS; do
  rate=$($PY -c "import pandas as pd; print(pd.read_csv('data/processed/$ds/labels.tsv', sep='\t').is_doublet.mean())")
  for pk in ${PKS:-0.005 0.01 0.05 0.09 0.2 0.3}; do
    f=$OUT/pred/df_${ds}_pk${pk}.tsv
    [[ -f $f ]] || DF_PK=$pk $R methods/doubletfinder_run.R data/processed/$ds $f $ds doubletfinder 0 $rate \
        > logs/sensitivity/df_${ds}_pk${pk}.log 2>&1
    $PY scripts/evaluate.py $f data/processed/$ds/labels.tsv $OUT/eval/df_${ds}_pk${pk}.tsv
    echo "$ds pK=$pk done"
  done
done
DATASETS="$DATASETS" TAG="${TAG:-}" $PY - <<'EOF'
import glob, os, re, pandas as pd
keep = set(os.environ["DATASETS"].split())
rows = []
for f in sorted(glob.glob("results/sensitivity/eval/df_*.tsv")):
    r = pd.read_csv(f, sep="\t").iloc[0]
    if r.dataset not in keep:
        continue
    rows.append({"dataset": r.dataset, "pK": float(re.search(r"pk([\d.]+)\.tsv", f).group(1)),
                 "auprc": r.auprc, "auroc": r.auroc})
d = pd.DataFrame(rows).sort_values(["dataset", "pK"])
d.to_csv(f"results/sensitivity/df_pk{os.environ['TAG']}.tsv", sep="\t", index=False)
print(d.pivot(index="pK", columns="dataset", values="auprc").round(3).to_string())
EOF
