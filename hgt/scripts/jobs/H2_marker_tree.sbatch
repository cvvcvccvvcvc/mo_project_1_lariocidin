#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-hgt-H2-marker
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
set -euo pipefail
: "${HGT_ROOT:?Isolated H2 root required}"
export PYTHONDONTWRITEBYTECODE=1
: "${HGT_BIN:?Set HGT_BIN to the pinned environment bin directory}"
export HGT_BIN
export HGT_SMALL_BIN="$HGT_ROOT/tools/annotation-trimming/bin"
panel=()
if [[ ${HGT_ACTINO:-0} == 1 ]]; then panel=(--actino); fi
"$HGT_BIN/python" "$HGT_ROOT/scripts/steps/h2_tree_job.py" --mode marker --index "${SLURM_ARRAY_TASK_ID:?}" "${panel[@]}"
