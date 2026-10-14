#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-HGT-H2-gCF
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
: "${HGT_ROOT:?Isolated H2 root required}"
export PYTHONDONTWRITEBYTECODE=1
: "${HGT_BIN:?Set HGT_BIN to the pinned environment bin directory}"
export HGT_BIN
panel=()
if [[ ${HGT_ACTINO:-0} == 1 ]]; then panel=(--actino); fi
"$HGT_BIN/python" "$HGT_ROOT/scripts/steps/h2_concordance.py" "${panel[@]}"
