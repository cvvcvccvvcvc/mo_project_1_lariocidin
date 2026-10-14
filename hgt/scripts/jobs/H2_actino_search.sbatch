#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-hgt-H2-actino
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:45:00
set -euo pipefail
: "${HGT_ROOT:?Isolated H2 root required}"
export PYTHONDONTWRITEBYTECODE=1
: "${HGT_BIN:?Set HGT_BIN to the pinned environment bin directory}"
export HGT_BIN
"$HGT_BIN/python" "$HGT_ROOT/scripts/steps/h2_actino_search.py"
