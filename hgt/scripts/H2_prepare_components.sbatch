#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-hgt-H2-gene-inputs
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
set -euo pipefail
: "${HGT_ROOT:?Isolated H2 root required}"
export PYTHONDONTWRITEBYTECODE=1
"$HGT_BIN/python" "$HGT_ROOT/scripts/h2_prepare_trees.py" --components
