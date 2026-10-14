#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-hgt-H2-marker-QC
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
: "${HGT_ROOT:?Isolated H2 root required}"
export PYTHONDONTWRITEBYTECODE=1
"$HGT_BIN/python" "$HGT_ROOT/scripts/steps/h2_markers.py"
