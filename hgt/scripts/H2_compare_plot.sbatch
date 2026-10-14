#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-HGT-H2-results
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
export MPLCONFIGDIR="$HGT_ROOT/tools/matplotlib-cache"
options=()
if [[ ${HGT_READY_ONLY:-0} == 1 ]]; then options=(--ready-only); fi
"$HGT_BIN/python" "$HGT_ROOT/scripts/h2_compare_trees.py" "${options[@]}"
"$HGT_BIN/python" "$HGT_ROOT/scripts/h2_review_relations.py" "${options[@]}"
"$HGT_BIN/python" "$HGT_ROOT/scripts/h2_plot_trees.py" "${options[@]}"
