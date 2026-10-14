#!/usr/bin/env bash
set -euo pipefail
: "${HGT_ROOT:?Isolated H2 root required}"
cd "$HGT_ROOT"
mkdir -p logs
export HGT_ROOT PYTHONDONTWRITEBYTECODE=1
: "${HGT_BIN:?Set HGT_BIN to the pinned environment bin directory}"
export HGT_BIN
search=$(sbatch --parsable --array=0-7%8 --export=ALL --output=logs/H2_actino_search_%A_%a.out scripts/H2_actino_search.sbatch)
collect=$(sbatch --parsable --dependency=afterok:"$search" --partition=main --qos=normal --cpus-per-task=1 --mem=8G --time=00:20:00 --job-name=HGT-H2-actino-collect --output=logs/H2_actino_collect_%j.out --export=ALL --wrap='"$HGT_BIN/python" "$HGT_ROOT/scripts/h2_markers.py" --actino')
markers=$(sbatch --parsable --dependency=afterok:"$collect" --array=0-137%8 --cpus-per-task=4 --export=ALL,HGT_ACTINO=1 --output=logs/H2_actino_marker_%A_%a.out scripts/H2_marker_tree.sbatch)
printf 'search=%s\ncollect=%s\nmarkers=%s\n' "$search" "$collect" "$markers" > config/H2_actino_jobs.txt
cat config/H2_actino_jobs.txt
