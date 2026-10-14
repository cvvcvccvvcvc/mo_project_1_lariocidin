#!/usr/bin/env bash
set -euo pipefail
: "${HGT_ROOT:?Isolated H2 root required}"
marker_array=${1:?Completed or running Bacteria74 marker array job ID required}
cd "$HGT_ROOT"
mkdir -p logs
export HGT_ROOT PYTHONDONTWRITEBYTECODE=1
: "${HGT_BIN:?Set HGT_BIN to the pinned environment bin directory}"
export HGT_BIN
prepare=$(sbatch --parsable --dependency=afterok:"$marker_array" --partition=main --qos=normal --cpus-per-task=1 --mem=8G --time=00:15:00 --job-name=HGT-H2-host-inputs --output=logs/H2_host_prepare_%j.out --export=ALL --wrap='"$HGT_BIN/python" "$HGT_ROOT/scripts/h2_prepare_trees.py" --hosts --fixed-host-models')
fit=$(sbatch --parsable --dependency=afterok:"$prepare" --array=0-2%3 --export=ALL,HGT_FIXED_HOST_MODELS=1 --output=logs/H2_host_%A_%a.out scripts/H2_host_tree.sbatch)
printf 'prepare=%s\nfit=%s\n' "$prepare" "$fit" > config/H2_host_jobs.txt
cat config/H2_host_jobs.txt
