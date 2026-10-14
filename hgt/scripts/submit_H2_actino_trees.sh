#!/usr/bin/env bash
set -euo pipefail
: "${HGT_ROOT:?Isolated H2 root required}"
markers=${1:?Actinobacteria138 marker array required}
cd "$HGT_ROOT"
mkdir -p logs
export HGT_ROOT HGT_ACTINO=1 PYTHONDONTWRITEBYTECODE=1
: "${HGT_BIN:?Set HGT_BIN to the pinned environment bin directory}"
export HGT_BIN
prepare=$(sbatch --parsable --dependency=afterok:"$markers" --partition=main --qos=normal --cpus-per-task=1 --mem=8G --time=00:20:00 --job-name=HGT-H2-actino-inputs --output=logs/H2_actino_prepare_%j.out --export=ALL --wrap='"$HGT_BIN/python" "$HGT_ROOT/scripts/h2_prepare_trees.py" --hosts --components --actino --fixed-host-models')
# Main noncanonical panel plus a separate uncertain-Nature diagnostic with relatives.
host=$(sbatch --parsable --dependency=afterok:"$prepare" --array=0 --export=ALL,HGT_FIXED_HOST_MODELS=1 --output=logs/H2_actino_host_%A_%a.out scripts/H2_host_tree.sbatch)
nature=$(sbatch --parsable --dependency=afterok:"$prepare" --array=1 --export=ALL,HGT_FIXED_HOST_MODELS=0 --output=logs/H2_Nature_host_%A_%a.out scripts/H2_host_tree.sbatch)
components=$(sbatch --parsable --dependency=afterok:"$prepare" --array=0-3%3 --export=ALL --output=logs/H2_actino_component_%A_%a.out scripts/H2_component_tree.sbatch)
printf 'prepare=%s\nhost=%s\nNature_host=%s\ncomponents=%s\n' "$prepare" "$host" "$nature" "$components" > config/H2_actino_tree_jobs.txt
cat config/H2_actino_tree_jobs.txt
