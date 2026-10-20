#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-H3-AU-array
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -euo pipefail

: "${HGT_AU_QUEUE:?Set frozen AU_queue.tsv path}"
: "${SLURM_ARRAY_TASK_ID:?Submit as an array}"
row="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$HGT_AU_QUEUE")"
test -n "$row"
IFS=$'\t' read -r index gcf_rank scope family split model <<< "$row"
test "$index" = "$SLURM_ARRAY_TASK_ID"
test -n "$model"
bash "$HGT_RESULT_ROOT/H3_screen_AU.sbatch" "$scope" "$family" "$split" "$model"
