#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-hgt-H2-native
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
set -euo pipefail
: "${HGT_ROOT:?Isolated H2 root required}"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1
: "${HGT_BIN:?Set HGT_BIN to the pinned environment bin directory}"
export HGT_BIN
export HGT_SMALL_BIN="$HGT_ROOT/tools/annotation-trimming/bin"
: "${HGT_GENOME_ROOT:?Set HGT_GENOME_ROOT to the frozen H1 genome directory}"
export HGT_GENOME_ROOT
"$HGT_BIN/python" "$HGT_ROOT/scripts/h2_native.py" --shard "${SLURM_ARRAY_TASK_ID:?}" --shards 16
