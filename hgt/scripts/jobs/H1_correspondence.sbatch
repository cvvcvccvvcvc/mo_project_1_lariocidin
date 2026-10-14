#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-hgt-H1
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
set -euo pipefail
: "${HGT_ROOT:?Set isolated HGT root}"
: "${HGT_PYTHON:?Set recorded read-only Python runtime}"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1
HGT_SHARD=${SLURM_ARRAY_TASK_ID:?}
"$HGT_PYTHON" -c 'import sys,Bio,numpy; print("runtime",sys.version.split()[0],Bio.__version__,numpy.__version__)'
"$HGT_PYTHON" "$HGT_ROOT/scripts/steps/audit_h1.py" --shard "$HGT_SHARD" --shards 8 \
    --output-dir "$HGT_ROOT/results/H1_shards/$HGT_SHARD"
