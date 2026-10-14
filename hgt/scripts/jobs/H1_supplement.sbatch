#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-hgt-H1-supplement
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
set -euo pipefail
: "${HGT_ROOT:?Set isolated HGT root}"
: "${HGT_PYTHON:?Set read-only Python runtime}"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1
"$HGT_PYTHON" "$HGT_ROOT/scripts/steps/audit_h1.py" --shard "${SLURM_ARRAY_TASK_ID:?}" --shards 8 \
    --inputs "$HGT_ROOT/results/tables/H1_supplement_assembly_inputs.tsv" \
    --metadata "$HGT_ROOT/inputs/H1_2026-10-13/ncbi/H1_supplement_candidate_metadata.json" \
    --output-dir "$HGT_ROOT/results/H1_supplement_shards/$SLURM_ARRAY_TASK_ID"
