#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-hgt-H1-4BL1
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -euo pipefail
: "${HGT_ROOT:?Set isolated HGT root}"
export PYTHONDONTWRITEBYTECODE=1
cd "$HGT_ROOT"
sha256sum -c config/H1_4BL1_inputs.sha256 > results/slurm/4BL1_input_hash_check.txt
"$HGT_BIN/python" scripts/audit_h1.py \
    --inputs results/tables/H1_Nature_4BL1_assembly_inputs.tsv \
    --metadata inputs/H1_2026-10-13/ncbi/H1_Nature_4BL1_candidate_metadata.json \
    --output-dir results/H1_4BL1_shards/0
