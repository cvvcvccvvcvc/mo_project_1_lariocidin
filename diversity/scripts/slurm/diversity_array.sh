#!/bin/bash
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --time=24:00:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
selection=()
if [[ -n "${DIVERSITY_SELECTION:-}" ]]; then selection=(--selection "$DIVERSITY_SELECTION"); fi
exec .env/bin/python scripts/steps/screen_diversity_proteomes.py --limit 0 "${selection[@]}" \
  --shard-count "${DIVERSITY_SHARDS:-8}" --shard-index "$SLURM_ARRAY_TASK_ID" \
  --cpu "$SLURM_CPUS_PER_TASK" --workers 1 --batch-size 100 \
  --output-name "${DIVERSITY_PREFIX:-diversity_species}_${SLURM_ARRAY_TASK_ID}"
