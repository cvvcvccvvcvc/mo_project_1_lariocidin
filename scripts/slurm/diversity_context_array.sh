#!/bin/bash
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --time=24:00:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
screen="${DIVERSITY_PREFIX:-diversity_species}_${SLURM_ARRAY_TASK_ID}"
.env/bin/python scripts/extract_diversity_contexts.py --screen-name "$screen"
exec .env/bin/python scripts/curate_diversity_contexts.py --screen-name "$screen" --cpu "$SLURM_CPUS_PER_TASK"
