#!/bin/bash
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --time=24:00:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
exec .env/bin/python scripts/steps/screen_diversity_proteomes.py "$@"
