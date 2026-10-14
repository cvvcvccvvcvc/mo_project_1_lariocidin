#!/bin/bash
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --time=00:20:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
read -r screen output < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" config/diversity_selection_review_tasks.tsv)
exec .env/bin/python scripts/curate_diversity_contexts.py --screen-name "$screen" \
    --reuse-evidence --output-name "$output" --protocol config/diversity_selection_review_protocol.json
