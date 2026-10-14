#!/bin/bash
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --time=00:30:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
.env/bin/python scripts/steps/collect_diversity_catalogue.py --output-prefix diversity_species_review
.env/bin/python scripts/steps/collect_diversity_catalogue.py --prefix diversity_extension --output-prefix diversity_extension_review \
    --selection inputs/selections/extension/diversity_extension_selection.tsv --refseq-metadata inputs/selections/extension/diversity_extension_refseq_metadata.tsv \
    --mgnify-metadata inputs/selections/extension/diversity_extension_mgnify_metadata.tsv
.env/bin/python scripts/steps/collect_diversity_catalogue.py --prefix diversity_controls --output-prefix diversity_controls_review --shards 1 \
    --selection inputs/selections/controls/diversity_controls_selection.tsv --refseq-metadata inputs/selections/controls/diversity_controls_refseq_metadata.tsv \
    --mgnify-metadata inputs/selections/extension/diversity_extension_mgnify_metadata.tsv
.env/bin/python scripts/steps/merge_diversity_catalogues.py --input-suffix review --output-prefix diversity_final
.env/bin/python scripts/steps/characterize_diversity_candidates.py
