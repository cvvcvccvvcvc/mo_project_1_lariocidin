#!/usr/bin/env bash
#SBATCH --job-name=diversity-baseline
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
bash scripts/run_reference.sh
bash scripts/run_literature_reference.sh
bash scripts/run_profile_calibration.sh
bash scripts/run_maize_search.sh
.env/bin/python scripts/steps/prepare_pseudomonadota_core.py \
  --features data/processed/reference/all_cds.tsv \
  --all-cds data/processed/reference/all_cds.faa \
  --selected data/processed/reference/lrcc_selected.tsv \
  --table-out results/tables/pseudomonadota_core_roles.tsv \
  --fasta-dir results/pseudomonadota_core
.env/bin/python scripts/steps/fetch_research_inputs.py --diversity
.env/bin/python scripts/steps/prepare_diversity_sources.py
.env/bin/python scripts/steps/annotate_diversity_domains.py
.env/bin/python scripts/steps/audit_diversity_baseline.py
.env/bin/python scripts/steps/resolve_diversity_ftp.py
.env/bin/python scripts/steps/prepare_diversity_published.py
.env/bin/python scripts/steps/prepare_diversity_role_reference.py
