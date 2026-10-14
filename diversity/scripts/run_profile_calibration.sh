#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
env_dir="$project_dir/.env"
processed_dir="$project_dir/data/processed"
table_dir="$project_dir/results/tables"
profile_dir="$project_dir/results/profile_calibration"

"$env_dir/bin/python" "$project_dir/scripts/steps/calibrate_family_profiles.py" \
  --cyclase-fasta "$processed_dir/literature/literature_unique_cyclases.faa" \
  --cyclase-table "$processed_dir/literature/literature_unique_cyclases.tsv" \
  --all-v-all "$processed_dir/literature/literature_cyclase_all_vs_all.tsv" \
  --mafft "$env_dir/bin/mafft" \
  --hmmbuild "$env_dir/bin/hmmbuild" \
  --hmmsearch "$env_dir/bin/hmmsearch" \
  --work-dir "$profile_dir" \
  --scores-out "$table_dir/literature_profile_cv_scores.tsv" \
  --folds-out "$table_dir/literature_profile_cv_folds.tsv" \
  --summary-out "$table_dir/literature_profile_calibration_summary.json"
