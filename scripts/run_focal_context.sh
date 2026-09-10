#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
raw_dir="$project_dir/data/raw/focal_context"

"$project_dir/scripts/fetch_focal_context.sh"
"$project_dir/.env/bin/python" "$project_dir/scripts/summarize_focal_context.py" \
  --config "$project_dir/config/focal_context_accessions.tsv" \
  --raw-dir "$raw_dir" \
  --selected "$project_dir/data/processed/lrcc_selected.tsv" \
  --core-roles "$project_dir/results/tables/pseudomonadota_core_roles.tsv" \
  --all-cds "$project_dir/data/processed/all_cds.faa" \
  --manifest-out "$project_dir/results/tables/focal_context_manifest.tsv" \
  --features-out "$project_dir/results/tables/focal_context_features.tsv"
