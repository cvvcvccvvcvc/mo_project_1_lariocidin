#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

"$project_dir/.env/bin/python" "$project_dir/scripts/fetch_research_inputs.py"

"$project_dir/.env/bin/python" "$project_dir/scripts/fetch_hgt_assemblies.py" \
  --selected "$project_dir/data/processed/lrcc_selected.tsv" \
  --raw-dir "$project_dir/data/raw/s5_ncbi" \
  --manifest "$project_dir/config/s5_download_manifest.json" \
  --output "$project_dir/config/hgt_assemblies.tsv"
