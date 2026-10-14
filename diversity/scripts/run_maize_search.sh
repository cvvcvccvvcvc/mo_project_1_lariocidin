#!/usr/bin/env bash
# Requires the retained literature profiles; pinned metadata are restored below.
set -euo pipefail
project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
python="$project_dir/.env/bin/python"
"$python" "$project_dir/scripts/steps/fetch_research_inputs.py"
"$python" "$project_dir/scripts/steps/fetch_s3_sources.py"
"$python" "$project_dir/scripts/steps/validate_final_profiles.py"
"$python" "$project_dir/scripts/steps/search_maize_catalogue.py"
"$python" "$project_dir/scripts/steps/maize_genomes.py"
"$python" "$project_dir/scripts/steps/search_maize_direct.py"
"$python" "$project_dir/scripts/steps/audit_maize_provenance.py"
bash "$project_dir/scripts/run_maize_context_review.sh"
"$python" "$project_dir/scripts/steps/plot_maize_checkpoint.py"
