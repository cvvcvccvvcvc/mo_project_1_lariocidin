#!/usr/bin/env bash
set -euo pipefail
project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
env_dir="$project_dir/.env"
processed="$project_dir/data/processed"
tables="$project_dir/results/tables"

"$env_dir/bin/python" "$project_dir/scripts/steps/extract_maize_contexts.py"
"$env_dir/bin/blastp" -query "$processed/maize/maize_direct_candidates.faa" \
  -subject "$processed/literature/literature_unique_cyclases.faa" -evalue 10 -max_target_seqs 307 \
  -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore qcovs' \
  -out "$tables/maize_direct_vs_literature.tsv"
"$env_dir/bin/blastp" -query "$processed/maize/maize_candidate_neighbors.faa" \
  -subject "$processed/literature/literature_all_cds.faa" -evalue 1e-5 -max_target_seqs 10 \
  -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore qcovs' \
  -out "$tables/maize_neighbors_vs_literature.tsv"
"$env_dir/bin/blastp" -task blastp-short -seg no -query "$processed/maize/maize_candidate_small_orfs.faa" \
  -subject "$processed/literature/literature_unique_precursors.faa" -evalue 10 -max_target_seqs 5 \
  -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore qcovs' \
  -out "$tables/maize_small_orfs_vs_precursors.tsv"
"$env_dir/bin/python" "$project_dir/scripts/steps/summarize_maize_candidates.py"
"$env_dir/bin/blastp" -query "$processed/maize/maize_direct_candidates.faa" \
  -subject "$project_dir/data/raw/s3/uniprot_reviewed_bacterial_asnb.faa" -evalue 10 -max_target_seqs 5 \
  -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore qcovs' \
  -out "$tables/maize_candidates_vs_asnb.tsv"
"$env_dir/bin/python" "$project_dir/scripts/steps/curate_maize_checkpoint.py"
