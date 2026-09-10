#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
env_dir="$project_dir/.env"
processed_dir="$project_dir/data/processed"
table_dir="$project_dir/results/tables"
raw_dir="$project_dir/data/raw/research_2026-09-10"

"$env_dir/bin/python" "$project_dir/scripts/prepare_literature_reference.py" \
  --archive "$raw_dir/triculamin_lasso_BGCs.zip" \
  --lar29-fasta "$processed_dir/lrcc_reference.faa" \
  --anchors "$project_dir/config/literature_anchors.tsv" \
  --processed-dir "$processed_dir" \
  --table-dir "$table_dir"

"$env_dir/bin/blastp" \
  -query "$processed_dir/literature_unresolved_bgc_cds.faa" \
  -subject "$processed_dir/literature_noncanonical_cyclases.faa" \
  -seg yes -evalue 10 -max_target_seqs 5 \
  -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore qcovs' \
  -out "$processed_dir/literature_unresolved_vs_noncanonical.tsv"

"$env_dir/bin/python" "$project_dir/scripts/rank_unresolved_literature.py" \
  --manifest "$processed_dir/literature_bgc_manifest.tsv" \
  --features "$processed_dir/literature_all_cds.tsv" \
  --blast "$processed_dir/literature_unresolved_vs_noncanonical.tsv" \
  --out "$table_dir/literature_unresolved_cyclase_ranking.tsv" \
  --candidates-out "$table_dir/literature_inferred_cyclase_candidates.tsv"

"$env_dir/bin/blastp" \
  -query "$processed_dir/literature_unique_cyclases.faa" \
  -subject "$processed_dir/literature_unique_cyclases.faa" \
  -seg yes -evalue 10 -max_target_seqs 307 \
  -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore qcovs qcovhsp' \
  -out "$processed_dir/literature_cyclase_all_vs_all.tsv"

"$env_dir/bin/python" "$project_dir/scripts/summarize_literature_similarity.py" \
  --cyclases "$processed_dir/literature_unique_cyclases.tsv" \
  --blast "$processed_dir/literature_cyclase_all_vs_all.tsv" \
  --nearest-out "$table_dir/literature_cyclase_nearest_by_group.tsv" \
  --summary-out "$table_dir/literature_similarity_summary.json"
