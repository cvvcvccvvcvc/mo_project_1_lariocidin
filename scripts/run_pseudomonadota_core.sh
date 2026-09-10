#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
env_dir="$project_dir/.env"
processed_dir="$project_dir/data/processed"
table_dir="$project_dir/results/tables"
core_dir="$processed_dir/pseudomonadota_core"

"$env_dir/bin/python" "$project_dir/scripts/prepare_pseudomonadota_core.py" \
  --features "$processed_dir/all_cds.tsv" \
  --all-cds "$processed_dir/all_cds.faa" \
  --selected "$processed_dir/lrcc_selected.tsv" \
  --table-out "$table_dir/pseudomonadota_core_roles.tsv" \
  --fasta-dir "$core_dir"

for role in A B1 B2 C D; do
  "$env_dir/bin/blastp" \
    -query "$core_dir/pseudomonadota_core_${role}.faa" \
    -subject "$core_dir/pseudomonadota_core_${role}.faa" \
    -seg yes -evalue 10 -max_target_seqs 14 \
    -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore qcovs' \
    -out "$core_dir/pseudomonadota_core_${role}.tsv"
done

"$env_dir/bin/python" "$project_dir/scripts/summarize_pseudomonadota_core.py" \
  --roles "$table_dir/pseudomonadota_core_roles.tsv" \
  --blast-dir "$core_dir" \
  --all-pairs-out "$table_dir/pseudomonadota_core_pairwise.tsv" \
  --focal-out "$table_dir/p_asym_core_affinity.tsv" \
  --summary-out "$table_dir/pseudomonadota_core_summary.json"
