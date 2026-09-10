#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
env_dir="$project_dir/.env"
processed_dir="$project_dir/data/processed"
alignment_dir="$project_dir/results/alignment"
iqtree_dir="$project_dir/results/iqtree"
log_dir="$project_dir/results/logs"
table_dir="$project_dir/results/tables"

mkdir -p "$alignment_dir" "$iqtree_dir" "$log_dir" "$table_dir"

"$env_dir/bin/python" "$project_dir/scripts/prepare_tree_sensitivity.py" \
  --fasta "$processed_dir/lrcc_reference.faa" \
  --metadata "$processed_dir/lrcc_selected.tsv" \
  --subset-fasta "$processed_dir/lrcc_pseudomonadota.faa"

"$env_dir/bin/mafft" --localpair --maxiterate 1000 --thread 4 \
  "$processed_dir/lrcc_pseudomonadota.faa" \
  > "$alignment_dir/lrcc_pseudomonadota.aligned.faa" \
  2> "$log_dir/mafft_pseudomonadota.log"

# Masking requires the completed alignment, so regenerate only the masked file.
"$env_dir/bin/python" "$project_dir/scripts/prepare_tree_sensitivity.py" \
  --fasta "$processed_dir/lrcc_reference.faa" \
  --metadata "$processed_dir/lrcc_selected.tsv" \
  --subset-fasta "$processed_dir/lrcc_pseudomonadota.faa" \
  --aligned "$alignment_dir/lrcc_pseudomonadota.aligned.faa" \
  --masked "$alignment_dir/lrcc_pseudomonadota.occupancy80.faa"

run_tree() {
  local name=$1
  local alignment=$2
  shift 2
  "$env_dir/bin/iqtree" \
    -s "$alignment" -m Q.PFAM+I+G4 -B 1000 --alrt 1000 \
    -T 2 --seed 20260910 -redo "$@" \
    --prefix "$iqtree_dir/$name" \
    > "$log_dir/$name.stdout.log" 2> "$log_dir/$name.stderr.log"
}

run_tree pseudomonadota_default "$alignment_dir/lrcc_pseudomonadota.aligned.faa" &
pid_default=$!
run_tree pseudomonadota_bnni "$alignment_dir/lrcc_pseudomonadota.aligned.faa" -bnni &
pid_bnni=$!
run_tree pseudomonadota_occupancy80_bnni "$alignment_dir/lrcc_pseudomonadota.occupancy80.faa" -bnni &
pid_masked=$!
wait "$pid_default" "$pid_bnni" "$pid_masked"

"$env_dir/bin/python" "$project_dir/scripts/summarize_tree_sensitivity.py" \
  --variant full_original "$iqtree_dir/reference_lrcc.treefile" "$alignment_dir/lrcc_with_outgroup.aligned.faa" \
  --variant canonical_original "$iqtree_dir/canonical_lrcc.treefile" "$alignment_dir/lrcc_canonical_with_outgroup.aligned.faa" \
  --variant pseudomonadota_default "$iqtree_dir/pseudomonadota_default.treefile" "$alignment_dir/lrcc_pseudomonadota.aligned.faa" \
  --variant pseudomonadota_bnni "$iqtree_dir/pseudomonadota_bnni.treefile" "$alignment_dir/lrcc_pseudomonadota.aligned.faa" \
  --variant pseudomonadota_occupancy80_bnni "$iqtree_dir/pseudomonadota_occupancy80_bnni.treefile" "$alignment_dir/lrcc_pseudomonadota.occupancy80.faa" \
  --out "$table_dir/focal_split_sensitivity.tsv"
