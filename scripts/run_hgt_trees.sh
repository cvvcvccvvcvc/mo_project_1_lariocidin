#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
python="$project_dir/.env/bin/python"
mafft="$project_dir/.env/bin/mafft"
iqtree="$project_dir/.env/bin/iqtree3"
marker_dir="$project_dir/data/processed/hgt_host_markers"
alignment_dir="$project_dir/results/alignments/hgt"
tree_dir="$project_dir/results/iqtree/hgt"
table_dir="$project_dir/results/tables"
core_dir="$project_dir/data/processed/pseudomonadota_core"

mkdir -p "$alignment_dir/host_markers" "$alignment_dir/bgc" "$tree_dir"

"$python" "$project_dir/scripts/build_host_markers.py" \
  --assemblies "$project_dir/config/hgt_assemblies.tsv" \
  --selected "$project_dir/data/processed/lrcc_selected.tsv" \
  --control-review "$project_dir/config/hgt_control_review.tsv" \
  --core-dir "$core_dir" \
  --output-dir "$marker_dir" \
  --inventory-out "$table_dir/host_marker_inventory.tsv" \
  --taxa-out "$table_dir/host_taxa.tsv" \
  --bgc-validation-out "$table_dir/carrier_bgc_assembly_validation.tsv"

for fasta in "$marker_dir"/*.faa; do
  marker=$(basename "$fasta" .faa)
  "$mafft" --auto --quiet "$fasta" > "$alignment_dir/host_markers/$marker.aln.faa"
  "$python" "$project_dir/scripts/trim_protein_alignment.py" \
    --input "$alignment_dir/host_markers/$marker.aln.faa" \
    --output "$alignment_dir/host_markers/$marker.trim.faa" \
    --min-occupancy 0.8
done

"$python" "$project_dir/scripts/concat_host_markers.py" \
  --alignment-dir "$alignment_dir/host_markers" \
  --taxa "$table_dir/host_taxa.tsv" \
  --inventory "$table_dir/host_marker_inventory.tsv" \
  --output "$alignment_dir/host_concat.trim.faa" \
  --boundaries "$table_dir/host_marker_boundaries.tsv"

"$iqtree" -s "$alignment_dir/host_concat.trim.faa" -st AA -m MFP \
  -B 1000 --alrt 1000 -T AUTO -seed 20260911 --redo \
  --prefix "$tree_dir/host"

for role in B1 B2 C D; do
  "$mafft" --auto --quiet "$core_dir/pseudomonadota_core_${role}.faa" \
    > "$alignment_dir/bgc/$role.aln.faa"
  "$python" "$project_dir/scripts/trim_protein_alignment.py" \
    --input "$alignment_dir/bgc/$role.aln.faa" \
    --output "$alignment_dir/bgc/$role.trim.faa" \
    --min-occupancy 0.8
  "$iqtree" -s "$alignment_dir/bgc/$role.trim.faa" -st AA -m MFP \
    -B 1000 --alrt 1000 -T AUTO -seed 20260911 --redo \
    --prefix "$tree_dir/$role.free"
done

"$python" "$project_dir/scripts/assess_free_bgc_trees.py" \
  --selected "$project_dir/data/processed/lrcc_selected.tsv" \
  --tree-dir "$tree_dir" \
  --constraint-out "$project_dir/config/hgt_order_constraint.nwk" \
  --output "$table_dir/hgt_free_tree_assessment.tsv"

while IFS=$'\t' read -r role compatible model; do
  if [[ "$role" == "role" || "$compatible" == "True" ]]; then
    continue
  fi
  "$iqtree" -s "$alignment_dir/bgc/$role.trim.faa" -st AA -m "$model" \
    -g "$project_dir/config/hgt_order_constraint.nwk" -T AUTO -seed 20260911 \
    --redo --prefix "$tree_dir/$role.constrained"
  "$python" "$project_dir/scripts/combine_candidate_trees.py" \
    --free "$tree_dir/$role.free.treefile" \
    --constrained "$tree_dir/$role.constrained.treefile" \
    --output "$tree_dir/$role.candidates.trees"
  "$iqtree" -s "$alignment_dir/bgc/$role.trim.faa" -st AA -m "$model" \
    -z "$tree_dir/$role.candidates.trees" -n 0 -zb 10000 -zw -au \
    -T AUTO -seed 20260911 --redo --prefix "$tree_dir/$role.au"
done < "$table_dir/hgt_free_tree_assessment.tsv"
