#!/usr/bin/env bash
# Rebuild S5 evidence and figures from the saved trees and pinned assemblies.
set -euo pipefail
project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$project_dir"
python="$project_dir/.env/bin/python"
tables="results/tables"
core="data/processed/pseudomonadota_core"

bash scripts/run_hgt_assembly_fetch.sh
"$python" scripts/summarize_hgt_trees.py \
  --selected data/processed/lrcc_selected.tsv --taxa "$tables/host_taxa.tsv" \
  --free-assessment "$tables/hgt_free_tree_assessment.tsv" \
  --alignment-dir results/alignments/hgt/bgc --tree-dir results/iqtree/hgt \
  --host-alignment results/alignments/hgt/host_concat.trim.faa \
  --bgc-validation "$tables/carrier_bgc_assembly_validation.tsv" \
  --tests-out "$tables/hgt_topology_tests.tsv" --pairs-out "$tables/host_control_pairs.tsv" \
  --host-summary-out "$tables/host_tree_summary.json"
"$python" scripts/fetch_hgt_biosamples.py \
  --assemblies config/hgt_assemblies.tsv --raw-dir data/raw/s5_ncbi \
  --manifest config/s5_download_manifest.json \
  --summary-out "$tables/hgt_biosamples.tsv" --attributes-out "$tables/hgt_biosample_attributes.tsv"
"$python" scripts/extract_hgt_control_proteomes.py \
  --assemblies config/hgt_assemblies.tsv \
  --fasta-out data/processed/hgt_control_proteomes.faa \
  --features-out data/processed/hgt_control_features.tsv
"$python" scripts/screen_hgt_controls.py \
  --fasta data/processed/hgt_control_proteomes.faa --features data/processed/hgt_control_features.tsv \
  --output-dir data/interim/hgt_control_search --hits-out "$tables/hgt_control_profile_hits.tsv"
"$python" scripts/curate_hgt_controls.py \
  --assemblies config/hgt_assemblies.tsv --review config/hgt_control_review.tsv \
  --features data/processed/hgt_control_features.tsv --proteins data/processed/hgt_control_proteomes.faa \
  --known-c "$core/pseudomonadota_core_C.faa" --core-dir "$core" \
  --homology-out "$tables/hgt_control_neighbor_homology.tsv" \
  --hits "$tables/hgt_control_profile_hits.tsv" --status-out "$tables/hgt_control_bgc_status.tsv" \
  --contexts-out "$tables/hgt_control_bgc_contexts.tsv"
"$python" scripts/audit_core_assembly_dna.py \
  --assemblies config/hgt_assemblies.tsv --validation "$tables/carrier_bgc_assembly_validation.tsv" \
  --core-dir "$core" --output "$tables/carrier_bgc_dna_validation.tsv"
"$python" scripts/plot_hgt_summary.py \
  --host-tree results/iqtree/hgt/host.treefile --taxa "$tables/host_taxa.tsv" \
  --assemblies config/hgt_assemblies.tsv --tests "$tables/hgt_topology_tests.tsv" \
  --control-status "$tables/hgt_control_bgc_status.tsv" --biosamples "$tables/hgt_biosamples.tsv" \
  --output results/figures/s5_hgt_summary.png
