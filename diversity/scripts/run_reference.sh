#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
env_dir="$project_dir/.env"
python_bin="$env_dir/bin/python"
reference_dir="$project_dir/data/raw/reference/LAR-BGC-42280aa5d43ef01f766fa207740dd9914166134c"
processed_dir="$project_dir/data/processed"
results_dir="$project_dir/results"

test -x "$python_bin"
"$project_dir/scripts/fetch_reference.sh"
mkdir -p "$processed_dir" "$results_dir/alignment" "$results_dir/iqtree" "$results_dir/logs" "$results_dir/tables" "$results_dir/figures"

"$python_bin" "$project_dir/scripts/steps/prepare_reference.py" \
  --gbk-dir "$reference_dir" \
  --metadata "$project_dir/inputs/selections/reference_metadata.tsv" \
  --m2-genome "$project_dir/data/raw/CP169648.1.gb" \
  --paeninodin-genome "$project_dir/data/raw/AHKH01000064.1.gb" \
  --out-dir "$processed_dir"

"$env_dir/bin/blastp" \
  -query "$processed_dir/all_cds.faa" \
  -subject "$processed_dir/lrcc_m2_seed.faa" \
  -seg yes -evalue 100 \
  -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore qcovs' \
  -out "$processed_dir/all_cds_vs_m2_lrcc.tsv"

"$python_bin" "$project_dir/scripts/steps/rank_lrcc_candidates.py" \
  --features "$processed_dir/all_cds.tsv" \
  --blast "$processed_dir/all_cds_vs_m2_lrcc.tsv" \
  --out "$results_dir/tables/lrcc_candidate_ranking.tsv" --top 4

"$python_bin" "$project_dir/scripts/steps/finalize_lrcc.py" \
  --manifest "$processed_dir/reference_manifest.tsv" \
  --features "$processed_dir/all_cds.tsv" \
  --all-cds "$processed_dir/all_cds.faa" \
  --ranking "$results_dir/tables/lrcc_candidate_ranking.tsv" \
  --selection "$project_dir/inputs/selections/lrcc_selection.tsv" \
  --outgroup "$processed_dir/paeninodin_outgroup.faa" \
  --out-dir "$processed_dir"

"$env_dir/bin/mafft" --localpair --maxiterate 1000 --thread 4 \
  "$processed_dir/lrcc_with_outgroup.faa" \
  > "$results_dir/alignment/lrcc_with_outgroup.aligned.faa" \
  2> "$results_dir/logs/mafft.log"
"$env_dir/bin/mafft" --localpair --maxiterate 1000 --thread 4 \
  "$processed_dir/lrcc_canonical_with_outgroup.faa" \
  > "$results_dir/alignment/lrcc_canonical_with_outgroup.aligned.faa" \
  2> "$results_dir/logs/mafft_canonical.log"

"$env_dir/bin/iqtree" \
  -s "$results_dir/alignment/lrcc_with_outgroup.aligned.faa" \
  -m MFP -B 1000 --alrt 1000 -T AUTO --seed 20260909 -redo \
  --prefix "$results_dir/iqtree/reference_lrcc"
"$env_dir/bin/iqtree" \
  -s "$results_dir/alignment/lrcc_canonical_with_outgroup.aligned.faa" \
  -m MFP -B 1000 --alrt 1000 -T AUTO --seed 20260909 -redo \
  --prefix "$results_dir/iqtree/canonical_lrcc"

MPLCONFIGDIR="$project_dir/.mplconfig" "$python_bin" "$project_dir/scripts/steps/build_report_assets.py" \
  --manifest "$processed_dir/reference_manifest.tsv" \
  --features "$processed_dir/all_cds.tsv" \
  --selected "$processed_dir/lrcc_selected.tsv" \
  --alignment "$results_dir/alignment/lrcc_with_outgroup.aligned.faa" \
  --tree "$results_dir/iqtree/reference_lrcc.treefile" \
  --canonical-alignment "$results_dir/alignment/lrcc_canonical_with_outgroup.aligned.faa" \
  --canonical-tree "$results_dir/iqtree/canonical_lrcc.treefile" \
  --table-dir "$results_dir/tables" \
  --figure-dir "$results_dir/figures"

"$python_bin" "$project_dir/scripts/steps/validate_reference.py" \
  --manifest "$processed_dir/reference_manifest.tsv" \
  --features "$processed_dir/all_cds.tsv" \
  --selected "$processed_dir/lrcc_selected.tsv" \
  --reference-fasta "$processed_dir/lrcc_reference.faa" \
  --alignment "$results_dir/alignment/lrcc_with_outgroup.aligned.faa" \
  --tree "$results_dir/iqtree/reference_lrcc.treefile"
"$python_bin" "$project_dir/scripts/steps/write_versions.py" --out "$results_dir/software_versions.tsv"

