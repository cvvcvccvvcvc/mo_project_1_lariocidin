#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
config="$project_dir/config/focal_context_accessions.tsv"
raw_dir="$project_dir/data/raw/focal_context"

mkdir -p "$raw_dir"
tail -n +2 "$config" | while IFS=$'\t' read -r sample_id accession protein_id database source_url; do
  output="$raw_dir/${sample_id}.gb"
  temporary="$output.download"
  curl --fail --location --retry 3 --silent --show-error \
    --user-agent 'lariocidin-course-research/1.0' \
    --get 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi' \
    --data-urlencode 'db=nuccore' \
    --data-urlencode "id=$accession" \
    --data-urlencode 'rettype=gbwithparts' \
    --data-urlencode 'retmode=text' \
    --output "$temporary"
  test -s "$temporary"
  grep -q '^LOCUS ' "$temporary"
  mv "$temporary" "$output"
  sleep 1
done
