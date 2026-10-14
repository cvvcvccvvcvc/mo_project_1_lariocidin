#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root"
for path in .env/bin/python protocols/diversity_protocol.json inputs/selections/species/diversity_refseq_species.tsv data/processed/discovery/diversity_refseq_metadata.tsv; do
  test -e "$path" || { echo "Missing prerequisite: $path" >&2; exit 1; }
done
mkdir -p results/logs results/tables data/processed/discovery/contexts
search=$(sbatch --parsable --array=0-7%8 --cpus-per-task=8 --mem=8G \
  --output=results/logs/diversity_screen_%A_%a.out scripts/slurm/diversity_array.sh)
contexts=$(sbatch --parsable --dependency="aftercorr:$search" --array=0-7%8 \
  --cpus-per-task=8 --mem=8G --output=results/logs/diversity_context_%A_%a.out \
  scripts/slurm/diversity_context_array.sh)
collect=$(sbatch --parsable --dependency="afterok:$contexts" --cpus-per-task=4 --mem=8G \
  --output=results/logs/diversity_species_collect_%j.out \
  --wrap='.env/bin/python scripts/steps/collect_diversity_catalogue.py')
printf 'search=%s\ncontexts=%s\ncollect=%s\n' "$search" "$contexts" "$collect"
