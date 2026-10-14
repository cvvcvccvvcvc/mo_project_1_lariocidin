#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root"
case "${1:-}" in
  extension)
    prefix=diversity_extension; shards=8
    selection=inputs/selections/diversity_extension_selection.tsv
    metadata=inputs/selections/diversity_extension_refseq_metadata.tsv
    mgnify=inputs/selections/diversity_extension_mgnify_metadata.tsv
    ;;
  controls)
    prefix=diversity_controls; shards=1
    selection=inputs/selections/diversity_controls_selection.tsv
    metadata=inputs/selections/diversity_controls_refseq_metadata.tsv
    mgnify=inputs/selections/diversity_extension_mgnify_metadata.tsv
    ;;
  *) echo 'Usage: submit_additional.sh extension|controls' >&2; exit 2 ;;
esac
for path in .env/bin/python "$selection" "$metadata" "$mgnify" data/processed/discovery/diversity_published_records.tsv; do
  test -e "$path" || { echo "Missing prerequisite: $path" >&2; exit 1; }
done
mkdir -p results/logs results/tables data/processed/discovery/contexts
export_args="ALL,DIVERSITY_SELECTION=$selection,DIVERSITY_PREFIX=$prefix,DIVERSITY_SHARDS=$shards"
search=$(sbatch --parsable --export="$export_args" --array="0-$((shards-1))%8" --cpus-per-task=8 --mem=8G \
  --output="results/logs/${prefix}_screen_%A_%a.out" scripts/slurm/diversity_array.sh)
contexts=$(sbatch --parsable --export="$export_args" --dependency="aftercorr:$search" --array="0-$((shards-1))%8" \
  --cpus-per-task=8 --mem=8G --output="results/logs/${prefix}_context_%A_%a.out" \
  scripts/slurm/diversity_context_array.sh)
collect=$(sbatch --parsable --export="$export_args" --dependency="afterok:$contexts" --cpus-per-task=4 --mem=8G \
  --output="results/logs/${prefix}_collect_%j.out" \
  --wrap=".env/bin/python scripts/steps/collect_diversity_catalogue.py --shards $shards --prefix $prefix --selection $selection --refseq-metadata $metadata --mgnify-metadata $mgnify")
printf 'search=%s\ncontexts=%s\ncollect=%s\n' "$search" "$contexts" "$collect"
