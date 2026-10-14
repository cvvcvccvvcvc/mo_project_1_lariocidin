#!/usr/bin/env bash
# Recurate the frozen context shards, then collect the final diversity catalogue.
set -euo pipefail
project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$project_dir"
test -x .env/bin/python || { echo 'Missing .env/bin/python' >&2; exit 1; }
test -f config/diversity_selection_review_protocol.json || { echo 'Missing review protocol' >&2; exit 1; }
tasks=config/diversity_selection_review_tasks.tsv
n=$(wc -l < "$tasks")
[[ "$n" -eq 17 ]] || { echo "Expected 17 frozen review shards; found $n" >&2; exit 1; }
while read -r screen output; do
  [[ -n "$screen" && -n "$output" ]] || { echo 'Invalid review task' >&2; exit 1; }
  test -f "data/processed/${screen}_contexts/validation.tsv" || {
    echo "Missing completed context: $screen" >&2; exit 1;
  }
done < "$tasks"
mkdir -p results/logs
review_job=$(sbatch --parsable --array="0-$((n-1))%8" --cpus-per-task=1 --mem=8G   --output=results/logs/diversity_review_%A_%a.out scripts/slurm/diversity_recurate.sh)
collect_job=$(sbatch --parsable --dependency="afterok:$review_job" --cpus-per-task=4 --mem=8G   --output=results/logs/diversity_collect_%j.out scripts/slurm/diversity_collect_review.sh)
printf 'review=%s\ncollect=%s\n' "$review_job" "$collect_job"
