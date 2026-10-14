#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root"
test -x .env/bin/python || { echo 'Create .env from environment.linux-64.lock first' >&2; exit 1; }
test -f inputs/archives/diversity_metadata_snapshots.zip || { echo 'Missing pinned diversity metadata archive' >&2; exit 1; }
sha256sum -c inputs/provenance/off_git_assets.sha256
mkdir -p results/logs results/tables data/processed/discovery
sbatch --parsable --cpus-per-task=16 --mem=64G --time=1-12:00:00 \
  --output=results/logs/baseline_%j.out scripts/slurm/diversity_baseline.sbatch
