#!/usr/bin/env bash
set -euo pipefail
HGT_ROOT=${HGT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
cd "$HGT_ROOT"
mkdir -p results/slurm
sha256sum -c config/H1/H1_slurm_inputs.sha256 > results/slurm/input_hash_check.txt
export HGT_ROOT
export HGT_PYTHON="${HGT_PYTHON:-${HGT_BIN:?Set HGT_BIN}/python}"
sbatch --parsable --array=0-7%8 --chdir="$HGT_ROOT" \
    --output="$HGT_ROOT/results/slurm/H1_%A_%a.out" \
    --error="$HGT_ROOT/results/slurm/H1_%A_%a.err" scripts/jobs/H1_correspondence.sbatch
