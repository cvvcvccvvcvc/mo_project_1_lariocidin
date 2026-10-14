#!/usr/bin/env bash
set -euo pipefail
HGT_ROOT=${HGT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
cd "$HGT_ROOT"
sha256sum -c config/H1_supplement_inputs.sha256 > results/slurm/supplement_input_hash_check.txt
export HGT_ROOT
export HGT_PYTHON="${HGT_PYTHON:-${HGT_BIN:?Set HGT_BIN}/python}"
sbatch --parsable --array=0-7%8 --chdir="$HGT_ROOT" \
    --output="$HGT_ROOT/results/slurm/H1_supplement_%A_%a.out" \
    --error="$HGT_ROOT/results/slurm/H1_supplement_%A_%a.err" scripts/H1_supplement.sbatch
