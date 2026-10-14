#!/usr/bin/env bash
set -euo pipefail
HGT_ROOT=${HGT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
cd "$HGT_ROOT"
sha256sum -c config/H2_compute_inputs.sha256 > results/slurm/native_input_hash_check.txt
export HGT_ROOT
sbatch --parsable --array=0-15%8 --output=results/slurm/H2_native_%A_%a.out \
    --error=results/slurm/H2_native_%A_%a.err scripts/H2_native.sbatch
