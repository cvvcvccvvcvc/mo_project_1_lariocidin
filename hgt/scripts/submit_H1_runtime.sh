#!/usr/bin/env bash
set -euo pipefail
HGT_ROOT=${HGT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
: "${HGT_BIN:?Set HGT_BIN to the pinned environment bin directory}"
export HGT_ROOT HGT_BIN
mkdir -p "$HGT_ROOT/results/slurm"
sbatch --parsable --partition=main --qos=normal \
    --cpus-per-task=1 --mem=2G --time=00:05:00 \
    --job-name=HGT-runtime \
    --output=$HGT_ROOT/results/slurm/runtime_%j.out \
    --wrap='PYTHONDONTWRITEBYTECODE=1 "$HGT_BIN/python" -c "import sys,Bio,numpy; print(sys.version.split()[0], Bio.__version__, numpy.__version__)"'
