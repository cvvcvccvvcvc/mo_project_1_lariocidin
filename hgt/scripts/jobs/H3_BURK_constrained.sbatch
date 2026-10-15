#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-H3-BURK-constraint
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail

: "${HGT_SOURCE_ROOT:?Set frozen H2 root}"
: "${HGT_RESULT_ROOT:?Set separate H3 output root}"
: "${HGT_BIN:?Set pinned IQ-TREE environment bin}"

output="$HGT_RESULT_ROOT/topology/BURK_C_LG_G4"
mkdir -p "$output"
"$HGT_BIN/iqtree" \
  -s "$HGT_SOURCE_ROOT/results/H2/component_trees/canonical_C/trimmed_alignment.faa" \
  -m LG+G4 -g "$HGT_RESULT_ROOT/constraints/BURK_C_host_split_constraint.nwk" \
  -seed 131017 -T "${SLURM_CPUS_PER_TASK:?}" --prefix "$output/constrained"
