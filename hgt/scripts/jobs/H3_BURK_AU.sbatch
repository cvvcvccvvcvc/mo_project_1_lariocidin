#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-H3-BURK-AU
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail

: "${HGT_SOURCE_ROOT:?Set frozen H2 root}"
: "${HGT_RESULT_ROOT:?Set separate H3 output root}"
: "${HGT_BIN:?Set pinned IQ-TREE environment bin}"

input="$HGT_SOURCE_ROOT/results/H2/component_trees/canonical_C"
output="$HGT_RESULT_ROOT/topology/BURK_C_LG_G4"
test -s "$input/ML.treefile"
test -s "$output/constrained.treefile"
cat "$input/ML.treefile" "$output/constrained.treefile" > "$output/candidate_trees.nwk"
"$HGT_BIN/iqtree" \
  -s "$input/trimmed_alignment.faa" -m LG+G4 \
  -z "$output/candidate_trees.nwk" -te "$input/ML.treefile" \
  -zb 10000 -au -seed 131018 -T "${SLURM_CPUS_PER_TASK:?}" \
  --prefix "$output/AU"
