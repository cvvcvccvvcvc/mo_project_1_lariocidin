#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-H3-screen-AU
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -euo pipefail

: "${HGT_SOURCE_ROOT:?Set frozen H2 root}"
: "${HGT_RESULT_ROOT:?Set separate H3 output root}"
: "${HGT_BIN:?Set pinned IQ-TREE environment bin}"

scope="${1:?Bacteria74 or Actinobacteria138}"
family="${2:?component family}"
split="${3:?host split ID}"
model="${4:?H2 component model}"
case "$scope" in
  Bacteria74) panel="" ;;
  Actinobacteria138) panel="actino/" ;;
  *) exit 2 ;;
esac

input="$HGT_SOURCE_ROOT/results/H2/${panel}component_trees/$family"
constraint="$HGT_RESULT_ROOT/global_screen/constraints/${scope}_${family}_${split}.nwk"
output="$HGT_RESULT_ROOT/global_screen/AU/${scope}_${family}_${split}"
test -s "$input/trimmed_alignment.faa"
test -s "$input/ML.treefile"
test -s "$constraint"
mkdir -p "$output"

"$HGT_BIN/iqtree" -s "$input/trimmed_alignment.faa" -m "$model" \
  -g "$constraint" -seed 20261020 -T "${SLURM_CPUS_PER_TASK:?}" \
  --prefix "$output/constrained"
test -s "$output/constrained.treefile"
cat "$input/ML.treefile" "$output/constrained.treefile" > "$output/candidate_trees.nwk"
"$HGT_BIN/iqtree" -s "$input/trimmed_alignment.faa" -m "$model" \
  -z "$output/candidate_trees.nwk" -te "$input/ML.treefile" \
  -zb 10000 -au -seed 20261021 -T "$SLURM_CPUS_PER_TASK" \
  --prefix "$output/AU"
