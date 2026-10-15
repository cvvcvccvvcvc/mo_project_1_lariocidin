#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-H3-host-sCF
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
set -euo pipefail

: "${HGT_SOURCE_ROOT:?Set frozen H2 root}"
: "${HGT_RESULT_ROOT:?Set separate H3 output root}"
: "${HGT_BIN:?Set pinned IQ-TREE environment bin}"

case "${SLURM_ARRAY_TASK_ID:?}" in
  0) panel=Bacteria74_canonical; input=results/H2/tree_inputs/canonical_hosts.faa; tree=results/H2/host_trees/canonical/ML.treefile ;;
  1) panel=Bacteria74_noncanonical; input=results/H2/tree_inputs/noncanonical_hosts.faa; tree=results/H2/host_trees/noncanonical/ML.treefile ;;
  2) panel=Actino138_noncanonical; input=results/H2/actino/tree_inputs/noncanonical_hosts.faa; tree=results/H2/actino/host_trees/noncanonical/ML.treefile ;;
  *) echo "Unknown panel index" >&2; exit 2 ;;
esac

output="$HGT_RESULT_ROOT/site_concordance/$panel"
mkdir -p "$output"
"$HGT_BIN/iqtree" -t "$HGT_SOURCE_ROOT/$tree" -s "$HGT_SOURCE_ROOT/$input" \
  --scf 100 -seed 131015 -T "${SLURM_CPUS_PER_TASK:?}" --prefix "$output/concordance"
