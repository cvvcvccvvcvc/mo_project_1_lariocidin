#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-hgt-H2-components
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:45:00
set -euo pipefail
: "${HGT_ROOT:?Isolated H2 root required}"
export PYTHONDONTWRITEBYTECODE=1
: "${HGT_GENOME_ROOT:?Set HGT_GENOME_ROOT to the frozen H1 genome directory}"
export HGT_GENOME_ROOT
: "${HGT_BIN:?Set HGT_BIN to the pinned environment bin directory}"
export HGT_BIN
"$HGT_BIN/hmmsearch" --cpu 2 --cut_ga \
    --tblout "$HGT_ROOT/results/H2/publication_domains.tbl" \
    --domtblout "$HGT_ROOT/results/H2/publication_domains.domtbl" \
    "$HGT_ROOT/inputs/H2_2026-10-13/references/BGC_domains.hmm" \
    "$HGT_ROOT/inputs/H2_2026-10-13/H1/proteins.faa" \
    > "$HGT_ROOT/results/H2/publication_domains.log"
"$HGT_BIN/python" "$HGT_ROOT/scripts/steps/h2_components.py"
