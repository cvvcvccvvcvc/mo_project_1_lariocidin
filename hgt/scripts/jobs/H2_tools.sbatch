#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-hgt-H2-tools
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
set -euo pipefail
: "${HGT_ROOT:?Isolated HGT root required}"
export MAMBA_ROOT_PREFIX="$HGT_ROOT/tools/mamba"
export PYTHONDONTWRITEBYTECODE=1
micromamba --no-rc create -y \
    -p "$HGT_ROOT/tools/annotation-trimming" -c conda-forge -c bioconda \
    prodigal=2.6.3 trimal=1.4.1
micromamba --no-rc list \
    -p "$HGT_ROOT/tools/annotation-trimming" --explicit > "$HGT_ROOT/config/H2/H2_tools_explicit.txt"
: "${HGT_BIN:?Set HGT_BIN to the pinned environment bin directory}"
export HGT_BIN
"$HGT_BIN/hmmsearch" -h | awk 'NR<=5'
"$HGT_BIN/mafft" --version
"$HGT_BIN/iqtree" --version
"$HGT_BIN/blastp" -version
"$HGT_BIN/makeblastdb" -version
"$HGT_ROOT/tools/annotation-trimming/bin/prodigal" -v
"$HGT_ROOT/tools/annotation-trimming/bin/trimal" --version
