#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-hgt-H1-B
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:15:00
set -euo pipefail
: "${HGT_ROOT:?Set isolated HGT root}"
: "${HGT_BIN:?Set HGT_BIN to the pinned environment bin directory}"
export HGT_BIN
cd "$HGT_ROOT"
sha256sum -c config/H1_B_inputs.sha256 > results/slurm/B_input_hash_check.txt
"$HGT_BIN/blastp" -version | head -2
mkdir -p results/H1_B_homology
"$HGT_BIN/makeblastdb" -in results/H1_sequences/B_author_references.faa -dbtype prot \
    -out results/H1_B_homology/references > results/H1_B_homology/makeblastdb.log
"$HGT_BIN/blastp" -query results/H1_sequences/canonical_component_queries.faa \
    -db results/H1_B_homology/references -evalue 1e-5 -qcov_hsp_perc 70 \
    -max_target_seqs 6 -num_threads 2 \
    -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen' \
    -out results/tables/H1_B_homology.tsv
