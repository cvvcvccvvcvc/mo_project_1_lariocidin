#!/usr/bin/env bash
set -euo pipefail
HGT_ROOT=$(cd "$(dirname "$0")/.." && pwd)
: "${HGT_BIN:?Set HGT_BIN to the pinned environment bin directory}"
HGT_SEQUENCES="$HGT_ROOT/results/H1_sequences"
"$HGT_BIN/makeblastdb" -in "$HGT_SEQUENCES/merrild_author_C.faa" -dbtype prot \
    -out "$HGT_SEQUENCES/author_C_db" > "$HGT_SEQUENCES/blast_build.log" 2>&1
"$HGT_BIN/blastp" -query "$HGT_SEQUENCES/proteins.faa" -db "$HGT_SEQUENCES/author_C_db" \
    -evalue 1e-10 -qcov_hsp_perc 70 -max_target_seqs 5 -num_threads 4 \
    -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen' \
    -out "$HGT_ROOT/results/tables/H1_C_homology.tsv" > "$HGT_SEQUENCES/blast_C.log" 2>&1
