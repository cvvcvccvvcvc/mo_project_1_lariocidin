#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-H3-H0163-reads
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail

: "${HGT_SOURCE_ROOT:?Set frozen H2 root}"
: "${HGT_RESULT_ROOT:?Set separate H3 output root}"

output="$HGT_RESULT_ROOT/read_validation/H0163_SRR7457462"
genome="$HGT_SOURCE_ROOT/results/H2/native/H0163/genome.fna"
mkdir -p "$output"

for mate in 1 2; do
  curl -fL --retry 3 --silent --show-error \
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR745/002/SRR7457462/SRR7457462_${mate}.fastq.gz" \
    -o "$output/R${mate}.fastq.gz"
done
(
  cd "$output"
  printf '%s  %s\n' '33bbb5e36c87ffee20b6c654eef83cc6' R1.fastq.gz \
    '2e6a28e06e1711faa287e814b153c701' R2.fastq.gz | md5sum -c -
)

bwa index -p "$output/H0163_index" "$genome"
printf 'NZ_QTQP01000053.1\t49900\t50550\n' > "$output/target.bed"
bwa mem -t 4 "$output/H0163_index" "$output/R1.fastq.gz" "$output/R2.fastq.gz" \
  | samtools view -b -h -L "$output/target.bed" -F 0x900 - \
  | samtools sort -@ 4 -m 1G -o "$output/target.sorted.bam" -
samtools index -@ 4 "$output/target.sorted.bam"
samtools flagstat "$output/target.sorted.bam" > "$output/target.flagstat.txt"
samtools mpileup -aa -q 20 -Q 20 -f "$genome" \
  -r NZ_QTQP01000053.1:50000-50452 "$output/target.sorted.bam" \
  > "$output/target.mpileup.tsv"
