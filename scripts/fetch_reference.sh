#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
raw_dir="$project_dir/data/raw"
archive="$raw_dir/LAR-BGC-42280aa.zip"
reference_dir="$raw_dir/reference/LAR-BGC-42280aa5d43ef01f766fa207740dd9914166134c"
mkdir -p "$raw_dir/reference"

if [[ ! -s "$archive" ]]; then
  curl -fsSL \
    https://codeload.github.com/jangrm1/LAR-BGC/zip/42280aa5d43ef01f766fa207740dd9914166134c \
    -o "$archive"
fi
if [[ ! -s "$raw_dir/CP169648.1.gb" ]]; then
  curl -fsSL \
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP169648.1&rettype=gbwithparts&retmode=text' \
    -o "$raw_dir/CP169648.1.gb"
fi
if [[ ! -s "$raw_dir/AHKH01000064.1.gb" ]]; then
  curl -fsSL \
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AHKH01000064.1&rettype=gbwithparts&retmode=text' \
    -o "$raw_dir/AHKH01000064.1.gb"
fi

(cd "$project_dir" && shasum -a 256 -c config/checksums.sha256)
if [[ ! -d "$reference_dir" ]]; then
  unzip -q "$archive" -d "$raw_dir/reference"
fi
test "$(find "$reference_dir" -maxdepth 1 -name '*.gbk' | wc -l | tr -d ' ')" = 29
echo "Reference inputs are present and checksums match."

