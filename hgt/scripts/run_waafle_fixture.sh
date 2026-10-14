#!/usr/bin/env bash
set -euo pipefail
HGT_ROOT=$(cd "$(dirname "$0")/.." && pwd)
HGT_PYTHON=${1:-"${HGT_BIN:?Set HGT_BIN or pass a Python executable}/python"}
HGT_TOOL="$HGT_ROOT/inputs/2026-10-13/method_sources/waafle_v1.0.0"
HGT_DEMO="$HGT_ROOT/inputs/2026-10-13/method_sources/waafle/demo"
HGT_OUTPUT="$HGT_ROOT/results/waafle_fixture_v1.0.0"
mkdir -p "$HGT_OUTPUT"
export PYTHONPATH="$HGT_TOOL${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONDONTWRITEBYTECODE=1
"$HGT_PYTHON" -m waafle.waafle_genecaller \
    "$HGT_DEMO/output/demo_contigs.blastout" \
    --gff "$HGT_OUTPUT/demo_contigs.gff" > "$HGT_OUTPUT/genecaller.log" 2>&1
"$HGT_PYTHON" -m waafle.waafle_orgscorer \
    "$HGT_DEMO/input/demo_contigs.fna" \
    "$HGT_DEMO/output/demo_contigs.blastout" \
    "$HGT_OUTPUT/demo_contigs.gff" \
    "$HGT_DEMO/input/demo_taxonomy.tsv" \
    --outdir "$HGT_OUTPUT" --quiet > "$HGT_OUTPUT/orgscorer.log" 2>&1
# ORF input can be checked against the author's archived coordinates.
# Scoring tables are measured outputs: their discrepancy is reported in the review.
cmp "$HGT_OUTPUT/demo_contigs.gff" "$HGT_TOOL/demo/output/demo_contigs.gff"
echo "WAAFLE v1.0.0 completed; gene coordinates match the author fixture. Scoring differs from the archived tables; see the review. Search and read QC were not run."
