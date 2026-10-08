#!/usr/bin/env python3
"""Run the fixed final LrcC profiles on the bounded host-tree controls."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from profile_search import search  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--hits-out", type=Path, required=True)
    args = parser.parse_args()

    features = pd.read_csv(args.features, sep="\t", keep_default_na=False)
    hits = search(args.fasta, args.output_dir, exhaustive=True)
    if not hits.empty:
        hits = hits.merge(features, left_on="target_id", right_on="feature_id", validate="many_to_one")
        hits = hits.sort_values(["sample_id", "sequence_score"], ascending=[True, False])
    hits.to_csv(args.hits_out, sep="\t", index=False)
    print(f"profile hits >=80 bits: {len(hits)} across {hits['sample_id'].nunique() if len(hits) else 0} controls")


if __name__ == "__main__":
    main()
