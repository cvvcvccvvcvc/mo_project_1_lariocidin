#!/usr/bin/env python3
"""Join BLASTP results to CDS metadata and retain the best candidates per BGC."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd


BLAST_FIELDS = ["feature_id", "seed_id", "pident", "alignment_length", "query_length", "seed_length", "evalue", "bitscore", "query_coverage"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--top", type=int, default=3)
    args = parser.parse_args()

    features = pd.read_csv(args.features, sep="\t", keep_default_na=False)
    blast = pd.read_csv(args.blast, sep="\t", names=BLAST_FIELDS)
    blast = (
        blast.sort_values(["feature_id", "bitscore"], ascending=[True, False])
        .drop_duplicates("feature_id", keep="first")
    )
    joined = features.merge(blast, on="feature_id", how="left")
    joined["bitscore"] = pd.to_numeric(joined["bitscore"], errors="coerce").fillna(-1)
    ranked = (
        joined.sort_values(["sample_id", "bitscore", "aa_length"], ascending=[True, False, False])
        .groupby("sample_id", sort=False)
        .head(args.top)
        .copy()
    )
    ranked["rank_within_bgc"] = ranked.groupby("sample_id").cumcount() + 1
    cols = [
        "sample_id", "rank_within_bgc", "feature_id", "cds_index", "aa_length",
        "label", "product", "locus_tag", "protein_id", "pident",
        "alignment_length", "query_coverage", "evalue", "bitscore",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    ranked[cols].to_csv(args.out, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"wrote {len(ranked)} ranked rows for {ranked['sample_id'].nunique()} BGC records")


if __name__ == "__main__":
    main()
