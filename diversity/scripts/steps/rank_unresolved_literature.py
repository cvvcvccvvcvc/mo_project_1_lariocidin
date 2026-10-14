#!/usr/bin/env python3
"""Rank CDS from published BGC records whose cyclase was not explicitly labelled."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd


BLAST_FIELDS = [
    "feature_id", "seed_id", "pident", "alignment_length", "query_length",
    "seed_length", "evalue", "bitscore", "query_coverage",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--candidates-out", type=Path, required=True)
    parser.add_argument("--top", type=int, default=4)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t", keep_default_na=False)
    unresolved_ids = set(manifest.loc[manifest["author_cyclase_count"] == 0, "bgc_id"])
    features = pd.read_csv(args.features, sep="\t", keep_default_na=False)
    features = features.loc[features["bgc_id"].isin(unresolved_ids)].copy()
    blast = pd.read_csv(args.blast, sep="\t", names=BLAST_FIELDS)
    best_hit = (
        blast.sort_values(["feature_id", "bitscore"], ascending=[True, False])
        .drop_duplicates("feature_id", keep="first")
    )
    ranked = features.merge(best_hit, on="feature_id", how="left")
    ranked["bitscore"] = pd.to_numeric(ranked["bitscore"], errors="coerce").fillna(-1)
    ranked = (
        ranked.sort_values(["bgc_id", "bitscore", "aa_length"], ascending=[True, False, False])
        .groupby("bgc_id", sort=False)
        .head(args.top)
        .copy()
    )
    ranked["rank_within_bgc"] = ranked.groupby("bgc_id").cumcount() + 1
    ranked["exact_explicit_cyclase_match"] = ranked["sequence_sha256"].map(
        lambda value: f"cyc_{value[:16]}" if value else ""
    ) == ranked["seed_id"]
    columns = [
        "bgc_id", "archive_member", "rank_within_bgc", "feature_id", "cds_index",
        "start_0based", "end_0based_exclusive", "aa_length", "label", "product",
        "protein_id", "seed_id", "exact_explicit_cyclase_match", "pident",
        "alignment_length", "query_coverage", "evalue", "bitscore",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    ranked[columns].to_csv(args.out, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    candidates = ranked.loc[ranked["rank_within_bgc"] == 1, columns].copy()
    candidates.insert(3, "annotation_status", "inferred_not_author_labelled")
    candidates.insert(
        4,
        "candidate_evidence",
        candidates.apply(
            lambda row: "exact_sequence_in_explicit_cyclase_set"
            if row["exact_explicit_cyclase_match"]
            else (
                "high_similarity_to_explicit_noncanonical_cyclase"
                if row["pident"] >= 80 and row["query_coverage"] >= 80
                else "requires_manual_review"
            ),
            axis=1,
        ),
    )
    candidates.to_csv(
        args.candidates_out, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL
    )
    assert ranked["bgc_id"].nunique() == len(unresolved_ids) == 4
    assert (candidates["candidate_evidence"] != "requires_manual_review").all()
    print(f"wrote {len(ranked)} ranked CDS rows for {len(unresolved_ids)} unresolved records")


if __name__ == "__main__":
    main()
