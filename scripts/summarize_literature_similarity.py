#!/usr/bin/env python3
"""Summarize leave-one-out BLASTP similarity within the published cyclase set.

This is a calibration diagnostic, not a trained family classifier.  It asks whether
the author-defined source groups are recoverable from cyclase sequence similarity
alone and records where genomic context will be required.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


BLAST_FIELDS = [
    "query_id",
    "subject_id",
    "pident",
    "alignment_length",
    "query_length",
    "subject_length",
    "evalue",
    "bitscore",
    "query_coverage_all_hsps",
    "query_coverage_hsp",
]


def family_class(source_group: str) -> str:
    if source_group in {"canonical", "noncanonical"}:
        return "lariocidin_triculamin_like"
    return "other_lasso_comparator"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cyclases", type=Path, required=True)
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--nearest-out", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, required=True)
    args = parser.parse_args()

    cyclases = pd.read_csv(args.cyclases, sep="\t", keep_default_na=False)
    if cyclases["source_groups"].str.contains(";").any():
        raise ValueError("Expected each unique sequence to occur in one source group")
    group_by_id = cyclases.set_index("cyclase_id")["source_groups"].to_dict()

    blast = pd.read_csv(args.blast, sep="\t", names=BLAST_FIELDS)
    for column in [
        "pident",
        "alignment_length",
        "query_length",
        "subject_length",
        "evalue",
        "bitscore",
        "query_coverage_all_hsps",
        "query_coverage_hsp",
    ]:
        blast[column] = pd.to_numeric(blast[column], errors="raise")
    blast["query_group"] = blast["query_id"].map(group_by_id)
    blast["subject_group"] = blast["subject_id"].map(group_by_id)
    if blast[["query_group", "subject_group"]].isna().any().any():
        raise ValueError("BLAST output contains an ID absent from the cyclase table")

    # BLAST may report more than one HSP for a sequence pair.  Keep the strongest
    # HSP while retaining qcovs, which represents coverage summed over all HSPs.
    pair_best = (
        blast.sort_values(
            ["query_id", "subject_id", "bitscore", "query_coverage_all_hsps"],
            ascending=[True, True, False, False],
        )
        .drop_duplicates(["query_id", "subject_id"], keep="first")
        .copy()
    )
    self_scores = (
        pair_best.loc[pair_best["query_id"] == pair_best["subject_id"]]
        .set_index("query_id")["bitscore"]
        .to_dict()
    )
    nonself = pair_best.loc[pair_best["query_id"] != pair_best["subject_id"]].copy()
    nonself["query_family_class"] = nonself["query_group"].map(family_class)
    nonself["subject_family_class"] = nonself["subject_group"].map(family_class)
    nonself["bitscore_fraction_of_self"] = nonself.apply(
        lambda row: row["bitscore"] / self_scores[row["query_id"]], axis=1
    )

    nearest = (
        nonself.sort_values(
            ["query_id", "subject_group", "bitscore", "query_coverage_all_hsps"],
            ascending=[True, True, False, False],
        )
        .drop_duplicates(["query_id", "subject_group"], keep="first")
        .sort_values(["query_id", "subject_group"])
    )
    nearest_columns = [
        "query_id",
        "query_group",
        "query_family_class",
        "subject_group",
        "subject_family_class",
        "subject_id",
        "pident",
        "alignment_length",
        "query_length",
        "subject_length",
        "query_coverage_all_hsps",
        "query_coverage_hsp",
        "evalue",
        "bitscore",
        "bitscore_fraction_of_self",
    ]
    args.nearest_out.parent.mkdir(parents=True, exist_ok=True)
    nearest[nearest_columns].to_csv(args.nearest_out, sep="\t", index=False)

    nearest_overall = (
        nonself.sort_values(
            ["query_id", "bitscore", "query_coverage_all_hsps"],
            ascending=[True, False, False],
        )
        .drop_duplicates("query_id", keep="first")
        .copy()
    )
    nearest_overall["same_source_group"] = (
        nearest_overall["query_group"] == nearest_overall["subject_group"]
    )
    nearest_overall["same_family_class"] = (
        nearest_overall["query_family_class"]
        == nearest_overall["subject_family_class"]
    )

    def confusion(column: str) -> dict[str, dict[str, int]]:
        table = pd.crosstab(nearest_overall["query_group"], nearest_overall[column])
        return {
            str(index): {str(col): int(value) for col, value in row.items()}
            for index, row in table.to_dict(orient="index").items()
        }

    positive = nearest_overall[
        nearest_overall["query_family_class"] == "lariocidin_triculamin_like"
    ]
    comparator = nearest_overall[
        nearest_overall["query_family_class"] == "other_lasso_comparator"
    ]
    summary = {
        "interpretation": (
            "Leave-one-out BLASTP calibration only; source labels derive from the "
            "published archive and are not independent experimental truth."
        ),
        "unique_sequences": int(len(cyclases)),
        "source_group_counts": {
            str(k): int(v) for k, v in cyclases["source_groups"].value_counts().items()
        },
        "cross_source_group_exact_sequence_duplicates": int(
            cyclases["source_groups"].str.contains(";").sum()
        ),
        "queries_with_nonself_hit": int(nearest_overall["query_id"].nunique()),
        "nearest_nonself_same_source_group_fraction": float(
            nearest_overall["same_source_group"].mean()
        ),
        "nearest_nonself_same_family_class_fraction": float(
            nearest_overall["same_family_class"].mean()
        ),
        "nearest_subject_group_by_query_group": confusion("subject_group"),
        "positive_queries_nearest_comparator": int(
            (positive["subject_family_class"] == "other_lasso_comparator").sum()
        ),
        "comparator_queries_nearest_positive": int(
            (comparator["subject_family_class"] == "lariocidin_triculamin_like").sum()
        ),
    }
    args.summary_out.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
