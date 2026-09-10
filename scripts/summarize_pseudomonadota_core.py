#!/usr/bin/env python3
"""Compare each core protein to same-order and cross-order homologues."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


FIELDS = [
    "query_id", "subject_id", "pident", "alignment_length", "query_length",
    "subject_length", "evalue", "bitscore", "query_coverage_all_hsps",
]
FOCAL_QUERY = "photorhabdus_asymbiotica_dsm15149"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--roles", type=Path, required=True)
    parser.add_argument("--blast-dir", type=Path, required=True)
    parser.add_argument("--all-pairs-out", type=Path, required=True)
    parser.add_argument("--focal-out", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, required=True)
    args = parser.parse_args()

    roles = pd.read_csv(args.roles, sep="\t", keep_default_na=False)
    group = roles.drop_duplicates("sample_id").set_index("sample_id")["taxonomic_group"]
    rows = []
    for role in ("A", "B1", "B2", "C", "D"):
        blast = pd.read_csv(args.blast_dir / f"pseudomonadota_core_{role}.tsv", sep="\t", names=FIELDS)
        blast = blast.loc[blast["query_id"] != blast["subject_id"]].copy()
        blast["role"] = role
        blast["query_group"] = blast["query_id"].map(group)
        blast["subject_group"] = blast["subject_id"].map(group)
        blast = (
            blast.sort_values(["query_id", "subject_id", "bitscore"], ascending=[True, True, False])
            .drop_duplicates(["query_id", "subject_id"], keep="first")
        )
        rows.append(blast)
    all_pairs = pd.concat(rows, ignore_index=True)
    args.all_pairs_out.parent.mkdir(parents=True, exist_ok=True)
    all_pairs.to_csv(args.all_pairs_out, sep="\t", index=False)

    focal = all_pairs.loc[all_pairs["query_id"] == FOCAL_QUERY].copy()
    focal = (
        focal.sort_values(["role", "subject_group", "bitscore"], ascending=[True, True, False])
        .drop_duplicates(["role", "subject_group"], keep="first")
        .sort_values(["role", "subject_group"])
    )
    focal.to_csv(args.focal_out, sep="\t", index=False)

    winners = (
        all_pairs.sort_values(["query_id", "role", "bitscore"], ascending=[True, True, False])
        .drop_duplicates(["query_id", "role"], keep="first")
    )
    focal_winners = winners.loc[winners["query_id"] == FOCAL_QUERY]
    summary = {
        "status": "descriptive_not_independent_HGT_test",
        "role_assignment": "narrow canonical-order/size heuristic; inspect role table",
        "loci": int(roles["sample_id"].nunique()),
        "complete_role_sets": int((roles.groupby("sample_id")["role"].nunique() == 5).sum()),
        "focal_query": FOCAL_QUERY,
        "focal_nearest_nonself_by_role": {
            row.role: {
                "sample_id": row.subject_id,
                "taxonomic_group": row.subject_group,
                "percent_identity": float(row.pident),
                "query_coverage_percent": float(row.query_coverage_all_hsps),
                "bitscore": float(row.bitscore),
            }
            for row in focal_winners.itertuples(index=False)
        },
        "focal_roles_whose_nearest_sequence_is_Burkholderiales": int(
            (focal_winners["subject_group"] == "Burkholderiales").sum()
        ),
        "limitations": [
            "Genes in one BGC are linked and are not independent statistical replicates.",
            "Nearest-sequence affinity is descriptive and does not distinguish HGT from vertical inheritance plus loss.",
            "The author-selected 14-locus set lacks many close host relatives.",
        ],
    }
    args.summary_out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
