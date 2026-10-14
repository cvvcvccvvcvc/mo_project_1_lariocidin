#!/usr/bin/env python3
"""Infer A/B1/B2/C/D roles for the 14 canonical Pseudomonadota BGCs.

The heuristic is deliberately narrow and is recorded in the output.  C comes from
the curated reference selection; A is the nearest preceding 20--60 aa CDS; B1/B2
are the final two intervening CDS; D is the first following 450--650 aa CDS on
the same strand.  CDS order is used because valid neighbouring annotations can
overlap.  This is suitable for these author-selected canonical loci, not a general
BGC caller.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


ROLES = ("A", "B1", "B2", "C", "D")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--all-cds", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--table-out", type=Path, required=True)
    parser.add_argument("--fasta-dir", type=Path, required=True)
    args = parser.parse_args()

    features = pd.read_csv(args.features, sep="\t", keep_default_na=False)
    selected = pd.read_csv(args.selected, sep="\t", keep_default_na=False)
    selected = selected.loc[selected["phylum"] == "Pseudomonadota"].copy()
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(args.all_cds, "fasta")}

    rows: list[dict[str, object]] = []
    for selection in selected.itertuples(index=False):
        locus = features.loc[features["sample_id"] == selection.sample_id].sort_values("cds_index")
        c_row = locus.loc[locus["feature_id"] == selection.feature_id]
        if len(c_row) != 1:
            raise ValueError(f"Expected one C for {selection.sample_id}")
        c_row = c_row.iloc[0]
        if int(c_row["strand"]) != 1:
            raise ValueError(f"Heuristic assumes forward canonical locus: {selection.sample_id}")

        upstream = locus.loc[
            (locus["cds_index"] < c_row["cds_index"])
            & locus["aa_length"].between(20, 60)
            & (locus["strand"] == c_row["strand"])
        ]
        if upstream.empty:
            raise ValueError(f"No precursor-sized upstream ORF for {selection.sample_id}")
        a_row = upstream.sort_values("cds_index").iloc[-1]
        intervening = locus.loc[
            (locus["cds_index"] > a_row["cds_index"])
            & (locus["cds_index"] < c_row["cds_index"])
            & (locus["strand"] == c_row["strand"])
        ].sort_values("cds_index")
        if len(intervening) < 2:
            raise ValueError(f"Fewer than two B candidates for {selection.sample_id}")
        b1_row, b2_row = intervening.iloc[-2], intervening.iloc[-1]
        downstream = locus.loc[
            (locus["cds_index"] > c_row["cds_index"])
            & locus["aa_length"].between(450, 650)
            & (locus["strand"] == c_row["strand"])
        ].sort_values("cds_index")
        if downstream.empty:
            raise ValueError(f"No D-sized downstream CDS for {selection.sample_id}")
        d_row = downstream.iloc[0]

        for role, feature_row, basis in [
            ("A", a_row, "nearest_preceding_20_60aa_cds_same_strand"),
            ("B1", b1_row, "penultimate_intervening_cds_after_A"),
            ("B2", b2_row, "last_intervening_cds_before_C"),
            ("C", c_row, "curated_lrcC_selection"),
            ("D", d_row, "first_following_450_650aa_cds_same_strand"),
        ]:
            feature_id = str(feature_row["feature_id"])
            rows.append(
                {
                    "sample_id": selection.sample_id,
                    "display_name": selection.display_name,
                    "taxonomic_group": selection.group,
                    "role": role,
                    "role_assignment_status": (
                        "curated_reference" if role == "C" else "narrow_heuristic_inference"
                    ),
                    "role_assignment_basis": basis,
                    "feature_id": feature_id,
                    "cds_index": int(feature_row["cds_index"]),
                    "start_0based": int(feature_row["start_0based"]),
                    "end_0based_exclusive": int(feature_row["end_0based_exclusive"]),
                    "aa_length": int(feature_row["aa_length"]),
                    "product": feature_row["product"],
                    "sequence": sequences[feature_id],
                }
            )

    table = pd.DataFrame(rows).sort_values(["sample_id", "role"])
    if len(table) != 14 * len(ROLES):
        raise ValueError("Expected 14 complete five-gene role assignments")
    args.table_out.parent.mkdir(parents=True, exist_ok=True)
    table.drop(columns="sequence").to_csv(args.table_out, sep="\t", index=False)
    args.fasta_dir.mkdir(parents=True, exist_ok=True)
    for role in ROLES:
        subset = table.loc[table["role"] == role]
        SeqIO.write(
            [
                SeqRecord(Seq(row.sequence), id=row.sample_id, description=f"role={role}")
                for row in subset.itertuples(index=False)
            ],
            args.fasta_dir / f"pseudomonadota_core_{role}.faa",
            "fasta",
        )
    print(f"assigned {len(table)} roles across {table['sample_id'].nunique()} loci")


if __name__ == "__main__":
    main()
