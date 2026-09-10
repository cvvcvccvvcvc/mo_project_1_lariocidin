#!/usr/bin/env python3
"""Validate an explicit LrcC selection and produce analysis-ready FASTA/table files."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--all-cds", type=Path, required=True)
    parser.add_argument("--ranking", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--outgroup", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t", keep_default_na=False)
    features = pd.read_csv(args.features, sep="\t", keep_default_na=False)
    ranking = pd.read_csv(args.ranking, sep="\t", keep_default_na=False)
    selection = pd.read_csv(args.selection, sep="\t", keep_default_na=False)

    if selection["sample_id"].duplicated().any() or selection["feature_id"].duplicated().any():
        raise SystemExit("selection must contain unique sample_id and feature_id values")
    if set(selection["sample_id"]) != set(manifest["sample_id"]):
        raise SystemExit("selection must contain exactly one row for every manifest sample")

    selected = selection.merge(features, on=["sample_id", "feature_id"], how="left", validate="one_to_one")
    if selected["file_name"].eq("").any() or selected["file_name"].isna().any():
        raise SystemExit("at least one selected feature does not exist")
    selected = selected.merge(
        manifest[["sample_id", "display_name", "accession", "phylum", "group"]],
        on="sample_id", how="left", validate="many_to_one",
    )
    best_scores = ranking.loc[ranking["rank_within_bgc"] == 1, [
        "sample_id", "feature_id", "pident", "alignment_length", "query_coverage", "evalue", "bitscore"
    ]].rename(columns={"feature_id": "top_blast_feature_id"})
    selected = selected.merge(best_scores, on="sample_id", how="left", validate="one_to_one")
    selected["is_top_blast_hit"] = selected["feature_id"] == selected["top_blast_feature_id"]

    seqs = SeqIO.index(str(args.all_cds), "fasta")
    records: list[SeqRecord] = []
    for row in selected.sort_values("sample_id").itertuples():
        if row.feature_id not in seqs:
            raise SystemExit(f"selected sequence missing from FASTA: {row.feature_id}")
        records.append(SeqRecord(seqs[row.feature_id].seq, id=row.sample_id, description=""))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, args.out_dir / "lrcc_reference.faa", "fasta")
    outgroup_records = list(SeqIO.parse(args.outgroup, "fasta"))
    SeqIO.write(records + outgroup_records, args.out_dir / "lrcc_with_outgroup.faa", "fasta")
    canonical_ids = set(selected.loc[selected["phylum"] != "Actinomycetota", "sample_id"])
    canonical_records = [record for record in records if record.id in canonical_ids]
    SeqIO.write(
        canonical_records + outgroup_records,
        args.out_dir / "lrcc_canonical_with_outgroup.faa",
        "fasta",
    )
    selected.to_csv(args.out_dir / "lrcc_selected.tsv", sep="\t", index=False)
    print(
        f"validated {len(records)} LrcC selections; exported full and "
        f"{len(canonical_records)}-sequence canonical subsets plus outgroup"
    )


if __name__ == "__main__":
    main()
