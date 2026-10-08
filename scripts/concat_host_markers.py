#!/usr/bin/env python3
"""Concatenate retained marker alignments and record their boundaries."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment-dir", type=Path, required=True)
    parser.add_argument("--taxa", type=Path, required=True)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--boundaries", type=Path, required=True)
    args = parser.parse_args()

    taxa_table = pd.read_csv(args.taxa, sep="\t", keep_default_na=False)
    taxa = taxa_table.loc[taxa_table["host_tree_status"] == "included", "sample_id"].tolist()
    paths = sorted(args.alignment_dir.glob("*.trim.faa"))
    inventory = pd.read_csv(args.inventory, sep="\t")
    retained = set(inventory.loc[inventory.marker_retained, "marker"])
    if {path.name.removesuffix(".trim.faa") for path in paths} != retained:
        raise ValueError("Missing or stale marker alignments relative to the extraction inventory")
    if not paths:
        raise ValueError("No trimmed marker alignments")
    concatenated = {taxon: "" for taxon in taxa}
    boundaries = []
    start = 1
    for path in paths:
        marker = path.name.removesuffix(".trim.faa")
        alignment = AlignIO.read(path, "fasta")
        sequences = {record.id: str(record.seq) for record in alignment}
        if len(sequences) != len(alignment):
            raise ValueError(f"Duplicate sequence IDs in {path}")
        unexpected = set(sequences) - set(taxa)
        if unexpected:
            raise ValueError(f"Unexpected taxa in {path}: {sorted(unexpected)}")
        length = alignment.get_alignment_length()
        for taxon in taxa:
            concatenated[taxon] += sequences.get(taxon, "-" * length)
        boundaries.append({"marker": marker, "start_1based": start, "end_1based": start + length - 1})
        start += length

    args.output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(
        [SeqRecord(Seq(concatenated[taxon]), id=taxon, description="") for taxon in taxa],
        args.output,
        "fasta",
    )
    pd.DataFrame(boundaries).to_csv(args.boundaries, sep="\t", index=False)
    print(f"concatenated {len(paths)} markers, {start - 1} columns, {len(taxa)} taxa")


if __name__ == "__main__":
    main()
