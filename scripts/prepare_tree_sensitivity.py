#!/usr/bin/env python3
"""Prepare the prespecified local and occupancy-masked LrcC alignments."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--phylum", default="Pseudomonadota")
    parser.add_argument("--subset-fasta", type=Path, required=True)
    parser.add_argument("--aligned", type=Path)
    parser.add_argument("--masked", type=Path)
    parser.add_argument("--occupancy", type=float, default=0.8)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep="\t", keep_default_na=False)
    selected = set(metadata.loc[metadata["phylum"] == args.phylum, "sample_id"])
    records = [record for record in SeqIO.parse(args.fasta, "fasta") if record.id in selected]
    if {record.id for record in records} != selected:
        raise ValueError("Metadata/FASTA mismatch for the selected phylum")
    args.subset_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, args.subset_fasta, "fasta")

    if args.aligned is None and args.masked is None:
        print(f"sequences={len(records)}")
        return
    if args.aligned is None or args.masked is None:
        raise ValueError("--aligned and --masked must be supplied together")

    alignment = AlignIO.read(args.aligned, "fasta")
    minimum_non_gap = math.ceil(args.occupancy * len(alignment))
    keep = [
        index
        for index in range(alignment.get_alignment_length())
        if sum(record.seq[index] not in "-.?" for record in alignment) >= minimum_non_gap
    ]
    masked_records = [
        SeqRecord(
            Seq("".join(str(record.seq[index]) for index in keep)),
            id=record.id,
            description="",
        )
        for record in alignment
    ]
    SeqIO.write(masked_records, args.masked, "fasta")
    print(
        f"sequences={len(records)} columns={alignment.get_alignment_length()} "
        f"masked_columns={len(keep)} minimum_non_gap={minimum_non_gap}"
    )


if __name__ == "__main__":
    main()
