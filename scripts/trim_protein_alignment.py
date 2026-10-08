#!/usr/bin/env python3
"""Remove alignment columns below a fixed non-gap occupancy threshold."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-occupancy", type=float, default=0.8)
    args = parser.parse_args()

    alignment = AlignIO.read(args.input, "fasta")
    keep = []
    for column in range(alignment.get_alignment_length()):
        occupied = sum(record.seq[column] not in "-.?" for record in alignment)
        if occupied / len(alignment) >= args.min_occupancy:
            keep.append(column)
    if not keep:
        raise ValueError(f"No columns retained from {args.input}")
    records = [
        SeqRecord(Seq("".join(str(record.seq[index]) for index in keep)), id=record.id, description="")
        for record in alignment
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    AlignIO.write(MultipleSeqAlignment(records), args.output, "fasta")
    print(args.input.stem, alignment.get_alignment_length(), len(keep))


if __name__ == "__main__":
    main()
