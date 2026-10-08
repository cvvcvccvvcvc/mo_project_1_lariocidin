#!/usr/bin/env python3
"""Extract translated CDS from the bounded same-species assembly controls."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemblies", type=Path, required=True)
    parser.add_argument("--fasta-out", type=Path, required=True)
    parser.add_argument("--features-out", type=Path, required=True)
    args = parser.parse_args()

    assemblies = pd.read_csv(args.assemblies, sep="\t", keep_default_na=False)
    assemblies = assemblies.loc[assemblies["genome_role"] == "same_species_control"]
    root = args.assemblies.resolve().parents[1]
    records = []
    rows = []
    for assembly in assemblies.itertuples(index=False):
        index = 0
        with gzip.open(root / assembly.local_gbff, "rt") as handle:
            for contig in SeqIO.parse(handle, "genbank"):
                for feature in contig.features:
                    if feature.type != "CDS" or "translation" not in feature.qualifiers:
                        continue
                    index += 1
                    feature_id = f"{assembly.sample_id}__cds{index:05d}"
                    sequence = feature.qualifiers["translation"][0]
                    records.append(SeqRecord(Seq(sequence), id=feature_id, description=""))
                    rows.append({
                        "feature_id": feature_id,
                        "sample_id": assembly.sample_id,
                        "assembly_accession": assembly.assembly_accession,
                        "contig_accession": contig.id,
                        "cds_index": index,
                        "start_0based": int(feature.location.start),
                        "end_0based_exclusive": int(feature.location.end),
                        "strand": int(feature.location.strand or 0),
                        "aa_length": len(sequence),
                        "gene": feature.qualifiers.get("gene", [""])[0],
                        "locus_tag": feature.qualifiers.get("locus_tag", [""])[0],
                        "protein_id": feature.qualifiers.get("protein_id", [""])[0],
                        "product": feature.qualifiers.get("product", [""])[0],
                    })
    args.fasta_out.parent.mkdir(parents=True, exist_ok=True)
    args.features_out.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, args.fasta_out, "fasta")
    pd.DataFrame(rows).to_csv(args.features_out, sep="\t", index=False)
    print(f"extracted {len(records)} CDS from {len(assemblies)} controls")


if __name__ == "__main__":
    main()
