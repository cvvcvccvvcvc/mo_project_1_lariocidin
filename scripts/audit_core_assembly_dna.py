#!/usr/bin/env python3
"""Locate reference BGC proteins without exact CDS matches in source assembly DNA."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable


ROLES = ("A", "B1", "B2", "C", "D")
START_CODONS = set(CodonTable.unambiguous_dna_by_id[11].start_codons)


def occurrences(haystack: str, needle: str):
    start = 0
    while (index := haystack.find(needle, start)) >= 0:
        yield index
        start = index + 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemblies", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--core-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = args.assemblies.resolve().parents[1]
    with args.assemblies.open(newline="") as handle:
        assemblies = {row["sample_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with args.validation.open(newline="") as handle:
        unmatched = [row for row in csv.DictReader(handle, delimiter="\t")
                     if row["exact_translation_matches"] == "0"]
    peptides = {}
    for role in ROLES:
        for record in SeqIO.parse(args.core_dir / f"pseudomonadota_core_{role}.faa", "fasta"):
            peptides[record.id, role] = str(record.seq).rstrip("*")

    by_sample = defaultdict(list)
    for row in unmatched:
        by_sample[row["sample_id"]].append(row)
    results = []
    for sample_id, rows in by_sample.items():
        hits = {(sample_id, row["role"]): [] for row in rows}
        assembly = assemblies[sample_id]
        with gzip.open(root / assembly["local_gbff"], "rt") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for strand, dna in (("+", record.seq), ("-", record.seq.reverse_complement())):
                    for frame in range(3):
                        length = (len(dna) - frame) // 3 * 3
                        amino_acids = str(dna[frame:frame + length].translate())
                        for key in hits:
                            query = peptides[key]
                            exact = list(occurrences(amino_acids, query))
                            # NCBI bacterial CDS translation renders an alternative
                            # initiation codon as M even when its ordinary translation
                            # is another residue (TTG translates to L internally).
                            alternative = []
                            if query.startswith("M"):
                                for pos in occurrences(amino_acids, query[1:]):
                                    if pos < 1:
                                        continue
                                    if amino_acids[pos - 1] == "M":
                                        continue  # already counted as a standard exact match
                                    codon = str(dna[frame + 3 * (pos - 1):frame + 3 * pos])
                                    if codon in START_CODONS:
                                        alternative.append((pos - 1, codon))
                            for pos, codon in [(p, "") for p in exact] + alternative:
                                begin = frame + 3 * pos
                                end = begin + 3 * len(query)
                                if strand == "-":
                                    begin, end = len(record) - end, len(record) - begin
                                hits[key].append({
                                    "sample_id": sample_id,
                                    "assembly_accession": assembly["assembly_accession"],
                                    "role": key[1],
                                    "dna_status": "exact_standard_translation" if not codon else "exact_with_alternative_start",
                                    "contig_accession": record.id,
                                    "start_0based": begin,
                                    "end_0based_exclusive": end,
                                    "strand": strand,
                                    "start_codon": codon,
                                })
        for row in rows:
            found = hits[sample_id, row["role"]]
            # A standard exact match takes precedence over an alternative start.
            found.sort(key=lambda hit: hit["dna_status"] != "exact_standard_translation")
            chosen = found[0] if found else {
                "sample_id": sample_id,
                "assembly_accession": assembly["assembly_accession"],
                "role": row["role"],
                "dna_status": "not_found",
                "contig_accession": "", "start_0based": "", "end_0based_exclusive": "",
                "strand": "", "start_codon": "",
            }
            chosen["matching_loci"] = len(found)
            results.append(chosen)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ("sample_id", "assembly_accession", "role", "dna_status", "contig_accession",
              "start_0based", "end_0based_exclusive", "strand", "start_codon", "matching_loci")
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(results, key=lambda row: (row["sample_id"], ROLES.index(row["role"]))))
    print(f"Audited {len(results)} role sequences with no exact CDS translation match")


if __name__ == "__main__":
    main()
