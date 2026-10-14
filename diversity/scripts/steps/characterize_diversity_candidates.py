#!/usr/bin/env python3
"""Compare unmatched BGC groups with literature-traced groups in the final catalogue."""

import csv
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results/tables"
PROCESSED = ROOT / "data/processed/discovery"
TRACED = "published_carrier_and_C_DNA_traced"
UNMATCHED = "unmatched_exported_core_literature_scope_limited"
BLAST_FIELDS = (
    "query_id", "subject_id", "identity_percent", "alignment_length", "query_length",
    "subject_length", "evalue", "bitscore", "query_start", "query_end",
    "subject_start", "subject_end",
)


def read_table(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def records(groups, sequences, key):
    output = []
    owners = {}
    for group in groups:
        identifiers = group[key].split(";")
        for identifier in identifiers:
            if not identifier:
                continue
            record_id = f"seq{len(output)}"
            original = sequences[identifier]
            output.append(SeqRecord(original.seq, id=record_id, description=""))
            owners[record_id] = group["merged_locus_id"]
    return output, owners


def nearest(query, subject, query_owner, subject_owner, short, workspace):
    query_path = workspace / "query.faa"
    subject_path = workspace / "subject.faa"
    hits_path = workspace / "hits.tsv"
    SeqIO.write(query, query_path, "fasta")
    SeqIO.write(subject, subject_path, "fasta")
    command = [
        str(ROOT / ".env/bin/blastp"), "-query", str(query_path),
        "-subject", str(subject_path), "-out", str(hits_path),
        "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore qstart qend sstart send",
        "-evalue", "100" if short else "1e-3",
        "-max_target_seqs", str(len(subject)),
    ]
    if short:
        command.extend(["-task", "blastp-short", "-seg", "no"])
    subprocess.run(command, check=True, capture_output=True)
    best = {}
    with hits_path.open(newline="") as handle:
        for cells in csv.reader(handle, delimiter="\t"):
            hit = dict(zip(BLAST_FIELDS, cells))
            query_coverage = 100 * (int(hit["query_end"]) - int(hit["query_start"]) + 1) / int(hit["query_length"])
            subject_coverage = 100 * (int(hit["subject_end"]) - int(hit["subject_start"]) + 1) / int(hit["subject_length"])
            if min(query_coverage, subject_coverage) < 80:
                continue
            group = query_owner[hit["query_id"]]
            score = (float(hit["bitscore"]), float(hit["identity_percent"]),
                     query_coverage + subject_coverage)
            if group not in best or score > best[group][0]:
                best[group] = (score, {
                    "nearest_group": subject_owner[hit["subject_id"]],
                    "identity_percent": hit["identity_percent"],
                    "query_coverage_percent": f"{query_coverage:.1f}",
                    "subject_coverage_percent": f"{subject_coverage:.1f}",
                })
    return {group: match for group, (_, match) in best.items()}


def main():
    catalogue = read_table(TABLES / "diversity_final_BGC_catalogue.tsv")
    traced = [row for row in catalogue if row["locus_relation_to_literature"] == TRACED]
    unmatched = [row for row in catalogue if row["locus_relation_to_literature"] == UNMATCHED]
    cyclases = SeqIO.to_dict(SeqIO.parse(PROCESSED / "diversity_final_C.faa", "fasta"))
    precursors = SeqIO.to_dict(SeqIO.parse(PROCESSED / "diversity_final_precursors.faa", "fasta"))
    known_genera = {row["organism"].split()[0] for row in traced}
    comparisons = {}
    with tempfile.TemporaryDirectory() as temporary:
        workspace = Path(temporary)
        for label, sequences, key, short in (
            ("C", cyclases, "C_query_id", False),
            ("A", precursors, "precursor_query_ids", True),
        ):
            query, query_owner = records(unmatched, sequences, key)
            subject, subject_owner = records(traced, sequences, key)
            comparisons[label] = nearest(query, subject, query_owner, subject_owner,
                                         short, workspace)
    fields = [
        "organism", "all_source", "family", "NCBI_species_taxid",
        "species_in_resolved_literature", "genus_in_traced", "C_nearest_traced_group",
        "C_identity_percent", "C_query_coverage_percent", "C_traced_coverage_percent",
        "A_nearest_traced_group", "A_identity_percent", "A_query_coverage_percent",
        "A_traced_coverage_percent", "merged_locus_id",
    ]
    destination = TABLES / "diversity_unmatched_characterization.tsv"
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in unmatched:
            group = row["merged_locus_id"]
            c = comparisons["C"].get(group, {})
            a = comparisons["A"].get(group, {})
            writer.writerow({
                "merged_locus_id": group,
                "organism": row["organism"],
                "all_source": row["all_source"],
                "family": row["family"],
                "NCBI_species_taxid": row["NCBI_species_taxid"],
                "species_in_resolved_literature": (
                    str(bool(row["published_same_resolved_species_records"]))
                    if row["NCBI_species_taxid"] else ""
                ),
                "genus_in_traced": str(row["organism"].split()[0] in known_genera),
                **{f"{label}_{column}": value for label, match in (("C", c), ("A", a))
                   for column, value in (
                       ("nearest_traced_group", match.get("nearest_group", "")),
                       ("identity_percent", match.get("identity_percent", "")),
                       ("query_coverage_percent", match.get("query_coverage_percent", "")),
                       ("traced_coverage_percent", match.get("subject_coverage_percent", "")),
                   )},
            })
    print(f"{len(unmatched)} unmatched groups compared with {len(traced)} traced groups: {destination}")


if __name__ == "__main__":
    main()
