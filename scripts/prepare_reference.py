#!/usr/bin/env python3
"""Normalize the author's GenBank records without guessing LrcC identities."""

from __future__ import annotations

import argparse
import csv
import hashlib
import warnings
from pathlib import Path

from Bio import BiopythonParserWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


FEATURE_FIELDS = [
    "sample_id", "feature_id", "file_name", "cds_index", "start_0based",
    "end_0based_exclusive", "strand", "aa_length", "label", "gene",
    "locus_tag", "protein_id", "product", "gene_kind", "sec_met_domain",
    "gene_functions", "note",
    "full_genome_start_0based", "full_genome_end_0based_exclusive",
    "full_genome_strand", "full_genome_locus_tag", "full_genome_protein_id",
    "full_genome_product",
]


def first(qualifiers: dict[str, list[str]], key: str) -> str:
    values = qualifiers.get(key, [])
    return values[0] if values else ""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gbk-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--m2-genome", type=Path, required=True)
    parser.add_argument("--paeninodin-genome", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    metadata = read_tsv(args.metadata)
    by_file = {row["file_name"]: row for row in metadata}
    gbk_files = sorted(args.gbk_dir.glob("*.gbk"))
    if set(by_file) != {path.name for path in gbk_files}:
        missing_metadata = sorted({path.name for path in gbk_files} - set(by_file))
        missing_files = sorted(set(by_file) - {path.name for path in gbk_files})
        raise SystemExit(f"metadata/file mismatch: metadata missing={missing_metadata}; files missing={missing_files}")

    manifest: list[dict[str, object]] = []
    features: list[dict[str, object]] = []
    fasta_records: list[SeqRecord] = []
    m2_seed: SeqRecord | None = None
    m2_genome = SeqIO.read(args.m2_genome, "genbank")

    for path in gbk_files:
        meta = by_file[path.name]
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", BiopythonParserWarning)
            records = list(SeqIO.parse(path, "genbank"))
        if len(records) != 1:
            raise SystemExit(f"expected one GenBank record in {path}, found {len(records)}")
        record = records[0]
        accessions = [value for value in record.annotations.get("accessions", []) if value != "."]
        accession = meta["accession_override"] or (accessions[0] if accessions else "")
        source = next((item for item in record.features if item.type == "source"), None)
        source_organism = first(source.qualifiers, "organism") if source else ""
        cds_features = [item for item in record.features if item.type == "CDS"]
        translation_feature_count = sum(
            1 for item in record.features if first(item.qualifiers, "translation")
        )
        full_match_orientation = ""
        full_match_start: int | str = ""
        full_match_end: int | str = ""
        if meta["sample_id"] == "paenibacillus_m2":
            full_text = str(m2_genome.seq)
            forward_start = full_text.find(str(record.seq))
            reverse_start = full_text.find(str(record.seq.reverse_complement()))
            if forward_start >= 0:
                full_match_orientation = "+"
                full_match_start = forward_start
            elif reverse_start >= 0:
                full_match_orientation = "-"
                full_match_start = reverse_start
            else:
                raise SystemExit("M2 BGC sequence does not occur in CP169648.1")
            full_match_end = int(full_match_start) + len(record)
        manifest.append({
            **meta,
            "embedded_record_id": record.id,
            "embedded_organism": record.annotations.get("organism", ""),
            "source_organism": source_organism,
            "accession": accession,
            "record_bp": len(record),
            "cds_count": len(cds_features),
            "translation_feature_count": translation_feature_count,
            "parser_warning_count": len(caught),
            "sha256": sha256(path),
            "relative_raw_path": str(Path("reference") / args.gbk_dir.name / path.name),
            "full_genome_match_orientation": full_match_orientation,
            "full_genome_start_0based": full_match_start,
            "full_genome_end_0based_exclusive": full_match_end,
        })

        for cds_index, feature in enumerate(cds_features, start=1):
            qualifiers = feature.qualifiers
            translation = first(qualifiers, "translation")
            if not translation:
                continue
            feature_id = f"{meta['sample_id']}__cds{cds_index:02d}"
            full_feature_start: int | str = ""
            full_feature_end: int | str = ""
            full_feature_strand: int | str = ""
            full_feature = None
            if meta["sample_id"] == "paenibacillus_m2":
                if full_match_orientation == "+":
                    full_feature_start = int(full_match_start) + int(feature.location.start)
                    full_feature_end = int(full_match_start) + int(feature.location.end)
                    full_feature_strand = int(feature.location.strand)
                else:
                    full_feature_start = int(full_match_start) + len(record) - int(feature.location.end)
                    full_feature_end = int(full_match_start) + len(record) - int(feature.location.start)
                    full_feature_strand = -int(feature.location.strand)
                overlapping = []
                for candidate in m2_genome.features:
                    if candidate.type != "CDS":
                        continue
                    overlap = min(full_feature_end, int(candidate.location.end)) - max(
                        full_feature_start, int(candidate.location.start)
                    )
                    if overlap <= 0:
                        continue
                    coverage = overlap / min(
                        full_feature_end - full_feature_start,
                        int(candidate.location.end) - int(candidate.location.start),
                    )
                    overlapping.append((coverage, candidate))
                if overlapping and max(overlapping, key=lambda item: item[0])[0] >= 0.8:
                    full_feature = max(overlapping, key=lambda item: item[0])[1]
            full_qualifiers = full_feature.qualifiers if full_feature is not None else {}
            row = {
                "sample_id": meta["sample_id"],
                "feature_id": feature_id,
                "file_name": path.name,
                "cds_index": cds_index,
                "start_0based": int(feature.location.start),
                "end_0based_exclusive": int(feature.location.end),
                "strand": feature.location.strand,
                "aa_length": len(translation),
                "label": first(qualifiers, "label"),
                "gene": first(qualifiers, "gene"),
                "locus_tag": first(qualifiers, "locus_tag"),
                "protein_id": first(qualifiers, "protein_id"),
                "product": first(qualifiers, "product"),
                "gene_kind": first(qualifiers, "gene_kind"),
                "sec_met_domain": " | ".join(qualifiers.get("sec_met_domain", [])),
                "gene_functions": " | ".join(qualifiers.get("gene_functions", [])),
                "note": " | ".join(qualifiers.get("note", [])),
                "full_genome_start_0based": full_feature_start,
                "full_genome_end_0based_exclusive": full_feature_end,
                "full_genome_strand": full_feature_strand,
                "full_genome_locus_tag": first(full_qualifiers, "locus_tag"),
                "full_genome_protein_id": first(full_qualifiers, "protein_id"),
                "full_genome_product": first(full_qualifiers, "product"),
            }
            features.append(row)
            seq_record = SeqRecord(Seq(translation), id=feature_id, description="")
            fasta_records.append(seq_record)
            if meta["sample_id"] == "paenibacillus_m2" and row["label"].lower() in {"larc", "lrcc"}:
                m2_seed = seq_record

    if m2_seed is None:
        raise SystemExit("could not identify the M2 LrcC seed from label LarC/LrcC")

    outgroup_genome = SeqIO.read(args.paeninodin_genome, "genbank")
    outgroup_feature = next(
        (
            item for item in outgroup_genome.features
            if item.type == "CDS" and first(item.qualifiers, "protein_id") == "EHQ60560.1"
        ),
        None,
    )
    if outgroup_feature is None:
        raise SystemExit("paeninodin cyclase EHQ60560.1 not found")
    outgroup = SeqRecord(
        Seq(outgroup_feature.qualifiers["translation"][0]),
        id="paeninodin_outgroup",
        description="EHQ60560.1 Paenibacillus dendritiformis C454",
    )

    manifest_fields = list(metadata[0]) + [
        "embedded_record_id", "embedded_organism", "source_organism", "accession",
        "record_bp", "cds_count", "translation_feature_count", "parser_warning_count",
        "sha256", "relative_raw_path",
        "full_genome_match_orientation", "full_genome_start_0based",
        "full_genome_end_0based_exclusive",
    ]
    write_tsv(args.out_dir / "reference_manifest.tsv", manifest, manifest_fields)
    write_tsv(args.out_dir / "all_cds.tsv", features, FEATURE_FIELDS)
    SeqIO.write(fasta_records, args.out_dir / "all_cds.faa", "fasta")
    SeqIO.write([m2_seed], args.out_dir / "lrcc_m2_seed.faa", "fasta")
    SeqIO.write([outgroup], args.out_dir / "paeninodin_outgroup.faa", "fasta")
    print(f"normalized {len(manifest)} BGC records and {len(features)} translated CDS")


if __name__ == "__main__":
    main()
