#!/usr/bin/env python3
"""Normalize the published triculamin-like BGC archive without editing raw data."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import warnings
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


GROUP_LABELS = {
    "canonical": "triculamin_like_canonical",
    "noncanonical": "triculamin_like_noncanonical",
    "Barrett et al review": "other_lasso_comparator",
}
GENBANK_SUFFIXES = (".gb", ".gbk", ".gbff")
KNOWN_REPAIR_MEMBER = "noncanonical/WP_199826149.1.gbk"
BROKEN_COMMENT = (
    "            project_name        :: SeaBioTech actinobacteria experimental_factor\r\n"
    "            :: not applicable"
)
REPAIRED_COMMENT = (
    "            project_name        :: SeaBioTech actinobacteria\r\n"
    "            experimental_factor :: not applicable"
)


def sha256(data: bytes | str) -> str:
    if isinstance(data, str):
        data = data.encode()
    return hashlib.sha256(data).hexdigest()


def joined(values: list[str] | None) -> str:
    return ";".join(values or [])


def slug(value: str) -> str:
    clean = re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")
    return clean[:48] or "record"


def is_author_cyclase(feature) -> bool:
    return any(value.startswith("C_") for value in feature.qualifiers.get("product", []))


def is_author_precursor(feature) -> bool:
    return any(
        value == "A" or value.startswith(("A_", "A1_"))
        for value in feature.qualifiers.get("label", [])
    )


def parse_member(member: str, raw: bytes):
    text = raw.decode("utf-8")
    repair = ""
    if member == KNOWN_REPAIR_MEMBER:
        if BROKEN_COMMENT not in text:
            raise ValueError(f"known repair target changed: {member}")
        text = text.replace(BROKEN_COMMENT, REPAIRED_COMMENT, 1)
        repair = "split_malformed_MIGS_project_name_and_experimental_factor"
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        records = list(SeqIO.parse(io.StringIO(text), "genbank"))
    if len(records) != 1:
        raise ValueError(f"expected one GenBank record in {member}, found {len(records)}")
    return records[0], repair, len(caught)


def write_fasta(records: list[SeqRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, path, "fasta")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--lar29-fasta", type=Path, required=True)
    parser.add_argument("--anchors", type=Path, required=True)
    parser.add_argument("--processed-dir", type=Path, required=True)
    parser.add_argument("--table-dir", type=Path, required=True)
    args = parser.parse_args()

    args.processed_dir.mkdir(parents=True, exist_ok=True)
    args.table_dir.mkdir(parents=True, exist_ok=True)

    lar29 = defaultdict(list)
    for record in SeqIO.parse(args.lar29_fasta, "fasta"):
        lar29[str(record.seq)].append(record.id)
    anchors = pd.read_csv(args.anchors, sep="\t", keep_default_na=False).set_index("archive_member")

    manifests: list[dict] = []
    cds_rows: list[dict] = []
    cyclase_rows: list[dict] = []
    precursor_rows: list[dict] = []
    all_cds_fasta: list[SeqRecord] = []
    unresolved_fasta: list[SeqRecord] = []
    cyclase_sequences: dict[str, str] = {}
    precursor_sequences: dict[str, str] = {}

    with zipfile.ZipFile(args.archive) as archive:
        members = sorted(
            member for member in archive.namelist()
            if member.lower().endswith(GENBANK_SUFFIXES) and not member.startswith("__MACOSX/")
        )
        for member in members:
            source_group = member.split("/", 1)[0]
            if source_group not in GROUP_LABELS:
                raise ValueError(f"unknown source group for {member}")
            raw = archive.read(member)
            record, repair, warning_count = parse_member(member, raw)
            bgc_id = f"lit_{slug(Path(member).stem)}_{sha256(member)[:10]}"
            source = next((feature for feature in record.features if feature.type == "source"), None)
            source_q = source.qualifiers if source else {}
            cds = [feature for feature in record.features if feature.type == "CDS"]
            author_cyclases = [feature for feature in cds if is_author_cyclase(feature)]
            author_precursors = [feature for feature in cds if is_author_precursor(feature)]
            accession = joined(record.annotations.get("accessions", [])) or record.id
            anchor = anchors.loc[member].to_dict() if member in anchors.index else {}

            manifests.append({
                "bgc_id": bgc_id,
                "archive_member": member,
                "source_group": source_group,
                "family_label": GROUP_LABELS[source_group],
                "record_id": record.id,
                "accession": accession,
                "organism": joined(source_q.get("organism")) or record.annotations.get("organism", ""),
                "strain": joined(source_q.get("strain")),
                "isolate": joined(source_q.get("isolate")),
                "host": joined(source_q.get("host")),
                "isolation_source": joined(source_q.get("isolation_source")),
                "geo_loc_name": joined(source_q.get("geo_loc_name")),
                "collection_date": joined(source_q.get("collection_date")),
                "record_bp": len(record.seq),
                "cds_count": len(cds),
                "author_cyclase_count": len(author_cyclases),
                "author_precursor_count": len(author_precursors),
                "product_name": anchor.get("product_name", ""),
                "evidence_status": anchor.get("evidence_status", "predicted_in_published_collection"),
                "evidence_source": anchor.get("evidence_source", "Merrild_et_al_2025_archive"),
                "parse_repair": repair,
                "parser_warning_count": warning_count,
                "member_sha256": sha256(raw),
                "record_sequence_sha256": sha256(str(record.seq).upper()),
            })

            explicit_cyclase_feature_ids = set()
            for index, feature in enumerate(cds, start=1):
                translation = joined(feature.qualifiers.get("translation")).replace(";", "").replace(" ", "").replace("\n", "").rstrip("*")
                feature_id = f"{bgc_id}__cds{index:03d}"
                start = int(feature.location.start)
                end = int(feature.location.end)
                row = {
                    "bgc_id": bgc_id,
                    "feature_id": feature_id,
                    "archive_member": member,
                    "source_group": source_group,
                    "family_label": GROUP_LABELS[source_group],
                    "cds_index": index,
                    "start_0based": start,
                    "end_0based_exclusive": end,
                    "strand": feature.location.strand or "",
                    "aa_length": len(translation),
                    "label": joined(feature.qualifiers.get("label")),
                    "gene": joined(feature.qualifiers.get("gene")),
                    "locus_tag": joined(feature.qualifiers.get("locus_tag")),
                    "protein_id": joined(feature.qualifiers.get("protein_id")),
                    "product": joined(feature.qualifiers.get("product")),
                    "author_role": "cyclase" if is_author_cyclase(feature) else ("precursor" if is_author_precursor(feature) else "context_other"),
                    "sequence_sha256": sha256(translation) if translation else "",
                }
                cds_rows.append(row)
                if translation:
                    all_cds_fasta.append(SeqRecord(Seq(translation), id=feature_id, description=""))
                if is_author_cyclase(feature):
                    explicit_cyclase_feature_ids.add(feature_id)
                    cyclase_id = f"cyc_{sha256(translation)[:16]}"
                    cyclase_sequences.setdefault(cyclase_id, translation)
                    cyclase_rows.append({
                        **row,
                        "cyclase_id": cyclase_id,
                        "exact_lar29_sample_ids": ";".join(sorted(lar29.get(translation, []))),
                        "exact_lar29": bool(lar29.get(translation)),
                        "anchor_product_name": anchor.get("product_name", ""),
                        "anchor_evidence_status": anchor.get("evidence_status", ""),
                        "anchor_evidence_source": anchor.get("evidence_source", ""),
                    })
                if is_author_precursor(feature):
                    precursor_id = f"pre_{sha256(translation)[:16]}"
                    precursor_sequences.setdefault(precursor_id, translation)
                    precursor_rows.append({**row, "precursor_id": precursor_id})

            if not explicit_cyclase_feature_ids:
                for row in cds_rows[-len(cds):]:
                    if row["aa_length"]:
                        sequence = next(record.seq for record in all_cds_fasta if record.id == row["feature_id"])
                        unresolved_fasta.append(SeqRecord(sequence, id=row["feature_id"], description=""))

    manifest = pd.DataFrame(manifests)
    cds_table = pd.DataFrame(cds_rows)
    cyclases = pd.DataFrame(cyclase_rows)
    precursors = pd.DataFrame(precursor_rows)
    manifest.to_csv(args.processed_dir / "literature_bgc_manifest.tsv", sep="\t", index=False)
    cds_table.to_csv(args.processed_dir / "literature_all_cds.tsv", sep="\t", index=False)
    cyclases.to_csv(args.processed_dir / "literature_cyclase_occurrences.tsv", sep="\t", index=False)
    precursors.to_csv(args.processed_dir / "literature_precursor_occurrences.tsv", sep="\t", index=False)
    unique_cyclases = (
        cyclases.groupby("cyclase_id", sort=True)
        .agg(
            aa_length=("aa_length", "first"),
            occurrence_count=("bgc_id", "size"),
            source_groups=("source_group", lambda values: ";".join(sorted(set(values)))),
            family_labels=("family_label", lambda values: ";".join(sorted(set(values)))),
            archive_members=("archive_member", lambda values: ";".join(sorted(set(values)))),
            protein_ids=("protein_id", lambda values: ";".join(sorted({value for value in values if value}))),
            exact_lar29_sample_ids=("exact_lar29_sample_ids", lambda values: ";".join(sorted({item for value in values for item in value.split(";") if item}))),
            anchor_products=("anchor_product_name", lambda values: ";".join(sorted({value for value in values if value}))),
        )
        .reset_index()
    )
    unique_cyclases.to_csv(args.processed_dir / "literature_unique_cyclases.tsv", sep="\t", index=False)
    write_fasta(all_cds_fasta, args.processed_dir / "literature_all_cds.faa")
    write_fasta(unresolved_fasta, args.processed_dir / "literature_unresolved_bgc_cds.faa")
    write_fasta(
        [SeqRecord(Seq(sequence), id=cyclase_id, description="") for cyclase_id, sequence in sorted(cyclase_sequences.items())],
        args.processed_dir / "literature_unique_cyclases.faa",
    )
    for source_group, file_stem in [
        ("canonical", "literature_canonical_cyclases"),
        ("noncanonical", "literature_noncanonical_cyclases"),
        ("Barrett et al review", "literature_other_lasso_cyclases"),
    ]:
        ids = set(cyclases.loc[cyclases["source_group"] == source_group, "cyclase_id"])
        write_fasta(
            [SeqRecord(Seq(cyclase_sequences[cyclase_id]), id=cyclase_id, description="") for cyclase_id in sorted(ids)],
            args.processed_dir / f"{file_stem}.faa",
        )
    write_fasta(
        [SeqRecord(Seq(sequence), id=precursor_id, description="") for precursor_id, sequence in sorted(precursor_sequences.items())],
        args.processed_dir / "literature_unique_precursors.faa",
    )

    group_counts = manifest.groupby("source_group").agg(
        bgc_records=("bgc_id", "size"),
        records_with_one_explicit_cyclase=("author_cyclase_count", lambda values: int((values == 1).sum())),
        records_without_explicit_cyclase=("author_cyclase_count", lambda values: int((values == 0).sum())),
    ).reset_index()
    group_counts.to_csv(args.table_dir / "literature_family_counts.tsv", sep="\t", index=False)
    summary = {
        "archive_members_parsed": len(manifest),
        "parse_repairs": int((manifest["parse_repair"] != "").sum()),
        "translated_cds": int((cds_table["aa_length"] > 0).sum()),
        "explicit_cyclase_occurrences": len(cyclases),
        "unique_explicit_cyclase_sequences": len(cyclase_sequences),
        "unique_cyclase_sequences_by_source_group": {
            source_group: int(cyclases.loc[cyclases["source_group"] == source_group, "cyclase_id"].nunique())
            for source_group in GROUP_LABELS
        },
        "explicit_precursor_occurrences": len(precursors),
        "unique_explicit_precursor_sequences": len(precursor_sequences),
        "records_without_explicit_cyclase": int((manifest["author_cyclase_count"] == 0).sum()),
        "lar29_sequences_with_exact_cyclase_overlap": len({sample for value in cyclases["exact_lar29_sample_ids"] for sample in value.split(";") if sample}),
        "source_group_records": {key: int(value) for key, value in Counter(manifest["source_group"]).items()},
    }
    (args.table_dir / "literature_reference_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    assert len(manifest) == 393
    assert len(cyclases) == 389
    assert (manifest["author_cyclase_count"] <= 1).all()
    assert summary["parse_repairs"] == 1
    assert summary["lar29_sequences_with_exact_cyclase_overlap"] == 25
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
