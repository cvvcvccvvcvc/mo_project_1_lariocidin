#!/usr/bin/env python3
"""Verify focal proteins against full NCBI records and describe 20 kb flanks."""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO


MOBILE_PATTERN = re.compile(
    r"transpos|integrase|recombinase|insertion sequence|mobile element|phage|prophage",
    re.IGNORECASE,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def joined(values) -> str:
    return ";".join(str(value).replace("\n", " ") for value in values or [])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--core-roles", type=Path, required=True)
    parser.add_argument("--all-cds", type=Path, required=True)
    parser.add_argument("--manifest-out", type=Path, required=True)
    parser.add_argument("--features-out", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=20_000)
    args = parser.parse_args()

    config = pd.read_csv(args.config, sep="\t", keep_default_na=False)
    selected = pd.read_csv(args.selected, sep="\t", keep_default_na=False).set_index("sample_id")
    core_roles = pd.read_csv(args.core_roles, sep="\t", keep_default_na=False)
    local_sequences = {record.id: str(record.seq) for record in SeqIO.parse(args.all_cds, "fasta")}
    roles_by_sample_sequence = {
        sample_id: {
            local_sequences[row.feature_id]: row.role
            for row in group.itertuples(index=False)
        }
        for sample_id, group in core_roles.groupby("sample_id")
    }
    manifest_rows = []
    feature_rows = []
    for request in config.itertuples(index=False):
        path = args.raw_dir / f"{request.sample_id}.gb"
        records = list(SeqIO.parse(path, "genbank"))
        if len(records) != 1:
            raise ValueError(f"Expected one NCBI record in {path}, found {len(records)}")
        record = records[0]
        focal = []
        for feature in record.features:
            if feature.type != "CDS":
                continue
            protein_ids = feature.qualifiers.get("protein_id", [])
            if request.focal_protein_id in protein_ids:
                focal.append(feature)
        if len(focal) != 1:
            raise ValueError(
                f"Expected one {request.focal_protein_id} in {record.id}, found {len(focal)}"
            )
        focal_feature = focal[0]
        focal_translation = joined(focal_feature.qualifiers.get("translation")).replace(" ", "")
        expected_length = int(selected.loc[request.sample_id, "aa_length"])
        exact_translation_length = len(focal_translation) == expected_length
        expected_c_feature = str(selected.loc[request.sample_id, "feature_id"])
        exact_translation = focal_translation == local_sequences[expected_c_feature]
        focal_start = int(focal_feature.location.start)
        focal_end = int(focal_feature.location.end)
        window_start = max(0, focal_start - args.flank)
        window_end = min(len(record.seq), focal_end + args.flank)
        source = next((feature for feature in record.features if feature.type == "source"), None)
        source_q = source.qualifiers if source else {}
        mobile_count = 0
        cds_count = 0
        record_feature_rows = []
        for feature in record.features:
            if feature.type != "CDS":
                continue
            start, end = int(feature.location.start), int(feature.location.end)
            if end < window_start or start > window_end:
                continue
            cds_count += 1
            product = joined(feature.qualifiers.get("product"))
            note = joined(feature.qualifiers.get("note"))
            translation = joined(feature.qualifiers.get("translation")).replace(" ", "")
            core_role = roles_by_sample_sequence[request.sample_id].get(translation, "")
            mobile_keyword = bool(MOBILE_PATTERN.search(f"{product} {note}"))
            mobile_count += int(mobile_keyword)
            record_feature_rows.append(
                {
                    "sample_id": request.sample_id,
                    "record_accession": record.id,
                    "focal_protein_id": request.focal_protein_id,
                    "start_0based": start,
                    "end_0based_exclusive": end,
                    "strand": feature.location.strand or "",
                    "distance_to_focal_bp": (
                        focal_start - end if end < focal_start else start - focal_end if start > focal_end else 0
                    ),
                    "locus_tag": joined(feature.qualifiers.get("locus_tag")),
                    "protein_id": joined(feature.qualifiers.get("protein_id")),
                    "product": product,
                    "exact_core_role": core_role,
                    "mobile_keyword": mobile_keyword,
                }
            )
        matched_core = [row for row in record_feature_rows if row["exact_core_role"]]
        matched_roles = {row["exact_core_role"] for row in matched_core}
        if not {"B1", "B2", "C", "D"} <= matched_roles:
            raise ValueError(f"Incomplete exact core mapping for {request.sample_id}: {matched_roles}")
        core_start = min(row["start_0based"] for row in matched_core)
        core_end = max(row["end_0based_exclusive"] for row in matched_core)
        for row in record_feature_rows:
            start, end = row["start_0based"], row["end_0based_exclusive"]
            row["distance_to_core_bgc_bp"] = (
                core_start - end if end < core_start else start - core_end if start > core_end else 0
            )
        feature_rows.extend(record_feature_rows)
        manifest_rows.append(
            {
                "sample_id": request.sample_id,
                "requested_nucleotide_accession": request.requested_nucleotide_accession,
                "resolved_record_accession": record.id,
                "record_length_bp": len(record.seq),
                "topology": record.annotations.get("topology", ""),
                "organism": record.annotations.get("organism", ""),
                "isolation_source": joined(source_q.get("isolation_source")),
                "geo_loc_name": joined(source_q.get("geo_loc_name")),
                "focal_protein_id": request.focal_protein_id,
                "focal_start_0based": focal_start,
                "focal_end_0based_exclusive": focal_end,
                "focal_strand": focal_feature.location.strand or "",
                "focal_translation_length": len(focal_translation),
                "reference_translation_length": expected_length,
                "translation_length_matches_reference": exact_translation_length,
                "translation_matches_reference": exact_translation,
                "exact_core_roles_found": ";".join(sorted(matched_roles)),
                "exact_core_role_count": len(matched_roles),
                "core_start_0based": core_start,
                "core_end_0based_exclusive": core_end,
                "left_context_bp_available": focal_start,
                "right_context_bp_available": len(record.seq) - focal_end,
                "window_start_0based": window_start,
                "window_end_0based_exclusive": window_end,
                "window_cds_count": cds_count,
                "mobile_keyword_cds_count": mobile_count,
                "raw_file_bytes": path.stat().st_size,
                "raw_file_sha256": sha256(path),
                "source_database": request.source_database,
                "source_url": request.source_url,
            }
        )
    manifest = pd.DataFrame(manifest_rows)
    features = pd.DataFrame(feature_rows).sort_values(["sample_id", "start_0based"])
    args.manifest_out.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.manifest_out, sep="\t", index=False)
    features.to_csv(args.features_out, sep="\t", index=False)
    if not manifest["translation_matches_reference"].all():
        raise ValueError("At least one focal translation differs from the selected reference")
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
