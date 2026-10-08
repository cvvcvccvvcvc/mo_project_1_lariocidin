#!/usr/bin/env python3
"""Extract the frozen ribosomal marker panel and audit BGC/assembly identity."""

from __future__ import annotations

import argparse
import gzip
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


MARKERS = (
    "rplB rplC rplD rplE rplF rplK rplL rplM rplN rplO rplP rplQ rplR "
    "rplS rplT rplV rplX rpsB rpsC rpsD rpsE rpsG rpsH rpsI rpsJ rpsK "
    "rpsL rpsM rpsN rpsO rpsP rpsQ rpsR rpsS"
).split()

L_NUMBERS = {
    "rplB": "L2", "rplC": "L3", "rplD": "L4", "rplE": "L5",
    "rplF": "L6", "rplK": "L11", "rplL": "L7/L12", "rplM": "L13",
    "rplN": "L14", "rplO": "L15", "rplP": "L16", "rplQ": "L17",
    "rplR": "L18", "rplS": "L19", "rplT": "L20", "rplV": "L22",
    "rplX": "L24",
}
S_NUMBERS = {
    "rpsB": "S2", "rpsC": "S3", "rpsD": "S4", "rpsE": "S5",
    "rpsG": "S7", "rpsH": "S8", "rpsI": "S9", "rpsJ": "S10",
    "rpsK": "S11", "rpsL": "S12", "rpsM": "S13", "rpsN": "S14",
    "rpsO": "S15", "rpsP": "S16", "rpsQ": "S17", "rpsR": "S18",
    "rpsS": "S19",
}
PRODUCT_TO_MARKER = {
    **{f"50S ribosomal protein {number}": marker for marker, number in L_NUMBERS.items()},
    **{f"30S ribosomal protein {number}": marker for marker, number in S_NUMBERS.items()},
    "50S ribosomal protein L12": "rplL",
}


def read_genome(path: Path) -> tuple[dict[str, list[tuple[str, str]]], set[str], int]:
    hits: dict[str, list[tuple[str, str]]] = defaultdict(list)
    translations: set[str] = set()
    cds_count = 0
    with gzip.open(path, "rt") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type != "CDS" or "translation" not in feature.qualifiers:
                    continue
                cds_count += 1
                sequence = feature.qualifiers["translation"][0]
                translations.add(sequence)
                gene = feature.qualifiers.get("gene", [""])[0]
                product = feature.qualifiers.get("product", [""])[0]
                if gene in MARKERS:
                    hits[gene].append((sequence, "gene_qualifier"))
                elif product in PRODUCT_TO_MARKER:
                    hits[PRODUCT_TO_MARKER[product]].append((sequence, "exact_product_fallback"))
    return hits, translations, cds_count


def read_bgc_sequences(core_dir: Path) -> dict[tuple[str, str], str]:
    sequences: dict[tuple[str, str], str] = {}
    for role in ("A", "B1", "B2", "C", "D"):
        for record in SeqIO.parse(core_dir / f"pseudomonadota_core_{role}.faa", "fasta"):
            sequences[(record.id, role)] = str(record.seq)
    return sequences


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemblies", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--control-review", type=Path, required=True)
    parser.add_argument("--core-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--inventory-out", type=Path, required=True)
    parser.add_argument("--taxa-out", type=Path, required=True)
    parser.add_argument("--bgc-validation-out", type=Path, required=True)
    args = parser.parse_args()

    assemblies = pd.read_csv(args.assemblies, sep="\t", keep_default_na=False)
    selected = pd.read_csv(args.selected, sep="\t", keep_default_na=False)
    control_review = pd.read_csv(args.control_review, sep="\t", keep_default_na=False)
    same_strain_controls = set(control_review.loc[
        control_review["strain_relationship"] == "same_strain_reassembly", "control_sample_id"
    ])
    groups = dict(zip(selected["sample_id"], selected["group"]))
    extracted: dict[str, dict[str, list[tuple[str, str]]]] = {}
    proteomes: dict[str, set[str]] = {}
    observed_cds: dict[str, int] = {}
    root = args.assemblies.resolve().parents[1]

    for row in assemblies.itertuples(index=False):
        hits, translations, cds_count = read_genome(root / row.local_gbff)
        extracted[row.sample_id] = hits
        proteomes[row.sample_id] = translations
        observed_cds[row.sample_id] = cds_count

    eligible_carriers = [
        row.sample_id for row in assemblies.itertuples(index=False)
        if row.genome_role == "carrier" and observed_cds[row.sample_id] > 0
    ]
    omitted_carriers = [
        row.sample_id for row in assemblies.itertuples(index=False)
        if row.genome_role == "carrier" and observed_cds[row.sample_id] == 0
    ]
    retained_markers = [
        marker for marker in MARKERS
        if all(len(extracted[sample].get(marker, [])) == 1 for sample in eligible_carriers)
    ]
    if len(retained_markers) < 10:
        retained_markers = [
            marker for marker in MARKERS
            if sum(len(extracted[sample].get(marker, [])) == 1 for sample in eligible_carriers)
            / len(eligible_carriers) >= 0.8
        ]
    if len(retained_markers) < 10:
        raise ValueError(f"Only {len(retained_markers)} markers pass the frozen rule")

    included_taxa: list[str] = []
    for row in assemblies.itertuples(index=False):
        available = sum(len(extracted[row.sample_id].get(marker, [])) == 1 for marker in retained_markers)
        if row.sample_id in same_strain_controls:
            continue
        if available / len(retained_markers) >= 0.8:
            included_taxa.append(row.sample_id)
        elif row.genome_role == "carrier" and observed_cds[row.sample_id] > 0:
            raise ValueError(f"Annotated carrier fails 80% marker coverage: {row.sample_id}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for marker in retained_markers:
        records = []
        for sample_id in included_taxa:
            values = extracted[sample_id].get(marker, [])
            if len(values) == 1:
                records.append(SeqRecord(Seq(values[0][0]), id=sample_id, description=""))
        SeqIO.write(records, args.output_dir / f"{marker}.faa", "fasta")

    inventory_rows = []
    for row in assemblies.itertuples(index=False):
        for marker in MARKERS:
            values = extracted[row.sample_id].get(marker, [])
            inventory_rows.append({
                "sample_id": row.sample_id,
                "genome_role": row.genome_role,
                "marker": marker,
                "copy_count": len(values),
                "extraction_basis": values[0][1] if len(values) == 1 else "",
                "marker_retained": marker in retained_markers,
                "taxon_included": row.sample_id in included_taxa,
            })
    pd.DataFrame(inventory_rows).to_csv(args.inventory_out, sep="\t", index=False)

    taxa_rows = []
    for row in assemblies.itertuples(index=False):
        carrier_id = row.sample_id.removeprefix("control__")
        taxa_rows.append({
            "sample_id": row.sample_id,
            "genome_role": row.genome_role,
            "carrier_context": carrier_id,
            "taxonomic_group": groups.get(carrier_id, ""),
            "organism": row.organism,
            "assembly_accession": row.assembly_accession,
            "markers_present": sum(
                len(extracted[row.sample_id].get(marker, [])) == 1
                for marker in retained_markers
            ),
            "markers_total": len(retained_markers),
            "host_tree_status": (
                "included" if row.sample_id in included_taxa else
                "same_strain_reassembly" if row.sample_id in same_strain_controls else
                "unannotated_or_low_coverage"
            ),
        })
    pd.DataFrame(taxa_rows).to_csv(args.taxa_out, sep="\t", index=False)

    bgc = read_bgc_sequences(args.core_dir)
    validation_rows = []
    for row in assemblies.loc[assemblies["genome_role"] == "carrier"].itertuples(index=False):
        for role in ("A", "B1", "B2", "C", "D"):
            sequence = bgc[(row.sample_id, role)]
            validation_rows.append({
                "sample_id": row.sample_id,
                "assembly_accession": row.assembly_accession,
                "role": role,
                "exact_translation_matches": int(sequence in proteomes[row.sample_id]),
                "assembly_has_translated_cds": bool(proteomes[row.sample_id]),
            })
    pd.DataFrame(validation_rows).to_csv(args.bgc_validation_out, sep="\t", index=False)

    print(
        f"retained {len(retained_markers)} markers; included {len(included_taxa)} taxa; "
        f"omitted unannotated carriers={','.join(omitted_carriers) or 'none'}"
    )


if __name__ == "__main__":
    main()
