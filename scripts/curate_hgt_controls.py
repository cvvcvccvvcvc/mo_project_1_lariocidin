#!/usr/bin/env python3
"""Apply an explicit local-context rule to profile hits in same-species controls."""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def processing_homology(proteins, core_dir, output):
    """Check B1/B2/D sequence homology; protein size alone does not identify a role."""
    refs = []
    for role in ("B1", "B2", "D"):
        for record in SeqIO.parse(core_dir / f"pseudomonadota_core_{role}.faa", "fasta"):
            refs.append(SeqRecord(record.seq, id=f"{role}|{record.id}", description=""))
    columns = ["feature_id", "reference_id", "identity_pct", "query_coverage_pct", "evalue", "bitscore"]
    with tempfile.TemporaryDirectory() as directory:
        directory = Path(directory)
        query, reference, hits = (directory / name for name in ("query.faa", "reference.faa", "hits.tsv"))
        SeqIO.write([SeqRecord(Seq(seq), id=key, description="") for key, seq in proteins.items()], query, "fasta")
        SeqIO.write(refs, reference, "fasta")
        if proteins:
            subprocess.run([
                str(Path(__file__).resolve().parents[1] / ".env/bin/blastp"),
                "-query", str(query), "-subject", str(reference), "-seg", "yes",
                "-evalue", "1e-5", "-max_target_seqs", str(len(refs)),
                "-outfmt", "6 qseqid sseqid pident qcovs evalue bitscore", "-out", str(hits),
            ], check=True)
        data = pd.read_csv(hits, sep="\t", names=columns) if proteins and hits.stat().st_size else pd.DataFrame(columns=columns)
    data["role"] = data.reference_id.str.split("|").str[0]
    data = data.sort_values("bitscore", ascending=False).drop_duplicates(["feature_id", "role"])
    data["passes_homology_rule"] = (data.identity_pct >= 30) & (data.query_coverage_pct >= 70)
    data.to_csv(output, sep="\t", index=False)
    return data.loc[data.passes_homology_rule]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemblies", type=Path, required=True)
    parser.add_argument("--review", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--proteins", type=Path, required=True)
    parser.add_argument("--known-c", type=Path, required=True)
    parser.add_argument("--core-dir", type=Path, required=True)
    parser.add_argument("--homology-out", type=Path, required=True)
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--status-out", type=Path, required=True)
    parser.add_argument("--contexts-out", type=Path, required=True)
    args = parser.parse_args()

    assemblies = pd.read_csv(args.assemblies, sep="\t", keep_default_na=False)
    assemblies = assemblies.loc[assemblies["genome_role"] == "same_species_control"]
    review = pd.read_csv(args.review, sep="\t", keep_default_na=False)
    features = pd.read_csv(args.features, sep="\t", keep_default_na=False)
    hits = pd.read_csv(args.hits, sep="\t", keep_default_na=False)
    proteins = {record.id: str(record.seq) for record in SeqIO.parse(args.proteins, "fasta")}
    known_c = {record.id: str(record.seq) for record in SeqIO.parse(args.known_c, "fasta")}
    precursor_sequences = {str(record.seq) for record in SeqIO.parse(args.core_dir / "pseudomonadota_core_A.faa", "fasta")}
    high_hits = hits.loc[(hits.profile == "canonical") & (hits.sequence_score >= 250) & (hits.model_coverage >= 0.70)]
    neighbor_ids = set()
    for hit in high_hits.itertuples(index=False):
        neighbor_ids.update(features.loc[
            (features.sample_id == hit.sample_id) & (features.contig_accession == hit.contig_accession)
            & (features.strand == hit.strand) & features.cds_index.between(hit.cds_index - 5, hit.cds_index + 5), "feature_id"
        ])
    homologs = processing_homology({key: proteins[key] for key in sorted(neighbor_ids)}, args.core_dir, args.homology_out)

    context_rows = []
    status_rows = []
    for assembly in assemblies.itertuples(index=False):
        sample_hits = hits.loc[
            (hits["sample_id"] == assembly.sample_id)
            & (hits["profile"] == "canonical")
            & (hits["sequence_score"] >= 250)
            & (hits["model_coverage"] >= 0.70)
        ]
        supported = []
        for hit in sample_hits.itertuples(index=False):
            neighborhood = features.loc[
                (features["sample_id"] == assembly.sample_id)
                & (features["contig_accession"] == hit.contig_accession)
                & features["cds_index"].between(hit.cds_index - 5, hit.cds_index + 5)
                & (features["strand"] == hit.strand)
            ]
            other = neighborhood.loc[neighborhood["feature_id"] != hit.target_id]
            precursor_count = int(other["aa_length"].between(20, 60).sum())
            processing_count = int(other["aa_length"].between(50, 200).sum())
            transporter_count = int(other["aa_length"].between(450, 700).sum())
            role_ids = {role: set(homologs.loc[
                (homologs.role == role) & homologs.feature_id.isin(other.feature_id), "feature_id"
            ]) for role in ("B1", "B2", "D")}
            precursor_ids = {key for key in other.feature_id if proteins[key] in precursor_sequences}
            # Distinct genes must support the roles; one short CDS cannot fill A and B1.
            context_supported = any(
                len({a, b1, b2, d}) == 4
                for a in precursor_ids for b1 in role_ids["B1"] for b2 in role_ids["B2"] for d in role_ids["D"]
            )
            exact_matches = [
                sample for sample, sequence in known_c.items()
                if sequence == proteins[hit.target_id]
            ]
            context_rows.append({
                "sample_id": assembly.sample_id,
                "assembly_accession": assembly.assembly_accession,
                "target_id": hit.target_id,
                "contig_accession": hit.contig_accession,
                "cds_index": hit.cds_index,
                "sequence_score": hit.sequence_score,
                "model_coverage": hit.model_coverage,
                "target_coverage": hit.target_coverage,
                "precursor_sized_same_strand_cds": precursor_count,
                "processing_sized_same_strand_cds": processing_count,
                "transporter_sized_same_strand_cds": transporter_count,
                "exact_known_c_matches": "|".join(exact_matches),
                "exact_known_precursor_ids": "|".join(sorted(precursor_ids)),
                **{f"{role}_homolog_ids": "|".join(sorted(ids)) for role, ids in role_ids.items()},
                "context_supported": context_supported,
            })
            supported.append(context_supported)

        relationship = review.loc[review["control_sample_id"] == assembly.sample_id]
        if len(relationship) != 1:
            raise ValueError(f"Expected one control review for {assembly.sample_id}")
        relationship = relationship.iloc[0]
        status_rows.append({
            "sample_id": assembly.sample_id,
            "assembly_accession": assembly.assembly_accession,
            "assembly_level": assembly.assembly_level,
            "carrier_source_label": relationship["carrier_source_label"],
            "control_source_label": relationship["control_source_label"],
            "strain_relationship": relationship["strain_relationship"],
            "independent_strain_label": relationship["strain_relationship"] != "same_strain_reassembly",
            "canonical_high_hits": len(sample_hits),
            "context_supported_bgc": any(supported),
            "status": "canonical_bgc_context_supported" if any(supported) else "target_bgc_not_detected",
            "absence_caveat": (
                "not an absence proof; depends on assembly and annotation"
                if not any(supported) else ""
            ),
        })

    statuses = pd.DataFrame(status_rows)
    contexts = pd.DataFrame(context_rows)
    statuses.to_csv(args.status_out, sep="\t", index=False)
    contexts.to_csv(args.contexts_out, sep="\t", index=False)
    print(statuses.groupby(["independent_strain_label", "status"]).size().to_string())


if __name__ == "__main__":
    main()
