#!/usr/bin/env python3
"""Set stable host aliases and experimentally/published annotated references."""
from collections import defaultdict
import csv
import json
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[2]
SNAP = ROOT / "inputs/H2_2026-10-13"
H1 = SNAP / "H1"
OUT = ROOT / "results/H2"


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    quality = read(H1 / "H1_genome_quality.tsv")
    assemblies = {r["genome_id"]: r for r in read(H1 / "H1_all_assembly_inputs.tsv")}
    grouped = defaultdict(list)
    for row in quality:
        grouped[row["sequence_content_sha256"]].append(row)
    representatives, aliases = [], []
    for number, (content_hash, rows) in enumerate(sorted(grouped.items()), 1):
        rows.sort(key=lambda r: (-int(r["translated_cds_count"]), r["full_assembly_input"] != "True",
                                 not r["sequence_provider_assembly"].startswith("GCF_"), r["genome_id"]))
        best = rows[0]
        host_id = f"H{number:04}"
        representatives.append({"host_id": host_id, **best, "path": assemblies[best["genome_id"]]["path"],
                                "represented_inputs": ";".join(r["genome_id"] for r in rows)})
        for row in rows:
            aliases.append({"genome_id": row["genome_id"], "host_id": host_id,
                            "representative_genome_id": best["genome_id"], "genome_content_sha256": content_hash,
                            "biosample": row["biosample"], "organism": row["organism"], "strain": row["strain"],
                            "reason": "identical_genomic_sequence_computational_alias; biological_sample_identity_not_assumed"})
    write(OUT / "host_inputs.tsv", representatives)
    write(OUT / "host_aliases.tsv", aliases)
    genes = read(H1 / "H1_publication_genes.tsv")
    refs = []
    for gene in genes:
        if gene["author_role"] in ["B1", "B2"]:
            refs.append({"gene_id": gene["gene_id"], "family": "canonical_" + gene["author_role"],
                         "basis": "published_Bgel_Cpal_Ptar_annotation"})
    for feature, family in [(2, "canonical_B1"), (3, "canonical_B2"), (4, "canonical_C"),
                            (5, "canonical_ABC_D1"), (6, "canonical_ABC_D2"), (7, "canonical_GNAT_E")]:
        refs.append({"gene_id": f"Nature_2025_91bc539987e9.f{feature}", "family": family,
                     "basis": "Nature_M2_reference; ABC_subfamily_not_automatically_orthologous"})
    for feature, family in [(1, "noncanonical_GNAT_neighbour"), (4, "noncanonical_C"), (5, "noncanonical_ABC_D")]:
        refs.append({"gene_id": f"Merrild_2025_879faa010cad.f{feature}", "family": family,
                     "basis": "Merrild_triculamin_reference; upstream_GNAT_is_context_not_confirmed_core_function"})
    write(OUT / "role_references.tsv", refs)
    wanted = {r["gene_id"] for r in refs}
    count = SeqIO.write((r for r in SeqIO.parse(H1 / "proteins.faa", "fasta") if r.id in wanted), OUT / "role_references.faa", "fasta")
    if count != len(wanted):
        raise ValueError("Reference sequence missing")
    print(json.dumps({"genome_inputs": len(quality), "independent_sequence_representatives": len(representatives), "role_references": count}))


if __name__ == "__main__":
    main()
