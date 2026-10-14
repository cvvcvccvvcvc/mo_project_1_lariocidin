#!/usr/bin/env python3
"""Audit annotated roles in selected reference genomes against published CDSs."""
import csv
import json
from collections import defaultdict

from Bio import SeqIO

from source_cache import ROOT


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    processed, tables = ROOT / "data/processed", ROOT / "results/tables"
    features = {r["feature_id"]: r for r in rows(processed / "reference/all_cds.tsv")}
    literature_features = rows(processed / "literature/literature_all_cds.tsv")
    literature_sequences = {r.id: str(r.seq) for r in SeqIO.parse(
        processed / "literature/literature_all_cds.faa", "fasta")}
    exact = defaultdict(list)
    for feature in literature_features:
        if feature["feature_id"] in literature_sequences:
            exact[literature_sequences[feature["feature_id"]]].append(feature)
    sequences = {r.id: str(r.seq) for r in SeqIO.parse(processed / "reference/all_cds.faa", "fasta")}
    domains = {r["protein_id"]: r for r in rows(tables / "diversity_domain_evidence.tsv")}

    evidence = []
    for role in rows(tables / "pseudomonadota_core_roles.tsv"):
        feature = features[role["feature_id"]]
        matches = exact[sequences[role["feature_id"]]]
        domain = domains[role["feature_id"]]
        evidence.append({
            "locus": role["sample_id"], "protein_id": role["feature_id"],
            "assigned_role": role["role"], "original_assignment": role["role_assignment_basis"],
            "annotation_product": feature["product"], "annotation_label": feature["label"],
            "existing_domain_annotations": feature["sec_met_domain"],
            "exact_published_CDS_matches": ";".join(x["feature_id"] for x in matches),
            "exact_match_labels": ";".join(sorted({x["label"] for x in matches if x["label"]})),
            "exact_match_products": ";".join(sorted({x["product"] for x in matches if x["product"]})),
            "function_status": "inferred_role_not_per_locus_experimental_validation",
            "measured_Pfam_panel_domains": domain["domains"],
            "measured_Pfam_panel_accessions": domain["domain_accessions"],
            "domain_panel_status": domain["status"],
        })
    with (tables / "diversity_role_evidence.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(evidence[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(evidence)

    summary = {
        "historical_cyclase_database": {
            group: sum(r["source_groups"] == group for r in rows(processed / "literature/literature_unique_cyclases.tsv"))
            for group in sorted({r["source_groups"] for r in rows(processed / "literature/literature_unique_cyclases.tsv")})
        },
        "historical_reviewed_AsnB": 5,
        "roles_audited": len(evidence),
        "published_records_compared": len(literature_features),
    }
    (tables / "diversity_baseline_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
