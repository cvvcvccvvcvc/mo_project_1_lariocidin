#!/usr/bin/env python3
"""Evidence-based checkpoint: explicit supported cores, provisional background.

Post-screen triage is not a trained classifier or a new sensitivity estimate.
AsnB/other-lasso affinity labels must not be turned into proven family exclusions.
"""
import hashlib
import json
import zipfile
from collections import defaultdict

import pandas as pd
from Bio import SeqIO

from maize_genomes import read_genome
from prepare_literature_reference import parse_member
from profile_search import ROOT
from summarize_maize_candidates import BLAST_COLUMNS


def seqhash(sequence):
    return hashlib.sha256(str(sequence).encode()).hexdigest()


def main():
    tables, processed = ROOT / "results/tables", ROOT / "data/processed"
    triage = pd.read_csv(tables / "maize_candidate_triage.tsv", sep="\t")
    asnb = pd.read_csv(tables / "maize_candidates_vs_asnb.tsv", sep="\t", names=BLAST_COLUMNS)
    best = asnb.sort_values("bitscore", ascending=False).drop_duplicates("query_id")
    best = best.rename(columns={c: "asnb_"+c for c in best.columns})
    triage = triage.merge(best, left_on="target_id", right_on="asnb_query_id", how="left", validate="one_to_one")
    triage["asnb_minus_lasso_bits"] = triage.asnb_bitscore - triage.bitscore
    triage["curation_status"] = "unresolved"
    strong_roles = pd.read_csv(ROOT / "inputs/selections/baseline/maize_curated_core_roles.tsv", sep="\t")
    supported_ids = set(strong_roles.candidate_id)
    for index, row in triage.iterrows():
        if row.target_id in supported_ids:
            status = "context_supported_canonical_culture_not_MAG"
            reason = "Near-identical published C plus contiguous B1/B2, precursor-like CDS and ABC; primary BioSample says cell culture."
        elif row.asnb_minus_lasso_bits > 0:
            status = "AsnB_favored_target_unconfirmed"
            reason = "Reviewed AsnB has higher pairwise BLAST score than every reference cyclase; affinity evidence, not experimentally verified function."
        elif row.source_groups == "Barrett et al review":
            status = "other_lasso_affinity_target_unconfirmed"
            reason = "Nearest published cyclase is an other-lasso comparator; target-family membership not established."
        else:
            status, reason = "unresolved", "Insufficient contextual/family evidence."
        triage.loc[index, "curation_status"] = status
        triage.loc[index, "curation_reason"] = reason
    triage = triage.drop(columns=["interpretation_status"])
    loci = pd.read_csv(tables / "maize_candidate_loci.tsv", sep="\t")
    triage = triage.merge(loci[["target_id", "C_left_edge_bp", "C_right_edge_bp"]], on="target_id", validate="one_to_one")
    triage.to_csv(tables / "maize_curated_candidates.tsv", sep="\t", index=False)

    literature_proteins = defaultdict(list)
    for record in SeqIO.parse(processed / "literature/literature_all_cds.faa", "fasta"):
        literature_proteins[seqhash(record.seq)].append(record.id)
    reference_C = SeqIO.to_dict(SeqIO.parse(processed / "literature/literature_unique_cyclases.faa", "fasta"))
    reference_manifest = pd.read_csv(processed / "literature/literature_bgc_manifest.tsv", sep="\t")
    archive_records = []
    with zipfile.ZipFile(ROOT / "data/raw/research_2026-09-10/triculamin_lasso_BGCs.zip") as archive:
        for member in reference_manifest.archive_member:
            record, _, _ = parse_member(member, archive.read(member))
            archive_records.append((member, str(record.seq).upper(), str(record.seq.reverse_complement()).upper()))
    core_rows, protein_rows = [], []
    for candidate, assignments in strong_roles.groupby("candidate_id"):
        features, contigs = read_genome(candidate.split("_")[0])
        lookup = {f["target_id"]: f for f in features}
        selected = [lookup[identifier] for identifier in assignments.protein_id]
        assert len({f["contig"] for f in selected}) == 1 and len({f["strand"] for f in selected}) == 1
        low = min(f["start_0based"] for f in selected)
        high = max(f["end_0based_exclusive"] for f in selected)
        dna = contigs[selected[0]["contig"]].seq[low:high]
        if selected[0]["strand"] == "-":
            dna = dna.reverse_complement()
        dna_text = str(dna).upper()
        exact_core_members = [member for member, forward, reverse in archive_records if dna_text in forward or dna_text in reverse]
        aa = lookup[candidate]["protein"]
        closest = triage.set_index("target_id").loc[candidate]
        assert closest.provenance_status == "culture_conflicts_with_MAG_label"
        assert closest.pident >= 99 and closest.query_coverage == 100
        ref_aa = str(reference_C[closest.subject_id].seq)
        substitutions = ";".join(f"{a}{i+1}{b}" for i, (a, b) in enumerate(zip(ref_aa, aa)) if a != b) if len(aa) == len(ref_aa) else "length_diff_requires_alignment"
        core_rows.append({"candidate_id": candidate, "core_start_0based": low, "core_end_0based_exclusive": high,
                          "core_bp": high-low, "core_sha256": seqhash(dna), "C_sha256": seqhash(aa),
                          "nearest_reference_C": closest.subject_id, "C_changes_vs_nearest": substitutions,
                          "core_exact_archive_matches": ";".join(exact_core_members),
                          "precursor_like_copies": int(assignments.role.eq("precursor_like").sum())})
        for assignment in assignments.itertuples():
            f = lookup[assignment.protein_id]
            protein_rows.append({"candidate_id": candidate, "protein_id": assignment.protein_id,
                                 "role": assignment.role, "aa_length": len(f["protein"]),
                                 "sequence_sha256": seqhash(f["protein"]), "sequence": f["protein"],
                                 "exact_literature_CDS_ids": ";".join(literature_proteins.get(seqhash(f["protein"]), []))})
    pd.DataFrame(core_rows).to_csv(tables / "maize_supported_core_comparison.tsv", sep="\t", index=False)
    pd.DataFrame(protein_rows).to_csv(tables / "maize_supported_core_proteins.tsv", sep="\t", index=False)

    # Catalogue-only IDs: the direct GFF search must not imply these were read.
    missing = pd.read_csv(tables / "maize_catalogue_GFF_differences.tsv", sep="\t")
    missing = missing.loc[missing.presence == "catalogue_only"].copy()
    missing_ids, mapping = set(missing.target_id), {}
    with (ROOT / "data/raw/s3/protein_catalogue-100.tsv").open() as handle:
        for line in handle:
            rep, member = line.rstrip().split("\t")
            if member in missing_ids:
                mapping[member] = rep
    missing["representative_id"] = missing.target_id.map(mapping)
    missing["exact_sequence_searched"] = missing.target_id == missing.representative_id
    missing.to_csv(tables / "maize_catalogue_only_audit.tsv", sep="\t", index=False)
    summary = {"status": "checkpoint_for_human_review_not_proof_of_MAG_absence",
               "curation_counts": triage.curation_status.value_counts().to_dict(),
               "provenance_counts": pd.read_csv(tables / "maize_provenance_audit.tsv", sep="\t").provenance_status.value_counts().to_dict(),
               "supported_canonical_MAG": int((triage.target_id.isin(supported_ids) & triage.provenance_status.eq("primary_metadata_supports_metagenomic_origin")).sum()),
               "supported_canonical_culture_loci": int((triage.target_id.isin(supported_ids) & triage.provenance_status.eq("culture_conflicts_with_MAG_label")).sum()),
               "supported_C_unique_sequences": len({r["C_sha256"] for r in core_rows}),
               "supported_defined_core_unique_DNA": len({r["core_sha256"] for r in core_rows}),
               "supported_cores_exact_in_literature_archive": sum(bool(r["core_exact_archive_matches"]) for r in core_rows),
               "catalogue_only_CDS_exact_sequences_searched": int(missing.exact_sequence_searched.sum()),
               "catalogue_only_CDS_representative_only": missing.loc[~missing.exact_sequence_searched, "target_id"].tolist(),
               "limitations": ["Low-score candidates retain provisional, not validated exclusion labels.",
                               "No target detection is not proof of absence in incomplete MAG/annotations.",
                               "A novel sequence relative to this reference is not necessarily new to all public databases.",
                               "Culture-origin hits cannot answer the MAG expansion question."]}
    (tables / "maize_checkpoint_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
