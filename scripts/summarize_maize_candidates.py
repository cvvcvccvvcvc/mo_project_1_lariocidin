#!/usr/bin/env python3
"""Join evidence for human review without turning a profile score into a BGC call."""
import hashlib
import json

import pandas as pd
from Bio import SeqIO

from profile_search import ROOT

BLAST_COLUMNS = ["query_id", "subject_id", "pident", "alignment_length", "query_length", "subject_length",
                 "evalue", "bitscore", "query_coverage"]


def main():
    tables = ROOT / "results/tables"
    reference = pd.read_csv(ROOT / "data/processed/literature_unique_cyclases.tsv", sep="\t")
    hits = pd.read_csv(tables / "maize_direct_vs_literature.tsv", sep="\t", names=BLAST_COLUMNS)
    best = hits.sort_values(["bitscore", "subject_id"], ascending=[False, True]).drop_duplicates("query_id")
    best = best.merge(reference[["cyclase_id", "source_groups", "archive_members"]],
                      left_on="subject_id", right_on="cyclase_id", validate="many_to_one")
    scores = pd.read_csv(tables / "maize_direct_profile_hits.tsv", sep="\t")
    scores = scores.sort_values("sequence_score", ascending=False).drop_duplicates("target_id")
    summary = scores.merge(best, left_on="target_id", right_on="query_id", how="left", validate="one_to_one")
    summary["Genome"] = summary.target_id.str.split("_").str[0]
    metadata = pd.read_csv(ROOT / "config/maize_s3_genomes.tsv", sep="\t")
    provenance = pd.read_csv(tables / "maize_provenance_audit.tsv", sep="\t")
    summary = summary.merge(metadata, on="Genome", validate="many_to_one").merge(
        provenance[["Genome", "provenance_status", "sample_type", "strain"]], on="Genome", validate="many_to_one")
    reference_hashes = {}
    for r in SeqIO.parse(ROOT / "data/processed/literature_unique_cyclases.faa", "fasta"):
        reference_hashes[hashlib.sha256(str(r.seq).encode()).hexdigest()] = r.id
    exact = {}
    for r in SeqIO.parse(ROOT / "data/processed/maize_direct_candidates.faa", "fasta"):
        exact[r.id] = reference_hashes.get(hashlib.sha256(str(r.seq).encode()).hexdigest(), "")
    summary["exact_literature_sequence_id"] = summary.target_id.map(exact)
    summary["interpretation_status"] = "context_review_pending"
    summary.to_csv(tables / "maize_candidate_triage.tsv", sep="\t", index=False)
    neighbors = pd.read_csv(tables / "maize_candidate_neighbors.tsv", sep="\t")
    neighbor_hits = pd.read_csv(tables / "maize_neighbors_vs_literature.tsv", sep="\t", names=BLAST_COLUMNS)
    neighbor_best = neighbor_hits.sort_values("bitscore", ascending=False).drop_duplicates("query_id")
    feature_ref = pd.read_csv(ROOT / "data/processed/literature_all_cds.tsv", sep="\t")
    annotated = neighbors.merge(neighbor_best, left_on="target_id", right_on="query_id", how="left", validate="many_to_one")
    annotated = annotated.merge(feature_ref[["feature_id", "archive_member", "source_group", "label", "author_role", "product"]].rename(
        columns={"product": "reference_product"}), left_on="subject_id", right_on="feature_id", how="left", validate="many_to_one")
    annotated.to_csv(tables / "maize_neighbor_reference_matches.tsv", sep="\t", index=False)
    print(json.dumps({"candidate_CDS": len(summary),
                      "nearest_reference_groups": summary.source_groups.value_counts().to_dict(),
                      "priority_by_provenance": {str(k): int(v) for k, v in summary.groupby(["priority", "provenance_status"]).size().items()},
                      "exact_literature_sequence_matches": int(summary.exact_literature_sequence_id.ne("").sum())}, indent=2))


if __name__ == "__main__":
    main()
