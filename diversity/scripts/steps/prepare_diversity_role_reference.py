#!/usr/bin/env python3
"""Expose annotation and published-position evidence for role homology.

Published B1/B2 positional assignments are hypotheses, not experimentally
verified functions. Their provenance is retained when a candidate matches.
"""
import re
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from profile_search import ROOT


def main():
    processed = ROOT / "data/processed"
    features = pd.read_csv(processed / "literature/literature_all_cds.tsv", sep="\t", keep_default_na=False)
    roles = defaultdict(set)
    basis = defaultdict(set)
    for r in features.to_dict("records"):
        text = (r["label"] + " " + r["product"]).lower()
        if r["author_role"] == "precursor":
            roles[r["feature_id"]].add("A")
            basis[r["feature_id"]].add("explicit_published_precursor_label")
        if r["author_role"] == "cyclase":
            roles[r["feature_id"]].add("C")
            basis[r["feature_id"]].add("explicit_published_cyclase_label")
        for role, match in [("B1", re.search(r"\bb1\b|\bpqqd\b|\brre\b", text)),
                            ("B2", re.search(r"\bb2\b|lasso peptide biosynthesis protein", text)),
                            ("ABC", "abc transporter" in text),
                            ("GNAT", "gnat" in text or "acetyltransferase" in text),
                            ("C11", "clostripain" in text)]:
            if match:
                roles[r["feature_id"]].add(role)
                basis[r["feature_id"]].add("published_label_or_annotation_not_function_experiment")
    for _, group in features[features.source_group == "canonical"].groupby("bgc_id"):
        for c in group[group.author_role == "cyclase"].to_dict("records"):
            ordered = group[group.strand == c["strand"]].sort_values("start_0based", ascending=c["strand"] == 1)
            records = ordered.to_dict("records")
            ci = next(i for i, r in enumerate(records) if r["feature_id"] == c["feature_id"])
            precursors = [i for i in range(ci) if records[i]["author_role"] == "precursor"]
            if not precursors:
                continue
            between = records[precursors[-1]+1:ci]
            if len(between) == 2 and 40 <= between[0]["aa_length"] <= 250 and 80 <= between[1]["aa_length"] <= 400:
                for role, r in zip(("B1", "B2"), between):
                    roles[r["feature_id"]].add(role)
                    basis[r["feature_id"]].add("position_in_published_canonical_A_B1_B2_C_hypothesis")
    selected = features[features.feature_id.isin(roles)].copy()
    selected["inferred_roles"] = selected.feature_id.map(lambda x: ";".join(sorted(roles[x])))
    selected["role_evidence_basis"] = selected.feature_id.map(lambda x: ";".join(sorted(basis[x])))
    selected.to_csv(processed / "discovery/diversity_role_reference.tsv", sep="\t", index=False)
    records = [r for r in SeqIO.parse(processed / "literature/literature_all_cds.faa", "fasta") if r.id in roles]
    SeqIO.write(records, processed / "discovery/diversity_role_reference.faa", "fasta")
    print(f"role-reference proteins: {len(records)}; positional hypotheses explicitly flagged")


if __name__ == "__main__":
    main()
