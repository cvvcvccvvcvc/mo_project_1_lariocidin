#!/usr/bin/env python3
"""Final-profile operational controls; not an independent accuracy benchmark."""
import hashlib
import json

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from profile_search import ROOT, search


def main():
    processed = ROOT / "data/processed"
    labels, sequences = [], []

    def add(identifier, sequence, group, parent="", source=""):
        sequences.append(SeqRecord(Seq(str(sequence)), id=identifier, description=""))
        labels.append({"target_id": identifier, "group": group, "parent_id": parent,
                       "source": source, "aa_length": len(sequence)})

    table = pd.read_csv(processed / "literature/literature_unique_cyclases.tsv", sep="\t").set_index("cyclase_id")
    for record in SeqIO.parse(processed / "literature/literature_unique_cyclases.faa", "fasta"):
        source_group = table.loc[record.id, "source_groups"]
        group = source_group if source_group in ("canonical", "noncanonical") else "other_lasso"
        add(record.id, record.seq, group, source="Merrild archive")
        if group in ("canonical", "noncanonical"):
            for fraction in (.5, .75):
                size = int(len(record) * fraction)
                for terminus, seq in (("N", record.seq[:size]), ("C", record.seq[-size:])):
                    add(f"{record.id}_{terminus}{int(fraction*100)}", seq,
                        f"{group}_fragment_{terminus}{int(fraction*100)}", record.id)
    for record in SeqIO.parse(ROOT / "data/raw/s3/uniprot_reviewed_bacterial_asnb.faa", "fasta"):
        add(record.id, record.seq, "reviewed_asnb", source=record.description)
    for path in sorted((ROOT / "data/raw/focal_context").glob("*.gb")):
        for record in SeqIO.parse(path, "genbank"):
            for index, feature in enumerate(record.features):
                if feature.type == "CDS" and "translation" in feature.qualifiers:
                    identifier = f"bg_{path.stem}_{index}"
                    add(identifier, feature.qualifiers["translation"][0], "mixed_genomic_background",
                        source="|".join([record.id, feature.qualifiers.get("protein_id", [""])[0],
                                         feature.qualifiers.get("product", [""])[0]]))
    assert len({record.id for record in sequences}) == len(sequences)
    fasta = processed / "reference/final_profile_controls.faa"
    SeqIO.write(sequences, fasta, "fasta")
    labels = pd.DataFrame(labels)
    labels.to_csv(processed / "reference/final_profile_control_labels.tsv", sep="\t", index=False)
    summaries, selected_by_mode = {}, {}
    for mode in ("default", "max"):
        scores = search(fasta, ROOT / f"results/hmmer/controls_{mode}", exhaustive=mode == "max")
        scores.merge(labels, on="target_id", validate="many_to_one").to_csv(
            ROOT / f"results/tables/final_profile_controls_{mode}.tsv", sep="\t", index=False)
        best = scores.sort_values("sequence_score", ascending=False).drop_duplicates("target_id")
        combined = labels.merge(best, on="target_id", how="left", validate="one_to_one")
        selected_by_mode[mode] = set(best.target_id)
        summaries[mode] = {
            group: {"n": len(rows), "detected": int(rows.sequence_score.notna().sum()),
                    "high_priority": int((rows.priority == "high").sum()),
                    "minimum_reported_score": None if rows.sequence_score.isna().all() else float(rows.sequence_score.min())}
            for group, rows in combined.groupby("group")
        }
    positive_ids = set(labels[labels.group.isin(["canonical", "noncanonical"])].target_id)
    summary = {"status": "operational_control_not_external_validation", "score_threshold": 80,
               "mode_summary": summaries,
               "missed_default_vs_max": sorted(selected_by_mode["max"] - selected_by_mode["default"]),
               "missed_full_length_positives": sorted(positive_ids - selected_by_mode["default"]),
               "control_fasta_sha256": hashlib.sha256(fasta.read_bytes()).hexdigest()}
    (ROOT / "results/tables/final_profile_control_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    assert not summary["missed_full_length_positives"], "Stop before MAG: full-length control lost"
    assert not (positive_ids & set(summary["missed_default_vs_max"])), "Stop: heuristic sensitivity issue"


if __name__ == "__main__":
    main()
