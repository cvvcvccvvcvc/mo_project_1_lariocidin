#!/usr/bin/env python3
"""Search every annotated CDS of the fixed 399-MAG sample, without clustering."""
import hashlib
import json
from collections import defaultdict

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

from maize_genomes import digest, read_genome
from profile_search import ROOT, search


def main():
    table = pd.read_csv(ROOT / "results/tables/maize_catalogue_coverage.tsv", sep="\t")
    scope = table[table.in_scope]
    catalogue_members = defaultdict(set)
    with (ROOT / "data/raw/s3/protein_catalogue-100.tsv").open() as handle:
        for line in handle:
            member = line.rstrip().split("\t")[1]
            catalogue_members[member.split("_")[0]].add(member)
    representative_hashes = {}
    with (ROOT / "data/raw/s3/protein_catalogue-100.faa").open() as handle:
        for header, seq in SimpleFastaParser(handle):
            representative_hashes[header.split()[0]] = hashlib.sha256(seq.encode()).digest()
    interim = ROOT / "data/interim/maize"
    interim.mkdir(parents=True, exist_ok=True)
    fasta = interim / "all_399_MAG_annotated_CDS.faa"
    coverage, mismatches, membership_differences = [], [], []
    with fasta.open("w") as handle:
        for i, row in enumerate(scope.to_dict("records"), 1):
            features, sequences = read_genome(row["Genome"])
            assert sum(len(s) for s in sequences.values()) == row["Length"], row["Genome"]
            compared, mismatch_count, internal_stops = 0, 0, 0
            identifiers = set()
            for feature in features:
                identifier, protein = feature["target_id"], feature["protein"]
                assert identifier not in identifiers and identifier.startswith(row["Genome"] + "_")
                identifiers.add(identifier)
                internal_stops += "*" in protein
                if identifier in representative_hashes:
                    compared += 1
                    if hashlib.sha256(protein.encode()).digest() != representative_hashes[identifier]:
                        mismatch_count += 1
                        mismatches.append(identifier)
                handle.write(f">{identifier}\n{protein}\n")
            expected_ids = catalogue_members[row["Genome"]]
            catalogue_only = expected_ids - identifiers
            gff_only = identifiers - expected_ids
            membership_differences.extend({"Genome": row["Genome"], "target_id": identifier,
                                           "presence": "catalogue_only"} for identifier in sorted(catalogue_only))
            membership_differences.extend({"Genome": row["Genome"], "target_id": identifier,
                                           "presence": "GFF_only"} for identifier in sorted(gff_only))
            coverage.append({"Genome": row["Genome"], "annotated_CDS": len(features),
                             "representatives_compared": compared, "translation_mismatches": mismatch_count,
                             "internal_stop_CDS": internal_stops, "catalogue_only": len(catalogue_only),
                             "GFF_only": len(gff_only)})
            if i % 50 == 0:
                print(f"Translated/audited MAG: {i}/399", flush=True)
    pd.DataFrame(coverage).to_csv(ROOT / "results/tables/maize_direct_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(membership_differences, columns=["Genome", "target_id", "presence"]).to_csv(
        ROOT / "results/tables/maize_catalogue_GFF_differences.tsv", sep="\t", index=False)
    (ROOT / "results/tables/maize_translation_mismatches.json").write_text(json.dumps(mismatches, indent=2) + "\n")
    assert not mismatches, f"Investigate translation before search: {len(mismatches)} catalogue mismatches"
    scores = search(fasta, ROOT / "results/hmmer/maize_direct")
    scores.to_csv(ROOT / "results/tables/maize_direct_profile_hits.tsv", sep="\t", index=False)
    hit_ids = set(scores.target_id)
    rep_members = pd.read_csv(ROOT / "results/tables/maize_catalogue_candidate_members.tsv", sep="\t")
    nominated = set(rep_members.loc[rep_members.in_scope, "member_id"])
    with fasta.open() as handle, (ROOT / "data/processed/maize_direct_candidates.faa").open("w") as output:
        for header, protein in SimpleFastaParser(handle):
            if header in hit_ids:
                output.write(f">{header}\n{protein}\n")
    summary = {"status": "all_annotated_CDS_searched_BGC_curation_pending", "MAG_searched": len(coverage),
               "CDS_searched": sum(r["annotated_CDS"] for r in coverage),
               "translated_representatives_checked": sum(r["representatives_compared"] for r in coverage),
               "translation_mismatches": len(mismatches), "candidate_proteins": len(hit_ids),
               "catalogue_members_absent_from_GFF": sum(r["catalogue_only"] for r in coverage),
               "GFF_CDS_absent_from_catalogue": sum(r["GFF_only"] for r in coverage),
               "candidate_MAG": len({i.split("_")[0] for i in hit_ids}),
               "new_vs_representative_screen": sorted(hit_ids - nominated),
               "nominated_but_not_direct_hit": sorted(nominated - hit_ids),
               "input_sha256": digest(fasta),
               "limitation": "Direct search covers GFF CDS; catalogue-only members remain represented by the separate cluster screen. Missing gene calls/unassembled DNA remain undetectable."}
    (ROOT / "results/tables/maize_direct_search_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
