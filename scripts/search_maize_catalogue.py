#!/usr/bin/env python3
"""Audit all protein-cluster memberships, then screen the fixed catalogue.

100% identity with 80% overlap is NOT exact sequence deduplication. A hit on a
representative nominates each member for direct sequence/context verification.
"""
import hashlib
import json
import shutil
import tarfile
from collections import Counter

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

from profile_search import ROOT, search


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    raw = ROOT / "data/raw/s3"
    output = ROOT / "results/tables"
    metadata = pd.read_csv(ROOT / "data/raw/research_2026-09-10/maize_genomes-all_metadata.tsv", sep="\t")
    metadata["in_scope"] = ((metadata.Genome_type == "MAG") & metadata.Lineage.str.startswith("d__Bacteria;")
                            & (metadata.Completeness >= 50) & (metadata.Contamination < 10))
    metadata["priority_quality"] = (metadata.Completeness >= 90) & (metadata.Contamination < 5)
    assert metadata.Genome.is_unique and len(metadata) == 431 and metadata.in_scope.sum() == 399
    metadata.to_csv(ROOT / "config/maize_s3_genomes.tsv", sep="\t", index=False)
    paths = [raw / f"protein_catalogue-100.{extension}" for extension in ("faa", "tsv")]
    with tarfile.open(raw / "protein_catalogue-100.tar.gz", "r:gz") as archive:
        assert set(archive.getnames()) == {p.name for p in paths}
        for path in paths:
            if not path.exists():
                with archive.extractfile(path.name) as source, path.open("wb") as target:
                    shutil.copyfileobj(source, target)
    fasta, membership = paths
    representative_ids, counts, members = set(), Counter(), set()
    with fasta.open() as handle:
        for header, sequence in SimpleFastaParser(handle):
            identifier = header.split()[0]
            assert identifier not in representative_ids and sequence
            representative_ids.add(identifier)
    cluster_ids, member_representatives = set(), Counter()
    with membership.open() as handle:
        for line in handle:
            representative, member = line.rstrip().split("\t")
            assert representative in representative_ids and member not in members
            members.add(member)
            cluster_ids.add(representative)
            genome = member.split("_")[0]
            counts[genome] += 1
            if representative == member:
                member_representatives[genome] += 1
    assert cluster_ids == representative_ids and set(counts) == set(metadata.Genome)
    coverage = metadata.copy()
    coverage["catalogue_members"] = coverage.Genome.map(counts)
    coverage["own_representatives"] = coverage.Genome.map(member_representatives).fillna(0).astype(int)
    coverage.to_csv(output / "maize_catalogue_coverage.tsv", sep="\t", index=False)
    controls = json.loads((output / "final_profile_control_summary.json").read_text())
    assert not controls["missed_full_length_positives"] and not controls["missed_default_vs_max"]
    scores = search(fasta, ROOT / "results/hmmer/maize_catalogue")
    scores.to_csv(output / "maize_catalogue_profile_hits.tsv", sep="\t", index=False)
    hits = set(scores.target_id)
    rows = []
    with membership.open() as handle:
        for line in handle:
            representative, member = line.rstrip().split("\t")
            if representative in hits:
                rows.append({"representative_id": representative, "member_id": member,
                             "Genome": member.split("_")[0]})
    occurrences = pd.DataFrame(rows, columns=["representative_id", "member_id", "Genome"]).merge(
        metadata, on="Genome", validate="many_to_one")
    occurrences.to_csv(output / "maize_catalogue_candidate_members.tsv", sep="\t", index=False)
    with fasta.open() as source, (ROOT / "data/processed/maize_candidate_representatives.faa").open("w") as target:
        for header, sequence in SimpleFastaParser(source):
            if header.split()[0] in hits:
                target.write(f">{header}\n{sequence}\n")
    summary = {
        "stage": "representative_screen_not_direct_all_member_search_or_BGC_confirmation",
        "catalogue_genomes": len(metadata), "scope_MAG": int(metadata.in_scope.sum()),
        "priority_quality_MAG": int((metadata.in_scope & metadata.priority_quality).sum()),
        "representatives_searched": len(representative_ids), "mapped_members": len(members),
        "scope_members": int(coverage.loc[coverage.in_scope, "catalogue_members"].sum()),
        "scope_MAG_covered_by_membership": int((coverage.in_scope & (coverage.catalogue_members > 0)).sum()),
        "hit_representatives": len(hits), "hit_members_all": len(occurrences),
        "hit_members_in_scope": int(occurrences.in_scope.sum()),
        "hit_MAG_in_scope": int(occurrences.loc[occurrences.in_scope, "Genome"].nunique()),
        "inputs": {str(p.relative_to(ROOT)): digest(p) for p in paths},
        "profiles": {g: digest(ROOT / f"data/processed/literature_profiles/final.{g}.hmm")
                     for g in ("canonical", "noncanonical")},
        "limitation": "Membership covers the MAGs; representatives are not every exact sequence. MGnify clustering permits 80% target coverage.",
    }
    (output / "maize_catalogue_search_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
