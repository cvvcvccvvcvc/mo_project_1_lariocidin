#!/usr/bin/env python3
"""Focused integrity checks for the reference prototype."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import AlignIO, Phylo, SeqIO


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--tree", type=Path, required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t", keep_default_na=False)
    features = pd.read_csv(args.features, sep="\t", keep_default_na=False)
    selected = pd.read_csv(args.selected, sep="\t", keep_default_na=False)
    fasta = list(SeqIO.parse(args.reference_fasta, "fasta"))
    alignment = AlignIO.read(args.alignment, "fasta")
    tree = Phylo.read(args.tree, "newick")

    assert len(manifest) == 29, f"expected 29 BGC records, found {len(manifest)}"
    assert manifest["sample_id"].is_unique
    assert set(selected["sample_id"]) == set(manifest["sample_id"])
    assert selected["feature_id"].is_unique and len(selected) == 29
    assert selected["is_top_blast_hit"].all(), "a manual selection differs from its BGC's top BLASTP hit"
    assert set(selected["feature_id"]) <= set(features["feature_id"])
    assert len(fasta) == 29 and {r.id for r in fasta} == set(manifest["sample_id"])
    assert len(alignment) == 30 and len({len(r.seq) for r in alignment}) == 1
    expected_tree_ids = {r.id for r in alignment}
    assert {tip.name for tip in tree.get_terminals()} == expected_tree_ids
    m2 = selected.loc[selected["sample_id"] == "paenibacillus_m2"].iloc[0]
    assert m2["label"] == "LarC" and int(m2["aa_length"]) == 668
    assert m2["full_genome_locus_tag"] == "ACETMN_30970"
    assert m2["full_genome_product"] == "asparagine synthase-related protein"
    m2_record = manifest.loc[manifest["sample_id"] == "paenibacillus_m2"].iloc[0]
    assert m2_record["full_genome_match_orientation"] == "-"
    assert int(m2_record["full_genome_start_0based"]) == 6808154
    translated_cds = features.groupby("sample_id").size()
    translation_features = (
        manifest.set_index("sample_id")["translation_feature_count"].reindex(translated_cds.index)
    )
    assert translation_features.ge(translated_cds).all()
    assert (translation_features > translated_cds).any(), (
        "expected at least one source record where non-CDS features also carry translations"
    )
    print("OK: 29 BGCs, 29 explicit LrcC selections, 30 aligned/tree tips, M2 seed and CDS-only extraction validated")


if __name__ == "__main__":
    main()
