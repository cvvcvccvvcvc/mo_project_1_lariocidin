#!/usr/bin/env python3
"""Summarize the prespecified host/BGC topology checks without reinterpreting support."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import pandas as pd
from Bio import AlignIO, Phylo
from assess_free_bgc_trees import model_from_report, validate_taxa


ROLES = ("B1", "B2", "C", "D")
FOCAL = "photorhabdus_asymbiotica_dsm15149"


def support_values(clade) -> tuple[float | None, float | None]:
    label = (clade.name or "") if clade is not None else ""
    match = re.fullmatch(r"([0-9.]+)/([0-9.]+)", label)
    return (float(match.group(1)), float(match.group(2))) if match else (None, None)


def split_clade(tree, group: set[str]):
    leaves = {leaf.name for leaf in tree.get_terminals()}
    if not group or not group < leaves:
        raise ValueError("Split must be a nonempty proper subset of tree taxa")
    complement = leaves - group
    for clade in tree.find_clades(order="preorder"):
        descendants = {leaf.name for leaf in clade.get_terminals()}
        if clade is not tree.root and descendants in (group, complement):
            return clade
    return None


def composition_failures(path: Path) -> int:
    match = re.search(r"\*\*\*\*\s+TOTAL.*?(\d+) sequences failed composition", path.read_text())
    if not match:
        raise ValueError(f"Composition summary missing from {path}")
    return int(match.group(1))


def au_constrained_row(path: Path) -> tuple[float, float, float]:
    # Parse the named columns in USER TREES, never an unrelated numbered table.
    section = path.read_text().split("USER TREES", 1)
    if len(section) != 2:
        raise ValueError(f"USER TREES section missing from {path}")
    columns, rows = None, {}
    for line in section[1].splitlines():
        fields = [field for field in line.split() if field not in ("+", "-")]
        if fields[:2] == ["Tree", "logL"]:
            columns = fields
        elif columns and fields and fields[0] in ("1", "2"):
            if len(fields) != len(columns) or fields[0] in rows:
                raise ValueError(f"Malformed or duplicate AU row in {path}")
            rows[fields[0]] = dict(zip(columns, map(float, fields)))
        elif len(rows) == 2:
            break
    if set(rows) != {"1", "2"} or "p-AU" not in rows["2"]:
        raise ValueError(f"Two candidate trees and p-AU required in {path}")
    logl, p = rows["2"]["logL"], rows["2"]["p-AU"]
    delta = rows["1"]["logL"] - logl
    if not all(map(math.isfinite, (logl, delta, p))) or not 0 <= p <= 1 or delta < -1e-4:
        raise ValueError(f"Invalid AU result or constrained search outperforms free search: {path}")
    return logl, rows["2"]["deltaL"], p


def apply_holm(rows: list[dict]) -> None:
    """Four prespecified genes; compatibility is an implicit p=1, not a measured test."""
    if len(rows) != len(ROLES) or {row["role"] for row in rows} != set(ROLES):
        raise ValueError("Expected exactly B1, B2, C, D")
    values = {row["role"]: row["constrained_au_p"] for row in rows}
    adjusted, running = {}, 0.0
    ordered = sorted(values, key=lambda role: values[role] if values[role] is not None else 1.0)
    for index, role in enumerate(ordered):
        value = values[role] if values[role] is not None else 1.0
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"Invalid p-value for {role}")
        running = max(running, min(1.0, (len(ROLES) - index) * value))
        adjusted[role] = running
    for row in rows:
        tested = values[row["role"]] is not None
        row["holm_adjusted_p"] = adjusted[row["role"]] if tested else None
        row["holm_reject_0_05"] = tested and adjusted[row["role"]] < 0.05
        row["h0_decision"] = (
            "constraint_rejected" if row["holm_reject_0_05"] else
            "constraint_not_rejected" if tested else "free_tree_compatible"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--taxa", type=Path, required=True)
    parser.add_argument("--free-assessment", type=Path, required=True)
    parser.add_argument("--alignment-dir", type=Path, required=True)
    parser.add_argument("--tree-dir", type=Path, required=True)
    parser.add_argument("--host-alignment", type=Path, required=True)
    parser.add_argument("--bgc-validation", type=Path, required=True)
    parser.add_argument("--tests-out", type=Path, required=True)
    parser.add_argument("--pairs-out", type=Path, required=True)
    parser.add_argument("--host-summary-out", type=Path, required=True)
    args = parser.parse_args()

    selected = pd.read_csv(args.selected, sep="\t", keep_default_na=False)
    selected = selected.loc[selected["phylum"] == "Pseudomonadota"]
    group_by_taxon = dict(zip(selected["sample_id"], selected["group"]))
    enterobacterales = set(selected.loc[selected["group"] == "Enterobacterales", "sample_id"])
    assessment = pd.read_csv(args.free_assessment, sep="\t", keep_default_na=False)
    validate_taxa(assessment.role, set(ROLES), "free-tree assessment")
    expected = set(selected.sample_id)

    rows = []
    for item in assessment.itertuples(index=False):
        role = item.role
        tree = Phylo.read(args.tree_dir / f"{role}.free.treefile", "newick")
        validate_taxa((leaf.name for leaf in tree.get_terminals()), expected, f"{role} free tree")
        if item.model != model_from_report(args.tree_dir / f"{role}.free.iqtree"):
            raise ValueError(f"Stale selected model for {role}")
        clade = split_clade(tree, enterobacterales)
        if not isinstance(item.order_split_present, bool) or item.order_split_present != (clade is not None):
            raise ValueError(f"Stale or invalid split assessment for {role}")
        sh_alrt, ufboot = support_values(clade) if clade else (None, None)
        leaves = {leaf.name for leaf in tree.get_terminals()}
        distance, nearest = min(
            (tree.distance(FOCAL, taxon), taxon) for taxon in leaves - {FOCAL}
        )
        alignment = AlignIO.read(args.alignment_dir / f"{role}.trim.faa", "fasta")
        validate_taxa((record.id for record in alignment), expected, f"{role} alignment")
        constrained_logl = delta_logl = p_au = None
        if not item.order_split_present:
            constrained = Phylo.read(args.tree_dir / f"{role}.constrained.treefile", "newick")
            validate_taxa((leaf.name for leaf in constrained.get_terminals()), expected, f"{role} constrained tree")
            if split_clade(constrained, enterobacterales) is None:
                raise ValueError(f"Constrained tree violates the order split: {role}")
            constrained_logl, delta_logl, p_au = au_constrained_row(
                args.tree_dir / f"{role}.au.iqtree"
            )
        rows.append({
            "role": role,
            "taxa": len(alignment),
            "trimmed_columns": alignment.get_alignment_length(),
            "model": item.model,
            "composition_failures": composition_failures(args.tree_dir / f"{role}.free.log"),
            "free_order_split_present": bool(item.order_split_present),
            "order_split_sh_alrt": sh_alrt,
            "order_split_ufboot": ufboot,
            "focal_nearest_taxon_by_tree_distance": nearest,
            "focal_nearest_taxonomic_group": group_by_taxon[nearest],
            "constrained_log_likelihood": constrained_logl,
            "delta_log_likelihood_vs_free": delta_logl,
            "constrained_au_p": p_au,
            "h0_decision": None,  # assigned only after the four-gene Holm correction
        })

    apply_holm(rows)
    pd.DataFrame(rows).to_csv(args.tests_out, sep="\t", index=False)

    taxa = pd.read_csv(args.taxa, sep="\t", keep_default_na=False)
    host_tree = Phylo.read(args.tree_dir / "host.treefile", "newick")
    host_alignment = AlignIO.read(args.host_alignment, "fasta")
    included = set(taxa.loc[taxa.host_tree_status == "included", "sample_id"])
    validate_taxa((leaf.name for leaf in host_tree.get_terminals()), included, "host tree")
    validate_taxa((record.id for record in host_alignment), included, "host alignment")
    host_sequences = {record.id: str(record.seq) for record in host_alignment}
    host_leaves = set(host_sequences)
    host_enterobacterales = set(
        taxa.loc[
            taxa["sample_id"].isin(host_leaves)
            & (taxa["taxonomic_group"] == "Enterobacterales"),
            "sample_id",
        ]
    )
    host_split = split_clade(host_tree, host_enterobacterales)
    host_sh, host_uf = support_values(host_split)
    pair_rows = []
    for row in taxa.loc[
        (taxa["genome_role"] == "same_species_control")
        & (taxa["host_tree_status"] == "included")
    ].itertuples(index=False):
        carrier = row.carrier_context
        mrca = host_tree.common_ancestor({"name": carrier}, {"name": row.sample_id})
        descendants = {leaf.name for leaf in mrca.get_terminals()}
        pair_clade = split_clade(host_tree, {carrier, row.sample_id})
        pair_sh, pair_uf = support_values(pair_clade)
        pair_rows.append({
            "carrier": carrier,
            "control": row.sample_id,
            "carrier_assembly": taxa.loc[taxa["sample_id"] == carrier, "assembly_accession"].iloc[0],
            "control_assembly": row.assembly_accession,
            "exact_concatenated_marker_match": host_sequences[carrier] == host_sequences[row.sample_id],
            "patristic_distance": host_tree.distance(carrier, row.sample_id),
            "mrca_leaf_count": len(descendants),
            "sister_pair": pair_clade is not None,
            "pair_sh_alrt": pair_sh,
            "pair_ufboot": pair_uf,
        })
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(args.pairs_out, sep="\t", index=False)

    validation = pd.read_csv(args.bgc_validation, sep="\t")
    summary = {
        "host_tree": {
            "taxa": len(host_alignment),
            "carrier_taxa": int(sum((taxa.genome_role == "carrier") & (taxa.host_tree_status == "included"))),
            "same_species_controls": int(sum(
                (taxa.genome_role == "same_species_control") & (taxa.host_tree_status == "included")
            )),
            "omitted_carriers": taxa.loc[
                (taxa.genome_role == "carrier") & (taxa.host_tree_status != "included"), "sample_id"
            ].tolist(),
            "markers": int(taxa["markers_total"].max()),
            "alignment_columns": host_alignment.get_alignment_length(),
            "composition_failures": composition_failures(args.tree_dir / "host.log"),
            "order_split_present": host_split is not None,
            "order_split_sh_alrt": host_sh,
            "order_split_ufboot": host_uf,
        },
        "same_species_controls": {
            "pairs": len(pairs),
            "sister_pairs": int(pairs["sister_pair"].sum()),
            "exact_marker_matches": int(pairs["exact_concatenated_marker_match"].sum()),
        },
        "assembly_translation_audit": {
            "exact_role_matches": int(validation["exact_translation_matches"].sum()),
            "tested_role_sequences": len(validation),
            "unannotated_role_sequences": int(sum(~validation["assembly_has_translated_cds"])),
        },
    }
    args.host_summary_out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
