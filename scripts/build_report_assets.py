#!/usr/bin/env python3
"""Create auditable reference summaries and static figures from finished analyses."""

from __future__ import annotations

import argparse
import copy
import json
import itertools
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from Bio import AlignIO, Phylo, SeqIO
from matplotlib.patches import FancyArrowPatch, Patch


COLORS = {
    "Bacillota": "#2c7fb8",
    "Pseudomonadota": "#d95f0e",
    "Actinomycetota": "#756bb1",
    "Outgroup": "#666666",
}

GENE_COLORS = {
    "precursor": "#e41a1c",
    "RRE/B1": "#ffbf00",
    "B2 protease": "#4daf4a",
    "LrcC": "#377eb8",
    "transport": "#984ea3",
    "other": "#bdbdbd",
}


def classify_gene(row: pd.Series, selected_ids: set[str]) -> str:
    text = " ".join(
        str(row.get(field, "")) for field in
        ["label", "gene", "product", "gene_kind", "sec_met_domain", "gene_functions"]
    ).lower()
    label = str(row.get("label", "")).lower()
    if row["feature_id"] in selected_ids:
        return "LrcC"
    if label in {"lara", "lrca"} or (
        int(row["aa_length"]) <= 120
        and ("predicted lassopeptide" in text or "lasso peptide precursor" in text)
    ):
        return "precursor"
    if label in {"larb1", "lrcb1"} or "pqqd" in text or "lasso_rre" in text:
        return "RRE/B1"
    if label in {"larb2", "lrcb2"} or "biosynthesis b2" in text or "pf13471" in text:
        return "B2 protease"
    if label.startswith(("lard", "lrcd")) or any(
        term in text for term in ["transporter", "efflux", "permease", "abc transport"]
    ):
        return "transport"
    return "other"


def parse_support(label: str | None) -> tuple[float | None, float | None]:
    if not label or "/" not in label:
        return None, None
    try:
        left, right = label.split("/", maxsplit=1)
        return float(left), float(right)
    except ValueError:
        return None, None


def plot_tree(tree_path: Path, manifest: pd.DataFrame, out_path: Path, title: str) -> None:
    tree = Phylo.read(tree_path, "newick")
    outgroup = next(tree.find_clades(name="paeninodin_outgroup"), None)
    if outgroup is not None:
        tree.root_with_outgroup(outgroup)
    mapping = manifest.set_index("sample_id").to_dict("index")

    def label(clade):
        if not clade.is_terminal():
            return None
        if clade.name == "paeninodin_outgroup":
            return "Paeninodin cyclase (outgroup)"
        return mapping.get(clade.name, {}).get("display_name", clade.name)

    def support_label(clade):
        sh, uf = parse_support(clade.name)
        if sh is None or uf is None or uf < 80:
            return None
        return f"{uf:g}"

    tip_count = len(tree.get_terminals())
    fig, ax = plt.subplots(figsize=(15, max(8, 0.39 * tip_count + 2)))
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=label,
        branch_labels=support_label,
        show_confidence=False,
    )
    display_to_phylum = {
        row.display_name: row.phylum for row in manifest.itertuples(index=False)
    }
    display_to_phylum["Paeninodin cyclase (outgroup)"] = "Outgroup"
    for text in ax.texts:
        clean = text.get_text().strip()
        if clean in display_to_phylum:
            text.set_color(COLORS[display_to_phylum[clean]])
            text.set_fontweight("semibold")
            text.set_fontsize(9)
        else:
            text.set_fontsize(7)
    ax.set_title(title + "\nInternal labels: ultrafast bootstrap ≥80", loc="left", fontsize=13)
    ax.set_xlabel("Amino-acid substitutions per site")
    ax.set_ylabel("")
    present = list(dict.fromkeys(display_to_phylum.values()))
    ax.legend(
        handles=[Patch(color=COLORS[name], label=name) for name in present],
        frameon=False,
        loc="upper right",
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    rooted_path = out_path.with_suffix(".rooted.treefile")
    Phylo.write(tree, rooted_path, "newick")


def plot_architectures(
    features: pd.DataFrame,
    manifest: pd.DataFrame,
    selected: pd.DataFrame,
    sample_ids: list[str],
    out_path: Path,
) -> None:
    selected_ids = set(selected["feature_id"])
    labels = manifest.set_index("sample_id")["display_name"].to_dict()
    fig, ax = plt.subplots(figsize=(14, 5.5))
    for y, sample_id in enumerate(reversed(sample_ids)):
        part = features.loc[features["sample_id"] == sample_id].sort_values("start_0based")
        for _, row in part.iterrows():
            category = classify_gene(row, selected_ids)
            start = row["start_0based"] / 1000
            end = row["end_0based_exclusive"] / 1000
            if int(row["strand"]) >= 0:
                arrow_start, arrow_end = start, end
            else:
                arrow_start, arrow_end = end, start
            patch = FancyArrowPatch(
                (arrow_start, y), (arrow_end, y),
                arrowstyle="Simple,tail_width=0.55,head_width=5,head_length=6",
                color=GENE_COLORS[category], linewidth=0.5, mutation_scale=3.3,
            )
            ax.add_patch(patch)
        ax.text(-0.15, y, labels[sample_id], ha="right", va="center", fontsize=9)
    max_kb = features.loc[features["sample_id"].isin(sample_ids), "end_0based_exclusive"].max() / 1000
    ax.set_xlim(-0.2, max_kb + 0.3)
    ax.set_ylim(-0.65, len(sample_ids) - 0.15)
    ax.set_yticks([])
    ax.set_xlabel("Coordinates within submitted BGC record (kb)")
    ax.set_title("Representative published BGC records (gene functions are annotation-based)", loc="left")
    ax.legend(
        handles=[Patch(color=color, label=name) for name, color in GENE_COLORS.items()],
        ncol=3, frameon=False, loc="upper right",
    )
    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--canonical-alignment", type=Path, required=True)
    parser.add_argument("--canonical-tree", type=Path, required=True)
    parser.add_argument("--table-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t", keep_default_na=False)
    features = pd.read_csv(args.features, sep="\t", keep_default_na=False)
    selected = pd.read_csv(args.selected, sep="\t", keep_default_na=False)
    selected_ids = set(selected["feature_id"])
    features["functional_category"] = features.apply(classify_gene, axis=1, selected_ids=selected_ids)
    category_counts = (
        features.groupby(["sample_id", "functional_category"]).size().unstack(fill_value=0).reset_index()
    )
    bgc_summary = manifest[[
        "sample_id", "display_name", "accession", "phylum", "group", "record_bp",
        "cds_count", "translation_feature_count", "parser_warning_count"
    ]].merge(
        selected[["sample_id", "feature_id", "aa_length", "confidence", "decision_basis"]]
        .rename(columns={"aa_length": "lrcc_aa_length"}),
        on="sample_id", validate="one_to_one",
    ).merge(category_counts, on="sample_id", how="left", validate="one_to_one")
    args.table_dir.mkdir(parents=True, exist_ok=True)
    bgc_summary.to_csv(args.table_dir / "reference_bgc_summary.tsv", sep="\t", index=False)

    alignment = AlignIO.read(args.alignment, "fasta")
    canonical_alignment = AlignIO.read(args.canonical_alignment, "fasta")
    alignment_rows = []
    seq_to_ids: dict[str, list[str]] = defaultdict(list)
    for record in alignment:
        sequence = str(record.seq)
        ungapped = sequence.replace("-", "")
        seq_to_ids[ungapped].append(record.id)
        alignment_rows.append({
            "sample_id": record.id,
            "alignment_columns": len(sequence),
            "ungapped_aa": len(ungapped),
            "gap_fraction": sequence.count("-") / len(sequence),
        })
    alignment_qc = pd.DataFrame(alignment_rows)
    alignment_qc.to_csv(args.table_dir / "alignment_qc.tsv", sep="\t", index=False)
    duplicate_groups = [ids for ids in seq_to_ids.values() if len(ids) > 1]

    reference_alignment = [record for record in alignment if record.id != "paeninodin_outgroup"]
    identity_rows = []
    for left, right in itertools.combinations(reference_alignment, 2):
        comparable = [
            (a, b) for a, b in zip(str(left.seq), str(right.seq))
            if a not in {"-", "X", "?"} and b not in {"-", "X", "?"}
        ]
        identical = sum(a == b for a, b in comparable)
        identity_rows.append({
            "sample_a": left.id,
            "sample_b": right.id,
            "comparable_sites": len(comparable),
            "identical_sites": identical,
            "percent_identity": 100 * identical / len(comparable) if comparable else float("nan"),
        })
    identity = pd.DataFrame(identity_rows)
    identity.to_csv(args.table_dir / "lrcc_pairwise_identity.tsv", sep="\t", index=False)

    tree = Phylo.read(args.tree, "newick")
    supports = [parse_support(clade.name) for clade in tree.get_nonterminals()]
    supports = [(sh, uf) for sh, uf in supports if sh is not None and uf is not None]
    rooted_for_table = copy.deepcopy(tree)
    outgroup = next(rooted_for_table.find_clades(name="paeninodin_outgroup"), None)
    if outgroup is not None:
        rooted_for_table.root_with_outgroup(outgroup)
    taxonomy = manifest.set_index("sample_id")[["display_name", "phylum", "group"]].to_dict("index")
    clade_rows = []
    for index, clade in enumerate(rooted_for_table.get_nonterminals(order="preorder"), start=1):
        tips = [tip.name for tip in clade.get_terminals() if tip.name != "paeninodin_outgroup"]
        if len(tips) < 2:
            continue
        sh, uf = parse_support(clade.name)
        phyla = sorted({taxonomy[tip]["phylum"] for tip in tips})
        groups = sorted({taxonomy[tip]["group"] for tip in tips})
        clade_rows.append({
            "node_id": f"N{index:03d}",
            "tip_count": len(tips),
            "sh_alrt": "" if sh is None else sh,
            "ufboot": "" if uf is None else uf,
            "phyla": ";".join(phyla),
            "taxonomic_groups": ";".join(groups),
            "mixed_phylum": len(phyla) > 1,
            "mixed_taxonomic_group": len(groups) > 1,
            "sample_ids": ";".join(tips),
            "display_names": ";".join(taxonomy[tip]["display_name"] for tip in tips),
        })
    pd.DataFrame(clade_rows).to_csv(args.table_dir / "tree_clades.tsv", sep="\t", index=False)

    full_for_comparison = copy.deepcopy(tree)
    canonical_tree = Phylo.read(args.canonical_tree, "newick")
    common_tips = {tip.name for tip in canonical_tree.get_terminals()}
    for candidate_tree in [full_for_comparison, canonical_tree]:
        outgroup = next(candidate_tree.find_clades(name="paeninodin_outgroup"), None)
        if outgroup is not None:
            candidate_tree.root_with_outgroup(outgroup)
    for tip in list(full_for_comparison.get_terminals()):
        if tip.name not in common_tips:
            full_for_comparison.prune(tip)

    def clade_sets(candidate_tree):
        total = len(candidate_tree.get_terminals())
        return {
            frozenset(tip.name for tip in clade.get_terminals())
            for clade in candidate_tree.get_nonterminals()
            if 1 < len(clade.get_terminals()) < total
        }

    full_clades = clade_sets(full_for_comparison)
    canonical_clades = clade_sets(canonical_tree)
    shared_clades = full_clades & canonical_clades
    summary = {
        "reference_bgc_count": int(len(manifest)),
        "translated_cds_count": int(len(features)),
        "selected_lrcc_count": int(len(selected)),
        "unique_lrcc_sequence_count_including_outgroup": int(len(seq_to_ids)),
        "duplicate_sequence_groups": duplicate_groups,
        "selection_confidence": {key: int(value) for key, value in Counter(selected["confidence"]).items()},
        "phylum_counts": {key: int(value) for key, value in Counter(selected["phylum"]).items()},
        "lrcc_aa_range": [int(selected["aa_length"].min()), int(selected["aa_length"].max())],
        "alignment_columns": int(alignment.get_alignment_length()),
        "canonical_alignment_sequence_count": len(canonical_alignment),
        "canonical_alignment_columns": int(canonical_alignment.get_alignment_length()),
        "reference_pairwise_identity_percent_range": [
            round(float(identity["percent_identity"].min()), 3),
            round(float(identity["percent_identity"].max()), 3),
        ],
        "reference_pairwise_identity_percent_median": round(float(identity["percent_identity"].median()), 3),
        "sequences_over_50pct_gaps": alignment_qc.loc[alignment_qc["gap_fraction"] > 0.5, "sample_id"].tolist(),
        "tree_tip_count": len(tree.get_terminals()),
        "internal_branches_with_both_supports_ge_80": sum(sh >= 80 and uf >= 80 for sh, uf in supports),
        "internal_branches_with_support_labels": len(supports),
        "canonical_sensitivity_shared_internal_clades": len(shared_clades),
        "canonical_sensitivity_full_internal_clades": len(full_clades),
        "canonical_sensitivity_rebuilt_internal_clades": len(canonical_clades),
        "canonical_sensitivity_symmetric_difference": len(full_clades ^ canonical_clades),
    }
    (args.table_dir / "reference_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    plot_tree(args.tree, manifest, args.figure_dir / "lrcc_full_tree.png", "LrcC phylogeny: all 29 published BGCs")
    canonical_manifest = manifest.loc[manifest["phylum"] != "Actinomycetota"].copy()
    plot_tree(
        args.canonical_tree,
        canonical_manifest,
        args.figure_dir / "lrcc_canonical_tree.png",
        "Sensitivity analysis: 25 LrcC without the short Actinomycetota candidates",
    )
    plot_architectures(
        features,
        manifest,
        selected,
        [
            "paenibacillus_m2",
            "bacillus_subtilis_b1166",
            "paraburkholderia_mimosarum_stm3621",
            "streptomyces_griseocarneus_atcc29818",
        ],
        args.figure_dir / "representative_bgc_architectures.png",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
