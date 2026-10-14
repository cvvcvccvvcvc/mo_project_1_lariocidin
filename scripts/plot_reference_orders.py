#!/usr/bin/env python3
"""Recolour the saved S2 tree using its frozen reference taxonomic groups."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
from Bio import Phylo
from matplotlib.patches import Patch


COLORS = {
    "Bacillales": "#277D64",
    "Burkholderiales": "#236B9E",
    "Enterobacterales": "#C35B23",
    "Streptomycetales": "#785A9D",
    "Actinomycetales": "#AF467C",
    "Outgroup": "#707784",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.metadata.open(newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    metadata = {row["sample_id"]: row for row in records}
    if len(metadata) != len(records):
        raise ValueError("Duplicate reference sample IDs")
    tree = Phylo.read(args.tree, "newick")
    leaves = [leaf.name for leaf in tree.get_terminals()]
    if len(set(leaves)) != len(leaves) or set(leaves) != set(metadata) | {"paeninodin_outgroup"}:
        raise ValueError("Saved tree and reference metadata have different leaves")
    if set(row["group"] for row in records) - COLORS.keys():
        raise ValueError("A reference group lacks a colour")

    # Use the already rooted tree verbatim: no rooting, sorting or branch edits.
    def label(clade):
        if not clade.is_terminal():
            return None
        if clade.name == "paeninodin_outgroup":
            return "Paeninodin cyclase (outgroup)"
        return metadata[clade.name]["display_name"]

    def support(clade):
        if not clade.name or "/" not in clade.name:
            return None
        _, uf = map(float, clade.name.split("/"))
        return f"{uf:g}" if uf >= 80 else None

    plt.rcParams.update({"font.family": "DejaVu Sans"})
    fig, ax = plt.subplots(figsize=(16, max(8, .39 * len(leaves) + 2)))
    Phylo.draw(tree, axes=ax, do_show=False, label_func=label,
               branch_labels=support, show_confidence=False)
    name_to_group = {row["display_name"]: row["group"] for row in records}
    name_to_group["Paeninodin cyclase (outgroup)"] = "Outgroup"
    focal = metadata["photorhabdus_asymbiotica_dsm15149"]["display_name"]
    for text in ax.texts:
        name = text.get_text().strip()
        if name in name_to_group:
            text.set_color(COLORS[name_to_group[name]])
            text.set_fontweight("semibold")
            text.set_fontsize(9.5)
            if name == focal:
                text.set_bbox({"facecolor": "#FFF0E6", "edgecolor": "none", "pad": 3})
        else:
            text.set_fontsize(7)
    counts = Counter(row["group"] for row in records)
    ax.legend(handles=[Patch(color=color, label=(
        f"{group} ({counts[group]})" if group != "Outgroup" else "Outgroup"
    )) for group, color in COLORS.items()], loc="upper right", frameon=False, fontsize=10)
    ax.set_title("03б  То же дерево C: окраска по порядкам / группам reference\n"
                 "29 BGC + outgroup; подписи узлов — UFBoot ≥80", loc="left", fontsize=15, pad=12)
    ax.set_xlabel("Оценённые аминокислотные замены на позицию")
    ax.set_ylabel("")
    fig.text(.055, .017,
             "Окраска по зафиксированным метаданным проекта. Bacillales и Actinomycetales — исходные обозначения.\n"
             "Топология, корень, длины ветвей и поддержки сохранены. Цветовой конфликт не устанавливает HGT.",
             fontsize=10, color="#5C6B78")
    fig.tight_layout(rect=(0, .058, 1, 1))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Recoloured {len(records)} reference tips in {len(counts)} frozen groups")


if __name__ == "__main__":
    main()
