#!/usr/bin/env python3
"""Create the compact review-C figure for the bounded HGT case study."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo


COLORS = {"Enterobacterales": "#D55E00", "Burkholderiales": "#0072B2"}


def short_name(name: str) -> str:
    words = name.split()
    return f"{words[0][0]}. {' '.join(words[1:])}" if len(words) > 1 else name


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--host-tree", type=Path, required=True)
    parser.add_argument("--taxa", type=Path, required=True)
    parser.add_argument("--assemblies", type=Path, required=True)
    parser.add_argument("--tests", type=Path, required=True)
    parser.add_argument("--control-status", type=Path, required=True)
    parser.add_argument("--biosamples", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    taxa = pd.read_csv(args.taxa, sep="\t", keep_default_na=False)
    taxa = taxa.loc[taxa["host_tree_status"] == "included"].copy()
    assemblies = pd.read_csv(args.assemblies, sep="\t", keep_default_na=False)
    assembly_by_id = assemblies.set_index("sample_id")
    labels = {}
    group_by_label = {}
    for row in taxa.itertuples(index=False):
        assembly = assembly_by_id.loc[row.sample_id]
        source_label = assembly["source_strain"] or assembly["source_isolate"] or assembly["assembly_accession"]
        symbol = "●" if row.genome_role == "carrier" else "○"
        source_label = source_label.split(" (")[0]
        label = f"{symbol} {short_name(row.organism)} {source_label[:16]}"
        labels[row.sample_id] = label
        group_by_label[label] = row.taxonomic_group

    fig, (ax_tree, ax_test) = plt.subplots(
        1, 2, figsize=(17, 9), gridspec_kw={"width_ratios": [1.55, 1]}
    )
    fig.subplots_adjust(left=0.055, right=0.97, top=0.88, bottom=0.13, wspace=0.13)

    tree = Phylo.read(args.host_tree, "newick")
    tree.root_at_midpoint()
    Phylo.draw(
        tree,
        axes=ax_tree,
        do_show=False,
        show_confidence=False,
        label_func=lambda clade: labels.get(clade.name, "") if clade.is_terminal() else "",
    )
    for text in ax_tree.texts:
        if text.get_text().strip() in group_by_label:
            text.set_color(COLORS[group_by_label[text.get_text().strip()]])
            text.set_fontsize(9)
    ax_tree.set_xlim(0, max(tree.depths().values()) * 1.75)
    ax_tree.set_title("A  Independent host tree (21 genomes, 34 markers)", loc="left", fontsize=12)
    ax_tree.set_xlabel("amino-acid substitutions per site")
    ax_tree.set_ylabel("")
    ax_tree.text(0.01, -0.08, "● published carrier    ○ different strain label",
                 transform=ax_tree.transAxes, fontsize=9)
    ax_tree.text(0.01, -0.12, "Orange: Enterobacterales    Blue: Burkholderiales",
                 transform=ax_tree.transAxes, fontsize=9)

    tests = pd.read_csv(args.tests, sep="\t", keep_default_na=False)
    ax_test.axis("off")
    ax_test.set_title("B  Prespecified BGC-gene tests", loc="left", fontsize=12)
    columns = [(0.02, "Gene"), (0.24, "Free split"), (0.55, "AU p"), (0.76, "Comp. fail")]
    for x, title in columns:
        ax_test.text(x, 0.87, title, transform=ax_test.transAxes, fontsize=10, weight="bold")
    ax_test.plot([0.02, 0.98], [0.83, 0.83], transform=ax_test.transAxes, color="#888888", lw=0.7)
    for index, row in enumerate(tests.itertuples(index=False)):
        y = 0.75 - index * 0.12
        values = [row.role, "yes" if row.free_order_split_present else "no",
                  f"{float(row.constrained_au_p):.3f}" if row.constrained_au_p != "" else "not run",
                  f"{row.composition_failures}/14"]
        for (x, _), value in zip(columns, values):
            ax_test.text(x, y, value, transform=ax_test.transAxes, fontsize=11,
                         color="#B22222" if x == 0.76 and int(row.composition_failures) else "#222222")
    ax_test.text(0.02, 0.20, "B1/B2/D: order constraint not rejected",
                 transform=ax_test.transAxes, fontsize=11, weight="bold")
    ax_test.text(0.02, 0.15, "C: free tree already order-compatible",
                 transform=ax_test.transAxes, fontsize=11)
    ax_test.text(0.02, 0.10, "Host order split: SH-aLRT/UFBoot 100/100",
                 transform=ax_test.transAxes, fontsize=10)

    status = pd.read_csv(args.control_status, sep="\t", keep_default_na=False)
    distinct = status.loc[status["strain_relationship"] != "same_strain_reassembly"]
    supported = int(distinct["context_supported_bgc"].sum())
    ax_test.text(0.02, 0.02,
                 f"Same-species controls: target BGC in {supported}/{len(distinct)} different strain labels",
                 transform=ax_test.transAxes, fontsize=10)

    fig.suptitle("Focal HGT case: no significant rejection of host-order split", fontsize=17)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
