#!/usr/bin/env python3
"""Static scientific tree figures from completed H2 estimates."""
import copy
import argparse
from collections import defaultdict
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Patch
from matplotlib.patches import FancyArrow
from Bio import Phylo

from prepare_h2 import OUT, read
from h2_compare_trees import branch_support

COLORS = {"Actinomycetota": "#a9473b", "Bacillota": "#526ca9", "Pseudomonadota": "#278473", "Cyanobacteriota": "#7b62a3", "Planctomycetota": "#a17c2b"}
qc = {r["host_id"]: r for r in read(OUT / "marker_genome_QC.tsv")}


def layout(tree):
    tips = tree.get_terminals()
    y = {node: i for i, node in enumerate(tips)}
    for node in tree.get_nonterminals(order="postorder"):
        y[node] = (y[node.clades[0]] + y[node.clades[-1]]) / 2
    depths = tree.depths()
    maximum = max(depths.values()) or 1
    return {node: value / maximum for node, value in depths.items()}, y


def draw(ax, tree, labels=True, mirror=False, compact=False):
    x, y = layout(tree)
    transform = (lambda value: 1 - value) if mirror else (lambda value: value)
    for node in tree.find_clades(order="preorder"):
        if node.clades:
            ax.plot([transform(x[node])] * 2, [y[node.clades[0]], y[node.clades[-1]]], color="#555555", lw=0.6)
        for child in node.clades:
            ax.plot([transform(x[node]), transform(x[child])], [y[child]] * 2, color="#555555", lw=0.6)
        if not node.is_terminal() and not compact:
            sh, uf = branch_support(node)
            if sh >= 80 and uf >= 95:
                ax.plot(transform(x[node]), y[node], "o", color="#222222", markersize=2)
        if node.is_terminal() and labels:
            host = node.name.split("_", 1)[0]
            meta = qc[host]
            suffix = node.name[len(host):]
            species = meta["taxonomy_name"]
            label = f"{host}{suffix} {species}"
            if not suffix:
                label += " " + meta["sequence_provider_assembly"]
            ax.text(transform(x[node]) + (-0.02 if mirror else 0.02), y[node], label,
                    fontsize=5.5 if compact else 6, va="center", ha="right" if mirror else "left", color=COLORS.get(meta["phylum"], "#555555"))
    ax.set_ylim(len(tree.get_terminals()), -1)
    ax.axis("off")
    ax.set_xlim((-1.4, 1.02) if mirror else (-0.02, 2.4))
    return {node.name: (transform(x[node]), y[node]) for node in tree.get_terminals()}


def plot_nature_context(directory):
    rows = read(OUT / "Nature_native_neighborhood.tsv")
    groups = defaultdict(list)
    for row in rows:
        groups[row["host_id"]].append(row)
    fig, axes = plt.subplots(len(groups), 1, figsize=(14, 9))
    for ax, (host, genes) in zip(axes, sorted(groups.items())):
        c_start = int(genes[0]["C_start_0based"])
        for i, gene in enumerate(sorted(genes, key=lambda x: int(x["gene_start_0based"]))):
            start, end = int(gene["gene_start_0based"]) - c_start, int(gene["gene_end_0based"]) - c_start
            domains = gene["domain_accessions"].split(";")
            label, color = "", "#c6c9ce"
            if int(gene["gene_start_0based"]) == int(gene["C_start_0based"]):
                label, color = "target C" + gene["protein_aa"], "#b53e38"
            elif "PF00733" in domains:
                label, color = "Asn-like " + gene["protein_aa"], "#547cba"
            elif "PF05402" in domains:
                label, color = "PqqD " + gene["protein_aa"], "#ce8734"
            elif "PF13471" in domains:
                label, color = "Transglut " + gene["protein_aa"], "#ce8734"
            elif "PF00005" in domains or "PF00664" in domains:
                label, color = "ABC " + gene["protein_aa"], "#3b8c78"
            elif any(x in domains for x in ["PF00583", "PF13302", "PF13508"]):
                label, color = "GNAT " + gene["protein_aa"], "#9977ac"
            elif "lasso" in gene["product"].lower() or "ripp" in gene["product"].lower():
                label, color = "precursor? " + gene["protein_aa"], "#444444"
            direction = int(gene["strand"])
            x, dx = (start, end - start) if direction == 1 else (end, start - end)
            ax.add_patch(FancyArrow(x, 0, dx, 0, width=0.22, head_width=0.4,
                head_length=min(120, abs(dx) * 0.3), length_includes_head=True, color=color))
            if label:
                ax.text((start + end) / 2, 0.38 + (i % 2) * 0.4, label, fontsize=7, rotation=35, ha="center")
        ax.axvline(0, color="#b53e38", lw=0.5, ls=":")
        ax.set_xlim(-12500, 14000); ax.set_ylim(-0.5, 2.1); ax.set_yticks([])
        ax.set_title(qc[host]["taxonomy_name"] + " · " + qc[host]["sequence_provider_assembly"], loc="left", fontsize=10)
        for spine in ["left", "right", "top"]:
            ax.spines[spine].set_visible(False)
        ax.set_xlabel("bp relative to the DNA-verified target C start", fontsize=8)
    fig.suptitle("Four Nature systems: native genomic context", fontsize=12)
    fig.text(0.02, 0.015, "Labels show native annotations or domain families; activity and exact BGC boundaries are unproven. Grey: other annotated genes. Nearby loci remain distinct.", fontsize=8)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    fig.savefig(directory / "Nature_native_context.pdf")
    fig.savefig(directory / "Nature_native_context.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--context-only", action="store_true")
    parser.add_argument("--ready-only", action="store_true")
    args = parser.parse_args()
    directory = OUT / "figures"
    directory.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})
    plot_nature_context(directory)
    if args.context_only:
        print(json.dumps({"figure": "Nature_native_context.pdf"}))
        return
    inventory = []
    for scope in [OUT, OUT / "actino"]:
        prefix = "Actino138_" if scope.name == "actino" else "Bacteria74_"
        for job in read(scope / "tree_jobs.tsv"):
            source = scope / ("host_trees" if job["kind"] == "host" else "component_trees") / job["name"] / "ML.treefile"
            if not source.exists() or (args.ready_only and not source.with_name("task_summary.json").exists()):
                continue
            tree = Phylo.read(source, "newick")
            tree.ladderize()
            fig, ax = plt.subplots(figsize=(15, max(5, len(tree.get_terminals()) * 0.17)))
            draw(ax, tree)
            diagnostic = job['name'] == 'C_domain_background_diagnostic'
            ax.set_title(prefix + job["name"] + "\n" +
                         ("Approximate ML domain diagnostic; 1000 SH-aLRT only. " if diagnostic else "") +
                         "Unrooted estimate; display orientation has no historical direction", fontsize=10, loc="left")
            handles = [Patch(color=c, label=p) for p, c in COLORS.items() if any(qc[x.name.split("_", 1)[0]]["phylum"] == p for x in tree.get_terminals())]
            if handles:
                fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.015, 0.035), fontsize=7, ncol=len(handles), frameon=False)
            caption = ("Domain diagnostic: no UFBoot, no combined-support dots; background function is unconfirmed. " if diagnostic else
                       "Black dot: SH-aLRT ≥80 and UFBoot ≥95. ")
            fig.text(0.02, 0.014, caption + "Branch lengths normalized within this figure; exact lengths/support in Newick.", fontsize=7)
            fig.tight_layout(rect=(0, 0.08, 1, 1))
            path = directory / (prefix + job["name"] + ".pdf")
            fig.savefig(path)
            if job["name"] in ["canonical", "noncanonical", "canonical_C", "noncanonical_C", "uncertain_Nature", "uncertain_Nature_system_C"]:
                fig.savefig(path.with_suffix(".png"), dpi=130)
            plt.close(fig)
            inventory.append({"figure": path.name, "source_tree": str(source.relative_to(OUT)), "display": "arbitrary_unrooted_orientation; no_donor_or_recipient_inference"})
        for lineage in ["canonical", "noncanonical"]:
            hp = scope / "host_trees" / lineage / "ML.treefile"
            gp = scope / "component_trees" / (lineage + "_C") / "ML.treefile"
            if not hp.exists() or not gp.exists() or (args.ready_only and not hp.with_name("task_summary.json").exists()):
                continue
            host, gene = Phylo.read(hp, "newick"), Phylo.read(gp, "newick")
            by_genus = defaultdict(list)
            for tip in gene.get_terminals():
                by_genus[qc[tip.name.split("_", 1)[0]]["genus"]].append(tip.name)
            # A small taxonomic overview; full trees/copies remain in separate figures.
            keep = {sorted(names)[0] for names in by_genus.values()}
            keep_hosts = {name.split("_", 1)[0] for name in keep}
            for leaf in list(gene.get_terminals()):
                if leaf.name not in keep:
                    gene.prune(leaf)
            for leaf in list(host.get_terminals()):
                if leaf.name not in keep_hosts:
                    host.prune(leaf)
            host.ladderize(); gene.ladderize()
            fig, (left, right) = plt.subplots(1, 2, figsize=(17, max(5, len(keep) * 0.25)))
            a = draw(left, host, compact=True)
            b = draw(right, gene, mirror=True, compact=True)
            for name, point in b.items():
                h = name.split("_", 1)[0]
                # Connect rows across the central gap, leaving taxon labels unobstructed.
                fig.add_artist(ConnectionPatch((2.4, a[h][1]), (-1.4, point[1]), "data", "data", axesA=left, axesB=right,
                    color=COLORS.get(qc[h]["phylum"], "#555555"), linewidth=0.5, alpha=0.45, zorder=1))
            left.set_title("Independent host markers", fontsize=10)
            right.set_title("Published BGC C family", fontsize=10)
            fig.suptitle(prefix + lineage + " — one deterministic representative per genus", fontsize=11)
            fig.text(0.02, 0.014, "Overview subset; support omitted after pruning. Crossings are topology differences, not HGT events. See full trees and shared-taxon comparisons.", fontsize=7)
            fig.tight_layout(rect=(0, 0.04, 1, 0.96))
            fig.subplots_adjust(wspace=0.45)
            path = directory / (prefix + lineage + "_C_overview.pdf")
            fig.savefig(path); fig.savefig(path.with_suffix(".png"), dpi=150)
            plt.close(fig)
            inventory.append({"figure": path.name, "shown_hosts": ";".join(sorted(keep_hosts)), "display": "one_lexicographic_C_tip_per_genus; overview_only; no_support_or_history_claim_after_pruning"})
    (directory / ("figure_inventory_partial.json" if args.ready_only else "figure_inventory.json")).write_text(json.dumps(inventory, indent=2) + "\n")
    print(json.dumps({"figures": len(inventory)}))


if __name__ == "__main__":
    main()
