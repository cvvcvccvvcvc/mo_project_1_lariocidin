#!/usr/bin/env python3
"""Static scientific figures from the completed catalogue and saved contexts."""
import json
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
from profile_search import ROOT


COLORS = {"A": "#d58b35", "B1": "#6d91b4", "B2": "#65a890", "C": "#a363a9"}


def main():
    tables = ROOT / "results/tables"
    primary = pd.read_csv(tables / "diversity_species_BGC_catalogue.tsv", sep="\t", keep_default_na=False)
    summary = json.loads((tables / "diversity_species_summary.json").read_text())
    out = ROOT / "results/figures"
    out.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5))
    statuses = summary["context_status_counts"]
    axes[0].barh(range(len(statuses)), list(statuses.values()), color=["#327990" if x == "family_bgc_candidate" else "#a8a6a0" for x in statuses])
    axes[0].set_yticks(range(len(statuses)), [x.replace("_", " ") for x in statuses])
    axes[0].invert_yaxis()
    axes[0].set_title("A. Context outcomes", loc="left")
    axes[0].set_xlabel("Source locus rows / unresolved nominations")
    for i, count in enumerate(statuses.values()):
        axes[0].text(count, i, f" {count}", va="center")
    for axis, metric, novel, title in [
        (axes[1], "C_variants", "C_amino_acid_variants_with_unrepresented_DNA", "B. Distinct C protein records"),
        (axes[2], "defined_core_variants", "defined_cores_not_exact_published", "C. Distinct defined DNA spans")]:
        families = summary["families"]
        names = list(families)
        new = [families[name][novel] for name in names]
        starts = [families[name]["C_start_annotation_differences_on_published_DNA"] if metric == "C_variants" else 0 for name in names]
        exact = [families[name][metric]-families[name][novel]-start for name,start in zip(names,starts)]
        axis.bar(names, exact, label="Exact in frozen published set", color="#a8a6a0")
        axis.bar(names, starts, bottom=exact, label="Changed start on known DNA", color="#d58b35")
        axis.bar(names, new, bottom=[a+b for a,b in zip(exact,starts)], label="Variant with unrepresented DNA", color="#327990")
        axis.set_title(title, loc="left")
        axis.set_ylabel("Sequence variants")
        for i, (first, second, third) in enumerate(zip(exact, starts, new)):
            axis.text(i, first+second+third, f"{first+second+third}", ha="center", va="bottom")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.65,.105), ncol=3, fontsize=8)
    fig.text(.025, .015, f"{summary['proteomes']:,} annotated proteomes; {summary['proteins']:,} proteins. Supported candidates remain functional predictions.\n"
             "Novelty is relative to 29 LAR + 393 Merrild exported records; unconfirmed homologs remain separate.", fontsize=9)
    fig.tight_layout(rect=(0,.22,1,1), w_pad=2)
    for suffix in ("png", "svg"):
        fig.savefig(out / f"diversity_species_summary.{suffix}", dpi=170)
    plt.close(fig)
    if primary.empty:
        return
    # Show the most divergent C in each source/family represented by the main
    # catalogue. These are descriptive examples, not a second selected endpoint.
    examples = primary.sort_values(["C_reference_identity_percent", "locus_id"]).groupby(["source", "family"], sort=True).head(1)
    features = pd.concat([pd.read_csv(p, sep="\t", keep_default_na=False) for p in
                         sorted((ROOT / "data/processed").glob("diversity_species_[0-7]_contexts/features.tsv"))], ignore_index=True)
    fig, axes = plt.subplots(len(examples), 1, figsize=(14, max(3, len(examples)*2.3)), squeeze=False)
    for axis, row in zip(axes[:,0], examples.to_dict("records")):
        roles = defaultdict(set)
        for gene, role in zip(row["chosen_core_gene_ids"].split(";"), row["architecture"].split("-")):
            roles[gene].update(role.split("+"))
        low = max(int(row["core_start"])-1800, int(row["window_start"]))
        high = min(int(row["core_end"])+1800, int(row["window_end"]))
        local = features[(features.locus_id == row["locus_id"]) & (features.start >= low) & (features.end <= high)]
        local = local[(local.origin == "annotated_CDS") | local.query_id.isin(roles)]
        for gene in local.to_dict("records"):
            labels = sorted(roles[gene["query_id"]])
            x = (gene["start"]-int(row["core_start"]))/1000
            length = (gene["end"]-gene["start"])/1000
            start, delta = (x, length) if gene["strand"] == 1 else (x+length, -length)
            color = COLORS[labels[0]] if labels else "#c9c9c9"
            axis.arrow(start, 0, delta, 0, width=.28, head_width=.40, head_length=min(.14, length*.3),
                       length_includes_head=True, facecolor=color, edgecolor="#555555", linewidth=.5)
            if labels:
                label = "+".join(labels) + ("*" if gene["origin"] != "annotated_CDS" else "")
                axis.text(x+length/2, .4, label, ha="center", fontsize=9)
        axis.set_xlim((low-int(row["core_start"]))/1000, (high-int(row["core_start"]))/1000)
        axis.set_ylim(-.55,.85)
        axis.set_yticks([])
        axis.set_title(f"{row['organism']} | {row['genome_id']} | {row['family']} | C identity {row['C_reference_identity_percent']:.1f}%", loc="left", fontsize=10)
        axis.set_xlabel("kb relative to the left edge of the defined core")
    fig.legend(handles=[Patch(facecolor=color, label=role) for role,color in COLORS.items()], loc="upper right", ncol=4)
    fig.text(.025,.015,"A/B1/B2/C are evidence-based role hypotheses. Grey: other annotated CDS. *: addressed local precursor ORF.\n"
             "Contexts are plotted in assembly orientation; examples minimize C identity within each source/family.",fontsize=9)
    fig.tight_layout(rect=(0,.09,1,.96),h_pad=2)
    for suffix in ("png", "svg"):
        fig.savefig(out / f"diversity_species_contexts.{suffix}",dpi=170)


if __name__ == "__main__":
    main()
