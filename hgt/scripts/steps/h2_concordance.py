#!/usr/bin/env python3
"""Standard IQ-TREE gCF on independent marker trees; no HGT score."""
import argparse
import json
import os
import subprocess
from pathlib import Path

from Bio import Phylo

from prepare_h2 import ROOT, OUT, read

parser = argparse.ArgumentParser()
parser.add_argument("--actino", action="store_true")
parser.add_argument("--dataset", help="Only the named completed host dataset")
args = parser.parse_args()
destination = OUT / "actino" if args.actino else OUT
models = read(destination / "marker_models.tsv")
for job in read(destination / "tree_jobs.tsv"):
    if job["kind"] != "host":
        continue
    if args.dataset and job['name'] != args.dataset:
        continue
    directory = destination / "host_trees" / job["name"]
    reference = directory / "ML.treefile"
    tree = Phylo.read(reference, "newick")
    taxa = {x.name for x in tree.get_terminals()}
    markers, skipped = [], []
    for model in models:
        source = destination / "marker_trees" / model["marker_id"]
        status = json.loads((source / "task_summary.json").read_text())
        if not (source / "ML_fastdiagnostic.treefile").exists():
            skipped.append({"marker": model["marker_id"], "reason": status["tree_status"]})
            continue
        gene = Phylo.read(source / "ML_fastdiagnostic.treefile", "newick")
        for leaf in list(gene.get_terminals()):
            if leaf.name not in taxa:
                gene.prune(leaf)
        if len(gene.get_terminals()) >= 4:
            # Pruning can leave a degree-two display root, which IQ-TREE
            # interprets as biological rooting and a virtual '__root__' taxon.
            if len(gene.root.clades) == 2:
                gene.root_with_outgroup(next(c for c in gene.root.clades if not c.is_terminal()))
            gene.rooted = False
            assert len(gene.root.clades) >= 3
            markers.append(gene)
    gene_path = directory / "independent_marker_trees.nwk"
    # gCF uses topology only. An inherited positive basal branch makes
    # IQ-TREE infer rooting even with a degree-three display root.
    Phylo.write(markers, gene_path, "newick", plain=True)
    command = [Path(os.environ["HGT_BIN"]) / "iqtree", "-t", reference, "--gcf", gene_path, "-nt", "1", "-seed", "131013", "--prefix", directory / "concordance"]
    with (directory / "concordance.stdout.log").open("w") as handle:
        subprocess.run([str(x) for x in command], stdout=handle, stderr=subprocess.STDOUT, check=True)
    (directory / "concordance_inputs.json").write_text(json.dumps({"gene_trees": len(markers), "skipped": skipped,
        "serialization": "topology-only unrooted degree>=3 display root after pruning; Biopython internal-node reorientation preserves unrooted topology; original lengths/support retained in source marker trees",
        "limitation": "IQ-TREE -fast gene trees; gCF describes topology agreement conditional on taxon coverage, not bootstrap, p-value or HGT probability"}, indent=2) + "\n")
print(json.dumps({"scope": "Actinobacteria138" if args.actino else "Bacteria74", "status": "gCF_completed"}))
