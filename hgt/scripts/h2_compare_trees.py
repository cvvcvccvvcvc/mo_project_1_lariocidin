#!/usr/bin/env python3
"""Descriptive unrooted split comparisons on shared single-copy hosts."""
from collections import Counter, defaultdict
import argparse
from itertools import combinations
import hashlib
import json
import re

from Bio import Phylo, SeqIO

from prepare_h2 import ROOT, OUT, read, write

# All H2 alignments have IQ-TREE's default 1e-6 branch-length floor.
# Include its tiny numerical neighbourhood in the unresolved-edge mask.
MIN_RESOLVED_EDGE_LENGTH = 1.01e-6


def branch_support(clade):
    parts = str(clade.name or "").split("/")
    if len(parts) == 2:
        return tuple(float(x) for x in parts)
    return (0.0, 0.0)


def splits(tree, shared, gene=False):
    result = {}
    for node in tree.find_clades():
        if node.is_terminal() or (node.branch_length or 0) <= MIN_RESOLVED_EDGE_LENGTH:
            continue
        side = frozenset((x.name.split("_", 1)[0] if gene else x.name) for x in node.get_terminals()) & shared
        other = shared - side
        if min(len(side), len(other)) < 2:
            continue
        key = min((side, other), key=lambda x: (len(x), tuple(sorted(x))))
        support = branch_support(node)
        # If projection contracts multiple original edges, keep conservative minima.
        if key in result:
            support = tuple(min(a, b) for a, b in zip(support, result[key]))
        result[key] = support
    return result


def compatible(a, b, total):
    return any(not x for x in [a & b, a & (total - b), (total - a) & b, (total - a) & (total - b)])


def identifier(side):
    return hashlib.sha256(";".join(sorted(side)).encode()).hexdigest()[:12]


def exact_group_edge(tree, group):
    """Descriptive unrooted separation, not a functional or HGT classifier."""
    total = {tip.name for tip in tree.get_terminals()}
    supports = []
    for node in tree.get_nonterminals():
        if (node.branch_length or 0) <= MIN_RESOLVED_EDGE_LENGTH:
            continue
        side = {tip.name for tip in node.get_terminals()}
        if side == group or total - side == group:
            supports.append(float(node.confidence or 0))
    return {'tips_in_group': len(group), 'other_tips': len(total - group),
            'exact_separating_edge_present': bool(supports), 'SH_aLRT': min(supports) if supports else None,
            'limit': 'approximate_ML_orthology_diagnostic; separation_does_not_establish_function_or_one_to_one_orthology'}


def compare(first, second, shared, gene_first=True, gene_second=False):
    a = splits(first, shared, gene_first)
    b = splits(second, shared, gene_second)
    strong_a = {s: u for s, u in a.items() if u[0] >= 80 and u[1] >= 95}
    strong_b = {s: u for s, u in b.items() if u[0] >= 80 and u[1] >= 95}
    conflicts = [(x, y) for x in strong_a for y in strong_b if not compatible(x, y, shared)]
    return {"shared_single_copy_hosts": len(shared), "first_nonzero_splits": len(a), "second_nonzero_splits": len(b),
            "common_splits": len(a.keys() & b.keys()), "RF_nonzero_splits": len(a.keys() ^ b.keys()),
            "strong_first_splits": len(strong_a), "strong_second_splits": len(strong_b),
            "strong_conflicting_split_pairs": len(conflicts)}, [(x, y, strong_a[x], strong_b[y]) for x, y in conflicts]


def main(ready_only=False):
    qc = {r["host_id"]: r for r in read(OUT / "marker_genome_QC.tsv")}
    tree_quality, comparisons, conflicts, family_pairs, family_conflicts, sensitivity, relatives = [], [], [], [], [], [], []
    for scope in [OUT, OUT / "actino"]:
        scope_name = "Actinobacteria138" if scope.name == "actino" else "Bacteria74"
        jobs = read(scope / "tree_jobs.tsv")
        host_trees, genes, c_trees = {}, {}, {}
        for job in jobs:
            directory = scope / ("host_trees" if job["kind"] == "host" else "component_trees") / job["name"]
            if ready_only and not (directory / "task_summary.json").exists():
                continue
            stats = json.loads((directory / "task_summary.json").read_text())
            report = (directory / "ML.iqtree").read_text() if (directory / "ML.iqtree").exists() else ""
            log = (directory / "iqtree.log").read_text() if (directory / "iqtree.log").exists() else ""
            model = re.search(r"Best-fit model according to BIC:\s*(.+)", report)
            composition = re.search(r"(\d+) sequences failed composition chi2 test", log)
            tree_quality.append({"scope": scope_name, "kind": job["kind"], "dataset": job["name"], "dataset_role": job["dataset_role"],
                **stats["filtered"], "raw_columns": stats.get("raw", {}).get("columns", "pretrimmed_partitions"),
                "retained_fraction": stats.get("retained_fraction", "pretrimmed_partitions"),
                "best_model": model.group(1) if model else "partitioned_or_no_tree", "composition_chi2_failures": int(composition.group(1)) if composition else 0,
                "tree_status": stats["tree_status"]})
            if not (directory / "ML.treefile").exists():
                continue
            tree = Phylo.read(directory / "ML.treefile", "newick")
            expected_tips = {r.id for r in SeqIO.parse(ROOT / job['input_path'], 'fasta')}
            observed_tips = [r.name for r in tree.get_terminals()]
            assert len(observed_tips) == len(expected_tips) == stats['filtered']['tips']
            assert set(observed_tips) == expected_tips, (scope_name, job['name'], 'tree/input tip mismatch')
            if stats['tree_status'].startswith('ML_complete'):
                with (directory / 'ML.ufboot').open() as handle:
                    assert sum(bool(line.strip()) for line in handle) == 1000, (scope_name, job['name'], 'bootstrap count')
            if job['kind'] == 'component' and (job['name'].endswith('_C') or job['name'].endswith('_C_MAG_diagnostic')):
                c_trees[job['name']] = tree
            if job["kind"] == "host":
                host_trees[job["name"]] = tree
            elif (job["name"].startswith(("canonical_", "noncanonical_")) or job["name"] == "uncertain_Nature_system_C") and not job["name"].endswith(("_diagnostic", "_noMAG")):
                copies = Counter(x.name.split("_", 1)[0] for x in tree.get_terminals())
                if job["name"] == "canonical_ABC_family":
                    comparisons.append({"scope": scope_name, "family": job["name"], "multicopy_hosts_excluded": sum(n > 1 for n in copies.values()),
                        "shared_single_copy_hosts": "not_an_orthology_criterion", "status": "D1_D2_paralogous_reference_family; retain_all_copy_tree; no_one_to_one_host_or_component_RF"})
                    continue
                genes[job["name"]] = (tree, {h for h, n in copies.items() if n == 1}, sum(n > 1 for n in copies.values()))
        for family, (gene, singles, multicopy) in genes.items():
            lineage = "uncertain_Nature" if family == "uncertain_Nature_system_C" else "noncanonical" if family.startswith("noncanonical_") else "canonical"
            if lineage not in host_trees:
                continue
            host = host_trees[lineage]
            shared = frozenset(singles & {x.name for x in host.get_terminals()})
            if len(shared) < 4:
                comparisons.append({"scope": scope_name, "family": family, "multicopy_hosts_excluded": multicopy,
                    "shared_single_copy_hosts": len(shared), "status": "too_few_shared_single_copy_hosts"})
                continue
            values, pairs = compare(gene, host, shared)
            comparisons.append({"scope": scope_name, "family": family, "multicopy_hosts_excluded": multicopy, **values,
                                "status": "descriptive_split_comparison; neither_pvalue_nor_HGT_event_test"})
            for x, y, gs, hs in pairs:
                conflicts.append({"scope": scope_name, "family": family, "gene_split": identifier(x), "host_split": identifier(y),
                    "shared_hosts": ";".join(sorted(shared)), "gene_side": ";".join(sorted(x)), "host_side": ";".join(sorted(y)),
                    "gene_side_genera": ";".join(sorted({qc[h]["genus"] for h in x})),
                    "host_side_genera": ";".join(sorted({qc[h]["genus"] for h in y})),
                    "gene_SH_aLRT": gs[0], "gene_UFBoot": gs[1], "host_SH_aLRT": hs[0], "host_UFBoot": hs[1],
                    "interpretation": "supported_topology_disagreement_requires_orthology_model_and_DTL_alternatives; not_an_event"})
        panels = read(scope / 'host_panels.tsv')
        for lineage in ['canonical', 'noncanonical']:
            if lineage not in host_trees:
                continue
            tree = host_trees[lineage]
            admitted = [r for r in panels if r['lineage'] == lineage and r['primary_admission'] == 'True']
            carriers = [r['host_id'] for r in admitted if r['role'] == 'observed_locus_carrier']
            for row in admitted:
                if row['role'] != 'taxonomic_relative':
                    continue
                host = row['host_id']
                distance, nearest = min((tree.distance(host, candidate), candidate) for candidate in carriers)
                relatives.append({'scope': scope_name, 'lineage': lineage, 'relative_host': host,
                    'closest_sampled_carrier': nearest, 'ML_patristic_distance': distance,
                    'relative_species_taxid': qc[host]['species_taxid'], 'carrier_species_taxid': qc[nearest]['species_taxid'],
                    'relative_genus': qc[host]['genus'], 'carrier_genus': qc[nearest]['genus'],
                    'limit': 'closest_within_frozen_study_panel; selection_preceded_BGC_conflict_analysis; neither_global_nearest_genome_nor_BGC_negative_status'})
        for (a, (ta, sa, ma)), (b, (tb, sb, mb)) in combinations(genes.items(), 2):
            if a.split("_", 1)[0] != b.split("_", 1)[0]:
                continue
            shared = frozenset(sa & sb)
            if len(shared) >= 4:
                values, pairs = compare(ta, tb, shared, True, True)
                family_pairs.append({"scope": scope_name, "first_family": a, "second_family": b, **values,
                    "interpretation": "comparison_of_single_copy_contextual_homologues; congruence_does_not_establish_whole_BGC_transfer"})
                for x, y, first_support, second_support in pairs:
                    family_conflicts.append({"scope": scope_name, "first_family": a, "second_family": b,
                        "shared_hosts": ";".join(sorted(shared)), "first_side": ";".join(sorted(x)), "second_side": ";".join(sorted(y)),
                        "first_SH_aLRT": first_support[0], "first_UFBoot": first_support[1],
                        "second_SH_aLRT": second_support[0], "second_UFBoot": second_support[1],
                        "interpretation": "supported_component_disagreement; no_event_or_function_inference"})
        for name, primary in host_trees.items():
            if name.endswith("_diagnostic"):
                original = name.removesuffix("_MAG_diagnostic")
                if original in host_trees:
                    other = host_trees[original]
                    shared = frozenset(x.name for x in other.get_terminals())
                    values, _ = compare(primary, other, shared, False, False)
                    sensitivity.append({"comparison": "MAG_addition_diagnostic", "dataset": original, **values})
        for name, diagnostic in c_trees.items():
            if name.endswith('_MAG_diagnostic'):
                original = name.removesuffix('_MAG_diagnostic')
                if original in c_trees:
                    other = c_trees[original]
                    shared = frozenset(x.name.split('_', 1)[0] for x in other.get_terminals())
                    values, _ = compare(diagnostic, other, shared, True, True)
                    sensitivity.append({'comparison': 'MAG_addition_C_diagnostic', 'dataset': original, **values})
    broad = Phylo.read(OUT / "host_trees/noncanonical/ML.treefile", "newick")
    actino_directory = OUT / "actino/host_trees/noncanonical"
    if not ready_only or (actino_directory / "task_summary.json").exists():
        actino = Phylo.read(actino_directory / "ML.treefile", "newick")
        shared = frozenset(x.name for x in broad.get_terminals()) & frozenset(x.name for x in actino.get_terminals())
        values, _ = compare(broad, actino, shared, False, False)
        sensitivity.append({"comparison": "74_versus_138_marker_host_trees_on_shared_taxa", "dataset": "noncanonical_Actinomycetota", **values})
    broad_c = Phylo.read(OUT / 'component_trees/noncanonical_C/ML.treefile', 'newick')
    actino_c = Phylo.read(OUT / 'actino/component_trees/noncanonical_C/ML.treefile', 'newick')
    counts = [Counter(x.name.split('_', 1)[0] for x in t.get_terminals()) for t in [broad_c, actino_c]]
    shared = frozenset(h for h, n in counts[0].items() if n == 1 and counts[1].get(h) == 1)
    values, _ = compare(broad_c, actino_c, shared, True, True)
    sensitivity.append({'comparison': 'C_tree_resampling_with_Actino138_admission_on_shared_single_copy_hosts',
                        'dataset': 'noncanonical_C', **values})
    destination = OUT / ("comparisons_partial" if ready_only else "comparisons")
    destination.mkdir(exist_ok=True)
    for name, rows in [("tree_quality", tree_quality), ("component_host", comparisons), ("supported_split_conflicts", conflicts),
                       ("component_pairs", family_pairs), ("component_pair_conflicts", family_conflicts),
                       ("host_sensitivity", sensitivity), ("relative_placement", relatives)]:
        if rows:
            # A skipped comparison may have fewer fields than a completed one.
            fields = list(dict.fromkeys(k for r in rows for k in r))
            write(destination / (name + ".tsv"), [{k: r.get(k, "") for k in fields} for r in rows])
    domain_tree = Phylo.read(OUT / 'component_trees/C_domain_background_diagnostic/ML.treefile', 'newick')
    aliases = read(OUT / 'C_domain_diagnostic_aliases.tsv')
    diagnostic = {}
    for family in ['canonical_C', 'noncanonical_C']:
        group = {r['diagnostic_tip'] for r in aliases if r['source'] == family}
        diagnostic['domain_' + family] = exact_group_edge(domain_tree, group)
    extra_tree = Phylo.read(OUT / 'orthology_diagnostic/ML.treefile', 'newick')
    group = {r.id for r in SeqIO.parse(OUT / 'components/noncanonical_C.faa', 'fasta')}
    diagnostic['noncanonical_targets_versus_four_extra_full_length_homologues'] = exact_group_edge(extra_tree, group)
    (OUT / 'orthology_diagnostic/topology_summary.json').write_text(json.dumps(diagnostic, indent=2) + '\n')
    print(json.dumps({"tree_datasets": len(tree_quality), "component_host_comparisons": len(comparisons),
                      "supported_conflicting_split_pairs": len(conflicts), "status": "descriptive_H2; no_HGT_event_called"}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ready-only", action="store_true", help="Compare completed fits only, writing separately from final H2 results")
    main(parser.parse_args().ready_only)
