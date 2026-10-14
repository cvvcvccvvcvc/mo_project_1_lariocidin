#!/usr/bin/env python3
"""Retrieve evidence for manually selected H2 relations, without calling events."""
import argparse
import csv
import json

from Bio import Phylo

from prepare_h2 import ROOT, OUT, read, write
from h2_compare_trees import MIN_RESOLVED_EDGE_LENGTH, branch_support


def edges(tree, group, shared, gene=False):
    for node in tree.get_nonterminals():
        if (node.branch_length or 0) <= MIN_RESOLVED_EDGE_LENGTH:
            continue
        side = {t.name.split('_', 1)[0] if gene else t.name for t in node.get_terminals()} & shared
        if side == group or shared - side == group:
            yield node


def main(ready_only=False):
    directory = OUT / ('comparisons_partial' if ready_only else 'comparisons')
    conflicts = read(directory / 'supported_split_conflicts.tsv')
    qc = {r['host_id']: r for r in read(OUT / 'marker_genome_QC.tsv')}
    leaves = read(OUT / 'component_leaf_map.tsv')
    queries = json.loads((ROOT / 'config/H2/H2_review_relation_queries.json').read_text())['candidates']
    catalog, details = [], []
    for query in queries:
        matches = [r for r in conflicts if all(r[k] == query[k] for k in ['scope', 'family', 'gene_split', 'host_split'])]
        if not matches:
            details.append({**query, 'status': 'selected_supported_conflict_not_present_in_current_comparison'})
            continue
        row = matches[0]
        scope = OUT / 'actino' if row['scope'] == 'Actinobacteria138' else OUT
        host_name = 'uncertain_Nature' if row['family'] == 'uncertain_Nature_system_C' else row['family'].split('_')[0]
        host_directory = scope / 'host_trees' / host_name
        with (host_directory / 'concordance.cf.stat').open() as handle:
            gcf = {r['ID']: r for r in csv.DictReader((x for x in handle if not x.startswith('#')), delimiter='\t')}
        host_tree = Phylo.read(host_directory / 'concordance.cf.branch', 'newick')
        shared, gene_group, host_group = [set(row[k].split(';')) for k in ['shared_hosts', 'gene_side', 'host_side']]
        host_evidence = []
        for node in edges(host_tree, host_group, shared):
            branch_id = str(int(node.confidence))
            host_evidence.append({k: gcf[branch_id][k] for k in ['ID', 'Label', 'gCF', 'gCF_N', 'gN', 'gDF1', 'gDF2', 'gDFP', 'Length']})
        families = {}
        for job in read(scope / 'tree_jobs.tsv'):
            if job['kind'] != 'component' or not job['name'].startswith(row['family'].split('_')[0] + '_') or job['name'] == 'canonical_ABC_family':
                continue
            tree = Phylo.read(scope / 'component_trees' / job['name'] / 'ML.treefile', 'newick')
            host_ids = [t.name.split('_', 1)[0] for t in tree.get_terminals()]
            if not gene_group <= set(host_ids) or any(host_ids.count(h) != 1 for h in gene_group):
                families[job['name']] = {'status': 'missing_or_multiple_query_copies'}
                continue
            support = [branch_support(n) for n in edges(tree, gene_group, set(host_ids), gene=True)]
            families[job['name']] = {'exact_full_tree_group_edge': bool(support), 'SH_UF': support}
        carriers = sorted(gene_group | host_group)
        identity = [{'host_id': h, 'taxonomy': qc[h]['taxonomy_name'], 'assembly': qc[h]['sequence_provider_assembly'],
                     'BioSample': qc[h]['biosample'], 'assembly_level': qc[h]['assembly_level'],
                     'C_loci': [{k: r[k] for k in ['tip_id', 'contig', 'start_0based', 'end_0based', 'C_contig_edge_bp', 'publication_locus_ids']}
                                for r in leaves if r['host_id'] == h and r['family'] == row['family']]} for h in carriers]
        evidence = {**query, 'status': 'H2_relation_candidate; no_DTL_AU_or_HGT_event_call', 'comparison': row,
                    'host_gCF': host_evidence, 'context_families': families, 'carrier_and_locus_identity': identity}
        details.append(evidence)
        catalog.append({'candidate_id': query['candidate_id'], 'scope': row['scope'], 'component': row['family'],
            'gene_side_hosts': row['gene_side'], 'host_side_hosts': row['host_side'],
            'gene_SH_UF': row['gene_SH_aLRT'] + '/' + row['gene_UFBoot'], 'host_SH_UF': row['host_SH_aLRT'] + '/' + row['host_UFBoot'],
            'host_gCF_N_over_gN': ';'.join(r['gCF_N'] + '/' + r['gN'] for r in host_evidence),
            'host_gCF_percent': ';'.join(r['gCF'] for r in host_evidence), 'interpretation': query['interpretation'],
            'DTL_model_support': 'not_estimated', 'AU_p_value': 'not_estimated', 'direction_and_route': 'unresolved',
            'assembly_alternatives': 'reads_and_independent_assembly_validation_not_completed'})
    write(directory / 'H2_relationship_catalog.tsv', catalog)
    (directory / 'H2_relationship_evidence.json').write_text(json.dumps(details, indent=2) + '\n')

    diagnostic = OUT / 'orthology_diagnostic'
    tree = Phylo.read(diagnostic / 'ML.treefile', 'newick')
    targets = {r['tip_id'] for r in leaves if r['family'] == 'noncanonical_C'}
    extra = [t for t in tree.get_terminals() if t.name not in targets]
    placement = []
    for tip in extra:
        distance, nearest = min((tree.distance(tip, candidate), candidate) for candidate in targets)
        placement.append({'extra_tip': tip.name, 'nearest_sampled_target': nearest, 'ML_patristic_distance': distance,
            'extra_taxonomy': qc[tip.name.split('_')[0]]['taxonomy_name'], 'target_taxonomy': qc[nearest.split('_')[0]]['taxonomy_name'],
            'limit': 'approximate_ML_only; nearest_within_fixed_sample; not_function_or_orthology_assignment'})
    write(diagnostic / 'extra_homologue_placement.tsv', placement)
    print(json.dumps({'curated_relations_with_evidence': len(catalog), 'extra_homologues_placed': len(placement)}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ready-only', action='store_true')
    main(parser.parse_args().ready_only)
