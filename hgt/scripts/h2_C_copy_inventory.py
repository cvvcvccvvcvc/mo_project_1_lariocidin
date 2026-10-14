#!/usr/bin/env python3
"""Inventory unrepresented full-length C-reference homologues without classifying them."""
from collections import defaultdict
import json
from prepare_h2 import OUT, read, write

refs = {r['gene_id']: r['family'] for r in read(OUT / 'role_references.tsv') if r['family'].endswith('_C')}
represented = defaultdict(list)
for r in read(OUT / 'component_leaf_map.tsv'):
    if r['family'].endswith('_C'):
        represented[r['host_id']].append(r)
rows = []
for host in read(OUT / 'host_inputs.tsv'):
    directory = OUT / 'native' / host['host_id']
    genes = {r['protein_id']: r for r in read(directory / 'proteins.tsv')}
    best = {}
    with (directory / 'role_homology.tsv').open() as handle:
        for line in handle:
            f = line.split()
            if f[1] not in refs:
                continue
            qcov, scov = (int(f[5]) - int(f[4]) + 1) / int(f[10]), (int(f[7]) - int(f[6]) + 1) / int(f[11])
            if min(qcov, scov) < 0.7:
                continue
            gene = genes[f[0]]
            if any(r['contig'] == gene['contig'] and int(r['start_0based']) < int(gene['end_0based'])
                   and int(r['end_0based']) > int(gene['start_0based']) for r in represented[host['host_id']]):
                continue
            key = f[0], refs[f[1]]
            if key not in best or float(f[9]) > float(best[key]['bitscore']):
                best[key] = {'host_id': host['host_id'], 'native_protein_id': f[0], 'reference_family': refs[f[1]],
                    'reference_gene': f[1], 'identity_percent': f[2], 'query_coverage': qcov, 'reference_coverage': scov,
                    'E_value': f[8], 'bitscore': f[9], 'contig': gene['contig'], 'start_0based': gene['start_0based'],
                    'end_0based': gene['end_0based'], 'strand': gene['strand'], 'product': gene['product'],
                    'partial_annotation': gene['partial_location'], 'overlaps_mapped_window': gene['overlaps_mapped_export_window'],
                    'status': 'unrepresented_full_length_reference_homologue; metabolic_or_other_RiPP_role_unresolved; not_added_to_primary_tree'}
    rows.extend(best.values())
if rows:
    write(OUT / 'unrepresented_C_reference_homologues.tsv', rows)
print(json.dumps({'unrepresented_homologue_rows': len(rows), 'hosts': len({r['host_id'] for r in rows}),
                  'limit': 'one represented target copy does not establish genomic one-to-one orthology'}))
