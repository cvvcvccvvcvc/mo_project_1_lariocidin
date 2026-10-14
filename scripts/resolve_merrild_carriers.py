#!/usr/bin/env python3
"""Resolve remaining published carrier metadata; do not extend the frozen search."""
import csv
import json
import re
from collections import defaultdict

import pandas as pd

from collect_diversity_catalogue import wgs_project
from fetch_hgt_assemblies import Cache, ROOT, assembly_summaries, result_uids


def main():
    table = ROOT / 'results/tables/diversity_merrild_coverage.tsv'
    data = pd.read_csv(table, sep='\t', keep_default_na=False)
    unresolved = data[data.initial_selection_status == 'carrier_not_resolved']
    cache = Cache(ROOT / 'data/raw/diversity_extension_2026-10-13',
                  ROOT / 'config/diversity_extension_download_manifest.json')
    resolved = defaultdict(set)
    # Assembly search can resolve GenBank-only samples omitted from RefSeq.
    samples = sorted(set(unresolved.biosample_from_export) - {''})
    for start in range(0, len(samples), 20):
        group = samples[start:start+20]
        result = cache.get_json(f'merrild_unresolved_samples_{start}', 'esearch.fcgi',
                 {'db': 'assembly', 'term': ' OR '.join(f'{x}[All Fields]' for x in group),
                  'retmode': 'json', 'retmax': 1000})
        if int(result['esearchresult']['count']) > 1000:
            raise ValueError('Unexpected unbounded published-sample mapping')
        for row in assembly_summaries(cache, f'merrild_unresolved_samples_{start}', result_uids(result)):
            for sample in group:
                if row.get('biosampleaccn') == sample:
                    resolved[sample].add(row['assemblyaccession'])
    accessions = sorted(set(x.removeprefix('NZ_') for x in unresolved.accession
                           if re.fullmatch(r'(?:NZ_)?[A-Z]{2,6}\d{5,}(?:\.\d+)?', x)))
    for accession in accessions:
        result = cache.get_json(f'merrild_nucleotide_{accession}', 'esearch.fcgi',
                 {'db': 'nuccore', 'term': accession+'[Accession]', 'retmode': 'json', 'retmax': 2})
        ids = result_uids(result)
        if not ids:
            continue
        linked = cache.get_json(f'merrild_nucleotide_{accession}_assembly', 'elink.fcgi',
                 {'dbfrom': 'nuccore', 'db': 'assembly', 'id': ','.join(ids), 'retmode': 'json'})
        assembly_ids = [str(x) for linkset in linked.get('linksets', [])
                        for group in linkset.get('linksetdbs', []) for x in group.get('links', [])]
        for row in assembly_summaries(cache, f'merrild_nucleotide_{accession}', assembly_ids):
            resolved[accession].add(row['assemblyaccession'])
    wanted_projects = {wgs_project(x) for x in unresolved.accession} - {''}
    with (ROOT/'data/raw/diversity_2026-10-12/refseq_bacterial_assembly_summary.txt').open() as handle:
        for line in handle:
            if line.lstrip('# ').startswith('assembly_accession\t'):
                fields = line.lstrip('# ').rstrip().split('\t')
                break
        projects = defaultdict(set)
        for values in csv.reader(handle, delimiter='\t'):
            row = dict(zip(fields, values))
            project = wgs_project(row.get('wgs_master', ''))
            if project in wanted_projects:
                projects[project].add(row['assembly_accession'])
        for accession in unresolved.accession:
            resolved[accession.removeprefix('NZ_')].update(projects[wgs_project(accession)])
    data['additional_primary_assembly_mapping'] = [
        ';'.join(sorted(resolved.get(r['biosample_from_export'], set()) |
                        resolved.get(r['accession'].removeprefix('NZ_'), set())))
        for r in data.to_dict('records')]
    data['additional_mapping_limit'] = data.additional_primary_assembly_mapping.map(
        lambda x: 'BioSample_can_contain_multiple_MAGs_or_versions_verify_exported_locus' if x else '')
    data.to_csv(table, sep='\t', index=False)
    print(json.dumps({'published_records_with_additional_mapping': int((data.additional_primary_assembly_mapping != '').sum()),
                      'assemblies': len(set(x for cell in data.additional_primary_assembly_mapping for x in cell.split(';') if x))}))


if __name__ == '__main__':
    main()
