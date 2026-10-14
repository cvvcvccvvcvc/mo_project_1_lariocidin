#!/usr/bin/env python3
"""Trace exact published C coding DNA through the fixed screening stages."""
import argparse
import json
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from extract_diversity_contexts import read_genome
from fetch_hgt_assemblies import ROOT
from screen_diversity_proteomes import fetch_batch


def read(path):
    return pd.read_csv(path, sep='\t', keep_default_na=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefixes', nargs='+', default=['diversity_species', 'diversity_extension', 'diversity_controls'])
    args = parser.parse_args()
    tables = ROOT/'results/tables'
    selected = pd.concat([read(path) for prefix in args.prefixes
                          for path in sorted((ROOT/'config').glob(prefix+'_[0-9]*_selection.tsv'))], ignore_index=True)
    if selected.genome_id.duplicated().any():
        raise ValueError('Repeated physical input IDs in coverage audit')
    selected = selected.set_index('genome_id')
    contexts = pd.concat([read(tables/f'{prefix}_all_contexts.tsv') for prefix in args.prefixes], ignore_index=True)
    coverage = read(tables/'diversity_merrild_coverage.tsv')
    coding = read(ROOT/'data/processed/diversity_merrild_C_coding.tsv')
    dna = {r.id: str(r.seq).upper() for r in SeqIO.parse(ROOT/'data/processed/diversity_merrild_C_coding.fna', 'fasta')}
    published_by_bgc = defaultdict(list)
    for row in coding.to_dict('records'):
        published_by_bgc[row['bgc_id']].append(row)
    by_genome = defaultdict(list)
    for row in coverage.to_dict('records'):
        ids = set(x for field in ['resolved_refseq_assemblies', 'additional_primary_assembly_mapping']
                  for x in row.get(field, '').split(';') if x)
        for genome in ids & set(selected.index):
            by_genome[genome].append(row)
    outputs = []
    # These are bounded published-carrier controls, including no-nomination cases.
    for genome, records in sorted(by_genome.items()):
        entry = selected.loc[genome]
        url = entry.get('context_url', '') or entry.url.replace('_protein.faa.gz', '_genomic.gbff.gz')
        request = {'genome_id': genome, 'source': entry.source, 'url': url,
                   'compression': 'gzip', 'file_suffix': '.gbff.gz'}
        download = fetch_batch([request], ROOT/'data/raw/diversity_contexts', workers=1)[0]
        with (ROOT/'config/diversity_merrild_context_downloads.jsonl').open('a') as ledger:
            ledger.write(json.dumps(download)+'\n')
        features, contigs = read_genome(ROOT/download['local_path'], entry.source)
        local = contexts[contexts.genome_id == genome]
        for published in records:
            refs = published_by_bgc[published['bgc_id']]
            if not refs:
                outputs.append({'bgc_id': published['bgc_id'], 'genome_id': genome,
                                'stage': 'export_has_no_explicit_C_to_trace', 'locus_id': '', 'reason': 'Unresolved author C annotation'})
                continue
            for ref in refs:
                sequence = dna[ref['feature_id']]
                positions = []
                for name, contig in contigs.items():
                    for strand, query in [(1, sequence), (-1, str(Seq(sequence).reverse_complement()))]:
                        start = str(contig.seq).upper().find(query)
                        while start >= 0:
                            positions.append((name, start, start+len(query), strand))
                            start = str(contig.seq).upper().find(query, start+1)
                if not positions:
                    outputs.append({'bgc_id': published['bgc_id'], 'published_C_feature': ref['feature_id'],
                                    'genome_id': genome, 'stage': 'published_C_DNA_not_exact_in_current_carrier',
                                    'locus_id': '', 'reason': 'Metadata carrier link cannot establish the exact exported locus in this assembly version'})
                for name, start, end, strand in positions:
                    annotated = [f for f in features if f['contig'] == name and f['strand'] == strand
                                 and max(0, min(end, f['end'])-max(start, f['start'])) >= .95*len(sequence)]
                    starts = pd.to_numeric(local.C_start, errors='coerce')
                    ends = pd.to_numeric(local.C_end, errors='coerce')
                    stop = ends if strand == 1 else starts
                    published_stop = end if strand == 1 else start
                    overlap = ends.clip(upper=end)-starts.clip(lower=start)
                    focal = local[(local.contig == name) & (pd.to_numeric(local.C_strand,errors='coerce') == strand)
                                  & ((stop-published_stop).abs() <= 3) & (overlap >= .7*len(sequence))]
                    stage = 'C_DNA_present_no_C_nomination' if focal.empty else 'context_supported_BGC' if (focal.status == 'family_bgc_candidate').any() else 'C_nominated_BGC_criteria_not_met'
                    outputs.append({'bgc_id': published['bgc_id'], 'published_C_feature': ref['feature_id'],
                                    'genome_id': genome, 'contig': name, 'C_start': start, 'C_end': end,
                                    'strand': strand, 'stage': stage, 'locus_id': ';'.join(focal.locus_id),
                                    'status': ';'.join(focal.status), 'reason': ';'.join(focal.interpretation_reason),
                                    'overlapping_annotation_protein_ids': ';'.join(f['source_protein_id'] for f in annotated),
                                    'context_source_sha256': download['sha256']})
    result = pd.DataFrame(outputs).fillna('')
    result.to_csv(tables/'diversity_merrild_screen_stages.tsv', sep='\t', index=False)
    print(json.dumps({'traced_rows': len(result), 'stages': result.stage.value_counts().to_dict()}))


if __name__ == '__main__':
    main()
