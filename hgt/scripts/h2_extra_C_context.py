#!/usr/bin/env python3
"""Check four unrepresented C homologues against their sampled locus neighbours."""
import gzip
import json
import os
from pathlib import Path
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from prepare_h2 import ROOT, H1, OUT, read, write

directory = OUT / 'orthology_diagnostic' / 'context_check'
directory.mkdir(exist_ok=True)
extras = {r['native_protein_id']: r for r in read(OUT / 'unrepresented_C_reference_homologues.tsv')}
leaves = {r['tip_id']: r for r in read(OUT / 'component_leaf_map.tsv') if r['family'] == 'noncanonical_C'}
hosts = {r['host_id']: r for r in read(OUT / 'host_inputs.tsv')}
genes = read(H1 / 'H1_publication_genes.tsv')
coding = {r.id: str(r.seq) for r in SeqIO.parse(H1 / 'coding_sequences.fna', 'fasta')}
proteins = {r.id: str(r.seq) for r in SeqIO.parse(H1 / 'proteins.faa', 'fasta')}
source_root = Path(os.environ['HGT_GENOME_ROOT'])
summary, precursor_hits, windows = [], [], []
for placement in read(OUT / 'orthology_diagnostic/extra_homologue_placement.tsv'):
    extra = extras[placement['extra_tip']]
    target = leaves[placement['nearest_sampled_target']]
    paths, regions = [], []
    for name, locus in [('extra', extra), ('target', target)]:
        host = locus['host_id']
        seq = next(r.seq for r in SeqIO.parse(OUT / 'native' / host / 'genome.fna', 'fasta') if r.id == locus['contig'])
        start, end = int(locus['start_0based']), int(locus['end_0based'])
        left, right = max(0, start - 5000), min(len(seq), end + 5000)
        region = seq[left:right]
        regions.append(region)
        path = directory / (placement['extra_tip'] + '_' + name + '.fna')
        SeqIO.write([SeqRecord(region, id=host + '_' + name, description='')], path, 'fasta')
        paths.append(path)
        windows.append({'extra_tip': placement['extra_tip'], 'role': name, 'host_id': host, 'contig': locus['contig'],
                        'window_start_0based': left, 'window_end_0based': right, 'left_flank_bp': start - left, 'right_flank_bp': right - end})
    blast_path = directory / (placement['extra_tip'] + '_BLASTN.tsv')
    subprocess.run([str(Path(os.environ['HGT_BIN']) / 'blastn'), '-task', 'blastn', '-query', str(paths[0]),
                    '-subject', str(paths[1]), '-evalue', '1e-10', '-outfmt',
                    '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen', '-out', str(blast_path)], check=True)
    hsps = [line.split('\t') for line in blast_path.read_text().splitlines()]
    longest = max(hsps, key=lambda h: int(h[3])) if hsps else None
    matches = []
    locus_ids = set(target['publication_locus_ids'].split(';'))
    for gene in genes:
        if gene['locus_id'] not in locus_ids or gene['author_role'] not in ['A', 'A1']:
            continue
        dna = coding[gene['gene_id']]
        for strand, query in [(1, dna), (-1, str(Seq(dna).reverse_complement()))]:
            position = str(regions[0]).find(query)
            if position >= 0:
                matches.append(gene['gene_id'])
                precursor_hits.append({'extra_tip': placement['extra_tip'], 'target_precursor': gene['gene_id'],
                    'window_offset_0based': position, 'strand': strand, 'coding_DNA': dna, 'published_AA': proteins[gene['gene_id']],
                    'evidence': 'exact_published_precursor_coding_DNA_in_extra_C_window; product_unconfirmed'})
    host = extra['host_id']
    native_aa = next(str(r.seq) for r in SeqIO.parse(OUT / 'native' / host / 'proteins.faa', 'fasta') if r.id == extra['native_protein_id'])
    path = source_root / hosts[host]['path']
    opening = gzip.open if path.suffix == '.gz' else open
    translation = 'native_feature_not_found'
    with opening(path, 'rt') as handle:
        for record in SeqIO.parse(handle, 'genbank'):
            if record.id != extra['contig']:
                continue
            feature = next(f for f in record.features if f.type == 'CDS' and int(f.location.start) == int(extra['start_0based'])
                           and int(f.location.end) == int(extra['end_0based']) and f.qualifiers.get('translation'))
            offset = int(feature.qualifiers.get('codon_start', ['1'])[0]) - 1
            dna = feature.extract(record.seq)[offset:]
            aa = str(dna.translate(table=int(feature.qualifiers.get('transl_table', ['11'])[0]))).rstrip('*')
            translation = 'exact_translation' if aa == native_aa else 'exact_with_bacterial_initiator_M' if aa[1:] == native_aa[1:] and native_aa.startswith('M') else 'translation_mismatch'
            break
    summary.append({'extra_tip': placement['extra_tip'], 'nearest_target_tip': placement['nearest_sampled_target'],
        'extra_assembly': hosts[host]['sequence_provider_assembly'], 'extra_BioSample': hosts[host]['biosample'],
        'target_assembly': hosts[target['host_id']]['sequence_provider_assembly'], 'target_BioSample': hosts[target['host_id']]['biosample'],
        'translation_check': translation, 'partial_annotation': extra['partial_annotation'], 'longest_DNA_HSP_bp': longest[3] if longest else 0,
        'longest_HSP_identity_percent': longest[2] if longest else '', 'HSP_query_coverage': int(longest[3]) / int(longest[12]) if longest else 0,
        'exact_target_precursor_DNA_matches': ';'.join(sorted(set(matches))),
        'limit': 'local_similarity_and_precursor_alleles_do_not_establish_activity_or_independent_isolate_confirmation; no_primary_input_change'})
write(directory / 'window_inventory.tsv', windows)
write(directory / 'summary.tsv', summary)
if precursor_hits:
    write(directory / 'exact_precursor_alleles.tsv', precursor_hits)
print(json.dumps({'native_C_checked': len(summary), 'exact_precursor_hits': len(precursor_hits)}))
