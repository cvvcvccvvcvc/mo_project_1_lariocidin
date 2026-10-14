#!/usr/bin/env python3
"""Standard pairwise DNA comparison of the two observed multicopy C hosts."""
from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from prepare_h2 import OUT, read, write

copies = defaultdict(list)
for row in read(OUT / 'component_leaf_map.tsv'):
    if row['family'] == 'noncanonical_C':
        copies[row['host_id']].append(row)
directory = OUT / 'copy_context'
directory.mkdir(exist_ok=True)
inventory = []
for host, genes in sorted(copies.items()):
    if len(genes) <= 1:
        continue
    assert len(genes) == 2, 'Review any new copy-number pattern before extending this pilot'
    genome = {r.id: r.seq for r in SeqIO.parse(OUT / 'native' / host / 'genome.fna', 'fasta')}
    paths = []
    for gene in genes:
        contig = genome[gene['contig']]
        start, end = int(gene['start_0based']), int(gene['end_0based'])
        left, right = max(0, start - 5000), min(len(contig), end + 5000)
        region = contig[left:right]
        if gene['strand'] == '-1':
            region = region.reverse_complement()
        path = directory / (gene['tip_id'] + '.fna')
        SeqIO.write([SeqRecord(region, id=gene['tip_id'], description='')], path, 'fasta')
        paths.append(path)
        canonical = min(str(contig), str(contig.reverse_complement()))
        inventory.append({'host_id': host, 'tip_id': gene['tip_id'], 'contig': gene['contig'],
            'contig_bp': len(contig), 'contig_sha256_strand_independent': hashlib.sha256(canonical.encode()).hexdigest(),
            'C_start_0based': start, 'C_end_0based': end, 'strand': gene['strand'],
            'window_start_0based': left, 'window_end_0based': right,
            'left_flank_available_bp': start - left, 'right_flank_available_bp': right - end,
            'limit': 'contig and flank similarity cannot distinguish biological duplication from assembly repetition without reads'})
    subprocess.run([str(Path(os.environ['HGT_BIN']) / 'blastn'), '-task', 'blastn', '-query', str(paths[0]),
                    '-subject', str(paths[1]), '-evalue', '1e-10', '-outfmt',
                    '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
                    '-out', str(directory / (host + '_pairwise_BLASTN.tsv'))], check=True)
write(directory / 'copy_window_inventory.tsv', inventory)
print(json.dumps({'multicopy_hosts': len(inventory) // 2, 'DNA_windows': len(inventory),
                  'method': 'BLASTN on C plus up to5kb of each native flank; all HSPs retained'}))
