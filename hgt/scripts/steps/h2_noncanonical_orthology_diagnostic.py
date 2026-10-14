#!/usr/bin/env python3
"""Small, separate C-family diagnostic for four unrepresented homologues."""
import json
import os
from pathlib import Path
import re
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from prepare_h2 import OUT, read, write
from h2_components import domains

directory = OUT / 'orthology_diagnostic'
directory.mkdir(exist_ok=True)
extra = [r for r in read(OUT / 'unrepresented_C_reference_homologues.tsv') if r['reference_family'] == 'noncanonical_C']
records = list(SeqIO.parse(OUT / 'components/noncanonical_C.faa', 'fasta'))
context = []
for candidate in extra:
    native = OUT / 'native' / candidate['host_id']
    domain_hits = domains(native / 'domains.domtbl')
    record = next(r for r in SeqIO.parse(native / 'proteins.faa', 'fasta') if r.id == candidate['native_protein_id'])
    records.append(SeqRecord(record.seq, id=record.id, description=''))
    for gene in read(native / 'proteins.tsv'):
        if gene['contig'] == candidate['contig'] and max(0, int(gene['start_0based']) - int(candidate['end_0based']),
                                                       int(candidate['start_0based']) - int(gene['end_0based'])) <= 5000:
            context.append({'candidate_protein': candidate['native_protein_id'], **gene,
                            'domain_accessions': ';'.join(sorted({r['accession'] for r in domain_hits.get(gene['protein_id'], [])}))})
input_path = directory / 'noncanonical_C_plus_unrepresented.faa'
SeqIO.write(records, input_path, 'fasta')
write(directory / 'four_candidate_neighborhoods.tsv', context)
bins = Path(os.environ['HGT_BIN'])
small = Path(os.environ['HGT_SMALL_BIN'])
alignment, trimmed = directory / 'raw_alignment.faa', directory / 'trimmed_alignment.faa'
with alignment.open('w') as output, (directory / 'mafft.log').open('w') as error:
    subprocess.run([str(bins / 'mafft'), '--auto', '--thread', '1', str(input_path)], stdout=output, stderr=error, check=True)
with (directory / 'trimal.log').open('w') as output:
    subprocess.run([str(small / 'trimal'), '-in', str(alignment), '-out', str(trimmed), '-gappyout'], stdout=output, check=True)
source = OUT / 'component_trees/noncanonical_C/ML.iqtree'
model = re.search(r'Best-fit model according to BIC:\s*(\S+)', source.read_text()).group(1)
command = [str(bins / 'iqtree'), '-s', str(trimmed), '-nt', '1', '-seed', '141014', '-m', model,
           '-fast', '-alrt', '1000', '--keep-ident', '--prefix', str(directory / 'ML')]
with (directory / 'iqtree.log').open('w') as output:
    subprocess.run(command, stdout=output, stderr=subprocess.STDOUT, check=True)
result = {'tips': len(records), 'extra_native_homologues': len(extra), 'model': model,
          'model_basis': 'previously selected primary noncanonical_C class; parameters reestimated; no new model comparison',
          'command': command, 'status': 'separate_approximate_ML_orthology_diagnostic; 1000SHaLRT_only; not_HGT_event_tree',
          'limits': 'unrepresented homologue function unconfirmed; one partial contig-edge gene; no change to primary trees'}
(directory / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
