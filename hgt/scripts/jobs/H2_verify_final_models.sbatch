#!/usr/bin/env bash
#SBATCH --job-name=lariocidin-HGT-model-input-check
#SBATCH --partition=main
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
set -euo pipefail
: "${HGT_ROOT:?Isolated H2 root required}"
export PYTHONDONTWRITEBYTECODE=1
"$HGT_BIN/python" - <<'PY'
import csv
import hashlib
import json
import os
from pathlib import Path
import re
from Bio import SeqIO

root = Path(os.environ['HGT_ROOT'])
def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))
def sequences(path):
    return {record.id: str(record.seq) for record in SeqIO.parse(path, 'fasta')}
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

checked = []
for scope in ('results/H2', 'results/H2/actino'):
    base = root / scope
    models = rows(base / 'host_marker_model_choices.tsv')
    markers = []
    for row in models:
        marker = row['marker_id']
        alignment = base / 'marker_trees' / marker / 'trimmed_alignment.faa'
        report = base / 'marker_trees' / marker / 'ML_fastdiagnostic.iqtree'
        sequence = sequences(alignment)
        report_text = report.read_text()
        selected = re.search(r'Best-fit model according to BIC:\s*(\S+)', report_text).group(1)
        assert selected == row['model'], (scope, marker, selected, row['model'])
        length = len(next(iter(sequence.values())))
        assert all(len(s) == length for s in sequence.values())
        input_name = re.search(r'Input file name:\s*(.+)', report_text).group(1).strip()
        dimensions = re.search(r'Input data: (\d+) sequences with (\d+) amino-acid sites', report_text)
        assert input_name == str(alignment), (scope, marker, 'report input path')
        assert tuple(map(int, dimensions.groups())) == (len(sequence), length), (scope, marker, 'report dimensions')
        assert alignment.stat().st_mtime <= report.stat().st_mtime, (scope, marker, 'alignment rewritten after fit')
        markers.append((marker, selected, length, sequence))
        checked.append({'scope': scope, 'marker': marker, 'model': selected,
                        'alignment_sha256': sha(alignment), 'model_report_sha256': sha(report),
                        'taxa_for_model_selection': len(sequence), 'columns': length,
                        'report_input_path_and_dimensions_verified': True, 'alignment_not_rewritten_after_report': True})
    for job in rows(base / 'tree_jobs.tsv'):
        if job['kind'] != 'host' or 'fixed_marker' not in job['partition_path']:
            continue
        matrix_path, partition_path = root / job['input_path'], root / job['partition_path']
        matrix, partition = sequences(matrix_path), partition_path.read_text()
        offset = 0
        for marker, selected, length, sequence in markers:
            assert f'charset {marker} = {offset + 1}-{offset + length};' in partition
            assert f'{selected}:{marker}' in partition
            assert all(s[offset:offset + length] == sequence.get(host, '-' * length)
                       for host, s in matrix.items()), (scope, job['name'], marker)
            offset += length
        assert all(len(s) == offset for s in matrix.values())
        assert len(re.findall(r'\bcharset\b', partition)) == len(markers)
        checked.append({'scope': scope, 'host_dataset': job['name'], 'taxa': len(matrix),
                        'partitions': len(markers), 'columns': offset,
                        'matrix_sha256': sha(matrix_path), 'partition_sha256': sha(partition_path)})
result = {'status': 'all_selected_models_and_final_alignment_blocks_match',
          'model_taxon_scope': 'full frozen marker background; final host subsets retain exactly the same columns',
          'parameter_scope': 'parameters reestimated in each final host ML fit; subset BIC optimality not asserted',
          'records': checked}
(root / 'config/H2/H2_final_alignment_model_correspondence.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'status': result['status'], 'markers': sum('marker' in r for r in checked),
                  'host_matrices': sum('host_dataset' in r for r in checked)}))
PY
