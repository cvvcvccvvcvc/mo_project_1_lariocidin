#!/usr/bin/env python3
"""Export a self-contained diversity snapshot without touching active HGT files."""
import argparse
import gzip
import hashlib
import json
import shutil
from datetime import datetime, timezone

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from extract_diversity_contexts import translated
from profile_search import ROOT


def read(path):
    return pd.read_csv(path, sep='\t', keep_default_na=False)


def sha(sequence):
    return hashlib.sha256(str(sequence).upper().encode()).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    args=parser.parse_args()
    output=ROOT/args.output
    output.mkdir(parents=True, exist_ok=False)
    tables,processed=ROOT/'results/tables',ROOT/'data/processed'
    catalogue=read(tables/'diversity_final_BGC_catalogue.tsv')
    occurrences=read(tables/'diversity_final_occurrences.tsv')
    merged=occurrences.set_index('locus_id').merged_locus_id.to_dict()
    metadata={};sequences={};dna_records={};gene_rows=[];core_records=[]
    fragments=read(tables/'diversity_C_fragment_audit.tsv')
    with gzip.open(output/'context_windows.fna.gz','wt') as windows_out:
        for prefix,n in [('diversity_species',8),('diversity_extension',8),('diversity_controls',1)]:
            contexts=read(tables/f'{prefix}_review_all_contexts.tsv').set_index('locus_id')
            for i in range(n):
                directory=processed/'discovery/contexts'/f'{prefix}_{i}_contexts'
                features=read(directory/'features.tsv')
                by_locus={locus:group for locus,group in features.groupby('locus_id')}
                proteins={r.id:r for file in ['neighbors.faa','precursor_queries.faa'] for r in SeqIO.parse(directory/file,'fasta')}
                for window in SeqIO.parse(directory/'windows.fna','fasta'):
                    row=contexts.loc[window.id].to_dict();metadata[window.id]={'locus_id':window.id,**row}
                    window.description=f'contig={row["contig"]} start_0based={row["window_start"]} end_exclusive={row["window_end"]}'
                    SeqIO.write(window,windows_out,'fasta')
                    assignments=dict(zip(row['chosen_core_gene_ids'].split(';'),row['architecture'].split('-')))
                    extra=fragments[fragments.locus_id==window.id]
                    for col in ['first_gene','second_gene']:
                        for query in extra[col]:assignments.setdefault(query,'C_fragment_partner')
                    selected=by_locus[window.id]
                    selected=selected[selected.query_id.isin(assignments)]
                    for f in selected.to_dict('records'):
                        query=f['query_id'];dna=window.seq[f['start']-row['window_start']:f['end']-row['window_start']]
                        if f['strand']==-1:dna=dna.reverse_complement()
                        protein=proteins[query];match=translated(dna)==str(protein.seq)
                        if row['status']=='family_bgc_candidate' and not match:raise ValueError(f'Primary gene DNA conflict: {query}')
                        sequences[query]=protein;dna_records[query]=SeqRecord(dna,id=query,description='coding_orientation_coordinates_in_genes_tsv')
                        gene_rows.append({'locus_id':window.id,'merged_locus_id':merged.get(window.id,''),
                            'genome_id':row['genome_id'],'assembly':row['assembly'],'biosample':row['biosample'],
                            'organism':row['organism'],'strain':row['strain'],'contig':row['contig'],
                            'locus_status':row['status'],'family':row['family'],'export_role':assignments[query],
                            **f,'coordinate_DNA_translation_match':match,'coding_DNA_sha256':sha(dna),
                            'status_basis':row['interpretation_reason'],'translation_caveat':'table11_coordinate_slice;partial_annotation_phase_requires_original_source_if_mismatch'})
                    if row['status']=='family_bgc_candidate':
                        core=window.seq[row['core_start']-row['window_start']:row['core_end']-row['window_start']]
                        if row['C_strand']==-1:core=core.reverse_complement()
                        assert sha(core)==row['core_sha256']
                        core_records.append(SeqRecord(core,id=window.id,description='defined_core_span_not_product_sequence'))
        history=read(ROOT/'inputs/selections/baseline/maize_curated_core_roles.tsv')
        history_coords=read(tables/'maize_candidate_neighbors.tsv')
        history_proteins=SeqIO.to_dict(SeqIO.parse(processed/'maize/maize_candidate_neighbors.faa','fasta'))
        historical=catalogue[catalogue.search_stage=='historical_Maize_S3'].set_index('locus_id')
        for window in SeqIO.parse(processed/'maize/maize_candidate_windows.gbk','genbank'):
            if window.id not in historical.index:continue
            row=historical.loc[window.id].to_dict();metadata[window.id]={'locus_id':window.id,**row}
            SeqIO.write(SeqRecord(window.seq,id=window.id,description=f'historical_window_start_0based={row["window_start"]}'),windows_out,'fasta')
            assigned=history[history.candidate_id==window.id].set_index('protein_id').role.to_dict()
            for f in history_coords[(history_coords.candidate_id==window.id)&history_coords.target_id.isin(assigned)].to_dict('records'):
                query=f['target_id'];strand=1 if f['strand']=='+' else -1
                dna=window.seq[f['start_0based']-int(row['window_start']):f['end_0based_exclusive']-int(row['window_start'])]
                if strand==-1:dna=dna.reverse_complement()
                protein=history_proteins[query];assert translated(dna)==str(protein.seq)
                sequences[query]=protein;dna_records[query]=SeqRecord(dna,id=query,description='historical_coding_orientation')
                gene_rows.append({'locus_id':window.id,'merged_locus_id':row['merged_locus_id'],'genome_id':row['genome_id'],
                    'assembly':row['assembly'],'biosample':row['biosample'],'organism':row['organism'],'strain':row['strain'],
                    'contig':row['contig'],'locus_status':row['status'],'family':row['family'],'export_role':assigned[query],
                    'query_id':query,'start':f['start_0based'],'end':f['end_0based_exclusive'],'strand':strand,
                    'origin':'historical_annotated_CDS','product':f['product'],'aa_length':len(protein.seq),
                    'sequence_sha256':sha(protein.seq),'coding_DNA_sha256':sha(dna),'coordinate_DNA_translation_match':True,
                    'status_basis':'preserved_historical_manual_review;precursor_like_labels_not_reclassified'})
            core=window.seq[int(row['core_start'])-int(row['window_start']):int(row['core_end'])-int(row['window_start'])]
            if row['C_strand']==-1:core=core.reverse_complement()
            assert sha(core)==row['core_sha256']
            core_records.append(SeqRecord(core,id=window.id,description='historical_core_includes_ABC_not_redefined'))
    pd.DataFrame(gene_rows).fillna('').to_csv(output/'genes.tsv',sep='\t',index=False)
    pd.DataFrame(metadata.values()).fillna('').to_csv(output/'locus_coordinates_and_statuses.tsv',sep='\t',index=False)
    SeqIO.write(sequences.values(),output/'genes.faa','fasta')
    SeqIO.write(dna_records.values(),output/'genes.fna','fasta')
    SeqIO.write(core_records,output/'source_locus_core_spans.fna','fasta')
    for path in [*tables.glob('diversity_final_*'),tables/'diversity_C_fragment_audit.tsv',tables/'diversity_merrild_sequence_audit.tsv']:
        shutil.copy2(path,output/path.name)
    references=output/'references';references.mkdir()
    for name in ['diversity_published_DNA.fna','diversity_published_C.faa','literature_unique_precursors.faa',
                 'literature_precursor_occurrences.tsv','literature_bgc_manifest.tsv','diversity_role_reference.tsv',
                 'diversity_M2_precursor_references.faa','diversity_precursor_core_sites.tsv']:
        shutil.copy2((processed/"discovery"/name) if name.startswith("diversity_") else (processed/"literature"/name),references/name)
    provenance=output/'provenance';provenance.mkdir()
    for path in [ROOT/'protocols/diversity_selection_review_protocol.json',ROOT/'protocols/diversity_protocol.json',
                 ROOT/'inputs/selections/review/diversity_carrier_identity_review.tsv',*ROOT.glob('inputs/provenance/shards/diversity_*_context_downloads.jsonl'),
                 ROOT/'inputs/provenance/diversity_extension_download_manifest.json',ROOT/'inputs/provenance/historical/diversity_selection_preservation_manifest.json']:
        shutil.copy2(path,provenance/path.name)
    code=output/'code';code.mkdir()
    for name in ['curate_diversity_contexts.py','collect_diversity_catalogue.py','merge_diversity_catalogues.py',
                 'diagnose_merrild_precursors.py','compare_diversity_catalogues.py','export_diversity_handoff.py']:
        shutil.copy2(ROOT/'scripts/steps'/name,code/name)
    source_review=read(tables/'diversity_final_source_review.tsv')
    for path in [ROOT/name for name in sorted(set(source_review.primary_metadata_file))]:
        shutil.copy2(path,provenance/path.name)
    (output/'README.md').write_text('''# Завершённый снимок diversity для HGT, 13.10.2026

**Непрохождение diversity-фильтра не означает отсутствия BGC. Неопределённые локусы нельзя использовать как отрицательные HGT-контроли.**

Снимок содержит итоговый каталог, все его представления и отдельный неопределённый слой. Активные файлы HGT не изменены. Выборку для HGT и ортологию нужно определять под собственный вопрос; предсказанная роль не устанавливает ортологию или продукт.

- `diversity_final_BGC_catalogue.tsv`: группы носитель/локус/вариант после устранения физических и документированных коллекционных повторов.
- `diversity_final_occurrences.tsv`: все представления, связанные через `merged_locus_id`.
- `locus_coordinates_and_statuses.tsv`, `genes.tsv`: исходные contig-координаты, цепи, носители, основания статусов; `genes.faa`/`genes.fna` используют `query_id`. DNA записана в направлении кодирования.
- `context_windows.fna.gz`: исходная ориентация contig; абсолютное начало окна дано в метаданных. Координаты везде 0-based, конец не включён.
- `source_locus_core_spans.fna`: определённые интервалы всех положительных представлений. Исторические Maize-интервалы включают ABC, современные интервалы определены выбранными A/B1/B2/C; границы нельзя приравнивать к границам продукта.
- Таблицы доказательств, сопоставление именно с M2 и литературный аудит сохранены отдельно. BLAST, домены, аннотация и положение имеют разные основания; функции не подтверждены для каждого локуса.
- `references/`, `provenance/`, `code/` и `manifest.json`: фиксированные сравнения, протоколы, URL/хеши исходных сборок, код и контрольные суммы. Полные host assemblies остаются доступными по исходным ledger, в снимок включены необходимые окна.

Неопределённые гены могут быть частичными или иметь неразрешённую фазу аннотации; это видно в проверке перевода координатной DNA. Исторические `precursor_like` сохранены с прежним основанием. Совпадение предшественника или сегмента M2 не доказывает производство лариоцидина.
''')
    files={str(p.relative_to(output)):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
           for p in sorted(output.rglob('*')) if p.is_file()}
    manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'primary_groups':len(catalogue),
              'source_locus_core_spans':len(core_records),'contexts_with_coordinates':len(metadata),
              'exported_gene_records':len(gene_rows),'coordinate_system':'0-based half-open, original contig strand +/-1',
              'warning':'Failure of the diversity filter does not establish absence of a BGC. Uncertain loci are not negative HGT controls.',
              'HGT_active_files':'not_modified','files':files}
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:v for k,v in manifest.items() if k!='files'}))


if __name__=='__main__':
    main()
