#!/usr/bin/env python3
"""Check exported A coding DNA for nominated literature loci rejected as BGCs."""
import json

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from extract_diversity_contexts import read_genome, translated
from fetch_hgt_assemblies import ROOT
from curate_diversity_contexts import read_blast


def read(path):
    return pd.read_csv(path, sep='\t', keep_default_na=False)


def audit_saved_C_integrity():
    """Flag adjacent CDSs matching successive portions of a shared C reference.

    High identity/query coverage nominate an integrity review, not a biochemical
    classifier. Raw hits and both gene coordinates remain available for review.
    """
    references=read(ROOT/'data/processed/literature_all_cds.tsv')
    c_references=set(references.loc[references.author_role=='cyclase','feature_id'])
    rows=[]
    for directory in sorted((ROOT/'data/processed/discovery/contexts').glob('diversity_*_contexts')):
        features=read(directory/'features.tsv')
        hits=read_blast(directory/'role_blast.tsv')
        hits=hits[hits.subject_id.isin(c_references)&(hits.pident>=90)&(hits.qcov>=90)]
        hits=hits[hits.subject_id.isin(hits[hits.query_id.isin(features.loc[features.is_C,'query_id'])].subject_id)]
        by_query={q:g.to_dict('records') for q,g in hits.groupby('query_id')}
        for locus,local in features[features.origin=='annotated_CDS'].groupby('locus_id'):
            c=local[local.is_C].iloc[0]
            ordered=local[local.strand==c.strand].sort_values('start',ascending=int(c.strand)==1).to_dict('records')
            i=next(i for i,f in enumerate(ordered) if f['query_id']==c.query_id)
            for j in [i-1,i+1]:
                if not 0<=j<len(ordered):continue
                first,second=(ordered[j],ordered[i]) if j<i else (ordered[i],ordered[j])
                gap=second['start']-first['end'] if c.strand==1 else first['start']-second['end']
                for a in by_query.get(first['query_id'],[]):
                    for b in by_query.get(second['query_id'],[]):
                        if a['subject_id']!=b['subject_id'] or a['send']>=b['sstart']:continue
                        rows.append({'locus_id':locus,'C_query_id':c.query_id,'first_gene':first['query_id'],
                            'second_gene':second['query_id'],'reference':a['subject_id'],'genomic_gap_nt':gap,
                            'first_ref_start':a['sstart'],'first_ref_end':a['send'],'second_ref_start':b['sstart'],
                            'second_ref_end':b['send'],'ref_length':a['slen'],'first_identity':a['pident'],
                            'second_identity':b['pident'],'first_length':first['aa_length'],'second_length':second['aa_length'],
                            'assessment':'consecutive_C_segments_intact_single_C_unresolved',
                            'limitation':'does_not_distinguish_pseudogene_assembly_error_or_functional_split_system'})
    result=pd.DataFrame(rows).drop_duplicates(['locus_id','first_gene','second_gene'])
    result.to_csv(ROOT/'results/tables/diversity_C_fragment_audit.tsv',sep='\t',index=False)
    print(json.dumps({'C_integrity_flags':len(result),'flagged_genomes':result.locus_id.str.split('__').str[0].nunique()}))


def main():
    tables = ROOT/'results/tables'
    stages = read(tables/'diversity_merrild_screen_stages.tsv')
    failed = stages[stages.stage == 'C_nominated_BGC_criteria_not_met']
    a = read(ROOT/'data/processed/literature_precursor_occurrences.tsv')
    published = {r.id:str(r.seq).upper() for r in SeqIO.parse(ROOT/'data/processed/discovery/diversity_published_DNA.fna','fasta')}
    paths = {};contexts={};features={};admitted={}
    for path in (ROOT/'inputs/provenance').glob('diversity_*_context_downloads.jsonl'):
        for line in path.read_text().splitlines():
            record=json.loads(line);paths[record['genome_id']]=ROOT/record['local_path']
    for prefix in ['diversity_species','diversity_extension','diversity_controls']:
        contexts.update(read(tables/f'{prefix}_all_contexts.tsv').set_index('locus_id').to_dict('index'))
        for directory in sorted((ROOT/'data/processed/discovery/contexts').glob(prefix+'_[0-9]*_contexts')):
            f=read(directory/'features.tsv');features.update({x:g for x,g in f.groupby('locus_id')})
            e=read(directory/'precursor_evidence.tsv');admitted.update({x:g for x,g in e.groupby('locus_id')})
    sequences={r.id:str(r.seq) for r in SeqIO.parse(ROOT/'data/processed/literature_unique_precursors.faa','fasta')}
    results=[];details=[]
    for genome,records in failed.groupby('genome_id'):
        annotations,contigs=read_genome(paths[genome],'RefSeq')
        for r in records.to_dict('records'):
            reference=a[a.bgc_id==r['bgc_id']]
            for locus in r['locus_id'].split(';'):
                context=contexts[locus];local=features[locus];e=admitted.get(locus,pd.DataFrame())
                if reference.empty:
                    results.append({'bgc_id':r['bgc_id'],'genome_id':genome,'locus_id':locus,'A_stage':'no_explicit_exported_A_to_trace'});continue
                for ref in reference.to_dict('records'):
                    dna=Seq(published[r['bgc_id']][int(ref['start_0based']):int(ref['end_0based_exclusive'])])
                    if int(ref['strand'])==-1:dna=dna.reverse_complement()
                    placements=[]
                    for name,contig in contigs.items():
                        for strand,query in [(1,str(dna)),(-1,str(dna.reverse_complement()))]:
                            start=str(contig.seq).upper().find(query)
                            while start>=0:
                                placements.append((name,start,start+len(query),strand));start=str(contig.seq).upper().find(query,start+1)
                    relevant=[x for x in placements if x[0]==context['contig'] and x[1]>=context['window_start'] and x[2]<=context['window_end']]
                    stage='exported_A_DNA_absent_from_current_assembly' if not placements else 'exported_A_DNA_outside_addressed_window' if not relevant else 'exported_A_DNA_in_window'
                    for name,start,end,strand in relevant:
                        nearby=local[(local.strand==strand)&(local.start>=start-3)&(local.end<=end+3)]
                        overlap=local[(local.origin=='annotated_CDS')&(local.aa_length>200)&(local.start<end)&(local.end>start)]
                        passing=set(nearby.query_id)&set(e.query_id if 'query_id' in e else [])
                        same_frame=[];other_frame=[]
                        for f in overlap.to_dict('records'):
                            frame=(start-int(f['start']))%3 if strand==1 else (int(f['end'])-end)%3
                            target=same_frame if strand==f['strand'] and frame==0 else other_frame
                            target.append(f['query_id'])
                        blockers=overlap.to_dict('records')
                        details.append({'bgc_id':r['bgc_id'],'genome_id':genome,'locus_id':locus,
                            'published_A_feature':ref['feature_id'],'A_start':start,'A_end':end,'A_strand':strand,
                            'A_sequence':sequences[ref['precursor_id']], 'A_translation':translated(dna),
                            'A_translation_matches_export':translated(dna)==sequences[ref['precursor_id']],
                            'A_start_codon':str(dna[:3]),'A_stop_codon':str(dna[-3:]),
                            'A_same_strand_as_C':strand==context['C_strand'],
                            'overlap_same_frame_CDS_ids':';'.join(same_frame),
                            'overlap_other_frame_CDS_ids':';'.join(other_frame),
                            'overlap_CDS_products':';'.join(f['product'] for f in blockers),
                            'overlap_CDS_lengths':';'.join(str(f['aa_length']) for f in blockers),
                            'overlap_focal_C':any(f['is_C'] for f in blockers),
                            'C_score':context['sequence_score'],'C_model_coverage':context['model_coverage'],
                            'C_aa_length':context['C_aa_length'],'C_partial_CDS':context['C_partial_CDS'],
                            'C_reference_identity_percent':context['C_reference_identity_percent'],
                            'C_reference_qcov':context['C_reference_qcov'],'C_reference_scov':context['C_reference_scov'],
                            'C_reference_groups':context['C_reference_groups'],
                            'short_feature_ids':';'.join(nearby.query_id),
                            'passing_short_feature_ids':';'.join(sorted(passing))})
                        results.append({'bgc_id':r['bgc_id'],'genome_id':genome,'locus_id':locus,'published_A_feature':ref['feature_id'],'A_stage':stage,
                            'same_strand_as_C':strand==context['C_strand'],'published_A_aa_length':ref['aa_length'],
                            'addressed_short_feature_count':len(nearby),'A_alignment_pass_count':len(passing),'overlap_large_annotated_CDS':len(overlap),
                            'C_strong':context['sequence_score']>=250 and context['model_coverage']>=.7,
                            'pipeline_reason':r['reason']})
                    if not relevant:
                        results.append({'bgc_id':r['bgc_id'],'genome_id':genome,'locus_id':locus,'published_A_feature':ref['feature_id'],'A_stage':stage,'pipeline_reason':r['reason']})
    result=pd.DataFrame(results).fillna('');result.to_csv(tables/'diversity_merrild_precursor_diagnostics.tsv',sep='\t',index=False)
    print(json.dumps({'rows':len(result),'A_stages':result.A_stage.value_counts().to_dict()}))
    detail=pd.DataFrame(details).fillna('')
    detail.to_csv(tables/'diversity_merrild_sequence_audit.tsv',sep='\t',index=False)
    print(json.dumps({'sequence_audit_rows':len(detail),'translation_matches':int(detail.A_translation_matches_export.sum()),
                      'same_frame_long_CDS_overlap':int(detail.overlap_same_frame_CDS_ids.ne('').sum()),
                      'other_frame_long_CDS_overlap':int(detail.overlap_other_frame_CDS_ids.ne('').sum())}))
    audit_saved_C_integrity()


if __name__=='__main__':
    main()
