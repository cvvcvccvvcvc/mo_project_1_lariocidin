#!/usr/bin/env python3
"""Preserve former calls and describe the final selection changes and literature trace."""
import json

import pandas as pd
from profile_search import ROOT


def read(path):
    return pd.read_csv(path, sep='\t', keep_default_na=False)


def main():
    tables=ROOT/'results/tables'
    before=pd.concat([read(tables/f'{p}_all_contexts.tsv') for p in ['diversity_species','diversity_extension','diversity_controls']],ignore_index=True)
    after=pd.concat([read(tables/f'{p}_review_all_contexts.tsv') for p in ['diversity_species','diversity_extension','diversity_controls']],ignore_index=True)
    old=before[['locus_id','status','family','precursor_query_ids','core_sha256']]
    joined=old.merge(after,on='locus_id',suffixes=('_previous',''),validate='one_to_one')
    assert len(before)==len(after)==len(joined)==8540
    changed=joined[(joined.status!=joined.status_previous)|(joined.family!=joined.family_previous)
                   |(joined.precursor_query_ids!=joined.precursor_query_ids_previous)|(joined.core_sha256!=joined.core_sha256_previous)]
    changed.to_csv(tables/'diversity_final_status_changes.tsv',sep='\t',index=False)
    lookup=after.set_index('locus_id').to_dict('index')
    stages=read(tables/'diversity_merrild_screen_stages.tsv')
    stages['previous_stage']=stages.stage
    for i,row in stages.iterrows():
        matched=[lookup[q] for q in row.locus_id.split(';') if q in lookup]
        if not matched:continue
        passing=[x for x in matched if x['status']=='family_bgc_candidate']
        stages.loc[i,'stage']='context_supported_BGC' if passing else 'C_nominated_BGC_criteria_not_met'
        stages.loc[i,'status']=';'.join(sorted({x['status'] for x in matched}))
        stages.loc[i,'reason']=';'.join(sorted({x['interpretation_reason'] for x in matched}))
    stages.to_csv(tables/'diversity_final_merrild_screen_stages.tsv',sep='\t',index=False)
    diagnostics=read(tables/'diversity_merrild_sequence_audit.tsv')
    diagnostics['final_status']=diagnostics.locus_id.map(lambda q:lookup[q]['status'])
    diagnostics['final_reason']=diagnostics.locus_id.map(lambda q:lookup[q]['interpretation_reason'])
    diagnostics['final_precursor_ids']=diagnostics.locus_id.map(lambda q:lookup[q]['precursor_query_ids'])
    diagnostics.to_csv(tables/'diversity_final_merrild_sequence_audit.tsv',sep='\t',index=False)
    final=read(tables/'diversity_final_BGC_catalogue.tsv')
    source_review=read(tables/'diversity_final_source_review.tsv').set_index('genome_id')
    for column in source_review.columns:
        final[column]=final.genome_id.map(source_review[column]).fillna('')
    m2=read(tables/'diversity_final_M2_precursor_similarity.tsv')
    for column,reference,flag in [('M2_full_exact_A_count','M2_full_precursor','exact_entire_reference'),
                                  ('M2_mature_exact_segment_A_count','M2_mature_LAR','contains_exact_reference_segment')]:
        counts=m2[(m2.comparison_reference==reference)&m2[flag]].groupby('merged_locus_id').query_id.nunique()
        final[column]=final.merged_locus_id.map(counts).fillna(0).astype(int)
    final.to_csv(tables/'diversity_final_BGC_catalogue.tsv',sep='\t',index=False)
    columns=['merged_locus_id','family','status','organism','strain','primary_strain','all_assembly','all_source','all_biosample',
             'architecture','C_aa_length','C_reference_identity_percent','C_domains','precursor_C_strand_agreement',
             'precursor_C_overlap_nt','C_difference_interpretation','carrier_relation_to_literature',
             'locus_relation_to_literature','M2_full_exact_A_count','M2_mature_exact_segment_A_count',
             'primary_origin_assessment','primary_metadata_caveat']
    final[columns].to_csv(tables/'diversity_final_catalogue_summary.tsv',sep='\t',index=False)
    previous=read(tables/'diversity_merged_BGC_catalogue.tsv')
    new=final[~final.merged_locus_id.isin(previous.merged_locus_id)]
    new.to_csv(tables/'diversity_final_added_loci.tsv',sep='\t',index=False)
    transitions=changed.groupby(['status_previous','status']).size().reset_index(name='contexts').to_dict('records')
    best=stages.assign(rank=stages.stage.map({'context_supported_BGC':0,'C_nominated_BGC_criteria_not_met':1,
        'published_C_DNA_not_exact_in_current_carrier':2,'export_has_no_explicit_C_to_trace':3}))
    best=best.sort_values('rank').drop_duplicates('bgc_id')
    summary={'contexts_reevaluated':len(after),'previous_status_counts':before.status.value_counts().to_dict(),
        'final_status_counts':after.status.value_counts().to_dict(),'changed_contexts':len(changed),'transitions':transitions,
        'previous_merged_groups':len(previous),'final_merged_groups':len(final),'additional_merged_groups':len(new),
        'old_merged_groups_retained':int(previous.merged_locus_id.isin(final.merged_locus_id).sum()),
        'literature_export_best_stage_counts':best.stage.value_counts().to_dict(),
        'new_locus_relation':new.locus_relation_to_literature.value_counts().to_dict(),
        'source_primary_counts':after[after.status=='family_bgc_candidate'].groupby(['source','family']).size().reset_index(name='source_locus_rows').to_dict('records'),
        'new_database_records':0,'full_protein_screen_repeated':False,'HGT':'unchanged'}
    (tables/'diversity_final_selection_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ['contexts_reevaluated','changed_contexts','final_merged_groups','additional_merged_groups','old_merged_groups_retained']}))


if __name__=='__main__':
    main()
