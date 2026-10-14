#!/usr/bin/env python3
"""Freeze one interbiome extension and bounded exact Merrild carrier controls."""
import csv,json,re,zipfile
from collections import defaultdict
import pandas as pd
from prepare_literature_reference import parse_member
from profile_search import ROOT
from collect_diversity_catalogue import wgs_project
from source_cache import Cache
from prepare_diversity_published import reference_metadata

def main():
    raw=ROOT/'data/raw/diversity_extension_2026-10-13'
    cache=Cache(raw,ROOT/'inputs/provenance/diversity_extension_download_manifest.json')
    for biome,version in [('marine','v2.0'),('marine-sediment','v1.0'),('human-gut','v3.0')]:
        base=f'https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/{biome}/{version}/'
        for name in [f'README_{version}.txt','genomes-all_metadata.tsv']+(['known_issues.txt'] if biome=='human-gut' else []):
            cache.download(biome+'_'+name.replace('.tsv','').replace('.txt',''),base+name,'.tsv' if name.endswith('.tsv') else '.txt')
    old=pd.concat([pd.read_csv(x,sep='\t',keep_default_na=False) for x in sorted((ROOT/'inputs/selections/shards').glob('diversity_species_*_selection.tsv'))])
    # Original assembly aliases and species-taxid are read from the frozen official summary.
    with (ROOT/'data/raw/diversity_2026-10-12/refseq_bacterial_assembly_summary.txt').open() as h:
     for line in h:
      if line.lstrip('# ').startswith('assembly_accession\t'):fields=line.lstrip('# ').rstrip().split('\t');break
     assemblies=[dict(zip(fields,v)) for v in csv.reader(h,delimiter='\t') if v]
    by_sample,by_wgs,by_strain=defaultdict(list),defaultdict(list),defaultdict(list)
    by_acc={}
    for r in assemblies:
     by_acc[r['assembly_accession']]=r
     by_sample[r['biosample']].append(r);by_wgs[wgs_project(r['wgs_master'])].append(r)
     by_strain[r['infraspecific_name'].removeprefix('strain=')].append(r)
    known=set(old.genome_id)|set(old.assembly)
    known|={by_acc[x]['gbrs_paired_asm'] for x in old.genome_id if x in by_acc}
    meta=[];selection=[];counts={};overlaps=[]
    for biome,version,label in [('marine','v2.0','Marine'),('marine-sediment','v1.0','Marine_Sediment'),('human-gut','v3.0','Human_Gut')]:
     t=pd.read_csv(raw/(biome+'_genomes-all_metadata.tsv'),sep='\t',keep_default_na=False)
     q=t[(t.Genome==t.Species_rep)&pd.to_numeric(t.Completeness).ge(50)&pd.to_numeric(t.Contamination).lt(10)&t.Lineage.str.startswith('d__Bacteria;')].copy()
     counts[label]={'version':version,'all_records':len(t),'eligible_species_reps':len(q),'source_types':q.Genome_type.value_counts().to_dict()}
     for r in q.to_dict('records'):
      acc=r['Genome'];url=f'https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/{biome}/{version}/species_catalogue/{acc.split(".")[0][:-2]}/{acc}/genome/{acc}.faa'+('.gz' if biome=='human-gut' else '')
      # Retain duplicate source memberships separately; one physical proteome is screened.
      if acc in known or r['Genome_accession'] in known:
       overlaps.append({'genome_id':acc,'source':label,'assembly':r['Genome_accession'],'basis':'exact_accession_previously_selected'});continue
      known.add(acc);known.add(r['Genome_accession']) if r['Genome_accession'] else None
      selection.append({'source':label,'genome_id':acc,'assembly':r['Genome_accession'],'biosample':r['Sample_accession'],'url':url,'compression':'gzip' if biome=='human-gut' else 'plain','context_url':r['FTP_download'].replace('ftp://','https://',1)})
      meta.append({**r,'source':label,'catalogue_version':version})
     pd.DataFrame(q).to_csv(ROOT/f'inputs/selections/diversity_extension_{biome}_species.tsv',sep='\t',index=False)
    published=pd.read_csv(ROOT/'data/processed/literature_bgc_manifest.tsv',sep='\t',keep_default_na=False)
    published=published[published.source_group.isin(['canonical','noncanonical'])]
    coverage=[];supplement={}
    with zipfile.ZipFile(ROOT/'data/raw/research_2026-09-10/triculamin_lasso_BGCs.zip') as archive:
     for r in published.to_dict('records'):
      record,*_=parse_member(r['archive_member'],archive.read(r['archive_member']))
      source=next((f.qualifiers for f in record.features if f.type=='source'),{})
      strain=r['strain'] or (re.search(r'NPDC\d+',r['archive_member']).group() if re.search(r'NPDC\d+',r['archive_member']) else '')
      sample=re.match(r'(SAM[END][A-Z]?\d+)',r['accession']);sample=sample.group() if sample else ''
      taxids=[x.split(':')[1] for x in source.get('db_xref',[]) if x.startswith('taxon:')]
      candidates=[];basis='unresolved_export_metadata'
      wp=wgs_project(r['accession'])
      if sample and by_sample.get(sample):candidates=by_sample[sample];basis='exact_BioSample'
      elif wp and by_wgs.get(wp):candidates=by_wgs[wp];basis='exact_WGS_project'
      elif strain:
       candidates=[x for x in by_strain.get(strain,[]) if x['taxid'] in taxids or strain.startswith('NPDC') or (r['organism']!='unspecified' and (x['organism_name'].split()[:2]==r['organism'].split()[:2]))]
       if candidates:basis='strain_label_and_taxon_or_unique_NPDC_label'
      candidate_ids=sorted({x['assembly_accession'] for x in candidates})
      exact_old=sorted(set(candidate_ids)&set(old.genome_id))
      sp={x['species_taxid'] for x in candidates};alt=sorted(x for x in old.genome_id if x in by_acc and by_acc[x]['species_taxid'] in sp and x not in candidate_ids)
      status='exact_published_carrier_selected' if exact_old else 'other_strain_same_NCBI_species_selected' if alt else 'resolved_published_carrier_absent' if candidates else 'carrier_not_resolved'
      coverage.append({'bgc_id':r['bgc_id'],'family':r['source_group'],'organism':r['organism'],'strain':strain,'accession':r['accession'],'biosample_from_export':sample,'resolution_basis':basis,'resolved_refseq_assemblies':';'.join(candidate_ids),'previous_exact_assemblies':';'.join(exact_old),'previous_same_species_other_strain':';'.join(alt),'initial_selection_status':status,'export_DNA_sha256':r['record_sequence_sha256']})
      # Resolve only bounded exact carriers, never all non-reference assemblies of a species.
      for x in candidates:
       acc=x['assembly_accession']
       if acc in known:continue
       supplement[acc]=x
    for acc,r in sorted(supplement.items()):
     base=r['ftp_path'].replace('ftp://','https://',1).rstrip('/');folder=base.rsplit('/',1)[-1]
     if base=='na' or r['protein_coding_gene_count']=='0':continue
     selection.append({'source':'RefSeq','genome_id':acc,'assembly':acc,'biosample':r['biosample'],'url':f'{base}/{folder}_protein.faa.gz','compression':'gzip','context_url':f'{base}/{folder}_genomic.gbff.gz'})
    pd.DataFrame(meta).to_csv(ROOT/'inputs/selections/diversity_extension_mgnify_metadata.tsv',sep='\t',index=False)
    pd.DataFrame(coverage).to_csv(ROOT/'results/tables/diversity_merrild_coverage.tsv',sep='\t',index=False)
    pd.DataFrame(list(supplement.values())).to_csv(ROOT/'inputs/selections/diversity_extension_refseq_assemblies.tsv',sep='\t',index=False)
    pd.DataFrame(selection).to_csv(ROOT/'inputs/selections/diversity_extension_selection.tsv',sep='\t',index=False)
    pd.DataFrame(overlaps,columns=['genome_id','source','assembly','basis']).to_csv(ROOT/'results/tables/diversity_extension_selection_overlaps.tsv',sep='\t',index=False)
    accessions=[r['genome_id'] for r in selection if r['source']=='RefSeq'];references=[]
    for i in range(0,len(accessions),50):
        url='https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'+','.join(accessions[i:i+50])+'/dataset_report?page_size=1000'
        path=cache.download(f'extension_refseq_metadata_{i//50:02d}',url,'.json')
        payload=json.loads(path.read_text())
        if payload.get('next_page_token'):raise ValueError('Unexpected metadata pagination')
        references.extend(reference_metadata(r) for r in payload['reports'])
    if set(r['assembly'] for r in references)!=set(accessions):raise ValueError('Incomplete exact RefSeq metadata')
    pd.DataFrame(references).to_csv(ROOT/'inputs/selections/diversity_extension_refseq_metadata.tsv',sep='\t',index=False)
    summary={'status':'fixed_selection_before_launch','sources':counts,'selected_proteomes':pd.DataFrame(selection).source.value_counts().to_dict(),'exact_accession_overlaps_omitted':len(overlaps),'merrild_initial_coverage':pd.DataFrame(coverage).initial_selection_status.value_counts().to_dict(),'criteria':'unchanged_diversity_protocol','stopping_rule':'one_pooled_extension_then_result_interpretation_no_all_nonreference_RefSeq','HGT':'unchanged'}
    (ROOT/'protocols/diversity_extension_scope.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
