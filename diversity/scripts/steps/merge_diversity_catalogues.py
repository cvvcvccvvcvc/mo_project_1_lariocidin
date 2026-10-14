#!/usr/bin/env python3
"""Merge completed stages, preserve evidence axes, and compare precursors to M2."""
import argparse
import hashlib
import json
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from curate_diversity_contexts import BLAST_FIELDS, blast
from extract_diversity_contexts import read_genome
from source_cache import ROOT


def read(path):
    return pd.read_csv(path, sep='\t', keep_default_na=False)


def sha(sequence):
    return hashlib.sha256(str(sequence).upper().encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-suffix', default='', help='Separately recurated stage tables, e.g. review')
    parser.add_argument('--output-prefix', default='diversity_merged')
    args = parser.parse_args()
    tables, processed = ROOT/'results/tables', ROOT/'data/processed/discovery'
    prefixes = ['diversity_species', 'diversity_extension', 'diversity_controls']
    catalogues, uncertain, genes = [], [], {}
    for prefix in prefixes:
        result_prefix = prefix+('_'+args.input_suffix if args.input_suffix else '')
        data = read(tables/f'{result_prefix}_BGC_catalogue.tsv')
        data['search_stage'] = 'species_RefSeq_Soil' if prefix == 'diversity_species' else 'final_interbiome_extension'
        data['technical_batch'] = prefix
        catalogues.append(data)
        u = read(tables/f'{result_prefix}_uncertain.tsv');u['technical_batch'] = prefix;uncertain.append(u)
        genes.update(SeqIO.to_dict(SeqIO.parse(processed/f'{result_prefix}_BGC_genes.faa', 'fasta')))
    # Historical S3 calls retain their original review basis and precursor-like labels.
    history = read(tables/'diversity_bacillus_novelty.tsv')
    roles = read(ROOT/'inputs/selections/maize_curated_core_roles.tsv')
    maize_loci = read(tables/'maize_candidate_loci.tsv').set_index('target_id')
    maize_neighbors = read(tables/'maize_candidate_neighbors.tsv')
    maize_genes = SeqIO.to_dict(SeqIO.parse(processed.parent/'maize_candidate_neighbors.faa', 'fasta'))
    historical = []
    for r in history.to_dict('records'):
        candidate = r['candidate_id'];assigned = roles[roles.candidate_id == candidate]
        a_ids = list(assigned[assigned.role == 'precursor_like'].protein_id)
        genes.update({x: maize_genes[x] for x in assigned.protein_id})
        locus = maize_loci.loc[candidate]
        coordinates = maize_neighbors[(maize_neighbors.candidate_id == candidate) & maize_neighbors.target_id.isin(assigned.protein_id)]
        historical.append({'locus_id': candidate, 'genome_id': candidate.rsplit('_', 1)[0], 'source': 'Maize_Rhizosphere_historical',
            'organism': 'Bacillus wiedmannii', 'strain': r['strain'], 'assembly': r['refseq_assembly'], 'biosample': r['biosample'],
            'assembly_contig_multiset_sha256': r['assembly_contig_multiset_sha256'], 'C_sha256': r['C_sha256'],
            'core_sha256': r['core_sha256'], 'family': 'canonical', 'status': 'historical_review_supported_BGC',
            'contig': locus['contig'], 'C_start': locus['start_0based'], 'C_end': locus['end_0based_exclusive'],
            'C_strand': 1 if locus['strand'] == '+' else -1, 'C_query_id': candidate,
            'window_start': locus['window_start_0based'], 'window_end': locus['window_end_0based_exclusive'],
            'core_start': int(coordinates.start_0based.min()), 'core_end': int(coordinates.end_0based_exclusive.max()),
            'chosen_core_gene_ids': ';'.join(coordinates.sort_values('start_0based', ascending=locus['strand']=='+').target_id),
            'precursor_query_ids': ';'.join(a_ids), 'precursor_sha256': ';'.join(sha(genes[x].seq) for x in a_ids),
            'precursor_origins': 'historical_annotated_precursor_like', 'architecture': 'historical_A_candidates-B1-B2-C',
            'search_stage': 'historical_Maize_S3', 'technical_batch': 'historical_Maize_S3',
            'published_exact_core_records': r['exact_core_matches_in_422_records'],
            'carrier_comparison': 'not_matched_in_exported_carrier_metadata',
            'interpretation_reason': 'Original manual context review retained; not retrospectively reclassified'})
    all_rows = pd.concat(catalogues+[pd.DataFrame(historical)], ignore_index=True).fillna('')
    # Only assemblies occurring under multiple source IDs require physical reconciliation.
    duplicate_assemblies = set(all_rows.groupby('assembly_contig_multiset_sha256').genome_id.nunique().loc[lambda x:x>1].index)
    raw_paths = {}
    for path in (ROOT/'inputs/provenance/shards').glob('diversity_*_context_downloads.jsonl'):
        for line in path.read_text().splitlines():
            record=json.loads(line);raw_paths[record['genome_id']]=(ROOT/record['local_path'], record.get('source', 'RefSeq' if record['local_path'].endswith('.gbff.gz') else 'MGnify'))
    genomes = {};keys = []
    for r in all_rows.to_dict('records'):
        assembly_hash = r['assembly_contig_multiset_sha256']
        if assembly_hash in duplicate_assemblies:
            genome = r['genome_id']
            if genome not in genomes:
                if genome not in raw_paths:
                    raise ValueError(f'Physical duplicate needs its original full genome: {genome}')
                genomes[genome] = read_genome(*raw_paths[genome])[1]
            seq = str(genomes[genome][r['contig']].seq).upper();rev=str(Seq(seq).reverse_complement())
            reverse = rev < seq
            start,end = int(r['C_start']),int(r['C_end']);strand=int(r['C_strand'])
            if reverse:start,end,strand=len(seq)-end,len(seq)-start,-strand
            stop = end if strand == 1 else start
            keys.append(f'{assembly_hash}|{sha(min(seq,rev))}|{strand}|{stop}')
        else:
            keys.append(f"{assembly_hash}|{r['contig']}|{r['C_strand']}|{r['C_end'] if int(r['C_strand'])==1 else r['C_start']}")
    assembly_keys = ['assembly_bgc_'+sha(x)[:16] for x in keys]
    all_rows['assembly_locus_group_id'] = assembly_keys
    identity_review = read(ROOT/'inputs/selections/diversity_carrier_identity_review.tsv').set_index('genome_id').to_dict('index')
    for index,r in enumerate(all_rows.to_dict('records')):
        culture = identity_review.get(r['genome_id'],{}).get('culture_identity_id','')
        if culture:
            keys[index] = culture+'|'+r['core_sha256']
    all_rows['documented_culture_identity'] = all_rows.genome_id.map(lambda x:identity_review.get(x,{}).get('culture_identity_id',''))
    all_rows['merged_locus_id'] = ['bgc_'+sha(x)[:16] for x in keys]
    all_rows.to_csv(tables/f'{args.output_prefix}_occurrences.tsv', sep='\t', index=False)
    unique = all_rows.drop_duplicates('merged_locus_id').copy()
    groups = all_rows.groupby('merged_locus_id')
    for column in ['source','search_stage','genome_id','assembly','biosample','technical_batch']:
        unique['all_'+column] = unique.merged_locus_id.map(groups[column].agg(lambda x:';'.join(sorted(set(x)-{''}))))
    coverage = read(tables/'diversity_merrild_coverage.tsv')
    stages = read(tables/'diversity_merrild_screen_stages.tsv')
    known_loci = defaultdict(set)
    for r in stages[stages.stage.isin(['context_supported_BGC','C_nominated_BGC_criteria_not_met'])].to_dict('records'):
        for locus in r['locus_id'].split(';'):known_loci[locus].add(r['bgc_id'])
    published_C = {sha(r.seq) for r in SeqIO.parse(processed/'diversity_published_C.faa','fasta')}
    published_A = {sha(r.seq) for path in [processed.parent/'literature_unique_precursors.faa',processed.parent/'all_cds.faa'] for r in SeqIO.parse(path,'fasta')}
    unique['exact_literature_carrier_locus_records'] = unique.merged_locus_id.map(lambda identifier:';'.join(sorted({record for locus in groups.get_group(identifier).locus_id for record in known_loci[locus]})))
    unique['locus_relation_to_literature'] = [
        'published_carrier_and_C_DNA_traced' if r['exact_literature_carrier_locus_records'] else
        'exact_exported_core_carrier_relation_unresolved' if r['published_exact_core_records'] else
        'unmatched_exported_core_literature_scope_limited' for r in unique.to_dict('records')]
    unique['C_exact_in_2025_exports'] = unique.C_sha256.isin(published_C)
    unique['A_exact_in_2025_exports'] = unique.precursor_sha256.map(lambda cell:all(x in published_A for x in cell.split(';') if x))
    latest = SeqIO.to_dict(SeqIO.parse(processed/'diversity_2026_gelatinamin_roles.faa','fasta'))
    gel_c = sha(latest['gelatinamin_2026__C'].seq);gel_a = sha(latest['gelatinamin_2026__A'].seq)
    unique['C_exact_gelatinamin_2026_plasmid'] = unique.C_sha256 == gel_c
    unique['A_exact_gelatinamin_2026_plasmid'] = unique.precursor_sha256.map(lambda x:gel_a in x.split(';'))
    tria = list(SeqIO.parse(processed/'diversity_2026_TriA_precursors.faa', 'fasta'))
    native_tria = {sha(r.seq) for r in tria if r.id == 'TriA'}
    engineered_tria = {sha(r.seq) for r in tria if r.id != 'TriA'}
    unique['A_exact_2026_native_TriA'] = unique.precursor_sha256.map(lambda x:bool(set(x.split(';')) & native_tria))
    unique['A_exact_2026_engineered_TriA'] = unique.precursor_sha256.map(lambda x:bool(set(x.split(';')) & engineered_tria))
    taxonomy = read(ROOT/'inputs/selections/diversity_carrier_taxonomy.tsv')
    taxa = {}
    for r in taxonomy.to_dict('records'):
        for acc in [r['assembly_accession'], r['gbrs_paired_asm']]:
            taxa[acc] = r
    published_species, published_strains, published_samples = defaultdict(set), defaultdict(set), defaultdict(set)
    for r in coverage.to_dict('records'):
        assemblies = {x for field in ['resolved_refseq_assemblies','additional_primary_assembly_mapping']
                      for x in r[field].split(';') if x}
        for assembly in assemblies:
            if assembly not in taxa:
                continue
            ref = taxa[assembly];species = ref['species_taxid']
            published_species[species].add(r['bgc_id'])
            strain = ref['infraspecific_name'].removeprefix('strain=')
            if strain and strain != 'na':published_strains[(species,strain)].add(r['bgc_id'])
            if ref['biosample'] != 'na':published_samples[ref['biosample']].add(r['bgc_id'])
    carrier_axes = []
    for r in unique.to_dict('records'):
        taxon = taxa.get(r['assembly'], {})
        species = taxon.get('species_taxid', '')
        same_species = published_species.get(species, set())
        same_strain = published_strains.get((species,r['strain']), set())
        same_sample = published_samples.get(r['biosample'], set())
        carrier_class = 'carrier_extension_unresolved'
        if r['exact_literature_carrier_locus_records']:
            carrier_class = 'known_published_carrier_and_locus'
        elif same_strain or same_sample:
            carrier_class = 'published_carrier_metadata_present_locus_relation_unresolved'
        elif same_species and r['strain']:
            carrier_class = 'additional_strain_relative_to_resolved_published_metadata'
        elif species and not same_species:
            carrier_class = 'additional_species_candidate_literature_coverage_incomplete'
        carrier_axes.append({'NCBI_species_taxid':species,'published_same_resolved_species_records':';'.join(sorted(same_species)),
            'published_same_strain_records':';'.join(sorted(same_strain)), 'published_same_BioSample_records':';'.join(sorted(same_sample)),
            'carrier_relation_to_literature':carrier_class})
    unique = pd.concat([unique.reset_index(drop=True),pd.DataFrame(carrier_axes)],axis=1)
    for index,r in unique.iterrows():
        review = identity_review.get(r['genome_id'],{})
        if review:
            unique.loc[index,'carrier_relation_to_literature'] = review['reviewed_carrier_relation']
            unique.loc[index,'carrier_review_source'] = review['evidence_source']
    unique['novelty_limit'] = 'sequence/export comparison; unregistered or unavailable literature loci remain unresolved; product untested'
    unique.to_csv(tables/f'{args.output_prefix}_BGC_catalogue.tsv', sep='\t', index=False)
    pd.concat(uncertain, ignore_index=True).fillna('').to_csv(tables/f'{args.output_prefix}_uncertain.tsv', sep='\t', index=False)
    core_evidence = []
    for prefix in prefixes:
        result_prefix = prefix+('_'+args.input_suffix if args.input_suffix else '')
        evidence = read(tables/f'{result_prefix}_core_gene_evidence.tsv')
        evidence['technical_batch'] = prefix
        core_evidence.append(evidence)
    combined_evidence = pd.concat(core_evidence, ignore_index=True).fillna('')
    locus_lookup = all_rows.set_index('locus_id').merged_locus_id.to_dict()
    combined_evidence['merged_locus_id'] = combined_evidence.locus_id.map(locus_lookup)
    combined_evidence.to_csv(tables/f'{args.output_prefix}_core_gene_evidence.tsv',sep='\t',index=False)
    # M2 mature sequence is taken from the experimental 9DFC polymer entity, not inferred for other systems.
    original = SeqIO.to_dict(SeqIO.parse(processed.parent/'all_cds.faa','fasta'))
    m2 = original['paenibacillus_m2__cds01']
    mature = next(r for r in SeqIO.parse(ROOT/'data/raw/diversity_extension_2026-10-13/m2_PDB9DFC_fasta.faa','fasta') if r.id.split('|')[0]=='9DFC_57')
    if str(m2.seq) != 'MRKVQVAPLKSYGAFLTVTLSKKSKPGDGKFGRGVKRG' or not str(m2.seq).endswith(str(mature.seq)):
        raise ValueError('M2 precursor/experimental mature peptide reference conflict')
    reference = processed/'diversity_M2_precursor_references.faa'
    SeqIO.write([SeqRecord(m2.seq,id='M2_full_precursor',description='CP169648.1_38aa'),SeqRecord(mature.seq,id='M2_mature_LAR',description='PDB9DFC_entity57_18aa_precursor_positions21_38')], reference,'fasta')
    precursor_ids = sorted({x for cell in unique.precursor_query_ids for x in cell.split(';') if x})
    fasta=processed/f'{args.output_prefix}_precursors.faa';SeqIO.write([genes[x] for x in precursor_ids],fasta,'fasta')
    blast(fasta,reference,processed/('diversity_M2_blast.tsv' if args.output_prefix == 'diversity_merged' else f'{args.output_prefix}_M2_blast.tsv'),cpu=4,short=True)
    alignments=pd.read_csv(processed/('diversity_M2_blast.tsv' if args.output_prefix == 'diversity_merged' else f'{args.output_prefix}_M2_blast.tsv'),sep='\t',names=BLAST_FIELDS)
    alignments=alignments.sort_values('bitscore',ascending=False).drop_duplicates(['query_id','subject_id'])
    evidence=[]
    for r in unique.to_dict('records'):
        for query in r['precursor_query_ids'].split(';'):
            sequence=str(genes[query].seq)
            for subject,refseq in [('M2_full_precursor',str(m2.seq)),('M2_mature_LAR',str(mature.seq))]:
                match=alignments[(alignments.query_id==query)&(alignments.subject_id==subject)]
                hit=match.iloc[0].to_dict() if len(match) else {}
                evidence.append({'merged_locus_id':r['merged_locus_id'],'locus_id':r['locus_id'],'source':r['source'],
                    'family':r['family'],'organism':r['organism'],'query_id':query,'sequence':sequence,'comparison_reference':subject,
                    'exact_entire_reference':sequence==refseq,'contains_exact_reference_segment':refseq in sequence,
                    'query_coverage_percent':100*(hit['qend']-hit['qstart']+1)/len(sequence) if hit else 0,
                    'M2_reference_coverage_percent':100*(hit['send']-hit['sstart']+1)/len(refseq) if hit else 0,**hit,
                    'interpretation':'sequence proximity only; other precursor cleavage, topology, activity and product remain untested'})
    pd.DataFrame(evidence).fillna('').to_csv(tables/('diversity_M2_precursor_similarity.tsv' if args.output_prefix == 'diversity_merged' else f'{args.output_prefix}_M2_precursor_similarity.tsv'),sep='\t',index=False)
    SeqIO.write([genes[x] for x in sorted({r['C_query_id'] for r in unique.to_dict('records')})],processed/f'{args.output_prefix}_C.faa','fasta')
    summary={'source_locus_rows':len(all_rows),'unique_assembly_locus_groups':len(set(assembly_keys)), 'unique_carrier_locus_variant_groups':len(unique), 'duplicate_source_rows':len(all_rows)-len(unique),
        'assembly_content_groups':all_rows.assembly_contig_multiset_sha256.nunique(), 'biosamples_in_occurrences':all_rows.biosample.nunique(), 'families':unique.family.value_counts().to_dict(),'stages':unique.search_stage.value_counts().to_dict(),
        'C_sequence_variants':unique.C_sha256.nunique(),'A_sequence_variants':len({x for cell in unique.precursor_sha256 for x in cell.split(';') if x}),
        'literature_relation':unique.locus_relation_to_literature.value_counts().to_dict(),'carrier_relation':unique.carrier_relation_to_literature.value_counts().to_dict(),'HGT':'unchanged'}
    (tables/f'{args.output_prefix}_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))


if __name__=='__main__':
    main()
