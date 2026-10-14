#!/usr/bin/env python3
"""Map frozen publication sequences to frozen genomes and report uncertainty."""
from collections import Counter, defaultdict
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import warnings

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq

from prepare_h1 import ROOT, H1, TABLES, SEQUENCES, read_tsv, write_tsv


def digest(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def bare(accession):
    return accession.removeprefix("NZ_").split(".")[0]


def occurrences(sequence, query):
    start = sequence.find(query)
    while start != -1:
        yield start
        start = sequence.find(query, start + 1)


def select_c_roles(loci, genes, out=TABLES):
    by_locus = defaultdict(list)
    for gene in genes.values():
        if gene["author_role"] == "C":
            by_locus[gene["locus_id"]].append(gene["gene_id"])
    hits = defaultdict(dict)
    with (TABLES / "H1_C_homology.tsv").open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            query, reference = fields[:2]
            if genes[query]["lineage"] != genes[reference]["lineage"]:
                continue
            qcov = (int(fields[5]) - int(fields[4]) + 1) / int(fields[10])
            scov = (int(fields[7]) - int(fields[6]) + 1) / int(fields[11])
            if qcov < 0.7 or scov < 0.7:
                continue
            row = {"gene_id": query, "reference_gene": reference, "identity_percent": float(fields[2]),
                   "bitscore": float(fields[9]), "query_coverage": qcov, "reference_coverage": scov}
            previous = hits[genes[query]["locus_id"]].get(query)
            if previous is None or row["bitscore"] > previous["bitscore"]:
                hits[genes[query]["locus_id"]][query] = row
    roles = []
    for locus_id, locus in loci.items():
        annotated = by_locus[locus_id]
        candidates = sorted(hits[locus_id].values(), key=lambda r: (-r["bitscore"], r["gene_id"]))
        selected = annotated[0] if len(annotated) == 1 else candidates[0]["gene_id"] if candidates else ""
        comparison = hits[locus_id].get(selected, {})
        roles.append({"locus_id": locus_id, "lineage": locus["lineage"], "C_gene_id": selected,
                      "role_basis": "author_annotation" if len(annotated) == 1 else "published_context_and_C_homology_candidate",
                      "other_C_homology_genes": ";".join(r["gene_id"] for r in candidates if r["gene_id"] != selected),
                      "reference_gene": comparison.get("reference_gene", ""),
                      "identity_percent": comparison.get("identity_percent", ""),
                      "query_coverage": comparison.get("query_coverage", ""),
                      "reference_coverage": comparison.get("reference_coverage", ""),
                      "functional_status": "published_prediction_or_reference; not_inferred_from_BLAST_alone"})
    write_tsv(out / "H1_C_role_evidence.tsv", roles)
    return {r["locus_id"]: r for r in roles}


def translation_status(dna, gene, protein):
    if not protein:
        return "no_published_translation"
    table = int(gene["translation_table"])
    sequence = dna[int(gene["codon_start"]) - 1:]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        translated = str(Seq(sequence).translate(table=table)).rstrip("*")
    if protein == translated:
        return "exact_translation"
    if protein.startswith("M") and protein[1:] == translated[1:] and sequence[:3] in CodonTable.unambiguous_dna_by_id[table].start_codons:
        return "exact_with_bacterial_initiator_M"
    return "translation_mismatch_requires_review"


def metadata_strain(row):
    if row.get("infraspecific_name"):
        return row["infraspecific_name"]
    values = row.get("biosource", {}).get("infraspecieslist", [])
    return ";".join(x.get("sub_value", x.get("sub_name", "")) for x in values if x.get("sub_type") == "strain")


def n50(lengths):
    half, total = sum(lengths) / 2, 0
    for size in sorted(lengths, reverse=True):
        total += size
        if total >= half:
            return size
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=TABLES)
    parser.add_argument("--inputs", type=Path, default=TABLES / "H1_assembly_inputs.tsv")
    parser.add_argument("--metadata", type=Path, default=H1 / "ncbi/assembly_candidate_metadata.json")
    args = parser.parse_args()
    if not 0 <= args.shard < args.shards:
        parser.error("Invalid shard")
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    loci = {r["locus_id"]: r for r in read_tsv(TABLES / "H1_publication_loci.tsv")}
    genes = {r["gene_id"]: r for r in read_tsv(TABLES / "H1_publication_genes.tsv")}
    by_locus = defaultdict(list)
    for gene in genes.values():
        by_locus[gene["locus_id"]].append(gene)
    locus_dna = {r.id: str(r.seq).upper() for r in SeqIO.parse(SEQUENCES / "loci.fna", "fasta")}
    coding = {r.id: str(r.seq).upper() for r in SeqIO.parse(SEQUENCES / "coding_sequences.fna", "fasta")}
    locus_reverse = {key: str(Seq(value).reverse_complement()) for key, value in locus_dna.items()}
    coding_reverse = {key: str(Seq(value).reverse_complement()) for key, value in coding.items()}
    proteins = {r.id: str(r.seq) for r in SeqIO.parse(SEQUENCES / "proteins.faa", "fasta")}
    roles = select_c_roles(loci, genes, out)
    inputs = read_tsv(args.inputs)[args.shard::args.shards]
    candidate_metadata = json.loads(args.metadata.read_text())
    direct_links = {r["locus_id"]: r for r in read_tsv(TABLES / "H1_nucleotide_assembly_links.tsv")}
    manifest = json.loads((ROOT / "config/H1/H1_ncbi_download_manifest.json").read_text())
    quality, mappings, gene_checks, errors = [], [], [], []
    for genome_index, entry in enumerate(inputs, 1):
        genome_id = entry["genome_id"]
        path = ROOT / entry["path"]
        source = candidate_metadata.get(entry["assembly_accession"], {})
        ftp = manifest.get("genome__" + entry["assembly_accession"], {}).get("url", "")
        actual_assembly = re.search(r"(GC[AF]_\d+\.\d+)", ftp)
        actual_assembly = actual_assembly.group(1) if actual_assembly else entry["assembly_accession"]
        opening = gzip.open if path.suffix == ".gz" else open
        try:
            with opening(path, "rt") as handle, warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                records = list(SeqIO.parse(handle, "genbank"))
            if not records:
                raise ValueError("No GenBank records")
        except Exception as error:
            errors.append({"genome_id": genome_id, "error": str(error)})
            continue
        record_sequences = {r.id: str(r.seq).upper() for r in records}
        lengths = [len(r) for r in records]
        annotations = defaultdict(list)
        protein_locations = defaultdict(list)
        strains, biosamples, taxa, organisms = set(), set(), set(), set()
        cds_count = translated_count = 0
        for record in records:
            organisms.add(record.annotations.get("organism", ""))
            biosamples.update(x.split(":", 1)[1] for x in record.dbxrefs if x.startswith("BioSample:"))
            for feature in record.features:
                if feature.type == "source":
                    strains.update(feature.qualifiers.get("strain", []))
                    taxa.update(x.split(":", 1)[1] for x in feature.qualifiers.get("db_xref", []) if x.startswith("taxon:"))
                if feature.type != "CDS":
                    continue
                cds_count += 1
                protein = "".join(feature.qualifiers.get("translation", [])).rstrip("*")
                if protein:
                    translated_count += 1
                item = (int(feature.location.start), int(feature.location.end), feature.location.strand,
                        ";".join(feature.qualifiers.get("protein_id", [])))
                annotations[(record.id, *item[:3])].append(item[3])
                if protein:
                    protein_locations[protein].append((record.id, *item))
        sequence_hash = digest(";".join(sorted(digest(min(s, str(Seq(s).reverse_complement()))) for s in record_sequences.values())))
        biosample = entry["biosample"] or source.get("biosampleaccn", "") or ";".join(sorted(biosamples))
        quality.append({"genome_id": genome_id, "requested_assembly_accession": entry["assembly_accession"],
                        "sequence_provider_assembly": actual_assembly, "input_scope": entry["input_scope"],
                        "biosample": biosample, "organism": entry["organism"] or ";".join(sorted(organisms)),
                        "taxid": entry["taxid"] or source.get("taxid", "") or ";".join(sorted(taxa)),
                        "species_taxid": source.get("species_taxid", source.get("speciestaxid", "")),
                        "strain": entry["source_strain"] or metadata_strain(source) or ";".join(sorted(strains)),
                        "assembly_level": entry["assembly_level"], "records": len(records), "genome_bp": sum(lengths),
                        "reported_coverage": source.get("coverage", "not_available"),
                        "N50_bp": n50(lengths), "ambiguous_bases": sum(sum(b not in "ACGT" for b in s) for s in record_sequences.values()),
                        "cds_count": cds_count, "translated_cds_count": translated_count, "parser_warning_count": len(caught),
                        "sequence_content_sha256": sequence_hash, "full_assembly_input": entry["input_scope"] != "single_reference_chromosome",
                        "marker_tree_input_status": "annotated_input_available; markers_not_yet_checked" if translated_count else "requires_standard_annotation",
                        "read_assembly_QC": "not_checked", "binning_or_contamination_QC": "not_established"})
        for locus_id, locus in loci.items():
            query = locus_dna[locus_id]
            reverse_query = locus_reverse[locus_id]
            found = []
            for record_id, sequence in record_sequences.items():
                found.extend((record_id, position, position + len(query), strand, "whole_export_exact")
                             for strand, needle in [(1, query), (-1, reverse_query)] for position in occurrences(sequence, needle))
            if not found:
                c_id = roles[locus_id]["C_gene_id"]
                if not c_id:
                    continue
                c_gene, c_dna = genes[c_id], coding[c_id]
                for record_id, sequence in record_sequences.items():
                    for strand, needle in [(1, c_dna), (-1, coding_reverse[c_id])]:
                        for position in occurrences(sequence, needle):
                            orientation = strand * int(c_gene["strand"])
                            offset = position - int(c_gene["start_0based"]) if orientation == 1 else position - (len(query) - int(c_gene["end_0based"]))
                            found.append((record_id, offset, offset + len(query), orientation, "C_DNA_exact_window_not_exact"))
                if not found:
                    for record_id, start, end, native_strand, annotation_id in protein_locations.get(proteins.get(c_id, ""), []):
                        orientation = native_strand * int(c_gene["strand"])
                        offset = start - int(c_gene["start_0based"]) if orientation == 1 else start - (len(query) - int(c_gene["end_0based"]))
                        found.append((record_id, offset, offset + len(query), orientation, "C_protein_exact_coding_DNA_or_window_differs"))
            for record_id, left, right, orientation, evidence in found:
                sequence = record_sequences[record_id]
                source_match = bare(locus["accession_evidence"]) == bare(record_id)
                linked = entry["assembly_accession"] in direct_links[locus_id]["linked_assemblies"].split(";") or actual_assembly in direct_links[locus_id]["linked_assemblies"].split(";")
                barcode = re.search(r"NPDC\d+", locus["archive_member"])
                barcode_match = bool(barcode and barcode.group() in set(re.findall(r"NPDC\d+", ";".join(strains) + ";" + metadata_strain(source) + ";" + entry["source_strain"])))
                primary = source_match or linked or barcode_match
                mapping_id = f"{locus_id}|{genome_id}|{record_id}|{left}|{orientation}"
                c_id = roles[locus_id]["C_gene_id"]
                core_start = core_end = ""
                exact_genes = checked_genes = 0
                for gene in by_locus[locus_id]:
                    start, end = int(gene["start_0based"]), int(gene["end_0based"])
                    if orientation == 1:
                        native_start, native_end = left + start, left + end
                    else:
                        native_start, native_end = left + len(query) - end, left + len(query) - start
                    gene_strand = orientation * int(gene["strand"])
                    if gene["gene_id"] == c_id:
                        core_start, core_end = native_start, native_end
                    if native_start < 0 or native_end > len(sequence):
                        status, translated, native_annotation = "outside_source_contig", "not_checked", ""
                    else:
                        native_parts = []
                        for part_start, part_end, part_strand in json.loads(gene["parts"]):
                            part_left, part_right = (left + part_start, left + part_end) if orientation == 1 else (left + len(query) - part_end, left + len(query) - part_start)
                            part_dna = sequence[part_left:part_right]
                            native_parts.append(part_dna if orientation * part_strand == 1 else str(Seq(part_dna).reverse_complement()))
                        native_dna = "".join(native_parts)
                        checked_genes += 1
                        status = "exact_coding_DNA" if native_dna == coding[gene["gene_id"]] else "coding_DNA_differs_at_projected_coordinates"
                        exact_genes += status == "exact_coding_DNA"
                        translated = translation_status(native_dna, gene, proteins.get(gene["gene_id"], ""))
                        native_annotation = ";".join(annotations.get((record_id, native_start, native_end, gene_strand), []))
                    protein_matches = [x for x in protein_locations.get(proteins.get(gene["gene_id"], ""), [])
                                       if x[0] == record_id and max(x[1], native_start) < min(x[2], native_end)]
                    gene_checks.append({"mapping_id": mapping_id, "locus_id": locus_id, "genome_id": genome_id,
                                        "gene_id": gene["gene_id"], "role_evidence": "C" if gene["gene_id"] == c_id else gene["author_role"],
                                        "contig": record_id, "projected_start_0based": native_start, "projected_end_0based": native_end,
                                        "strand": gene_strand, "DNA_status": status, "translation_status": translated,
                                        "exact_coordinate_annotation_ids": native_annotation,
                                        "overlapping_identical_protein_annotation_ids": ";".join(x[4] for x in protein_matches),
                                        "identical_protein_annotation_locations": ";".join(f"{x[0]}:{x[1]}-{x[2]}:{x[3]}:{x[4]}" for x in protein_matches),
                                        "published_partial_location": gene["partial_location"]})
                c_protein = proteins.get(c_id, "")
                mappings.append({"mapping_id": mapping_id, "locus_id": locus_id, "lineage": locus["lineage"],
                                 "genome_id": genome_id, "assembly_accession": entry["assembly_accession"],
                                 "sequence_provider_assembly": actual_assembly, "biosample": biosample, "contig": record_id,
                                 "export_start_0based": left, "export_end_0based": right, "orientation": orientation,
                                 "C_start_0based": core_start, "C_end_0based": core_end,
                                 "sequence_evidence": evidence, "source_nucleotide_base_match": source_match,
                                 "nucleotide_assembly_link": linked,
                                 "published_NPDC_strain_match": barcode_match,
                                 "carrier_assignment": "primary_accession_or_link_supported" if primary else "sequence_occurrence_only; published_carrier_not_assigned",
                                 "CDS_checked": checked_genes, "CDS_exact_DNA": exact_genes,
                                 "export_window_within_contig": left >= 0 and right <= len(sequence),
                                 "C_distance_to_contig_edge_bp": min(core_start, len(sequence) - core_end) if core_start != "" else "",
                                 "identical_C_protein_annotation_copies_in_genome": len(protein_locations.get(c_protein, [])),
                                 "C_copies_limitation": "exact_protein_only; divergent_paralogs_not_counted"})
        if genome_index % 25 == 0:
            print(f"Genomes audited {genome_index}/{len(inputs)}; sequence mappings {len(mappings)}", flush=True)
    write_tsv(out / "H1_genome_quality.tsv", quality)
    if mappings:
        write_tsv(out / "H1_locus_genome_correspondence.tsv", mappings)
        write_tsv(out / "H1_gene_genome_correspondence.tsv", gene_checks)
    (out / "H1_genome_parse_errors.json").write_text(json.dumps(errors, indent=2) + "\n")
    by_publication = defaultdict(list)
    for row in mappings:
        by_publication[row["locus_id"]].append(row)
    catalogue = []
    for locus_id, locus in loci.items():
        matches = by_publication[locus_id]
        primary = [r for r in matches if r["carrier_assignment"] == "primary_accession_or_link_supported"]
        full = [r for r in primary if r["sequence_evidence"] == "whole_export_exact"]
        catalogue.append({**locus, "C_gene_id": roles[locus_id]["C_gene_id"], "C_role_basis": roles[locus_id]["role_basis"],
                          "primary_genome_inputs": ";".join(sorted(set(r["genome_id"] for r in primary))),
                          "other_sequence_occurrences": ";".join(sorted(set(r["genome_id"] for r in matches if r not in primary))),
                          "mapping_status": "whole_export_and_carrier_supported" if full else "C_and_carrier_supported; context_requires_review" if primary else "sequence_occurrence_without_primary_carrier" if matches else "not_resolved_in_frozen_genomes",
                          "full_export_primary_inputs": len(full), "core_context_review": "required_before_component_orthology_or_HGT",
                          "assembly_read_binning_verdict": "no_raw_read_validation; absence_of_artifact_not_established"})
    write_tsv(out / "H1_carrier_catalogue.tsv", catalogue)
    summary = {"publication_targets": len(loci), "genome_inputs": len(inputs), "parsed_genomes": len(quality),
               "parse_errors": len(errors), "locus_sequence_mappings": len(mappings), "gene_mappings": len(gene_checks),
               "carrier_statuses": dict(Counter(r["mapping_status"] for r in catalogue)),
               "statuses_by_lineage": {line: dict(Counter(r["mapping_status"] for r in catalogue if r["lineage"] == line)) for line in ["canonical", "noncanonical"]},
               "gene_DNA_statuses": dict(Counter(r["DNA_status"] for r in gene_checks)),
               "gene_translation_statuses": dict(Counter(r["translation_status"] for r in gene_checks))}
    (out / "H1_correspondence_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
