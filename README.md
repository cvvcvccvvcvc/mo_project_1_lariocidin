# Lariocidin-like BGC project

## Current state

Stages S0–S3 are complete; S4 candidate triage is ready for review B. The fixed MGnify search produced no context-supported target BGC among records whose primary metadata support metagenomic origin. Three supported canonical BGC belong to cell cultures mislabeled or ambiguously classified as MAG by the catalogue and are nearly identical to published systems. No formal HGT test has yet been performed.

Start with the [S3 analytical checkpoint](reviews/2026-09-11-s3-analysis.md), then use the [research chronology](reviews/research-log.md) for decision history. The active scope is in [`../../docs/lariocidin-plan.md`](../../docs/lariocidin-plan.md); durable source facts and limitations are in [`../../docs/lariocidin-evidence.md`](../../docs/lariocidin-evidence.md). Older packages remain evidence for their stages, not the current protocol.

## Reference-stage contract

Question for the complete project: how does known lariocidin-like BGC diversity expand in public genomes and MAGs, and which suitable cases have concordant evidence compatible with HGT?

The current work tests the diversity part on one fixed catalogue but has not yet completed the focal HGT test. Its unit of observation is a BGC occurrence in an assembly. A unique protein sequence and an organism/assembly occurrence are retained as separate concepts. The 29 records are an author-curated reference set, not 29 experimentally confirmed lariocidin producers and not a random ecological sample.

The operational family rule remains: LrcC homology must be supported by coherent BGC context; an isolated profile hit is insufficient. Canonical and noncanonical candidates remain separate. S3 low-score labels are provisional triage, not experimentally validated exclusions.

## Reproduce the reference prototype

From this directory:

```bash
micromamba create -y -p .env -f environment.yml
scripts/run_reference.sh
```

The run downloads only the fixed reference archive and two accessioned NCBI records when absent, checks their SHA-256 values, rebuilds tables and alignments, runs IQ-TREE model selection plus 1,000 SH-aLRT and 1,000 ultrafast bootstrap replicates, regenerates figures, and validates the result. Exact versions from the completed run are in [`results/software_versions.tsv`](results/software_versions.tsv). The two tree searches took about seven wall-clock minutes in total on the machine used for this run; environment creation is separate.

## Important inputs and decisions

- [`config/sources.tsv`](config/sources.tsv): source URLs, revisions/accessions, access date, and checksums.
- [`config/reference_metadata.tsv`](config/reference_metadata.tsv): stable sample IDs, display names, assembly/contig accessions, and broad taxonomic groups. Some names/accessions repair incomplete exported metadata; the basis is recorded per row.
- [`config/lrcc_selection.tsv`](config/lrcc_selection.tsv): the explicit one-per-BGC LrcC choices and confidence. This file, not an invisible heuristic, owns the reference selection.
- [`data/README.md`](data/README.md): storage and coordinate conventions.

The BLASTP comparison to M2 is a within-record ranking aid, not a validated family classifier. Every chosen feature is currently the highest-scoring CDS in its submitted BGC, but eight choices remain marked `moderate` because their source annotation or sequence architecture is incomplete/distant.

## Main outputs

| Output | Meaning |
|---|---|
| [`data/processed/reference_manifest.tsv`](data/processed/reference_manifest.tsv) | One row per published BGC, with provenance and parsing diagnostics |
| [`data/processed/all_cds.tsv`](data/processed/all_cds.tsv) | Only translated CDS features, with zero-based half-open BGC coordinates |
| [`data/processed/lrcc_selected.tsv`](data/processed/lrcc_selected.tsv) | Final LrcC selections joined to taxonomy, evidence, and M2 similarity |
| [`results/tables/lrcc_candidate_ranking.tsv`](results/tables/lrcc_candidate_ranking.tsv) | Top four within-BGC candidates by M2-seeded BLASTP |
| [`results/tables/reference_bgc_summary.tsv`](results/tables/reference_bgc_summary.tsv) | BGC, LrcC length/confidence, and annotation-derived gene-category summary |
| [`results/tables/alignment_qc.tsv`](results/tables/alignment_qc.tsv) | Gap fraction and ungapped length for every aligned sequence |
| [`results/tables/lrcc_pairwise_identity.tsv`](results/tables/lrcc_pairwise_identity.tsv) | Pairwise identity over positions where both sequences have residues |
| [`results/tables/tree_clades.tsv`](results/tables/tree_clades.tsv) | Descendants, taxonomy, and SH-aLRT/UFBoot for each rooted internal clade |
| [`results/iqtree/reference_lrcc.treefile`](results/iqtree/reference_lrcc.treefile) | Full unrooted ML tree; internal labels are SH-aLRT/UFBoot |
| [`results/iqtree/canonical_lrcc.treefile`](results/iqtree/canonical_lrcc.treefile) | Sensitivity tree excluding four short Actinomycetota candidates |
| [`results/figures/lrcc_full_tree.png`](results/figures/lrcc_full_tree.png) | Full tree rooted on the paeninodin cyclase for display |
| [`results/figures/lrcc_canonical_tree.png`](results/figures/lrcc_canonical_tree.png) | Sensitivity tree rooted on the same outgroup |
| [`results/figures/representative_bgc_architectures.png`](results/figures/representative_bgc_architectures.png) | Four annotation-based example architectures |
| [`data/processed/literature_bgc_manifest.tsv`](data/processed/literature_bgc_manifest.tsv) | 393 published BGC records with source group, metadata, parse status, and evidence status |
| [`results/tables/literature_profile_calibration_summary.json`](results/tables/literature_profile_calibration_summary.json) | Cluster-held-out internal benchmark of canonical/noncanonical HMMs |
| [`results/tables/focal_split_sensitivity.tsv`](results/tables/focal_split_sensitivity.tsv) | Presence and support of the focal unrooted split across five tree variants |
| [`results/tables/pseudomonadota_core_roles.tsv`](results/tables/pseudomonadota_core_roles.tsv) | Explicit A–D assignments and inference basis for 14 Pseudomonadota loci |
| [`results/tables/focal_context_manifest.tsv`](results/tables/focal_context_manifest.tsv) | Versioned NCBI accessions, exact C checks, available flanks, and checksums |
| [`results/tables/maize_curated_candidates.tsv`](results/tables/maize_curated_candidates.tsv) | All 50 S3 hits with sequence, provenance, context status, and stated reason |
| [`results/tables/maize_checkpoint_summary.json`](results/tables/maize_checkpoint_summary.json) | Machine-readable S3 denominators, result, and limitations |
| [`results/figures/s3_search_checkpoint.png`](results/figures/s3_search_checkpoint.png) | Review-B summary of provenance and candidate triage |

The rooted files next to the PNGs are visualization products. Rooting does not turn branch support into HGT probability or establish the direction of evolution.

## Known limitations

- The original source uses malformed/stripped GenBank fields; parser warnings and translation/CDS discrepancies are retained. Missing colored gene categories in the architecture figure often mean missing annotations, not biological absence.
- Twenty-one of thirty aligned sequences fail IQ-TREE's composition chi-square check in the full analysis. This is a model/data warning, not a failed program run.
- The four Actinomycetota candidates are only 376–442 aa and form a long, distinct part of the alignment; one has more than 50% gaps. The sensitivity analysis therefore excludes them.
- The paeninodin sequence is the outgroup used in the published analysis, but its branch is long. Interpret deep rooting cautiously.
- The authors' MEGA `.mtsx` session confirms their 29 taxa, WAG model and 100 bootstraps, but is not a portable Newick tree. Our modern analysis uses MAFFT and IQ-TREE with the BIC-selected `Q.PFAM+I+G4` model; it is a transparent reanalysis, not a claim of exact split-by-split reproduction.
- Broad phylum/order labels are working metadata for visualization. They do not replace an independent host species tree.

## Reproduce S2a

After creating the environment, run in order:

```bash
scripts/run_literature_reference.sh
scripts/run_profile_calibration.sh
scripts/run_tree_sensitivity.sh
scripts/run_pseudomonadota_core.sh
scripts/run_focal_context.sh
```

These runs normalize the fixed literature archive, build/test profiles, repeat the focal tree under prespecified variants, compare A–D, and retrieve four versioned NCBI contexts. Raw downloads remain outside version control; source manifests and derived results are retained.

## Reproduce S3

After S2a products exist, run:

```bash
scripts/run_maize_search.sh
```

This verifies/downloads the pinned MGnify and BioSample/ENA inputs, reruns profile controls, searches catalogue representatives and every GFF-annotated CDS, audits provenance, extracts candidate contexts, and regenerates the review checkpoint. Raw/cache inputs require roughly 2 GB and remain outside Git; URLs, retrieval metadata, and SHA-256 values are retained in `config/`. The run does not perform the HGT stage.

## Research records

Keep dated interpretation and audit notes in `reviews/`; maintain current decisions in the shared plan rather than extending this README into a second protocol. Additional downloaded research inputs have separate `config/research_sources_*.tsv` manifests with URLs, access dates, versions, sizes, and SHA-256. These files do not change the reference pipeline's original input contract. Store small derived tables in `results/tables/` and link their derivation from the corresponding audit. Local Git history begins with a baseline after S2a; the journal, not that first commit, reconstructs earlier decisions.
