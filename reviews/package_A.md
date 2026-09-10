# Review package A: published 29-BGC reference

> Interpretation update, 10 September 2026: see the [analytical report](2026-09-10-analysis.md) and [audit](2026-09-10-audit.md). The counts below describe the original run. UFBoot ≥80 is a plotting/counting cutoff, not the usual ≥95 support guideline; the focal mixed clade remains exploratory. Original figures and calculations are preserved.

## What is ready

This is the first human interpretation checkpoint. It contains measured reference-stage results, not a claim that HGT or new lariocidin producers have been found.

- 29 published BGC records, 238 translated CDS, one explicit LrcC candidate per BGC.
- All 29 LrcC sequences are unique in this set. Selection confidence: 21 `high`, 8 `moderate`.
- Taxonomic composition: 11 Bacillota, 14 Pseudomonadota, 4 Actinomycetota.
- LrcC length range: 376–674 aa. Actinomycetota form the short group (376–442 aa); Bacillota are 659–674 aa; most Pseudomonadota are 638–653 aa, with one 575-aa candidate.
- Full alignment: 30 sequences including the outgroup, 790 columns. `Streptomyces griseocarneus` has more than 50% gaps. Pairwise identity among the 29 reference candidates ranges from 13.486% to 94.578% (median 27.680%) over pairwise comparable sites.
- IQ-TREE selected `Q.PFAM+I+G4` by BIC for both analyses. In the full tree, 20 of 27 labelled internal branches have both SH-aLRT and UFBoot at least 80.
- Removing four short Actinomycetota candidates and rebuilding the alignment/tree preserves 21 of 24 rooted internal clades among the remaining tips; three clades differ in each reconstruction. This supports many local relationships but does not make the whole topology invariant.

## Figures to inspect

1. [Full LrcC tree](../results/figures/lrcc_full_tree.png). Tip colours are phyla. Display rooting uses the paeninodin cyclase. Node numbers shown in the PNG are UFBoot values ≥80; both SH-aLRT/UFBoot values remain in the Newick and [`tree_clades.tsv`](../results/tables/tree_clades.tsv).
2. [Sensitivity tree](../results/figures/lrcc_canonical_tree.png), rebuilt without the four short Actinomycetota candidates.
3. [Representative BGC architectures](../results/figures/representative_bgc_architectures.png). Colours reflect available annotations. A grey/missing precursor or RRE must not yet be interpreted as biological absence.

The detailed tables are [`reference_bgc_summary.tsv`](../results/tables/reference_bgc_summary.tsv), [`lrcc_selected.tsv`](../data/processed/lrcc_selected.tsv), [`lrcc_pairwise_identity.tsv`](../results/tables/lrcc_pairwise_identity.tsv), and [`alignment_qc.tsv`](../results/tables/alignment_qc.tsv).

## Concrete patterns present in the data

These are observations to check, not final explanations.

1. At the broadest level, Bacillota, Pseudomonadota and Actinomycetota occupy separate regions of the full LrcC tree. There is no obvious inter-phylum mixing in this small author-selected set.
2. `Paenibacillus sp. M2` and `Paenibacillus sp. FSL R5-0765` form a 100/100 SH-aLRT/UFBoot pair and share 92.216% amino-acid identity. `Bacillus glycinifermentans` and `Bacillus subtilis` also form a 100/100 pair and share 94.578% identity.
3. The named Paenibacillus sequences do not form one exclusive group: different Paenibacillus pairs are separated by Bacillus, Brevibacillus or Lysinibacillus branches. A real host tree is required before treating this as evolutionary discordance.
4. `Photorhabdus asymbiotica DSM 15149` forms an exact four-tip clade with `Burkholderia sp. IMCC1007`, `Ralstonia pickettii DTP0602` and `Paraburkholderia mimosarum STM 3621`: support is 95.7/92 in the full tree and 94.3/86 in the sensitivity tree. This is a taxonomically mixed Burkholderiales/Enterobacterales pattern worth checking later.
5. Other Photorhabdus and Xenorhabdus sequences are interleaved. For example, `Photorhabdus sp. RW14` and `Xenorhabdus bovienii MC_235` share 86.364% identity. This does not by itself distinguish HGT from ordinary gene/species-tree differences.
6. The four short Actinomycetota sequences group together with 100/100 support at the base of that group, but some internal placement is weak and the long branches/gaps make direct comparison to full-length systems uncertain. `Streptomyces sp. B15` and `S. lunaelactis MM31` form a 100/100 pair.
7. The architecture examples show that M2 is unusually easy to label. Many other exported records retain fewer functional annotations; absence of a coloured precursor/RRE/B2 in the figure is not yet evidence of a noncanonical biological architecture.

## Questions for the user's independent analysis

1. Which relationships look biologically interesting to you before reading an HGT interpretation?
2. Does the `P. asymbiotica` + three Burkholderiales clade deserve first priority for a host-tree/context check? What alternative explanations would you test first?
3. Do you consider the Photorhabdus/Xenorhabdus interleaving a separate candidate, or part of normal uncertainty within related Enterobacterales?
4. Should the short Actinomycetota systems be treated as a separate search branch rather than mixed with the full-length LrcC family?
5. Which 2–4 reference sequences would you choose as seeds for the expanded search, and why?

Please record your reading without changing thresholds or exclusions yet. After that, the next iteration can freeze the search definition and decide whether to analyse the Actinomycetota branch separately.
