# Gao et al. 2021 - Candida auris piggyBac (PRJNA679970)

piggyBac mutant pools of *C. auris*; Gao J et al. (2021) Nat Microbiol 6:842-851, doi:10.1038/s41564-021-00915-x
(the study that identified the lncRNA DINOR).

Six runs, deposited as `Cau01`-`Cau06` with sample titles `Cau0h11`, `Cau0h12`, `Cau0h21`, `Cau0h22`, `Cau0h51`,
`Cau0h52`. **What the deposit does not say, anywhere:** what those six are. The SRA study abstract gives only
"piggyBac transposon-mediated mutagenesis and genome-wide genetic screens in C. auris"; the experiment titles repeat
the sample names; the BioSample records add nothing; the run descriptions say "Candida auris PB sample Cau0N". The
paper (Gao et al. 2021, Nat Microbiol) is subscription-only, has no PubMed Central version and no preprint, so its
Methods could not be read when this dataset was written. **It should be read before these libraries are used for
anything beyond QC**, and this README updated.

Two things follow from that, and both are handled without guessing:

1. **Selection.** Every run is treated as an *unselected* pool (`library_role` parent, `condition` parent): QC,
   density and essentiality only, no condition comparison. The `0h` in the sample titles is consistent with time
   zero, but that is inference, not evidence. If the paper shows any of the six were selected under a stress, they
   move to `treated` and this dataset gains a comparison.
2. **Grouping.** The names group as {1,2,5} x {1,2}, which is equally consistent with three pools sequenced twice
   (so: three libraries, to be merged with `ytab merge-runs`) and with six independent pools (six libraries). This
   changes the per-library density that R26 H1 and H2 use, so the rule is fixed in advance rather than chosen after
   the numbers are seen: **if the two runs of a pair share insertion sites at Jaccard >= 0.9, they are runs of one
   library and are merged; below that they are counted as independent libraries.** 0.9 is the value that separates
   the known technical pairs from the known independent pools among the 174 libraries already analysed. Whichever
   way it falls is reported.

Read structure: 150 bp paired end, R1 carries the junction; it begins inside the piggyBac terminal repeat
(`CATGCGTCAATTTTACGCAGACTATCTTTCTAGGG`), then the TTAA target site and the genome (`piggybac` / `pb_tr_end`; 63% of
the first 2,000 reads of SRR13106162 carry it, 68% of the trimmed reads map at MAPQ >= 20).

Reference: `resources/species/auris` (GCF_002759435.1 Cand_auris_B8441_V3, clade I strain B8441; 5,424 protein-coding
genes). The isolate behind these libraries is not stated in the deposit; the mapping rate is the check.
