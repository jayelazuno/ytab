# Gao et al. 2021 - Candida auris piggyBac (PRJNA679970)

piggyBac insertional libraries of *C. auris*; Gao J et al. (2021) Nat Microbiol 6:842-851,
doi:10.1038/s41564-021-00915-x (the study that identified the lncRNA DINOR).

**The design comes from the paper, not the deposit.** The deposit names the runs only `Cau0h11`-`Cau0h52` and
describes none of them; the Methods state it. Everything below is from the paper.

- **Three independent PB insertional libraries**, each sequenced twice: the six runs are three libraries
  (`lib1`, `lib2`, `lib5`, after the deposit labels), and the two runs of each are merged with `ytab merge-runs`
  before any analysis. Counting them as six libraries would halve the apparent density of each.
- **Unselected.** All three are the mutant libraries themselves, with no drug or stress selection, so this dataset
  is QC, density and essentiality only. The genetic-interaction screen in the same paper was a colony screen with
  inverse PCR, not sequenced pools, so none of its material is here.
- **Strain CauW08**, derived from CauW49, a leucine auxotroph of **CBS10913 (clade II)**, the first *C. auris*
  isolate from Japan.
- **Library preparation** (PBISeq): DNA from 4-5 x 10^8 transposed cells, Nextera tagmentation to 200-600 bp, then
  two PCR rounds with a PB-specific primer and the Tn5 i7 adaptor. That is why Nextera read-through is trimmed.

Read structure: 150 bp paired end, R1 carries the junction, `piggybac` / `pb_tr_end`. The paper's own analysis
keeps "reads that have a prefix matching the terminus of the PB transposon (TTTCTAGGG)" and removes it, which is
the sequence this profile trims.

**Reference: `auris_B11220`** (GCF_003013715.1, strain B11220, clade II), which matches the documented clade of
CBS10913. The authors mapped to the clade I B8441 assembly (GCA_002759435.2) instead, so their published
coordinates are B8441-based; on the first 2,000 reads of SRR13106162 the two references gave 65% (B11220) and 64%
(B8441) of trimmed reads at MAPQ >= 20, so the choice barely changes yield and is made on the strain, not the rate.
B8441 remains available as species `auris`.

**Published numbers to reproduce** (R26, descriptive): an average of 77 unique PB insertion sites per 10 kb;
insertions in 4,503 of 5,417 genes (83%) across the three libraries, 3,568 of them in all three, leaving 914 genes
(17%) without an insertion. Note the authors counted insertions at one position in opposite orientations as two
independent insertions, whereas `ytab hit-files` unites adjacent positions, so site counts are not expected to
match exactly; the gene-level numbers are the comparable ones.
