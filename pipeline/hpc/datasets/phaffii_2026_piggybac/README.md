# Yao et al. 2026 - Komagataella phaffii piggyBac (PRJNA1376236)

piggyBac insertion library of *K. phaffii* (Pichia pastoris) and two sorted sub-populations; Yao X et al. (2026)
J Fungi 12:136, doi:10.3390/jof12020136, "Genome-Wide Landscape of Position Effect Variegation in Pichia pastoris"
(open access).

**The design comes from the paper.** The deposit labels the three runs `G15P`, `G15P-P4` and `G15P-P67` and its
experiment titles say only "sequence of PB transposon in Pichia pastoris". The Methods state what they are:

- `G15P` - the unsorted PB insertion library, made in strain P.pg0015 (GS115 carrying a PBL-eGFP-BleoR-PBR
  cassette) after methanol-induced transposition. Role parent.
- `G15P-P4` - the population sorted for **low** eGFP expression.
- `G15P-P67` - the population sorted for **high** eGFP expression.

They are not passages 4 and 67, which the names suggest.

**This is a reporter-expression sort, not a fitness selection.** The two sorted populations differ by how strongly a
reporter is expressed at the insertion site (position-effect variegation), not by growth. They are therefore used
only as a descriptive contrast and are kept out of every fitness hypothesis; only the unsorted `G15P` library feeds
density and essentiality.

Read structure: 250 bp paired end, R1 carries the junction (`piggybac` / `pb_tr_end`; 94% of the first 2,000 reads
of SRR36343507 begin with the terminal repeat). **Only 33% of the trimmed reads mapped at MAPQ >= 20**, the lowest
of the piggyBac sets, so its QC is read before the library is used for anything.

Reference: `resources/species/phaffii` (GCF_000027005.1 ASM2700v1, GS115; 5,040 protein-coding genes).
