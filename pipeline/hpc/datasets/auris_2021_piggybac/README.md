# Gao et al. 2021 - Candida auris piggyBac (PRJNA679970)

piggyBac mutant pools of *C. auris*; Gao J et al. (2021) Nat Microbiol 6:842-851, doi:10.1038/s41564-021-00915-x
(the study that identified the lncRNA DINOR).

Six runs in three pools of two (`h1`, `h2`, `h5`). The deposit gives no treatment or time metadata, so every run is
analysed as an unselected pool: essentiality and QC only, no condition comparison. The pool structure must be
confirmed against the paper before any claim that the pools are independent.

Read structure: 150 bp paired end, R1 carries the junction; it begins inside the piggyBac terminal repeat
(`CATGCGTCAATTTTACGCAGACTATCTTTCTAGGG`), then the TTAA target site and the genome (`piggybac` / `pb_tr_end`; 63% of
the first 2,000 reads of SRR13106162 carry it, 68% of the trimmed reads map at MAPQ >= 20).

Reference: `resources/species/auris` (GCF_002759435.1 Cand_auris_B8441_V3, clade I strain B8441; 5,424 protein-coding
genes). The isolate behind these libraries is not stated in the deposit; the mapping rate is the check.
