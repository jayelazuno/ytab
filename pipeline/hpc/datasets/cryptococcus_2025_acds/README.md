# Billmyre et al. 2025 - Cryptococcus neoformans Ds (PRJNA1134100)

Saturated mutagenesis of *C. neoformans* H99 (strain yBB119) with a Ds element mobilised by maize Ac transposase;
Billmyre RB et al. (2025) PLoS Biol 23:e3003184, doi:10.1371/journal.pbio.3003184. The first basidiomycete in YTAB.

One pool (yBB119): the library at time zero (`T0`, role parent), then two days in DMSO (`DMSO`, control) or at the
fluconazole IC50 (`fluconazole`, treated). Each culture was split into a CviAII and an MseI digest and sequenced
separately: those are runs of one library, merged with `ytab merge-runs` before any analysis, not replicates.

Read structure: 151 bp paired end, R1 carries the junction; it begins inside the Ds terminal repeat
(`...CGGTACGGGATTTTCCCATCCTACTTTCATCCCTG`) and runs into the genome (`ac_ds` / `ds_tir_restriction`; 99% of the first
2,000 reads of SRR29762408 carry it, 71% of the trimmed reads map at MAPQ >= 20).

Reference: `resources/species/neoformans` (GCF_000149245.1 CNA3, H99; 14 chromosomes + mitochondrion, 6,975
protein-coding genes).
