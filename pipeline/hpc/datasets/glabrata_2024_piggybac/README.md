# Chow et al. 2024 - Candida glabrata piggyBac, fluconazole (PRJNA967833)

piggyBac insertion mutants of *C. glabrata* (library from strain YWCg03); Chow EWL et al. (2024) Mol Microbiol,
doi:10.1111/mmi.15229.

The value of this dataset to YTAB is that it is **piggyBac in a species YTAB already has deep Hermes coverage of**, so
transposon chemistry can be compared with the host species held constant (R26 H1).

Three replicates (`R1`, `R2`, `R3`), each at 0 (control), 60 and 120 ug/mL fluconazole. The BioSample records give
only `R1-0`-style labels; the doses come from the **SRA experiment titles** (`<n> g/mL fluconazole treatment
replicate <r>`, where `g/mL` is ug/mL). Three controls support the z-score, resampling and the dose series.

Read structure: 150 bp paired end, R1 carries the junction (`piggybac` / `pb_tr_end`; 99% of the first 2,000 reads of
SRR24456510 begin with the terminal repeat, 79% of the trimmed reads map at MAPQ >= 20 to the BG2 reference).

Reference: `glabrata_CBS138`. YWCg03's background is not stated in the deposit; the first QC step compares the
mapping rate against both the CBS138 and the BG2 references and the choice is recorded before any analysis.
