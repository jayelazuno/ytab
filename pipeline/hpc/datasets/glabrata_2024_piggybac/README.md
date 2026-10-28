# Chow et al. 2024 - Candida glabrata piggyBac (PRJNA967833)

piggyBac insertion mutants of *C. glabrata* (library from strain YWCg03); Chow EWL et al. (2024) Mol Microbiol,
doi:10.1111/mmi.15229.

The value of this dataset to YTAB is that it is **piggyBac in a species YTAB already has deep Hermes coverage of**, so
transposon chemistry can be compared with the host species held constant (R26 H1).

Three pools (`R1`, `R2`, `R3`), each sampled at deposit levels 0, 60 and 120. The deposit gives neither the treatment
nor the units; `condition` therefore records the deposit labels (`level_0` as control, `level_60`, `level_120`) and
they must be resolved from the paper before any biological statement. Pool structure supports resampling, the
z-score (3 controls) and the multi-condition series.

Read structure: 150 bp paired end, R1 carries the junction (`piggybac` / `pb_tr_end`; 99% of the first 2,000 reads of
SRR24456510 begin with the terminal repeat, 79% of the trimmed reads map at MAPQ >= 20 to the BG2 reference).

Reference: `glabrata_CBS138`. YWCg03's background is not stated in the deposit; the first QC step compares the
mapping rate against both the CBS138 and the BG2 references and the choice is recorded before any analysis.
