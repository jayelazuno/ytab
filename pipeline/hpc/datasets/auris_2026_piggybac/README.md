# Song et al. 2026 - Candida auris piggyBac, fluconazole series (PRJNA1273510)

piggyBac pools of *C. auris* grown without drug and at three fluconazole concentrations; Song Y et al. (2026)
Nat Microbiol, doi:10.1038/s41564-026-02270-1 (vacuolar calcium pump, fluconazole efflux and resistance evolution).

Three replicates (`n1`, `n2`, `n3`), each at 0 (control), 30, 60 and 120 ug/mL fluconazole. The doses come from the
**SRA experiment titles** (`<n> g/mL fluconazole treatment replicate <r>`, where `g/mL` is ug/mL); the BioSample
records give only the `Cau30n1`-style labels. Whether the three replicates are independent pools or cultures split
from one is not stated, so each is treated as its own pool, which is the conservative reading: a treated library is
compared only with the control of its own replicate. Three controls make the z-score, resampling and the
multi-condition dose test all available.

Read structure as `auris_2021_piggybac` (`piggybac` / `pb_tr_end`; 70% of the first 2,000 reads of SRR33878342 carry
the terminal repeat, 86% of the trimmed reads map at MAPQ >= 20).

Reference: `auris_B11220` (GCF_003013715.1, clade II). The deposit does not name the strain, but Yabing Song is an
author of Gao et al. 2021, whose libraries are documented as CBS10913 (clade II), and this deposit repeats that
study's abstract wording and PB system; on the first 2,000 reads of SRR33878342 the two references gave 86%
(B11220) and 85% (B8441) of trimmed reads at MAPQ >= 20. The strain should be confirmed from the paper, and the
mapping rates of the full libraries against both references are recorded in QC before the data are used.
