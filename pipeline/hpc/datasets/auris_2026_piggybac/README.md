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

Reference: `resources/species/auris` (GCF_002759435.1, clade I strain B8441).
