# Song et al. 2026 - Candida auris piggyBac, fluconazole series (PRJNA1273510)

piggyBac pools of *C. auris* grown without drug and at three fluconazole levels; Song Y et al. (2026) Nat Microbiol,
doi:10.1038/s41564-026-02270-1 (vacuolar calcium pump, fluconazole efflux and resistance evolution).

Three pools (`n1`, `n2`, `n3`), each with an untreated control (deposit `Cau0n*`) and three fluconazole levels
(`Cau30n*`, `Cau60n*`, `Cau120n*`). The deposit labels the levels 30, 60 and 120 without units; the units must be
read from the paper before any dose statement. Three controls in three pools make the z-score, resampling and the
multi-condition dose test all available.

Read structure as `auris_2021_piggybac` (`piggybac` / `pb_tr_end`; 70% of the first 2,000 reads of SRR33878342 carry
the terminal repeat, 86% of the trimmed reads map at MAPQ >= 20).

Reference: `resources/species/auris` (GCF_002759435.1, clade I strain B8441).
