# Recommendations: library density, normalization and tests

These recommendations come from reprocessing 214 published yeast Tn-seq libraries (nine species and strain
backgrounds; Hermes, Ac/Ds and piggyBac) with
YTAB and testing each method against criteria fixed in advance. Each one names the evidence behind it: a pre-registered
test (R-number; its design, checksum and result are in the benchmark's folder under `output/`) or an analysis of
the reprocessed data. The reprocessed projects
(`Gale_2020`, `Segal_2018_albicans`, `Michel_2017_SATAY`, `Pombe_Hermes`) open in this app, so every example below can
be looked at tab by tab.

## Before the screen: how dense must a library be?

Density is measured as **insertion sites per coding kb** (median gene). The QC table reports it for every library as
`sites_per_coding_kb`, which is not the same as its `insertions_per_kb`: that one divides by the whole reference
rather than by coding sequence, so it reads several times higher and must not be compared with the thresholds below.
Because an essentiality run is fitted on a strain's unselected libraries together, the figure that decides whether a
call is worth making is their pooled density, and `ytab qc` writes that per strain to `essentiality_readiness.tsv`
(the Essentiality tab warns from it when a project falls short).

| analysis | recommendation | evidence |
|---|---|---|
| Essentiality without labels (HMM) | Pool the parent libraries to at least **≈ 36** sites per coding kb. Hermes libraries were within 0.02 AUC of their full density at 11–18; Ac/Ds libraries needed their full 15–36 (no plateau, so 36 is a lower bound). | thinning four species, R16 |
| Essentiality with labels (species-trained forest or hybrid) | About half: **≈ 18**. | R16 |
| Condition comparisons | Per library, at least the density of the published screens: **≈ 4.4** with two replicates per condition (Hermes, Gale fluconazole) and **≈ 12** with one library per condition (Ac/Ds, SATAY rapamycin). A pool with half the mutants lost a fifth of the published hits and half of all calls: denser is better. | thinning by pool, R23 |
| Is the HMM fit usable? | It must converge (`fit.json`: `converged`). Fits converged in every species at ≥ 18 pooled sites per coding kb; sparser fits often stopped at the iteration cap. | R16 |

The screen planner (Screen Planning tab) gives the number of colonies needed for a target fraction of genes hit.

### What density each transposon reaches

Medians over the 93 unselected libraries of the reprocessed studies, so these are the densities published screens
actually achieved rather than a theoretical ceiling. Read them as the starting point when planning a screen in a
species where only one chemistry has been tried.

| transposon | libraries | median sites per coding kb | range | species reprocessed |
|---|---|---|---|---|
| Ac/Ds | 28 | **19.86** | 0.49–47.87 | *S. cerevisiae* 36.30, *C. neoformans* 18.67, *C. albicans* 1.02 |
| Hermes | 55 | **6.46** | 0.00–46.56 | *C. glabrata* CBS138 9.55, *C. glabrata* BG2 7.43, *S. pombe* 6.57, *Y. lipolytica* 1.09 |
| piggyBac | 10 | **1.32** | 0.00–2.03 | *C. glabrata* CBS138 1.97, *C. auris* 1.24, *K. phaffii* 0.00 |

Two things follow for planning. **piggyBac libraries are the sparsest of the three**, and the gap is not only an
artefact of its TTAA restriction: no piggyBac library in these studies reached a quarter of the ≈ 36 that label-free
essentiality needs, and the densest reached 2.03. A piggyBac screen aimed at essentiality therefore needs either far
deeper pooling than published practice or labelled training data (the hybrid route, ≈ 18). **Density varies more
within a transposon than between species**, so a single published library is weak evidence about what a chemistry
can reach in a new organism.

Note that piggyBac densities are measured after off-motif positions are dropped, because piggyBac inserts only at
TTAA; Hermes and Ac/Ds declare no target site and are not filtered. The piggyBac ceiling is correspondingly lower
(8.25 insertable sites per coding kb in *C. glabrata*, against every position for the hAT elements), so as a share
of what the chemistry can occupy, the best piggyBac libraries sit near 31% — not far below where Hermes libraries
sit against their own ceiling. The *K. phaffii* libraries hold too few positions (15,008 in all three) for a
density estimate and are reported as such rather than as a measured 0.00.

## Experimental design

- **Grow control and treated cultures from the same pool.** Their shared mutants carry most of the test's power:
  thinning each library independently, which mimics libraries from different pools, called almost nothing at half
  density (R22), whereas thinning the pool kept 14.5 of 18 published hits (R23).
- **Replicates:** the z-score needs at least two control libraries per background. With one library per condition
  (as in SATAY), resampling still works but with less power.
- **Match the comparison to where selection acts.** If selection acts before the first time point (for example, a
  reporter selected when the donor plasmid is removed), compare strains or conditions at the same time point, not time
  points (S. pombe heterochromatin screen; R17).

## QC finding → normalization

| QC finding (Quality Control tab) | normalization | evidence |
|---|---|---|
| Libraries differ only in depth; no flag | Any scalar method (TTR, total reads). BGC leaves such libraries unchanged. | clean replicate pair, rank correlation 0.998 under every method |
| Jackpot: one or a few sites hold a large share of reads (top-site share or jackpot reads flagged) | **TTR or BGC, not total reads.** In S. pombe, a site with 42% of a library's reads shifted it from its replicate by 0.63 log2 with total reads, 0.04 with TTR and 0.01 with BGC. | replicate pairs |
| Heavy-tailed or skewed counts (top-1% share or tail index flagged) | **BGC.** It raised the rank agreement between replicates in every biased pair (0.67 → 0.73 in C. albicans, 0.86 → 0.89 in SATAY). | replicate pairs |
| Skewness above 30–50 only | Not a reason to drop a Hermes or Ac/Ds library: the rule comes from Himar1, and 42 of 50 libraries exceed it. | 50-library QC |
| Low mapping rate | Check the reference, adapters and contamination; normalization cannot fix it. | six S. pombe libraries at 44–52% |
| Donor-locus hotspot (Ac/Ds donor at ADE2; local re-insertion near a donor) | Mask the donor locus and marker genes before essentiality calls (done for the public projects). | SATAY, C. albicans |

Every flag in the QC table carries the matching normalization hint.

## Which test?

| question | method | evidence |
|---|---|---|
| Treated vs control | **Resampling with BGC**, or the **control-conditioned Gumbel test with TTR**: both kept false discoveries near the nominal rate on Hermes-calibrated simulations (aBGC, the adaptive form, performs within 0.01 AUC). The z-score is for ranking only (false discovery proportion 0.20). Never total reads under one-sided selection: it turned rapamycin's effect into 504 apparently depleted genes, against 188 with BGC. | R13, R17, R27 |
| Several conditions or mutant backgrounds | Multi-condition ANOVA with all pairwise comparisons, BGC (aBGC passes the same gates). | R11, R28 |
| Strain × treatment interactions | Genetic-interaction test. | R12a |
| Essentiality, no labels | The label-free HMM with default settings: windows of 50–100 bp give the same ranking, 200 bp is worse; the essential run is at least 3 windows or a fifth of the gene. | R14, R19–R21 |
| Essentiality, labels available | A species-trained random forest, or the hybrid (the forest with the HMM's posteriors as features), which matched the forest at full density and was better in sparse libraries. Both are 0.04–0.06 AUC better than the label-free HMM, and even 50 labelled genes beat it. | R15, R16 |

## What has not worked

- TRANSIT's U-test and ANOVA hold their error rate on Hermes data but find almost nothing (8 and 3 of 19 published
  fluconazole genes); TRANSIT's HMM calls no gene essential.
- Permutation and quasi-likelihood variants that were anticonservative on Hermes-like simulations were not adopted
  (R1, R2, R4, R5).
- A time-course genetic-interaction test where selection preceded the first time point called no gene (R17).

Each pre-registered test's design (`PREREGISTRATION.md`, with its SHA-256 in `PREREGISTRATION.sha256`) and outcome
(`RESULT.md`) are in its benchmark folder under `output/`.
