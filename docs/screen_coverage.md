# Screen Coverage

Screen Coverage separates observed library representation from colony-screen
recovery. It does not change other YTAB analyses and does not recommend a final
experimental screen size.

## Accepted representation

The mapping reference contains 5,220 protein-coding genes after excluding the
14 annotated pseudogenes. The accepted parent-sequencing union contains 169,099
distinct genomic insertion coordinates, 49,034 unambiguous coding coordinates,
and 4,379 protein-coding genes with at least one observed coding coordinate.
The remaining 841 reference protein-coding genes have no observed coding
insertion. Thus observed reference representation is 4,379/5,220 = 83.8889%.
This is a library-representation ceiling, not colony-screen recovery.

The accepted four-pool inventory is retained from the existing audit:

| Pool | Genomic coordinates (inventory) | Eligible coding coordinates C_j | Represented genes N_j |
|---|---:|---:|---:|
| 1 | 83,391 | 21,841 | 3,677 |
| 2 | 112,214 | 30,949 | 4,027 |
| 3 | 72,500 | 19,072 | 3,446 |
| 4 | 70,958 | 18,371 | 3,469 |

Pools 1–4 are independently constructed insertion libraries on the same
parental genetic background. Their overlap is overlap in insertion-coordinate
and gene content, not evidence of shared ancestry. Independent construction
does not imply equal insertion complexity. Repeated parent sequencing
observations are deduplicated within each pool, and only then combined into the
global coordinate union. Treated samples never enter the parent catalogue.

The 4,379 represented-gene denominator describes recovery relative to the
observed parent-union represented-gene set. Whole-reference recovery uses
5,220 as denominator and cannot exceed the observed representation ceiling.
Neither quantity should be called unqualified “genome coverage.”

## Primary scientific models

1. **Theoretical Poisson / Equal-Gene Benchmark.** For N=4,379,
   `C(K)=1-exp(-K/N)` and `K_q=ceil(-N ln(1-q))`. The 95% rule is
   `K_0.95=2.9957N≈3N`. This total-library benchmark assumes every gene has the
   same recovery probability; it does not allocate picks among pools or model
   insertion architecture or abundance.
2. **Hermes Empirical Monte Carlo — Equal Insertion.** This uses gene-hit
   occupancy over each independent pool's full observed insertion-coordinate
   catalogue `P_j`. `K` is total insertion-positive colony/mutant picks.
   Coding coordinates can recover protein-coding genes; noncoding/intergenic/
   unassigned insertion coordinates still consume a pick but do not recover a
   protein-coding gene.
3. **Hermes Empirical Monte Carlo — Observed-Abundance Sensitivity.** It uses
   the same full `P_j` sampling universe, replacing equal coordinate weights
   with parent/control read-abundance weights. Each sequencing observation is
   normalized by its own usable read depth; normalized coordinate abundance is
   averaged within pool and then normalized across all coordinates in `P_j`.
   Reads are a sensitivity proxy, not measured viable-mutant frequencies.

The analytical equal-gene Poisson curve is the benchmark; the primary empirical
simulation is site-to-gene occupancy from the observed catalogue, not a claim
to reproduce Baym et al.'s organism-specific experiment. The target selector
sets a target **expected recovery fraction** (the exact mean crossing); the
Monte Carlo interval shows replicate sampling variability and is not a
probability guarantee that a specified fraction will be reached.

For every empirical strategy, samples are drawn separately from each pool's
full catalogue and global recovery uses the union of recovered coordinates and
genes after sampling. A coordinate observed in multiple independent pools remains one unique genomic coordinate for allele counting,
although pool provenance can retain separate pool-derived occurrences.
Different coordinates in the same gene are distinct insertion-coordinate
alleles. Repeated recovery of the same coordinate is not another allele.

## Allocation strategies

Let `J` be the number of independent pools, `P_j` the full observed genomic
insertion-coordinate catalogue, `M_j=|P_j|`, `C_j⊂P_j` the coding/gene-assigned
subset, `N_j` the pool's represented genes, `m_gj` the coding coordinates
assigned to gene g, `K_j` the insertion-positive picks from pool j, and
`K=Σ_j K_j` the modeled budget. For Equal Insertion, each coordinate in `P_j`
has probability `1/M_j`, so a pick hits gene g with probability `m_gj/M_j`.
For abundance sensitivity, coordinate weights sum to one over all `P_j`; the
gene probability is the sum of weights for that gene's coding coordinates.
The implementation is parameterized by
J and builds controls and tables from the selected pool metadata.

**Complexity-normalized allocation** is the primary empirical strategy:

`w_j=M_j/Σ_l M_l`, `K_j≈K w_j`, and `λ_j=K_j/M_j`.

For the accepted full genomic-coordinate counts (83,391; 112,214; 72,500;
70,958; total 339,063), the calculated weights are approximately 24.60%,
33.10%, 21.38%, and 20.93%. Coding counts `C_j` remain reported as the gene-
recovery subset; coding fraction `C_j/M_j` describes how much of each pool's
observed insertion catalogue maps to eligible protein-coding genes.

Integer counts use Hamilton/largest-remainder apportionment: floor each quota,
then assign the remaining picks in descending fractional-remainder order;
pool order breaks exact ties. This approximately equalizes insertion-positive
picks per observed genomic coordinate. The weights are allocation weights only, not
abundance estimates or physical mixing proportions.

**Balanced saturation-aware allocation** uses each pool's own observed gene
probabilities and represented-gene repertoire `N_j`. For
`S_j(k)=E[G_j(k)]/N_j`, it maximizes `min_j S_j(K_j)` subject to integer
`K_j≥1` for every selected pool and `Σ_jK_j=K`. A binary search finds the
largest feasible common saturation threshold; each pool receives the minimum
integer count reaching that threshold. Any residual budget is apportioned by
largest remainder proportional to `M_j` (pool order breaks ties). This is a
max-min pool-saturation objective, not a global-gene-count objective.

**Equal-pool allocation** (`K_j≈K/J`) remains a transparent baseline. It
equalizes colony counts, not insertion depth or gene-saturation progress.
Manual counts are supported and are always summarized as `K=Σ_jK_j`. The
previous optimizer maximizing immediate expected global ≥1-gene recovery is
retained only under advanced allocation comparisons; it is not a default or
experimental recommendation.

## Mathematics shown in the UI

The Math & Equations tab displays and explains:

- Equal-gene Poisson: `λ=K/N`, `P(X≥1)=1-exp(-K/N)`,
  `E[G]=N[1-exp(-K/N)]`, and `K_q=-N ln(1-q)`; 95% gives the 3N rule.
- Full-catalogue complexity weights and normalized insertion-positive depth
  `K_j/M_j`.
- Whole-library gene recovery `1-(1-m_gj/M_j)^{K_j}` under equal insertion;
  abundance-weighted gene probability sums normalized weights over coding
  coordinates while retaining noncoding coordinates in the denominator.
- Global gene recovery after separate sampling,
  `P(g recovered)=1-Π_j P(g missed in pool j)`; simulation still uses the
  actual full insertion-coordinate catalogues and globally deduplicates them.
- The precise max-min saturation objective, equal allocation, the narrower
  immediate-gene optimizer, and the two distinct recovery denominators.

## Plots and planner

Screen Planner shows the theoretical benchmark and empirical curves under
complexity-normalized, balanced saturation-aware, and equal-pool allocations;
observed-abundance results remain visibly identified as sensitivity. The
practical planner reports expected represented-union genes, percent of 4,379,
percent of the 5,220 reference genes, missed represented genes, simulation
intervals, and recovery of genes with at least one, two, or three distinct
insertion coordinates. It shows the selected per-pool allocation and the
selected pool's expected recovery.

Pool Contribution / QC includes pool-specific gene saturation (denominator
`N_j`), insertion-coordinate saturation, normalized sampling depth, allocation
across K, and diminishing marginal gene recovery. The allocation view includes
equal, complexity-normalized and balanced strategies. Library QC, Parent
Sequencing, Statistical Details, and Methods & References preserve diagnostics
and explain that parent-sequencing discovery is distinct from colony-screen
saturation.

## Simulation and limitations

Equal-insertion and abundance-weighted draws are pool-stratified; pool
catalogues are never flattened before sampling. Both models draw total
insertion-positive picks from full `P_j`; noncoding/intergenic/unassigned
coordinates consume K but do not recover protein-coding genes. Gene length and
Hermes target-site preference are not fitted separately because realized
effects are reflected in the observed catalogue and `m_gj`. Coordinate metrics
globally deduplicate matching genomic coordinates. Fixed seeds make
simulation reproducible. Simulation intervals
describe sampling variability conditional on the model, not all biological or
technical uncertainty.

Parent sequencing is an observed catalogue, not proof of complete physical
library composition; nondetection does not establish absence. Sequencing reads
are imperfect proxies for viable-colony abundance and may reflect PCR, DNA
extraction, mapping, growth, depth and stochastic sampling. Hermes insertion
preferences and gene length can create unequal representation. A CDS insertion
does not guarantee complete loss of function. Individual colonies may carry
multiple insertions and linkage is not resolved by this model. These
limitations should accompany interpretation of every simulation interval.
CDS length, per-gene coding-coordinate counts (`m_gj`), coordinates per kb, and
number of represented pools are exposed as interpretive diagnostics; length and
Hermes +2/+7 target-site preference are not separately fitted because their
realized effects are reflected in the observed catalogue. Bottleneck estimates
remain QC context, not a separate model layer; current estimates are around
three or fewer reads in most pools, with some slightly higher. The primary empirical
planner does not retain a coding-conditioned-only K diagnostic: its K always
counts insertion-positive picks from full `P_j`.

## Validation

Focused tests cover arbitrary J, largest-remainder allocations, full-catalogue
probabilities and normalized depth, noncoding pick consumption, balanced saturation,
pool-specific denominators, global coordinate deduplication, manual allocation
totals, fixed seeds, accepted inventory, and the 5,220/4,379 denominator distinction. The R/Shiny
smoke test exercises the planner UI, plot/table outputs, dynamic pool controls,
worker execution, and isolation of cosmetic controls from scientific runs.

References: Baym et al. 2016, *Nature Communications* 7:13270,
doi:10.1038/ncomms13270 (single-library gene-hit occupancy precedent); Gale et
al. 2020, *G3* 10:3859–3870, doi:10.1534/g3.120.401595 (*C. glabrata* Hermes
insertion catalogue and biases). Mahmutovic et al. 2020,
*Computational and Structural Biotechnology Journal* 18:791–804,
doi:10.1016/j.csbj.2020.03.021, studies growth, fitness and bottleneck effects;
that process model is **not** implemented here. Conte et al. 2019 concerns
ISY100 delivery and is not a basis for the Hermes Monte Carlo.
