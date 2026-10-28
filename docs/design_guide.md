# Choosing the analysis for your experimental design

Each analysis in YTAB answers one question and needs a minimum design to answer it. The table gives each analysis's
requirements as the code enforces them, and which libraries to give it. `ytab design-check --project-config
P/config/project.yaml` applies these rules to your own sample sheet: it lists every analysis your design supports, with
the libraries to use, and says why the others are not possible.

Terms: a **pool** is one transposon induction (one set of mutants); a **library** is what was sequenced from one culture
of a pool; **parent** libraries are the pool itself, **control** libraries are untreated cultures grown alongside the
treated ones, and **treated** libraries are cultures under the condition. Set them in the sample sheet columns
`library_role` (parent, control, treated), `condition`, `background` (the strain) and `pool`.

## Requirements

| Question | Analysis | Minimum design | Libraries to use | Command |
|---|---|---|---|---|
| Which genes are essential? | Essentiality, Hermes window HMM (label-free) | ≥ 1 unselected library of the strain; ≥ ≈ 36 insertion sites per coding kb for calls close to full density (Hermes; Ac/Ds libraries needed their full density) | parent libraries; pool independent inductions to reach the density | `ytab essentiality --libraries ...` |
| | Hybrid (HMM + forest) | as above, ≈ 18 sites per coding kb, and curated labels for the species | parents, labels | `ytab essentiality --labels ...` |
| Which genes change under a treatment? | Resampling (BGC) | ≥ 1 control and ≥ 1 treated library of the strain and condition, **from the same pool** | treated cultures and the untreated (or parent) cultures of their pool | `ytab comparison-design --controls all`, then `ytab treated-vs-control --method resampling` |
| | Local-abundance z-score (Gale et al. 2020) | as resampling, and ≥ 2 controls in the strain: its noise model compares two controls | as above | `ytab fitness` |
| | Control-conditioned Gumbel | ≥ 2 controls; depletion only | as above | `ytab treated-vs-control --method gumbel` |
| How does the response change across doses, times or several conditions? | Multi-condition (permutation ANOVA + all-pairs) | ≥ 2 conditions in each pool, one of them the reference (untreated or time 0) | every condition grown from the pool | `ytab multicondition-project --reference-condition ...` |
| Do two strains differ? | Multi-condition, strains as conditions | the strains in the same condition; `--independent-libraries` when each strain is its own library (no shared mutants) | one library per strain or more | `ytab multicondition-project --independent-libraries` |
| Does a mutation change the response to a treatment? | Genetic interactions (strain × treatment) | a reference strain (e.g. wild type) and a query strain (e.g. a deletion), **each with ≥ 2 control/treated pairs in the same condition**: at least 2 × 2 × 2 libraries. A pair is a control and a treated culture of one pool (two pools, or replicate cultures of one pool) | control and treated of every pool of both strains | `ytab genetic-interactions --reference-strain WT --query-strain mutant` |
| How similar are my libraries? | Library correlation | ≥ 2 libraries | all | `ytab library-correlation` |

Before any of these: **runs of one library** (lanes, sequencing runs, restriction digests of one culture) are not
replicates. Merge them with `ytab merge-runs`; `design-check` flags runs that share strain, pool, condition and culture
replicate.

## Selection markers leave blind spots

The scheme that enriches a pool can make genes look essential that are not, because the mutants that lost them were
removed before sequencing. These are not statistical artefacts and no normalization repairs them; the genes must be
masked before essentiality is read.

| Scheme | Blind spot | Seen in |
|---|---|---|
| Hermes launched from a plasmid carrying wild-type *RPL28*, enriched with nourseothricin + **cycloheximide** (prototrophic *C. glabrata*, Zaeske et al. 2026) | **PDR1** and **CDR1**, both required for cycloheximide resistance, are called essential | every CHX-enriched library, in three strain backgrounds; neither gene is called essential in any FOA-enriched library of the same species |
| Hermes launched from a *URA3* plasmid in a *ura3* auxotroph, enriched with nourseothricin + **5-FOA** (Gale et al. 2020 and after) | no gene blind spot, but every mutant is *ura3* and so is outcompeted in nutrient-limited or host environments | the uracil auxotrophy is carried by the whole pool |
| Ac/Ds SATAY from a donor plasmid | the donor and marker loci pile up reads (up to 69% of a library at one site) | every SATAY library (R25) |

Before trusting an essentiality call, ask what the pool was selected on. `ytab design-check` reports the design; the
marker is in the protocol, not in the data.

## Design rules that follow from the reanalysis

- **Grow control and treated cultures from the same pool.** The shared mutants carry the test's power; controls from
  another pool share no mutants, and only much denser libraries detect effects then.
- **Compare strains or conditions at the same time point,** not time points: selection acts before the first one.
- **One control per drug series is enough for resampling, not for the z-score or GI.** A second control (or a
  replicate control culture) opens the z-score and, with a second treated culture per strain, genetic interactions.
- **Density before depth.** Most libraries already carry ≥ 10 reads per site; more sequencing adds few sites. Pool
  independent inductions instead.

## Examples from the public datasets

| Study | Design | Supports | Not possible, and why |
|---|---|---|---|
| Lee et al. 2020 (S. pombe) | WT and otr1R::ura4, 2 pools each, T0 control and T80 | resampling, z-score, multi-condition, strain comparison, **GI** (2 pairs per strain) | — |
| Gale et al. 2020 (C. glabrata) | 3 parents; one pool split into 2 no-drug and 2 fluconazole cultures | essentiality, resampling, z-score, multi-condition | GI (one strain) |
| Nickels et al. 2024, Nickels & Cunningham 2025 (C. glabrata) | one no-drug control per drug, a dose series per strain | resampling per dose, multi-condition per series, strain comparison | z-score in strains with one control; GI (one pair per strain) |
| Billmyre et al. 2022 (S. pombe) | one pool; 2 libraries each of glycerol, mitotic, short and long meiotic outgrowth; each library sequenced twice | merge runs; resampling, z-score (2 mitotic controls), multi-condition | GI (one strain) |
| Michel et al. 2017 (S. cerevisiae SATAY) | 2 wild-type libraries, rapamycin from one; mutant libraries made separately | resampling (rapamycin), strain comparison of the mutants (`--independent-libraries`) | z-score (one control); GI (no treatment of the mutants) |
| Kingma et al. 2025 (S. cerevisiae SATAY) | wild type and nrp1Δ libraries, each digested and sequenced several times | merge runs; essentiality, strain comparison | treatment comparisons (none) |

Every study's libraries and the analyses run on them are listed in *Experimental designs of the public datasets*
(`ytab datasets --find TEXT` from the command line).
