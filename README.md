# YTAB — Yeast Transposon Analysis Browser

YTAB takes transposon-sequencing (Tn-seq) data from raw reads to insertion maps, quality control, normalization,
condition comparisons and essentiality calls, in a browser interface backed by a command-line pipeline. It is built
for the transposons used in yeasts and yeast-like fungi — **Hermes**, **Ac/Ds** (including SATAY) and **piggyBac** —
rather than the TA-site systems most Tn-seq software assumes, and ships reference resources for *C. glabrata*,
*S. cerevisiae*, *S. pombe*, *C. albicans*, *C. neoformans*, *C. auris*, *K. phaffii* and *Y. lipolytica*.

Its statistics are adapted from **TRANSIT** (DeJesus et al. 2015), built for Himar1 TA-site data, and from the
**Hermes Tn-seq literature** (Gale et al. 2020; Segal et al. 2018; Michel et al. 2017). Each adaptation names its
source and states what was kept and what was changed — see [docs/references.md](docs/references.md).

Every analysis in the app is a YTAB command. The app only reads the tables the commands write, so anything seen on
screen can be reproduced, and the same analysis gives the same result whether it is started from the command line or
from the interface.

## Try the app

With [Docker](https://docs.docker.com/get-started/get-docker/), nothing else to install:

```bash
git clone --depth 1 https://github.com/jayelazuno/ytab.git
cd ytab
docker run --rm -p 3838:3838 -v "$PWD:/opt/ytab" ghcr.io/jayelazuno/ytab-app:v0.1
```

Then open <http://localhost:3838>. The clone carries demo projects, so there is data to explore immediately.
`docs/testing.md` is the short guide for anyone evaluating the interface.

## Install locally

```bash
./scripts/local/ytab_setup_env.sh          # creates the ytab-local environment from environment.local.yml
mamba activate ytab-local
./scripts/local/ytab_launch_app.sh         # or: --project-config output/projects/<PROJECT>/config/project.yaml
```

`ytab_setup_env.sh` uses `mamba` when available, installs it into an existing conda base, or bootstraps a repo-local
micromamba. Add `--launch-app` to set up and start in one step. The app binds to `127.0.0.1:3838`; another host or
port has to be asked for explicitly. Run R and Python from the same `ytab-local` environment, since the app calls
Python for several actions. Longer notes: [docs/setup_local.md](docs/setup_local.md).

## What the app does

| Tab | |
|---|---|
| **Dataset** | projects, samples, the comparison design |
| **Run** | preprocessing from FASTQ: trimming, mapping, hit files, gene summaries, with resumable cached stages |
| **Quality Control** | a flagged per-library table — reads and mapping (including reads from the transposon donor), insertions, genes, saturation — with the plot behind every metric; per-library statistics after TRANSIT's `tnseq_stats` (DeJesus et al. 2015), complexity and jackpot measures after Gale et al. (2020) |
| **Normalization** | every site-normalization method, before-and-after views, and the optional top-site treatment |
| **Condition Comparisons** | treated vs control, multi-condition dose series and genetic interactions |
| **Essentiality** | the label-free window HMM, and the hybrid model where curated labels exist |
| **Gene Explorer** | insertion profiles, domains, the genome browser, GO enrichment and cross-species views |
| **Screen Planning** | how large a screen has to be, and how to read coverage |

Preprocessing runs one sample at a time, can be cancelled and resumed, and skips completed stages. FASTQs stay where
they are and are never uploaded.

## Methods and defaults

- **Normalization.** BGC (beta-geometric correction) is the default: it corrects the heavy-tailed "jackpot" sites
  typical of these transposons and ranked best in the pre-registered benchmarks. TTR, total reads, nonzero mean,
  adaptive BGC, quantile, ZINB and MidLC are all available and measured. BGC and adaptive BGC are adapted from
  DeJesus & Ioerger (2016); the scalar and distribution methods from TRANSIT (DeJesus et al. 2015); MidLC from the
  half-saturation measure of Gale et al. (2020).
  See [docs/normalization_methods.md](docs/normalization_methods.md).
- **Treated vs control.** Site-restricted **resampling** (BGC), adapted from TRANSIT's resampling test
  (DeJesus et al. 2015), and the **control-conditioned Gumbel** test (TTR, depletion only), adapted from the Gumbel
  method of DeJesus et al. (2013), both report p-values with a Benjamini–Hochberg false discovery rate. The
  **local-abundance z-score** of Gale et al. (2020) is kept for ranking and for continuity with published screens; it
  has no error control.
  See [docs/methods/06_condition_comparisons.md](docs/methods/06_condition_comparisons.md).
- **Several conditions.** A permutation ANOVA with all pairwise comparisons for dose series and mutant panels,
  adapted from TRANSIT's multi-condition tests (DeJesus et al. 2015), and a genetic-interaction test for
  strain × treatment designs, adapted from TRANSIT GI (DeJesus et al. 2017).
- **Essentiality.** A window HMM that needs no training labels, adapted from the TRANSIT HMM
  (DeJesus & Ioerger 2013), with a hybrid model when a species has curated ones, whose gene features follow
  Segal et al. (2018).
- **Which to use when.** [docs/methods/00_recommendations.md](docs/methods/00_recommendations.md) maps a QC symptom
  or a design to the method that fits it, with the evidence for each.

Every adapted method names its source in its own documentation and code, down to the inspected TRANSIT commit;
[docs/references.md](docs/references.md) is the full table of what came from where, with the complete citations.
Claims rest on pre-registered benchmarks; the registry, including the tests that failed, is in `docs/manuscript/`.

## Command line

```bash
python scripts/local/ytab.py modules          # what is available
python scripts/local/ytab.py <module> --help  # one analysis
python scripts/local/ytab.py run --workflow WORKFLOW.yaml
```

Every module and option: [docs/cli_reference.md](docs/cli_reference.md), with worked examples in
[docs/command_line.md](docs/command_line.md). To inspect resumable state without running anything:

```bash
python scripts/local/ytab_project_status.py --project-config <project.yaml> --show-next
```

## On a cluster

YTAB is meant for ordinary desktop use; a cluster is optional. The [Argon stage jobs and submitter](pipeline/hpc/README.md)
are SGE scripts for one institutional installation rather than a general scheduler interface, and they call the same
CLI as the app, so the portable commands can be wrapped for any scheduler.

## Repository layout

- `app/` — the Shiny interface
- `src/ytab/` — the analysis library (mapping, QC, normalization, fitness, essentiality)
- `scripts/local/` — the command-line entry points
- `pipeline/hpc/` — cluster job scripts and public-dataset definitions
- `resources/` — references, annotations, transposon and donor sequences, curated labels
- `containers/` — the app image
- `docs/` — methods, CLI reference, dataset catalogue and release notes
- `tests/` — unit and integration tests
- `output/` — project outputs; the demo projects live here

## Key references

The methods YTAB adapts, with the full list and what changed in each case in
[docs/references.md](docs/references.md):

- DeJesus MA, Ambadipudi C, Baker R, Sassetti C, Ioerger TR. TRANSIT — a software tool for Himar1 TnSeq analysis.
  *PLoS Comput Biol* 2015;11(10):e1004401. (Normalization, resampling, multi-condition tests, Tn-seq statistics.)
- DeJesus MA, Ioerger TR. Normalization of transposon-mutant library sequencing datasets to improve identification of
  conditionally essential genes. *J Bioinform Comput Biol* 2016;14(3):1642004. (BGC and adaptive BGC.)
- DeJesus MA, Ioerger TR. A Hidden Markov Model for identifying essential and growth-defect regions in bacterial
  genomes from transposon insertion sequencing data. *BMC Bioinformatics* 2013;14:303. (Window HMM.)
- DeJesus MA, Zhang YJ, Sassetti CM, Rubin EJ, Sacchettini JC, Ioerger TR. Bayesian analysis of gene essentiality
  based on sequencing of transposon insertion libraries. *Bioinformatics* 2013;29(6):695–703. (Gumbel.)
- DeJesus MA, Nambi S, Smith CM, Baker RE, Sassetti CM, Ioerger TR. Statistical analysis of genetic interactions in
  Tn-Seq data. *Nucleic Acids Res* 2017;45(11):e93. (Genetic interactions.)
- Gale AN, Sakhawala RM, Levitan A, Sharan R, Berman J, Timp W, Cunningham KW. Identification of essential genes and
  fluconazole susceptibility genes in *Candida glabrata* by profiling Hermes transposon insertions. *G3*
  2020;10:3859–3870. (Local-abundance z-score, MidLC complexity and depth correction.)
- Segal ES, Gritsenko V, Levitan A, Yadav B, Dror N, Steenwyk JL, et al., Berman J. Gene essentiality analyzed by in
  vivo transposon mutagenesis and machine learning in a stable haploid isolate of *Candida albicans*. *mBio*
  2018;9:e02048-18. (Gene features for the hybrid essentiality model.)
- Michel AH, Hatakeyama R, Kimmig P, Arter M, Peter M, Matos J, De Virgilio C, Kornmann B. Functional mapping of yeast
  genomes by saturated transposition. *eLife* 2017;6:e23570. (SATAY design and donor locus.)

## Licence

YTAB is released under the GNU General Public License v3.0 (see `LICENSE`), the licence of TRANSIT
(DeJesus et al. 2015), from which several of its methods are adapted.
