# Release notes

## 2026-11-10: one dataset ships, and a launcher to start it

**Downloading and starting YTAB**
- **`YTAB.command`** (macOS) and **`YTAB.bat`** (Windows): double-click to start. They find a container runtime,
  download the image the first time, pick a free port, start the app and open the browser; closing the window
  stops YTAB. Docker Desktop, Podman and Colima all work, and the launcher uses whichever is installed.
- YTAB is downloaded from the releases page as a zip, about 0.8 GB, so git is not needed.
- **Your own reads are reachable.** *Start from FASTQ › Browse local folders* opens on **My files**, your home
  folder, shared with YTAB read-only. `YTAB_FILES` shares a narrower folder instead.

**What travels with the app**
- One worked dataset: **Gale et al. 2020**, the screen the tutorial follows, with 200,000 read pairs per library
  so the mapping stage can be run on first launch. Its QC table carries the donor-plasmid metric: the induction
  pools are 11% and 18% donor reads against about 1.2% in the treated and control cultures.
- The other published datasets are no longer carried in the repository. They are archived per dataset and fetched
  when wanted, which takes a clone from about 16 GB to 2.4 GB.
- References and bowtie2 indexes for eight species, and ten library designs across Hermes, Ac/Ds and piggyBac, so
  a library of your own can be analysed from reads to calls.

**Naming**
- `ytab <command>` works as the documentation writes it: put `scripts/local` on PATH, or call
  `python scripts/local/ytab.py <command>`.
- An analysis records itself by what it does rather than by one transposon: *site-restricted resampling*,
  *window HMM*, *control-conditioned Gumbel*, *site-paired genetic interactions* and the rest. The algorithm
  revisions behind them are unchanged, so existing results still match.

**Donor-plasmid reads**
- A library merged from several sequencing runs is counted across all of them, so merged libraries report a donor
  share like any other.
- Counting a deeply sequenced library no longer needs tens of gigabytes of memory.

**Validation**
- Python tests: 374 passed.
- Verified from a clean checkout, in the published image: the landing page, the demo dataset, all eight tabs, the
  directory chooser reaching the user's own files, and the launcher on a busy port.

## 2026-11-09: donor-plasmid reads, building a dataset, and the help centre

**Quality control**
- **Reads from the transposon donor (%)** and **donor reads placed on the genome (%)**: the share of a library's
  reads that come from the donor plasmid or an unexcised donor locus rather than from a genomic insertion. The reads
  are counted against the donor sequence with the transposon masked out, and only alignments of at least 30 bases
  count, so a read that merely carries the transposon end is not counted. Under *Reads and Mapping*, with their own
  plot group.
- `ytab donor-reads` reads the FASTQs directly: counting them needs no remapping to the genome.
- Donors for five constructs ship under `resources/transposons/donors/`, built from their published parts and
  checked against the real pCU-MET3-Hermes sequence.
- **Build QC table** picks up metrics added since a stored table was made.

**Building a dataset from FASTQ**
- The scan now lists the FASTQ files with checkboxes: tick the ones this dataset should analyse, with
  *Select all files* / *Clear selection* and a running count. Unticked files are still recorded as samples with
  `include=no`, so they can be taken up later from Run > Samples.
- The same choice on the command line: `ytab_init_local_project.py --include-fastq FILE` (repeatable) or
  `--include-fastq-list PATH`. The command the app shows is the command that ran.
- Mapping runs one sample at a time, at the chosen thread count, and can be stopped and resumed; the dataset card
  says so.

**Read cleanup is named for what it does**
- `hermes_cleanup` is now `detect_ends`, `--hermes-preprocess` is `--detect-ends`, `--hermes-reference` is
  `--transposon-reference`, `--hermes-overlap` is `--end-overlap`, and the report is written as
  `<sample>.read_cleanup.json`. The method was never specific to one transposon: it builds its adapters from the
  element's own terminal sequence and a catalogue of the primers used.
- A library design can name its `primer_catalogue` and `reference`, so an element other than Hermes can use detected
  ends; the designs that declare the sequence their construct leaves (Ac/Ds, piggyBac) are unchanged.
- Nothing has to be re-run: the former names are still read, and the QC table accepts either report filename.

**The app**
- Dialogs display as dialogs. Bootstrap is loaded for the whole app, not only the workspace, so the directory
  chooser, Methods, Datasets, Search and every confirmation open over the page instead of unstyled at the bottom of
  it. The landing page keeps the proportions it was designed with.
- **Maximize / Restore** on the Methods and Datasets dialogs.
- Help pages scroll inside the dialog, with the contents list fixed beside them.
- **Equations render.** The methods pages are converted with pandoc to MathML, which the browser draws itself: no
  JavaScript, no network, so the help centre still works offline.
- Public datasets opens at a readable width, one row per study.
- Buttons that are not ready are disabled, rather than looking enabled and refusing the click.
- Streaming a PCA no longer reports a missing function.

**Documentation**
- The tutorial follows the current tabs and names, including Multi-condition and Genetic Interactions, the donor
  metric and picking FASTQ files.
- *Read selection and transposon cleanup* states the two cleanup paths and lists the piggyBac design.
- *Window HMM and hybrid*: the model is not specific to one transposon. `v2` uses the Hermes consensus as its
  sequence covariate; `v3` learns the target-site model from the libraries' own insertion sites, and essentiality
  has been run this way on Ac/Ds libraries.
- The bench protocol for building a Hermes pool has left Methods: YTAB analyses sequenced libraries.
- The README describes the release as it is, with the adapted methods and their references.

**Testing the app**
- An image for testers: `docker run --rm -p 3838:3838 -v "$PWD:/opt/ytab" ghcr.io/jayelazuno/ytab-app:v0.1`. It
  carries the environment only and mounts the clone, so a `git pull` is enough to pick up changes.
- `docs/testing.md` is the guide, including how to mount a FASTQ folder the container would not otherwise see.

**Validation**
- Python tests: 374 passed.
- App smoke suite: 56 passed. The 15 failures all fail without this release's changes as well, which in turn fix
  five that were failing (app initialization, the essentiality and library-correlation tab UIs, run mapping and
  treated-vs-control agreement).
- Verified in the container against the demo dataset: the landing page, all eight tabs, the QC table, the directory
  chooser, both dialogs and the rendered equations.

## 2026-10-23: four species, two transposons

**Mapping**
- Transposon profiles: Hermes (`splinkerette_pcr`, `end_ligated_adapter`, `terminal_10bp`) and Ac/Ds
  (`ds_primer_tail`, `satay_restriction`), chosen per project or per sample (`--transposon`, `--library-design`).
  Hermes projects map exactly as before.
- `MapFastq.py --trim-3p-seq` (repeatable).
- Hit files no longer fail when a chromosome with genes has no insertions.

**Data**
- Public datasets reprocessed from raw reads: Segal et al. 2018 (C. albicans, Ac/Ds), S. pombe Hermes (Guo et al.
  2013, Grech et al. 2019, Lee et al. 2020) and Michel et al. 2017 (S. cerevisiae SATAY, all seven libraries).
  Sources, trimming and references: *Public datasets* in Methods.
- Argon: one download job per dataset (ENA, MD5-checked, resumable) and one submission per project.

**Essentiality**
- `ytab hermes-hmm --sequence-covariate learned` (hermes_hmm_v3): the HMM learns the transposon's target preference
  from the libraries' own insertion sites, so it works for Ac/Ds as well as Hermes. Pre-registered check: it
  rediscovers Hermes' nTnnnnAn from Gale's data alone and is not worse than the motif model (AUC 0.926 vs 0.922).

**Species resources**
- One layout for every species (`reference_genome`, `annotation`, `orthologs`, `essentiality_labels`, `tnseq_data`);
  S. pombe and S. cerevisiae references indexed.

## 2026-10-23: app and pipeline release

**Structure: the pipelines write, the app streams.** Every analysis writes the same tables, whatever the
normalization, library or species, and checks them before it finishes. The app only displays them, so a run started
from the command line appears in the app.

**Normalization**
- Site Normalization runs every method at once and shows each with its own view of the libraries before and after.
- BGC (beta-geometric correction) is the default.

**Condition comparisons**
- Three methods:
  - local-abundance z-score: ranking; replicates averaged per condition, as in Gale et al. 2020;
  - resampling: false-discovery control, both directions;
  - control-conditioned Gumbel: false-discovery control, depletion only.
- Validated in a pre-registered benchmark and on Gale's fluconazole data.
- Top hits follow published practice: called genes with at least a 2-fold change, ranked by |z|.
- Results say which normalization and call rule made them, and whether they match the current Run settings.
- Display options remember whether they are open, and plot labels never overlap.

**Essentiality**
- The Hermes window HMM needs no training labels.
- The hybrid model can be added when curated labels exist.
- Results show the calls, the HMM–hybrid agreement and chromosome tracks.

**Gene Explorer**
- An embedded genome browser built from pipeline files, plus **Open in IGV** and **Download for IGV**.
- GO enrichment works from the new top hits.
- Across Species opens on the current project.

**Screen Planning**
- A coverage planner.
- Mutant identification from Sanger reads, with a bench protocol for nested arbitrary PCR.

**Help**
- A Methods menu, a tutorial (Gale 2020), these release notes, and search across them and the tabs.

**Gene names**
- The common name first, then the gene ID. C. glabrata uses S. cerevisiae ortholog names; S. cerevisiae, C. albicans
  and other species use their own names (C. albicans names now come from CGD's annotation).

**Projects**
- No project depends on an external drive.
- Saved gene queries are limited to the 500 most recent per project.

**Validation**
- Python tests: 353 passed.
- App smoke suite: 58 passed (14 failures that predate this release).
- All stored results reproduce exactly from their recorded commands.
