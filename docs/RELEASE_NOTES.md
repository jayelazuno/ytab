# Release notes

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
