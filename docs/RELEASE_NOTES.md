# Release notes

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

**Projects**
- No project depends on an external drive.
- Saved gene queries are limited to the 500 most recent per project.

**Validation**
- Python tests: 353 passed.
- App smoke suite: 58 passed (14 failures that predate this release).
- All stored results reproduce exactly from their recorded commands.
