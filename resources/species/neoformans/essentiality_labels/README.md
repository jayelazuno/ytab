# Curated C. neoformans essentiality labels

Built by `scripts/local/ytab_curate_neoformans_labels.py` from direct *C. neoformans* deletion evidence and from
*S. cerevisiae* (SGD) viability through 1:1 orthologs only. Rules are in the script docstring; sources, retrieval
dates and SHA-256 hashes are in `sources.json`; counts in `summary.json`.

| tier | evidence | essential | nonessential |
|---|---|---|---|
| gold | a viable deletion strain exists in the Madhani 2024 collection | — | 4,268 |
| silver | *S. cerevisiae* null-mutant viability (SGD), 1:1 orthogroup | 709 | 590 |

8,334 annotated genes; 5,567 labelled, 2,767 unlabelled.

- **Unlabelled (2,767 genes):** no 1:1 *S. cerevisiae* ortholog in Billmyre's S4 Table (which lists orthologs for
  2,862 genes, 2,325 of them 1:1), or a 1:1 ortholog whose SGD viability is not curated, and no deletion strain.
- **Asymmetry, stated plainly:** there is no gold *essential* source for this species. Every essential label is
  silver, inferred from the *S. cerevisiae* ortholog; the nonessential labels are mostly gold. An evaluation using
  these labels is therefore more trustworthy about what is **not** essential than about what is.
- **Absence from the deletion collection is not evidence of essentiality.** A gene may simply not have been
  attempted, so `in_madhani_collection = false` yields no label — only `true` does.
- **Conflicts (41 genes):** a viable deletion strain exists although the *S. cerevisiae* ortholog is essential. The
  gold label wins and the silver one is kept in `conflict`. These are real species differences, as for *C. glabrata*.
- **`billmyre2025_comparator`:** the ESS/NESS classification of Billmyre et al. 2025. **These are never labels.**
  Their own S1 Table description gives the method — "the mean from 0-1 of essentiality prediction score from 100
  random forest model runs" — fitted to insertion features of the same libraries YTAB models are run on, so using
  them as truth would compare two pipelines on shared data. This is the rule the *C. glabrata* labels already apply
  to Gale et al. 2020 (`resources/species/glabrata/essentiality_labels/README.md`). Agreement with the curated
  labels is 91.9% (653 of 703 essential, 4,411 of 4,807 nonessential).

## Source

Billmyre RB, Sun S, Heitman J et al. (2025) Landscape of essential growth and fluconazole-resistance genes in the
human fungal pathogen *Cryptococcus neoformans*. PLoS Biol 23:e3003184. doi:10.1371/journal.pbio.3003184.
S1 and S4 Tables are in `raw/billmyre_2025_plosbiol/`.
