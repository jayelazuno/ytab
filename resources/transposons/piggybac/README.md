# piggyBac

piggyBac inserts only at TTAA and duplicates the tetramer, so its libraries differ from the hAT elements YTAB was
built for (Hermes, Ac/Ds), which insert almost anywhere:

- **Opportunities are discrete.** A gene's insertable positions are its TTAA sites, not its length. Density is still
  reported as insertion sites per coding kb so that it can be compared with Hermes and Ac/Ds, but a piggyBac library
  cannot reach the density of a hAT library of the same depth. Read the density recommendations
  (`docs/methods/00_recommendations.md`) with that in mind.
- **The target site is in the read.** After the terminal repeat is trimmed the read begins at the first base of the
  TTAA, so an insertion's two orientations give positions 3 bp apart. `ytab hit-files` unites adjacent positions
  within its minimal distance, as it does for the 8-bp duplications of Hermes and Ac/Ds.
- **The HMM needs no change.** The window HMM learns its target preference from the library's own sites
  (`--sequence-covariate learned`), which is how it already handles Ac/Ds; it is not given the TTAA motif.

**Counting convention.** Gao et al. 2021 counted insertions at one position in opposite orientations as two
independent insertions. `ytab hit-files` instead unites adjacent positions, as it does for the 8-bp duplications of
Hermes and Ac/Ds, so YTAB's site counts are not expected to match a piggyBac paper's exactly; gene-level counts are
the comparable quantity. The same paper's analysis keeps reads "that have a prefix matching the terminus of the PB
transposon (TTTCTAGGG)" and removes it, which is what `pb_tr_end` trims.

`profile.json` holds the read structure (`pb_tr_end`), which was checked against the first reads of one run of each
published piggyBac dataset; `used_by` lists them. The terminal-repeat sequence trimmed here is the 3' end of the
element as it appears in the read, ending at the TTAA junction.
