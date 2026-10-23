# Ac/Ds library designs

`profile.json` describes how reads from Ac/Ds (Ds, MiniDs) libraries are trimmed before mapping. A project selects it
with `mapping: {transposon: ac_ds, library_design: <design>}` in its config, or per sample with `transposon` and
`library_design` columns in the sample sheet.

- `ds_primer_tail`: the Ds primer and end used by Segal et al. (2018), copied from the authors' pipeline
  (github.com/berman-lab/transposon-pipeline, `MapFastq.py`: `TnPrimerAndTail`, `AdapterSeq`).
- `satay_restriction`: SATAY (Michel et al. 2017). The sequences after the DpnII (GATC) and NlaIII (CATG) sites were
  read from the wild-type libraries (ERR1533147, ERR1533148); they are the circularised MiniDs sequence that follows
  the genomic fragment.

Both designs report the insertion at the first genomic base of the read, as for Hermes.
