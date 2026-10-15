# Hermes reference for mutant identification

`L34807.1.fasta` is the 2,749 bp accessioned Hermes sequence downloaded from
https://www.ebi.ac.uk/ena/browser/api/fasta/L34807.1
(https://www.ncbi.nlm.nih.gov/nuccore/L34807.1).

Downloaded FASTA SHA-256:
`3ea00548a600ac382dc3f454a5778cd1c9e45bf4fe5327b6c68bf024224230ab`.

The HR-2 primer in the supplied method ends at base 2,716, leaving a 33 bp
Hermes tail. Both terminal orientations are used for clipping. The wild-type
reference does not describe the engineered marker cassette; users may supply
their complete construct FASTA. Runs record hashes, clipping evidence and tool
versions. No download is performed during a UI run.

`primers.json` records the library preparation and sequencing oligos in Gale et al.
2020 Supplemental Table S3 (doi:10.25387/g3.12445373), their purpose, and the
distinction between HR2 PCR and the custom Illumina sequencing primer. It also
contains common Illumina adaptor sequences, explicitly named separately from the
Gale splinkerette oligos. Phosphate and phosphorothioate notation is omitted from
the machine-readable nucleotide strings; these are chemical modifications, not
sequenced bases. Indices are retained as metadata, not standalone trimming anchors.
See `docs/hermes_alignment.md` for automatic trimming and declared-mate selection.
