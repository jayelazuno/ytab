# C. glabrata CBS138 essentiality labels

The same curated labels as `resources/species/glabrata/essentiality_labels/labels.tsv`, re-keyed on the CBS138
systematic identifier (`gene_id`, the `CAGL*` ORF names) so that projects mapped to the CBS138 reference can use
them. The BG2 identifier each row came from is kept as `bg2_gene_id` (not `cglab_id`, which the hybrid would read as a second key column); nothing else is changed, and no label was added,
dropped or re-decided here.

Rebuild it from the BG2 table whenever that table changes; the hybrid model keys on the first column, so a project
must be given the labels file for the strain it was mapped to.
