"""Single-condition essentiality: whether a gene tolerates insertions at all.

Every method here reads one condition's libraries, with no treatment to compare against:

  hermes_hmm   window HMM over each gene's sites, label-free
  Classifier   random forest on gene features, for the hybrid model when curated labels exist
  model_menu   one row per gene with every available model's score and call
  contract     the files and fields an Essentiality run always writes

Methods that compare conditions live in `ytab.fitness`, even when they borrow a statistic from an essentiality
method: a run of insertion-free sites in one library says a gene is essential, while a run of sites the treatment
depleted says only that the treatment acted there.
"""
