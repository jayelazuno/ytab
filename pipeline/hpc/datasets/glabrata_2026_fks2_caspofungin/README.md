# C. glabrata FKS2-F659del Hermes, caspofungin dose series (PRJNA1454099)

Hermes Tn-seq of *C. glabrata* BG2 carrying **FKS2-F659del**, grown without drug and at eight caspofungin
concentrations (250 to 32,000 ng/mL). Johns Hopkins deposit, "FKS2-F659del Caspofungin Resistance and
Heteroresistance"; no publication at the time of writing.

The BioSample records label these runs only `Sample1`-`Sample9`; the design is in the **SRA experiment titles**
(`Tn-seq of C. glabrata Strain BG2 FKS2-F659del: <dose> caspofungin`), which is where the dose of each run comes
from. Library preparation, from the deposit: sonication, end repair, indexed splinkerette ligation, bead size
selection, PCR of the insertion sites.

One pool, one no-drug control and eight doses, so resampling and the multi-condition dose series are available; the
z-score is not (one control). Its value is the pairing with the BG2 **wild-type** caspofungin series of Nickels &
Cunningham 2025 (PRJNA1247003, already in YTAB): the same drug in a resistant mutant and in the wild type, which is
a strain x treatment comparison (`ytab multicondition-project --independent-libraries`). It is not a genetic-
interaction design: that needs >= 2 control/treated pairs per strain and each strain here has one.

Reference: `resources/species/glabrata` (BG2).
