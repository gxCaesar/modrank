# Data card

No data were generated for this study and none are redistributed here. Everything below is a public
release by another group, and this card records what was used, in which form, and what was done to
it.

## The cohort

| | |
|---|---|
| cohort | TCGA-BLCA, muscle-invasive and non-muscle-invasive bladder urothelial carcinoma |
| cases scored | 359 |
| events | 113 |
| comparable pairs | 24,219 |
| endpoint | disease-specific survival, from the TCGA Pan-Cancer Clinical Data Resource |
| splits | the released five-fold case-ID splits distributed with SurvPath, used unmodified |
| evaluation | out of fold, five seeds, concordance index |

The case set is not ours to choose: it is whichever cases the released splits contain, which is what
makes the comparison against the published entries a comparison at all. A case absent from those
files is absent here.

## The three inputs

**Histology.** One 768-dimensional slide-level embedding per case, from the public precomputed TITAN
release. No slide is read at run time and no whole-slide image is in this repository. One
open-access tile is fetched by a figure builder, from the NCI Genomic Data Commons tile service, for
the purpose of showing a reader what the modality is.

**Transcriptome.** The 275 pathway means of SurvPath's `combine` pathway definition, computed from
the expression matrix that release distributes, using that release's own definition file. The
pathway definitions were not re-derived.

**Clinical.** Age, sex and pathologic tumour stage. After amendment A1 every clinical value is read
from DIMAF's released split files rather than from a query of our own. That amendment exists because
the frozen run's clinical block took stage from a cached query that selected the first diagnosis
record per case, which is the wrong tumour for a case with several diagnoses. It disagreed with the
incumbent's own file on 26 of 359 cases and had no stage at all for 35, understating stage
systematically. The amendment moved the primary result down, from 0.7260 to 0.7212, and is recorded
in the frozen protocol rather than applied silently.

## Provenance and licensing

Each release is cited in the manuscript's Data availability statement and is obtained from its
authors, not from us. The licences are theirs. Redistributing another group's release inside this
repository would be the wrong call even where the licence permits it, so the code takes the location
of the inputs through an environment variable and reads them where the user put them.

## What this data cannot support

**One cohort.** There is no external validation set here. Every number in this repository is
out-of-fold within TCGA-BLCA, and a concordance index measured that way answers a narrower question
than a clinician would ask.

**113 events.** That is the quantity that governs how finely anything can be resolved. It is why the
selection correction matters, why a 0.0048 spread across seeds is not the relevant uncertainty, and
why a parameter count of 24.7 million against 113 events is a statement about the sample size rather
than about deep learning.

**Retrospective, and not consecutive.** TCGA is an assembled research cohort. Nothing here
establishes how the score would behave on a consecutive clinical series, and no claim in the
manuscript extends to one.
