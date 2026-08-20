# Model card

## What it is

ModRank is three ridge Cox models and an average. One model is fitted per modality, each on that
modality alone, and the three risk scores are converted to within-fold ranks and averaged with equal
weight.

| | |
|---|---|
| slide block | ridge Cox on 768 TITAN dimensions |
| omics block | ridge Cox on 275 SurvPath pathway means |
| clinical block | ridge Cox on age, sex and pathologic stage, six columns after coding |
| combination | equal-weight average of within-fold ranks |
| coefficients | 1,049 |
| parameters fitted in the combination | 0 |

The combination has no weights, no gate, no attention and no temperature. That is the design, not an
ablation of one: a weighted fusion was in the candidate slate and the equal-weight rank average is
what survived it.

For reference, and taken verbatim from each trainer's own parameter file rather than reconstructed
by us, SurvPath reports 24,702,532 parameters and PIBD 26,930,702, against 1,049 here and 113 events
in the cohort.

## Intended use

Method-development research on this benchmark, and reproduction of the numbers in the accompanying
manuscript. It is a reference arm: something a benchmark entry should be asked to beat before its
architecture is credited with the result.

## Out of scope

**Not for clinical use, in any form.** It has never been evaluated prospectively, it has never seen
a consecutive series, and it has no regulatory status. Do not use it to inform the care of a
patient.

**Not validated outside TCGA-BLCA.** Every number here is out of fold within one cohort. Nothing
establishes transfer to another institution, another scanner, another assay platform, or another
disease.

## How it was evaluated

Five-fold case-ID splits released by SurvPath, five seeds, concordance index on disease-specific
survival, out of fold. Primary result 0.7212, standard deviation across seeds 0.0048, on 359 cases
and 113 events.

The confirmatory run was executed against a protocol frozen beforehand, and the protocol, its
amendment, the run command, the environment record and the run log are all in this repository.

## The result this model card exists to state plainly

**The pre-registered selection correction was mis-specified, and by a correct one the reported
result does not clear its bar.**

The frozen rule set the target at the incumbent's 0.679 plus a selection-inflation term computed as
one standard deviation times the square root of twice the log of the candidate count. The standard
deviation it used, 0.0073, is the spread of a *paired between-arm difference*. The maximum-of-N
construction needs the spread of *a single candidate's own absolute score*. Those are different
quantities and the second is much larger.

An outcome-permuted null measures the second directly: permuting the time-and-event pair across
cases, leaving folds, features and estimator untouched, over 200 permutations. It gives a
single-candidate standard deviation of 0.0397, five times the assumed value.

| bar | value | 0.7212 |
|---|---|---|
| as frozen, sigma 0.0073 at N=211 | 0.7029 | clears |
| Gaussian maximum-of-N at the measured sigma | 0.8087 | fails |
| empirical 95th percentile of the permuted maximum | 0.7716 | fails |

The frozen arm ranks 11th of the 206 comparable candidates the campaign scored, with 8 strictly
above it and 3 tied.

The manuscript reports this. It is the reason the contribution is stated as a benchmark finding
about clinical baselines rather than as a state-of-the-art claim, and `verify_release.py` recomputes
both corrected bars and asserts the failure rather than leaving a reader to take it on trust.

## Other limitations worth stating

**The lead is partly the encoder.** Holding the method fixed and swapping the slide block across
seven public encoders moves the full method between 0.6704 and 0.7225. The slide arm alone spans
0.1372 and the full method spans 0.0521, so the rank average absorbs most but not all of the
encoder's contribution, and five of the seven still land above the incumbent's 0.679.

**Neither input-parity margin survives Holm correction.** That is reported in the manuscript rather
than reframed.

**Stage is doing a great deal of the work**, which is the point of the paper's benchmark finding and
also a limitation of the model: a reader should expect the clinical block alone to be strong, and it
is, at roughly 0.665.

## Fairness and subgroups

Per-stage subgroup results are in the confirmatory run's results. No analysis by race, ethnicity or
socioeconomic variable was performed, and TCGA's composition would not support a reliable one at 113
events. Absence of a subgroup analysis here should not be read as evidence of equal performance
across subgroups.
