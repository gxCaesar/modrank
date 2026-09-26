# 20260912-five-cohort-intervals: what the five-study arms look like with intervals

The benchmark's five TCGA studies had been run with all three modalities since the S5 campaign, and
every one of those numbers was a point estimate. An external review on 2026-09-12 recommended desk
rejection partly because the method has no independent validation with all three modalities, and the
answer to that charge could not be written on point estimates. This run adds the paired case-level
bootstrap that every other interval in the manuscript uses (6,000 resamples of cases, comparable
pairs rebuilt inside each replicate, seed 0, identical to `s6_confirmatory.boot`).

**Known answers first.** The script reproduces `stage5.json`'s `OURS_wsi_omics_age_sex_stage` and
`control_same_arm_with_GRADE_instead` for all five cohorts to within 5e-5 before it writes anything.
It did: `known_answer_drift` is absent from the output. That is worth recording for a second reason.
The original five-study run used Python 3.8.20 with NumPy 1.24.4 and SciPy 1.10.1; this one used
Python 3.9.6 with NumPy 2.0.2 and SciPy 1.13.1, and the point estimates are identical. A library
change of that size moving nothing is a stronger reproducibility statement than the rerun that was
done to restore the inputs.

## Observed

2,241 cases and 411 events across the five studies.

| cohort | n | events | OURS − corrected clinical arm | OURS − same arm with grade |
|---|---|---|---|---|
| blca | 359 | 113 | **+0.0774 [+0.0284, +0.1264] p=0.001** | +0.0132 [−0.0136, +0.0396] p=0.318 |
| brca | 871 | 60 | **+0.1208 [+0.0525, +0.1904] p<0.001** | +0.0283 [−0.0155, +0.0723] p=0.202 |
| coadread | 298 | 37 | +0.0236 [−0.0686, +0.1197] p=0.640 | **+0.0911 [+0.0451, +0.1449] p<0.001** |
| hnsc | 394 | 117 | **+0.0726 [+0.0206, +0.1231] p=0.007** | +0.0164 [−0.0181, +0.0496] p=0.334 |
| stad | 319 | 84 | +0.0495 [−0.0127, +0.1122] p=0.113 | +0.0184 [−0.0162, +0.0529] p=0.299 |

**Three of five separate from the corrected clinical reference; two do not.** The two that do not
are the two smallest event counts after BRCA, and COADREAD's 37 events give an interval 0.19 wide.

## What this does NOT support, and the reason each is written here rather than discovered later

- **Not "three modalities beat two".** `ours_minus_wsi_omics` covers zero in every cohort, and is
  negative in BRCA (−0.0211 [−0.0727, +0.0296]). What the clinical block adds to the two molecular
  arms is not resolvable at these event counts.
- **Not "the swap raises the clinical arm in all five studies" as a statistical claim.** The point
  estimates are positive in all five, which is what the manuscript says, but with intervals only
  blca (+0.0563, p=0.036) and coadread (+0.2212, p<0.001) exclude zero. BRCA is p=0.068, HNSC 0.226,
  STAD 0.164. The manuscript sentence needs that qualification.
- **Not "beats the best published value".** These arms clear SurvPath's released-fold value in every
  cohort, and that is the only comparison the released folds support. Against the best value under
  any protocol they lose badly outside bladder: APL reports 0.794 on BRCA and DSCASurv 0.832 on
  COADREAD, against our 0.6892 and 0.7141. The output file carries both reference values per cohort
  so the narrow claim cannot be widened by accident.

## Not run

- Per-cohort noise floors. The bladder floor (24 re-partitions, paired SD 0.0073) has no counterpart
  in the other four, so each cohort's margin is read against its bootstrap interval only.
- Seed averaging. These arms use `fitapply(..., seed=0)` once, where the bladder primary averages
  five seeds. The BLCA value here is 0.6992 against the primary's 0.7212, and that gap is the
  clinical block's source plus the seed count, not a different method.
