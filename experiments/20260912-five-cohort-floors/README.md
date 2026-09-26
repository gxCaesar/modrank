# 20260912-five-cohort-floors: what each study's own re-partitioning says

The five-study validation read every margin against its bootstrap interval and nothing else. The
bladder result in the main text is held to a second and stricter standard, twice the standard
deviation of the same paired difference under re-partitioning, because the two answer different
questions: an interval says how far the estimate moves if these patients are resampled, a
re-partition says how far it moves if somebody re-runs the study. The other four studies had no such
number, so "how much of that margin is partition luck" had no answer for them.

Twenty-four re-partitions per study, built by `s5_step0_headroom.reseeded_folds` with rng 1000 + r
exactly as the bladder floor was, refitting the slide, transcriptome and clinical arms inside each
one. Partition 0 is the released folds and reproduced every committed point estimate before a single
re-partition was drawn; `known_answer_drift` is absent from the output. 796 seconds.

## Observed

| study | n | events | margin | floor (SD) | bar (2 SD) | clears | partitions positive | ModRank range |
|---|---|---|---|---|---|---|---|---|
| bladder | 359 | 113 | +0.0774 | 0.0117 | 0.0233 | **yes** | 24/24 | 0.6728–0.7178 |
| breast | 871 | 60 | +0.1206 | 0.0237 | 0.0473 | **yes** | 24/24 | 0.6764–0.7757 |
| colorectal | 298 | 37 | +0.0243 | 0.0135 | 0.0270 | **no** | 24/24 | 0.7509–0.8063 |
| head and neck | 394 | 117 | +0.0729 | 0.0132 | 0.0263 | **yes** | 24/24 | 0.6144–0.6620 |
| stomach | 319 | 84 | +0.0491 | 0.0162 | 0.0323 | **yes** | 24/24 | 0.6259–0.6804 |

**The paired difference is positive in every one of the 24 partitions in all five studies.** Four
margins exceed twice their own re-partition noise.

## The one place the two criteria disagree, and why it is reported rather than resolved

**Stomach clears this bar and its bootstrap interval covers zero** (+0.0491 [−0.0127, +0.1122]).
Colorectal fails both. So across the five studies: three pass both criteria, one passes only the
re-partition bar, one passes neither.

That is not a contradiction to be argued away. Resampling patients and re-partitioning folds are
different perturbations, and a margin can be robust to one and not the other — a stable ordering that
depends on which patients happen to be in the cohort looks exactly like this. Reporting only the
criterion stomach passes would turn a three-of-five result into a four-of-five one by choosing the
question after seeing the answer.

## The bladder row is not the main text's 0.0073

It is 0.0117 here. The main text measures that floor with the incumbent's clinical block and five
seeds; this protocol uses the benchmark's released clinical file and one seed, which is noisier. The
point of measuring bladder under both is that it gives the four non-bladder floors something already
in the paper to be read against.

## Not run

Per-study floors for the grade control and for the two-modality arm. The floor that matters for the
claim is the one on ModRank minus the clinical arm, and adding two more would have tripled the
runtime for quantities no sentence rests on.
