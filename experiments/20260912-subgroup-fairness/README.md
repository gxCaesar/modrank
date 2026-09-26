# 20260912-subgroup-fairness: the two TRIPOD rows that said "Not addressed"

TRIPOD+AI item 3c (health inequalities) and item 14 (fairness) were answered "Not addressed" in this
work's checklist. They were answerable from data already on disk: the benchmark's released metadata
carries sex and age for all 359 patients, and the seed-averaged out-of-fold scores are the ones the
main result is computed from. **Nothing is refitted here.** The run recomputes concordance inside
each stratum from committed vectors, and refuses to write unless the pooled values still equal
0.7237 for ModRank and 0.6644 for the clinical arm. They did, with no drift, so the strata are
partitions of the numbers the manuscript reports.

## Observed

| stratum | n | events | clinical | slide | transcriptome | ModRank | ModRank − clinical |
|---|---|---|---|---|---|---|---|
| men | 269 | 82 | 0.7015 | 0.6613 | 0.6637 | 0.7406 | +0.0396 [−0.0143, +0.0936] |
| women | 90 | 31 | **0.5652** | 0.6298 | 0.6631 | **0.6616** | +0.0963 [−0.0031, +0.1964] |
| age ≤ 68 | 180 | 58 | 0.6578 | 0.6817 | 0.6308 | 0.7285 | +0.0699 [+0.0072, +0.1335] |
| age > 68 | 179 | 55 | 0.6749 | 0.6272 | 0.6962 | 0.7213 | +0.0464 [−0.0309, +0.1189] |

Between strata: ModRank men − women **+0.0811 [−0.0231, +0.1900]**, p=0.132. Age below − above the
median **+0.0065 [−0.0849, +0.0996]**, p=0.886.

## Interpretation, and the two readings that do not survive

Every arm scores lower in women, and the clinical arm most of all. **Neither tempting conclusion is
supported.** ModRank's sex gap covers zero. And the difference between the clinical arm's sex gap and
ModRank's, which would say the method is less sex-dependent than stage, is +0.0555 [−0.0531,
+0.1695], p=0.326 — also covering zero. Ninety women carrying 31 events is why: the interval on any
statement about them is about 0.2 wide. Event rates are similar between the sexes (0.305 against
0.344), so this is not an artefact of how many events each group contributes.

One observation is worth recording because it runs with the paper's main finding rather than against
it: in women, where the clinical arm is weakest at 0.5652, the slide arm reaches 0.6298 and the
transcriptome arm 0.6631, both above it. That is what the paper's thesis predicts for a stratum
where the staging variable carries least. It is a point estimate in 90 patients.

## Race: obtained rather than assumed unavailable, and the cohort still cannot carry it

Race is absent from the benchmark's released files, which is a statement about the released subset
and not about the study, so it was requested from the Genomic Data Commons for all 359 patients
(`analysis/s30_fetch_gdc_race.py`; the patient-level map lands in the gitignored inputs tree and the
counts in `results/gdc-demographics-counts.json`).

| recorded race | n | events | comparable pairs | reportable |
|---|---|---|---|---|
| white | 285 | 93 | 16,454 | **yes** |
| asian | 41 | 5 | 50 | no, below the 100-pair floor |
| black or African American | 20 | 9 | 108 | no, below the 25-patient floor |
| not reported | 8 | 5 | 18 | no |
| unknown | 5 | 1 | 2 | no |

**Only one stratum clears the floors, so there is no comparison to make.** Inside the white stratum
the clinical arm reaches 0.6539, the slide arm 0.6505, the transcriptome arm 0.6544 and ModRank
0.7122, a difference over the clinical arm of +0.0588 [+0.0043, +0.1132].

Getting a race comparison out of this cohort would mean abandoning the floors in order to obtain
one, on five events in the second-largest group. That is a fact about the benchmark as much as about
this analysis: a public resource a subfield measures itself on does not support the subgroup check
that subfield's own reporting standard asks for.

Insurance status and socioeconomic position are recorded nowhere in these data and were never
available. Contributing site is a third axis and is reported separately, as whole sites held out in
the site-grouped cross-validation in the main text.
