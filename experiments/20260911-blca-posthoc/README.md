# 20260911-blca-posthoc: post-freeze analyses on TCGA-BLCA

These analyses ran after the 2026-08-17 freeze and are descriptive. They leave the frozen primary
(0.7212) unchanged. Code: `analysis/s19_unified_fusion_and_added_value.py` and
`analysis/s19b_modality_atlas_amended.py`. Command: `run.sh` (s19). s19b was run by hand with the
same interpreter; its command is below. Environment: `env.txt` (Mac, Python 3.8.20, NumPy 1.24.4,
SciPy 1.10.1: the frozen run's stack).

## Observed

**Known answers.** All eight reproduced exactly before anything was reported: 0.5666, 0.6638,
0.6856, 0.7225, 0.7212, 0.7237, 0.6897 and 0.6989. s19b first reproduced the development atlas on
its own pre-amendment block (0.2130, 0.1656, 0.5012, 0.6120) and only then swapped the clinical
block.

**Fusion on identical inputs.** Same cases, folds, seeds 0–4, estimator:

| construction | per-seed mean ± SD | canonical | paired vs ModRank (6,000 case resamples) |
|---|---|---|---|
| ModRank, equal-weight rank average | 0.7212 ± 0.0048 | 0.7237 | — |
| concatenated ridge Cox, 1,049 columns, grid to 32,768 | 0.6864 ± 0.0063 | 0.6898 | +0.0340 [+0.0014, +0.0674], p=0.042 |
| stacked, weights learned on inner out-of-fold scores | 0.7106 ± 0.0044 | 0.7129 | +0.0108 [−0.0052, +0.0269], p=0.185 |

**Evidence-gated variant.** Every arm cleared its within-training permutation null in all 25
fold-seed cells, so the variant is bit-identical to the frozen method on this cohort (0.7212). Inner
out-of-fold C by arm: slide 0.568–0.683, omics 0.570–0.712, clinical 0.602–0.678. The null q95 ran
0.549–0.560.

**Added value, and how much a grade-based reference inflates it.** Canonical vectors; weak
reference age+sex+grade C=0.5665, stage reference age+sex+stage C=0.6644.

| construction | added over grade | added over stage | D (inflation) |
|---|---|---|---|
| ModRank | +0.1204 [0.0684, 0.1708] | **+0.0594 [0.0113, 0.1063]** | +0.0610 [0.0284, 0.0955], p=0.0003 |
| SurvPath + clinical | +0.0670 [0.0249, 0.1085] | +0.0253 [−0.0139, 0.0634] | +0.0417 [0.0196, 0.0668], p<0.001 |
| PIBD + clinical, best-val | +0.0841 [0.0408, 0.1290] | +0.0345 [−0.0046, 0.0745] | +0.0496 [0.0238, 0.0769], p<0.001 |
| PIBD + clinical, final epoch | +0.0405 [−0.0099, 0.0886] | −0.0004 [−0.0432, 0.0394] | +0.0409 [0.0158, 0.0682], p=0.001 |

The p of 0.0 that the script writes for SurvPath and best-val PIBD means no replicate crossed zero
among 6,000, so it is reported as p<0.001. ModRank's added value over stage is the Holm table's "vs
clinical alone" comparison (p 0.0163 raw, q 0.0652 after Holm), so it is not significant after
correction.

**Modality paragraph, amended block, seed 0.**

| quantity | amended | pre-amendment value it replaces |
|---|---|---|
| Spearman slide vs clinical | 0.2445 | 0.213 |
| Spearman incumbent vs clinical | 0.1924 | 0.166 |
| tightest 5% of pairs | 1,213 pairs | 1,281 |
| clinical C on those pairs | 0.4959 | 0.5012 |
| slide C on those pairs | 0.6220 | 0.6120 |

**Resplits (s22).** The noise floor's own 24 random five-fold partitions, estimator seed 0.
Known answer on the released folds: 0.7225, 0.691 and 0.7132, all exact.

| difference | mean ± SD | positive in |
|---|---|---|
| ModRank − concatenated | +0.0464 ± 0.0158 | 24/24 |
| ModRank − stacked | +0.0111 ± 0.0065 | 24/24 (min +0.0036) |
| ModRank − clinical | +0.0549 ± 0.0083 | 24/24 |
| clinical stage − grade | +0.0919 ± 0.0114 | 24/24 |
| D, added-value inflation | +0.0628 ± 0.0086 | 24/24 |

Across the resplits ModRank scores 0.7158 ± 0.0089, with a range of 0.698–0.732.

**Site-grouped five-fold CV (s22).** 33 tissue source sites, assigned whole, so no site is on both
sides of any fold. Pooled results:

| construction | C-index |
|---|---|
| ModRank | 0.7210 |
| stacked | 0.6855 |
| omics | 0.6730 |
| clinical | 0.6568 |
| concatenated | 0.6487 |
| slide | 0.6312 |
| grade | 0.5636 |

ModRank is best in 4 of 5 site folds. In the fifth, stacked scores 0.7571 against ModRank's 0.7545.

**Calibration and decision curve (s23).** Fitted inside the training folds, seed 0; the known
answers 0.7225 and 0.6638 were reproduced first.

| | ModRank | clinical |
|---|---|---|
| calibration slope | 0.770 | 0.849 |
| 24-month mean predicted vs observed risk | 0.342 vs 0.341 | 0.337 vs 0.341 |
| IPA at 12 / 24 / 36 months | 0.098 / 0.152 / 0.105 | 0.047 / 0.079 / 0.105 |
| net benefit at 24 months, thresholds 0.30 / 0.40 / 0.50 | 0.144 / 0.136 / 0.071 | 0.108 / 0.090 / 0.023 |

Treat-all at the same thresholds: 0.059, −0.098, −0.318. At risk at 12, 24 and 36 months: 261, 127
and 78.

## Not run

- Bootstrap intervals on the IPA and net-benefit differences.
- The five-study version of D, which needs the four non-bladder input files restored (Day 5).

## s19b command

```
D=data/frozen-inputs; OMP_NUM_THREADS=4 $PYTHON -W ignore \
  analysis/s19b_modality_atlas_amended.py --root $D --titan $D/TCGA_TITAN_features.pkl \
  --sp-dir $D/survpath_preds --dimaf-dir $D/dimaf_splits \
  --out experiments/20260911-blca-posthoc/results/modality-atlas-amended.json
```

## s19c: the double-tied probe (added 2026-09-11, while binding the NC build's numbers)

The manuscript's "direct test" (0.6396 and 0.5604 on the pairs the transcriptome arm separates
least; 0.6212 on the 3,893 pairs neither the transcriptome arm nor the clinical model separates) had
no committed emitter. The values lived only in `development/analysis-completeness-audit.md`.

- **Reproduced exactly**, after five known answers (24,219 pairs; seed-0 slide 0.6596,
  transcriptome 0.6510, amended clinical 0.6638, development-era clinical 0.6856).
- **The note's "stage" was the clinical model** (age, sex and stage from the incumbent's split
  files, the amended block). Read literally as "same stage group", the restriction gives 3,183 pairs
  and 0.6030 (amended stage) or 2,335 and 0.6066 (development-era). The prose now names the
  clinical model. The stage-group reading is kept as a sensitivity row.
- **Rounding.** The committed per-case vectors are rounded to five decimals, which creates
  cross-fold collisions between within-fold percentiles. A gap-quantile restriction is sensitive to
  them: the rounded clinical vector gives 3,892 pairs and 0.6209. The script therefore refits the
  clinical arm at full precision.

```
R=$(pwd); D=$R/data/frozen-inputs
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 $PYTHON -W ignore \
  $R/analysis/s19c_double_tied_probe.py \
  --vectors $R/experiments/20260911-blca-posthoc/results/percase-canonical-vectors.json \
  --root $D --titan $D/TCGA_TITAN_features.pkl --sp-dir $D/survpath_preds \
  --dimaf-dir $D/dimaf_splits --tnm $D/tnm.json \
  --out $R/experiments/20260911-blca-posthoc/results/double-tied-probe.json
```

Output (stderr, verbatim):
```
amended block: 3893 pairs, slide 0.6212, transcriptome 0.5543, clinical 0.5634 | stage-group sensitivity (amended) 3183 pairs, slide 0.6030
```

## s7b: a committed emitter for why-grade-fails.json (added 2026-09-11)

`experiments/20260817-blca-confirm/results/why-grade-fails.json` (the 94.4%, 0.311 and 0.665 of the
mechanism, and the five-study entropy table) had no committed producer. `analysis/s7b_why_grade_fails.py`
recomputes every field from the benchmark's released clinical and metadata files and the five-study
results, and refuses to write unless the result equals the committed file. It did, on the first run.
The committed file is left byte-identical; the reproduced copy is not stored.

```
R=$(pwd)
/usr/bin/python3 $R/analysis/s7b_why_grade_fails.py --meta-dir $R/data/frozen-inputs/meta \
  --clin-dir $R/data/frozen-inputs/clin \
  --stage5 $R/experiments/20260911-five-study-repro/results/stage5.json \
  --committed $R/experiments/20260817-blca-confirm/results/why-grade-fails.json --out <scratch>
```
Output (stderr, verbatim):
```
reproduced: 5 studies, grade entropy {'blca': 0.311, 'hnsc': 0.643, 'stad': 0.623}, Pearson -0.976 over 3 studies with a grade
```

## s23b: intervals for the utility differences (added 2026-09-11)

s23 reported each model's IPA and net benefit, not the uncertainty of the difference. s23b rebuilds
s23's predictions line for line, refuses unless every committed s23 value reproduces (it did), then
resamples patients 6,000 times with the predictions fixed. Output (stderr, verbatim):
```
s23 reproduced: concordances, slopes, Brier and IPA at 12, 24, 36 months
IPA difference at 12 months: +0.0512 [-0.0003, +0.0993] p=0.0520
IPA difference at 24 months: +0.0724 [-0.0054, +0.1502] p=0.0703
IPA difference at 36 months: -0.0013 [-0.0897, +0.0842] p=0.9940
net benefit difference at 20%: +0.0090 [-0.0235, +0.0413] p=0.5860
net benefit difference at 30%: +0.0364 [-0.0038, +0.0768] p=0.0790
net benefit difference at 40%: +0.0454 [+0.0006, +0.0915] p=0.0467
net benefit difference at 50%: +0.0482 [-0.0075, +0.1035] p=0.0910
net benefit difference at 60%: +0.0250 [-0.0279, +0.0766] p=0.3383
```
**Reading.** Every IPA interval includes zero, and only the 40% net-benefit interval excludes it
(uncorrected for five thresholds). The manuscript's "improves calibrated two-year prediction over the
clinical model" was therefore stronger than the evidence, and now reads as a gain this cohort cannot
resolve.
