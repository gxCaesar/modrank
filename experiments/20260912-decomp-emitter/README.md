# 20260912-decomp-emitter: the last panel whose numbers had no producer

`development/s5-results/decomp.json` held the four clinical constructions that the missing-variable
figure's panel e drew. It was written on 2026-08-17 by a script that was never committed, which made
it the one file in this work whose values reached a display item with no way to recompute them. The
figure's Source Data disclosed that; disclosure is not a fix.

This rebuilds the constructions from committed code, using the recipe the amendment used, imported
rather than restated: the released folds, ridge Cox with the protocol's alpha grid at seed 0,
within-fold percentiles, and the equal-weight rank average of slide, transcriptome and clinical.

## Observed

| row | what | clinical alone | with both modalities | old file | agrees |
|---|---|---|---|---|---|
| A | frozen block: age, sex, one-hot T, N, M, stage from the cached GDC query (23 columns) | 0.6856 | 0.7291 | 0.6863 / 0.7285 | **no** |
| B | amended block: age, sex, one-hot stage from the incumbent's split files (6 columns) | 0.6638 | 0.7225 | 0.6638 / 0.7225 | yes |
| D | the same with grade in place of stage (4 columns) | 0.5666 | 0.6856 | 0.5666 / 0.6856 | yes |

**Two of the four rows were unprovenanced, not one.** That is the finding.

B and D are the gates here, because their pairs also appear in the amendment record, which has its
own producer. Both reproduce exactly, so the recipe is the one the paper used.

**A does not match, and the rebuild is the value that agrees with everything else.** 0.6856 is also
the frozen run's own clinical arm and the step-0 baseline for age, sex, T, N, M and stage. The old
file's 0.6863 matches neither, so that row came from some other handling of the same variables that
was never written down. The figure now draws the rebuilt pair.

**C was not attempted at all.** Its label is "corrected GDC, right diagnosis", and the rule that
chose one diagnosis per case for all 359 patients was never recorded. The amendment queried the
archive directly for *four* of the disagreeing cases, to diagnose the defect rather than to build a
corrected block. Re-querying would let me invent a selection rule and check whether my invention
happened to reproduce 0.6992 and 0.7335, which would establish nothing whichever way it came out.
The row is dropped from the figure rather than reconstructed by guess.

## What changed in the paper

Panel e goes from four constructions to three, and two drawn values move by less than 0.001. In
exchange every number in that panel has a script behind it, and `decomp.json` is no longer read by
any builder. The file stays in the tree as the record of what the figure used to draw.
