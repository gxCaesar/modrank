# Amendment A2: Breslow risk sets and average ranks, every post-freeze result regenerated

Protocol: `protocol.md` (written before any result was regenerated), following the author's decision
of 27 September 2026. The values before the amendment are listed beside the regenerated ones in
`results/before-after.json`.

## What changed in the code

- `analysis/blca_common.py`: `cox_objective`/`cox_fit` use the full Breslow risk set of each tied
  time from a stable sort; `pct` and `ranks` give tied values their average rank. One environment
  variable, `BLCA_A2` (unset or 1 = corrected, 0 = original), switches every primitive, and every
  `BLCA_A2=0` path is the old code byte for byte. `tests/test_a2_survival_primitives.py` (9 tests)
  checks the objective against an O(n^2) Breslow brute force, row-order invariance, equality with
  the old objective when no times tie, the gradient against finite differences, order-free ranks,
  the score test, and the tuned-stacking selection rule.
- Local rank helpers and `s5_omics_arm.score_test` follow the same switch (list in `protocol.md`).
- Known-answer pins that hold pre-A2 values are checked exactly under `BLCA_A2=0` and only as a
  gross bound (`A2_SANITY` = 0.01) under A2; the exact check is the legacy reproduction below.
- Protocol item 3, the added comparator: `blca_common.stack_weights_tuned` chooses the stacking
  ridge from {0.1, 1, 10, 100} by three-fold cross-validation over the training fold's inner
  out-of-fold percentiles (selection as in `fitapply`). s19 and s32 compute it under A2 only and
  add fields; no existing field changes (checked: s19's A2 output with and without it are equal in
  every other field, and s19, s32 still reproduce their committed files under `BLCA_A2=0`).
- Two reported numbers had no committed producer and gained one, each first shown to reproduce
  the pre-A2 value from the pre-A2 inputs: the clinical arm's sex gap minus ModRank's (s28,
  +0.0555 [-0.0531, +0.1695] reproduced exactly) and the count of retraining-fold pairs with a
  constant risk (s33).
- `analysis/s37_amendment_a2_before_after.py` writes `results/before-after.json`, the headline
  quantities before (read from the last commit whose results were all pre-A2) and after. SI Table 15 is bound to it cell by cell.

## Observed

Step 1, every producer in scope with `BLCA_A2=0` in the environment of its committed file
(`legacy-checks/`): zero differences for s6, s8, s9, s10, s11, s12, s14, s16, s19 (conda "ml"),
s19b, s19c, s22 (ml), s23, s23b, s27 (values; three metadata keys added), s27 with the second
encoder, s28, s29, s31, s32, s33, s34, s35 (blca, five, geo), s36 (analyze, summary), stage5,
s7b, s7 interpretation, reviewer gaps, survival metrics, SurvPath multiseed, and on sysu s24, s24b.
s13 reproduces every committed field and adds one the script gained after the file was written.
s7_ablation differs in its timing fields and in one value that is not an A2 effect (below).

Step 2, the same producers with `BLCA_A2=1`, in dependency order, reading upstream A2 files. The
headline moved from 0.7212 to 0.7214. `results/before-after.json` lists 49 quantities. The
statements that changed:

- external cohorts where the stage-based reference's interval excludes zero: 4 -> 5 of 6
  (GSE19915, gap +0.129 -> +0.161, whose clinical models predict from a few categories);
- DeepMISL returned one risk for the whole validation fold in 13 of 25 retraining-fold pairs.
  Row-order ranks had passed patient order into its concordance; with shared ranks its D is
  +0.005 (was +0.034) and its excess over the no-signal benchmark is -0.0204 [-0.0399, -0.0021].
  Three architectures, not four, have D resolved from zero;
- GSE48075's permutation-screened variant: 0.6326 -> 0.6170, now below the stage block (0.6334);
- selection bar 0.7716 -> 0.7708; the primary still does not clear it.

Tuned-ridge stacking: 0.7147 +/- 0.0048 on TCGA-BLCA, ModRank ahead by +0.0055 [-0.0071, +0.0180];
none of its five five-study differences separates (30 fusion comparisons, none separating).

## Interpretation

The two defects moved every number, mostly in the fourth decimal. Where predictions tie heavily
(DeepMISL, categorical clinical models in external cohorts) the movement was larger and in each
case the A2 value is the correct one: the old values depended on row order.

## Found on the way, not A2 effects

- `s7_ablation_and_generalisation`: brca `best_published_verified_folds` was 0.759 in the committed
  file; the script's 2026-09-11 correction (DIMAF ships bladder splits only) gives 0.736 and had
  never been written back. Supplementary Figure 1 recomputed it locally; it now asserts agreement.
- `geo-gated.json` (s24b), reported in SI Note 5, was missing from the dependency inventory because
  no gate read it; it was regenerated with the GEO chain.
- The SI table parser in `paper/check_numbers.py` matched from the first tabular in the document,
  so a row could be found in an earlier table; fixed, and the pre-A2 tree still passes with it.

## Not run

- The development campaign (pre-freeze `s5_*`, the iteration ledger, the regime evidence, the
  component slate) and the frozen confirmatory run (0.7260) are reported as recorded, per the
  protocol's scope.
- `run_jobs.py` is the driver used; it expects the scratch layout it created (`A2_SCRATCH`), the
  legacy interpreter in `A2_LEGACY_PYTHON` and the current one in `A2_PYTHON`.
