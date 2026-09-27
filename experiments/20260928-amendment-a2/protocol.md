# Amendment A2: the Cox tie handling and the rank ties, corrected everywhere

Written 2026-09-27, before any result is regenerated. PI decision H4-2026-09-27-amendment-a2, after
cold panel round 2 (gpt_reviews/2026-09-27-nc-coldpanel-r2/meta_review.md, V1).

## What is wrong

1. `analysis/blca_common.cox_fit` builds each event's risk set as the suffix of the time-sorted rows,
   so subjects that share an event's time and sort before it are left out, and the objective depends
   on the (unstable) sort order. The manuscript states the Breslow partial likelihood. TCGA-BLCA has
   36 tied time values, 24 with events, covering 75 patients.
2. `analysis/blca_common.pct` breaks exact ties in NumPy's unspecified argsort order, which differs
   between NumPy builds (experiments/20260927-field-inflation/README.md, finding 2).

## What changes, and nothing else

1. `cox_fit`: every event at time t uses the full risk set {j : t_j >= t} (Breslow), computed from a
   stable sort. Tested before use: equal to an O(n^2) brute force, invariant to row permutation,
   identical to the old objective when no times tie, gradient equal to finite differences
   (gpt_reviews/2026-09-27-nc-coldpanel-r2/breslow-check/). The optimiser, its options, the penalty
   grids, the folds, seeds and everything else are unchanged.
2. `pct`: tied values receive their average rank, (average rank - 1) / (n - 1), computed without any
   sort-order dependence. Every local re-implementation of rank-based percentiles or Spearman ranks in
   code that feeds a reported number is switched to the same rule.
3. One comparator is added, not substituted: stacking whose ridge is chosen from {0.1, 1, 10, 100} by
   three-fold cross-validation over the training fold's inner out-of-fold percentiles, beside the
   existing stacking at ridge 1.0 (cold panel M5).

The switch is one environment variable read by `blca_common`, `BLCA_A2` (unset or 1 = corrected,
0 = original). Development-phase records (the iteration ledger, pre-freeze `s5_*` runs) are history
and are not regenerated; any development number the manuscript reports as evidence is regenerated.

## Procedure, per producer script, in dependency order

1. Run it with `BLCA_A2=0`. It must reproduce its committed output in every field but runtime (or its
   own known-answer check must pass). A script that fails here stops the amendment: the pipeline, not
   the functions, would be the cause of any later difference.
2. Run it with `BLCA_A2=1` into the same path. Known-answer checks that pin pre-A2 values are replaced,
   for this run only, by step 1's reproduction and by upstream A2 files where the script reads them.
3. The pre-A2 value of every bound quantity is kept by git (the commit before regeneration is recorded
   in the amendment's JSON) and tabulated beside its A2 value in the Supplementary Information.

## Reporting

Every reported number takes its A2 value. A conclusion that changes (an interval that starts or stops
excluding zero, a count such as "two of four", an ordering the text states) is rewritten and listed in
the amendment note. No threshold is moved. The amendment is described in Methods as a correction of
an implementation defect found in review, with its measured effect.

## Scope, fixed 2026-09-27 after the dependency inventory and before any A2 result was written

The inventory (every result file read by the number gate, the claim emitter, the Source Data builder,
the release verifier and the thirteen figure builders, traced to its producer) sorts the results into
three kinds.

1. Regenerated under A2: every post-freeze result reported as evidence whose producer is committed.
   The TCGA-BLCA chain (amendment A1's construction, the unified fusion and added-value analysis and
   its per-case vectors, the reporting dump, calibration and utility, re-partitions and held-out
   sites, the modality atlas, ablations, biology, competitor parity, the cold-panel round-1
   measurements, case interpretation, reviewer gaps, survival metrics, the SurvPath multiseed summary,
   the figure data, the selection null and encoder parity, the clinical-block decomposition), the
   five-study chain (stage5 in both of its committed copies, intervals, the second encoder, floors,
   fusion, subgroups), the GEO chain (the external evaluation, its per-case rerun), and the
   2026-09-27 analyses A to D and C (s33 to s36). The Cox-derived stage-minus-grade values that
   why-grade-fails.json carries without a producer are taken from the regenerated stage5.
2. Reported as recorded, with the original tie handling, and labelled so in the manuscript: the
   development campaign before the freeze (the pre-freeze s5_* results, the regime evidence and
   component slate curated from them, the iteration ledger and the frozen protocol's constants) and the
   frozen confirmatory run with its code snapshot (0.7260). These record the path and the decisions as
   they were taken; several have no committed producer or command, and some use their own solvers. The
   chain reported in the manuscript becomes frozen run 0.7260, amendment A1 0.7212, amendment A2 the
   regenerated value.
3. Not Cox- or rank-derived and unchanged: parameter counts, the published-values table, cohort
   counts, gate-0 availability counts, identifiers.

Local rank-correlation helpers in producers of kind 1 (s8_biology, s7_interpret_calibrate_cases,
s5_omics_arm, s5_step1_atlas as imported by s19b) and s5_omics_arm.score_test, which had the same
tied-risk-set defect as cox_fit, follow the same switch; under BLCA_A2=0 each is byte-for-byte the old
code path. Development scripts with their own solvers are not changed.
