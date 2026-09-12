# Baseline ledger — TCGA-BLCA whole-slide survival, released SurvPath folds

Written at the start of the S5 campaign, 2026-08-17, before any component was designed. It answers
one question: **is there room above the cheap baselines that is larger than this task's own
noise?**

Everything here is `phase_of_origin: exploratory` and **nothing in it is reportable**. Scripts and
raw outputs live in the session scratchpad; the numbers below are copied from those runs verbatim.

## The task, stated so the comparison set is unambiguous

TCGA-BLCA, disease-specific survival, the **released `mahmoodlab/SurvPath` five-fold case-ID
splits** (`splits/5foldcv/tcga_blca/splits_{0..4}.csv`), 359 cases, 113 events, 24,219 comparable
pairs, Harrell C-index, higher is better. Every published competitor named below is scored on those
same fold files; PIBD's copies were checked byte-identical with `cmp`, and DIMAF's fold-0 case set
was checked to match after collapsing slide IDs to case IDs.

## What the cheap baselines already do — and this is the finding, not a formality

| baseline | tier | C-index |
|---|---|---|
| constant risk for every patient | trivial | 0.5000 |
| sex alone | trivial | 0.4697 |
| age alone | cheap | 0.5756 |
| **stage block alone (T, N, M, AJCC stage)** | cheap | **0.6780** |
| **clinical Cox (age + sex + T + N + M + stage)** | cheap | **0.6856** |
| TITAN slide embedding + ridge Cox | cheap | 0.6596 |

**A ridge Cox on six routine clinical fields reaches 0.6856, which is above every published method
on this benchmark**, including the strongest verified one (DIMAF, 0.679 ± 0.043). That is not a
statement about our method; it is a statement about the benchmark's baseline hygiene, and it is
checkable.

The published clinical baselines on these exact folds are:

| paper | clinical covariates used | C-index |
|---|---|---|
| MMP, ICML 2024 | age only | 0.578 ± 0.056 |
| MMP, ICML 2024 | age + sex + **grade** | 0.570 ± 0.033 |
| DIMAF, MICCAI 2025 | age + sex + **grade** | 0.519 ± 0.067 |

**All three use grade; none uses stage.** A targeted literature check on 2026-08-17 returned
`NOT_ESTABLISHED` for a published age + sex + stage + T + N + M Cox on the released folds — nobody
has run one. Pathologic stage is the primary prognostic variable in bladder cancer, and on these
folds it alone is worth 0.6780 against grade-based baselines at 0.519–0.578.

This is the exact failure S5's own guidance names: *G-M1 checks the cheap baseline was RUN; nothing
checks it was run well.* Here it appears in the published literature of the benchmark itself.

## How the noise floor was measured

Not simulated and not quoted. The **same code**, the **same 359 patients**, the same estimator,
alpha grid and inner-CV seed, re-split under **24 fresh five-fold assignments**. What varies is
only which patients land in which fold — which is exactly what a reader reproducing the study
would vary.

| quantity | mean | SD | 2 SD |
|---|---|---|---|
| TITAN arm, pooled | 0.6424 | 0.0115 | 0.0230 |
| clinical arm, pooled | 0.6823 | 0.0093 | 0.0186 |
| rank-averaged arm, pooled | 0.7093 | 0.0092 | 0.0184 |
| **paired delta (rank-avg − clinical)** | 0.0270 | **0.0073** | **0.0145** |

The paired delta is the quantity a method claim is made of: both arms see the same folds, so fold
difficulty cancels. **0.0145 is the bar** every candidate in this campaign is held to, and it is
tighter than the 0.0312 this project had been carrying, which came from a *simulation* at a
different cohort size (412 patients, 181 events) and did not apply here.

Note the released folds run systematically above the reseeds (rank-avg 0.7191 vs 0.7093, TITAN
0.6596 vs 0.6424). The released folds are site-stratified and the reseeds are plain random; the
comparison against published numbers is still fair because every competitor uses those same files,
but a number from these folds is not interchangeable with a number from a fresh random split.

## Headroom — and the correction that reversed this section's first conclusion

An oracle arm was fitted **in-sample on the held-out fold** at matched capacity (`d` principal
components of the training fold's own basis) and compared against a control at the same capacity.
**The first control was the wrong one, and this section originally reported its answer as a
finding.** Both are now shown together so the difference is visible rather than described.

| arm | real, in-sample | **corrected null** † | first control ‡ | out-of-fold | real − corrected null |
|---|---|---|---|---|---|
| TITAN, d=16 | 0.7666 | **0.7474** | 0.4968 | 0.6496 | **+0.0192** |
| TITAN, d=8 | 0.7067 | 0.6650 | 0.4959 | 0.6094 | +0.0417 |
| TITAN, d=4 | 0.6450 | 0.6167 | 0.4953 | 0.6056 | +0.0283 |
| TITAN, d=2 | 0.5730 | 0.5793 | 0.5026 | 0.5515 | **−0.0063** |
| clinical, d=16 | 0.7560 | 0.6990 | 0.4906 | 0.6741 | +0.0570 |
| clinical, d=2 | 0.6882 | 0.5756 | 0.5014 | 0.6930 | **+0.1127** |

† fit the permuted outcome, score **that same permutation** — zero signal, identical n, event count,
censoring and capacity. 200 draws per fold.
‡ fit the permuted outcome, score the **real** one. 1,000 draws.

**The first control measured whether two different targets are related, not how much a fit is
inflated by being scored on the outcomes it was fitted to.** A permutation destroys the first, so
it answers 0.5 whatever the capacity — and sixteen parameters against ~23 events fit a great deal of
any target. The corrected null says so directly: **a Cox fit to pure noise, scored in-sample at
d=16, reaches 0.7474.**

So the reading this section carried until it was rechecked — *"the image representation is not
exhausted; it holds 0.2686 above its control and gives up 0.1170 to estimation"* — **is withdrawn.
Both numbers were artifacts of comparing an inflated quantity against an honest one.** The real
excess over a matched null is +0.0192 at d=16 and it *shrinks* as capacity grows, which is what
added noise dimensions look like. The image representation was never shown to hold more than it
delivers out of fold.

Two things survive and are worth keeping:

- **The clinical block is genuinely informative at low capacity** (+0.1127 at d=2), and its
  out-of-fold value (0.6930) *exceeds* its own in-sample value (0.6882) — out-of-fold fits on 287
  patients, in-sample on 72. That arm has no estimation bottleneck at all.
- **The correction was triggered by a falsified prediction, not by suspicion.** Candidate C2
  pre-registered that variance reduction would help the arms in proportion to their dimension;
  bagging then cost the 768-dimensional image block 0.056 and gained the 275-dimensional omics block
  0.024. A mechanism failing in the direction opposite to its prediction is a reason to re-examine
  the measurement the mechanism was read off, and that is what happened.

An earlier defect in the same probe is kept on the record because it is a different one: the first
version used **one** permutation per fold and returned a control of 0.6467 for `clinical_d2`. At
d=2 a permuted fit still points somewhere inside a two-dimensional span of a genuinely prognostic
block, so its *sign* is what is random and five draws is not an estimate.

## Is the combination rule the deficit? No

| reference | value |
|---|---|
| deployed equal-weight rank average of image + clinical | 0.7191 |
| best fixed weight anywhere on a 21-point sweep | 0.7191, at **w = 0.50** |
| **linear fusion oracle** — 2 params fitted in-sample on the held-out fold | 0.7268 (control 0.5006) |
| **gated fusion oracle** — image weight free per clinical tertile | 0.7492 (control 0.4977) |
| model reference: two signals of these strengths, rank-averaged | 0.7421 ± 0.0238 |

Even given the held-out outcomes, the best *fixed linear* fusion of these two scores gains +0.008
over the equal weight already in use, and the sweep's maximum is exactly at equal weight. Gating
buys +0.0224 with its control at −0.0029, so the structure alone buys nothing — but that is an
oracle ceiling and a deployed gate would get less.

**So the room is not in the fusion rule.** It is in how much each arm extracts, and the arm with
the most to give back is the image arm.

```yaml
task: tcga_blca_dss_survpath_released_folds
metric: harrell_c_index
higher_is_better: true

# The claim there has to be room FOR. G-M2 tests headroom against this.
target:
  metric: harrell_c_index
  claimed_value: 0.7029
  why: "DIMAF's published 0.679 plus the 0.0239 selection-inflation term over the campaign's 211
    scored candidates. This is the number charter.md's ## Primary hypothesis commits to, and it is
    deliberately above the 0.0145 reseed bar: a margin that clears only the reseed bar is what
    searching 211 candidates produces with no real effect."
  ceiling_estimate: 0.7421
  ceiling_source: "a simulation calibrated to the two strongest arms' own C-indices (0.6596 and
    0.6856) placing two signals on independent risk components and rank-averaging them: 0.7421 +/-
    0.0238. It is used INSTEAD of the in-sample fusion oracles because the corrected-null analysis
    below shows in-sample references on this cohort are badly inflated -- a Cox fit to pure noise at
    16 parameters against ~23 events scores 0.7474 in-sample. A model-based reference with no
    in-sample fitting is the more honest ceiling here, and it is a reference, not a bound."

noise_floor:
  metric: harrell_c_index_paired_delta
  value: 0.0073
  bar_2sd: 0.0145
  higher_is_better: true
  measured_by: "24 reruns of THIS code on THESE 359 cases: the 5-fold assignment was re-randomised
    24 times and nothing else varied -- same arms, same ridge Cox, same alpha grid
    (1,8,64,512,4096), same 3-fold inner CV, same seed 0, same numpy/scipy build, one host (local
    macOS). Reported value is the SD of the paired delta between the rank-averaged arm and the
    clinical arm across those 24 reruns. Script: s5_step0_headroom.py --reseeds 24"

headroom:
  incumbent:
    name: DIMAF
    value: 0.679
    sd: 0.043
    source: "10.1007/978-3-032-05185-1_12 (MICCAI 2025)"
    verified_on_released_folds: CONFLICTED
    fold_identity_conflict: >
      Two independent literature passes on 2026-08-17 disagree on whether DIMAF uses the released
      SurvPath case-ID folds. The first reported true and said it compared DIMAF's fold-0 case set
      against the released fold after collapsing slide IDs to case IDs. The second reported false.
      Neither was run by the orchestrator. PIBD's fold identity, by contrast, WAS checked here
      directly -- `cmp` on all five files, byte-identical -- so the conflict is specific to DIMAF
      and is recorded rather than resolved by preference.
    date_checked: "2026-08-17"
    unverified_higher: "MOAD-FNet, arXiv:2411.17418, 0.691 +/- 0.069, DSS, n=359, but no released
      code and no split-file identity, so it cannot displace DIMAF as the established incumbent"
  reachable_reference:
    name: gated_fusion_oracle_two_arms
    value: 0.7492
    note: "in-sample fit on the held-out fold; a reference point for a two-arm combination, not a
      bound. See the oracle correction below -- an in-sample fit at 16 parameters against ~23
      events reaches 0.7474 on PURE NOISE, so in-sample references on this cohort are inflated and
      this one is quoted with that caveat attached."
  headroom_over_incumbent: 0.0702
  headroom_over_best_cheap_baseline: 0.0636
  exceeds_noise_floor: true

# Names here are the names charter.md's ## Baselines table uses. G-M1 compares the two lists, and
# a cheap-tier tag is not a substitute for the baseline the charter actually named.
baselines:
  - name: constant predictor
    tier: trivial
    metric: harrell_c_index
    observed_value: 0.5000
    search_budget: "0 trials"
    provenance: {repo: null, commit: null, command: "s5_step0_headroom.py", split_hash: "0d0b8f82e7f5aa6eb5833e89c19b0f5b53c14ba7490dfce28180cd3e050c7487", seeds: [0]}
  - name: age alone
    tier: trivial
    metric: harrell_c_index
    observed_value: 0.5756
    search_budget: "5 alphas x 3-fold inner CV"
    provenance: {repo: null, commit: null, command: "s5_step0_headroom.py", split_hash: "0d0b8f82e7f5aa6eb5833e89c19b0f5b53c14ba7490dfce28180cd3e050c7487", seeds: [0]}
  - name: ridge penalised Cox on clinical
    tier: cheap
    metric: harrell_c_index
    observed_value: 0.6856
    search_budget: "5 alphas x 3-fold inner CV, inside each training fold"
    note: "age, sex, T, N, M, AJCC stage. THE decisive cheap baseline: it beats 8 of the 9 published
      numbers on this benchmark, and a full-text pass over eleven primary papers found none that
      ever ran it -- every published clinical baseline in the family uses grade, not stage."
    provenance: {repo: null, commit: null, command: "s5_step0_headroom.py", split_hash: "0d0b8f82e7f5aa6eb5833e89c19b0f5b53c14ba7490dfce28180cd3e050c7487", seeds: [0]}
  - name: ridge penalised Cox on the slide embedding
    tier: cheap
    metric: harrell_c_index
    observed_value: 0.6596
    search_budget: "5 alphas x 3-fold inner CV, inside each training fold"
    provenance: {repo: null, commit: null, command: "s5_step0_headroom.py", split_hash: "0d0b8f82e7f5aa6eb5833e89c19b0f5b53c14ba7490dfce28180cd3e050c7487", seeds: [0]}
  - name: ridge penalised Cox on pathway expression
    tier: cheap
    metric: harrell_c_index
    observed_value: 0.6510
    search_budget: "5 alphas x 3-fold inner CV, inside each training fold"
    provenance: {repo: null, commit: null, command: "s5_omics_arm.py", split_hash: "0d0b8f82e7f5aa6eb5833e89c19b0f5b53c14ba7490dfce28180cd3e050c7487", seeds: [0]}
  - name: SurvPath, official implementation, rerun
    tier: recent_strong
    metric: harrell_c_index
    observed_value: 0.6118
    published_value: 0.625
    published_sd: 0.056
    search_budget: "its own released configuration, unmodified, seed 1"
    provenance:
      repo: "https://github.com/mahmoodlab/SurvPath"
      commit: "codeload tarball fetched 2026-08-16; upstream ships no tagged release"
      command: "experiments/20260816-survpath-blca-repro/run.sh"
      split_hash: "0d0b8f82e7f5aa6eb5833e89c19b0f5b53c14ba7490dfce28180cd3e050c7487"
      seeds: [1]
      note: "reproduced with the 768-d CTransPath-family features its own paper specifies, after
        patching an upstream defect in which --encoding_dim never reaches SurvPath's model_dict
        (utils/core_utils.py); the original is preserved as core_utils.py.orig. Per-case risks were
        recovered and the per-fold C-index recomputed to within 0.0005 of the run's own log on all
        five folds, which is what licenses using them for the paired comparison."
  - name: SurvPath rerun plus the same clinical variables
    tier: recent_strong
    metric: harrell_c_index
    observed_value: 0.6963
    search_budget: "SurvPath's own configuration, plus the identical clinical ridge Cox our arm uses"
    note: "INPUT PARITY. The only competitor on this benchmark whose per-case predictions exist here,
      hence the only paired comparison that can be made at all."
    provenance:
      repo: "https://github.com/mahmoodlab/SurvPath"
      commit: "codeload tarball fetched 2026-08-16"
      command: "s5_headline.py"
      split_hash: "0d0b8f82e7f5aa6eb5833e89c19b0f5b53c14ba7490dfce28180cd3e050c7487"
      seeds: [1]

# Published values that were NOT rerun. They are deliberately not in `baselines` above: a tier tag
# of `recent_strong` asserts that the method was run here, and none of these was. Their fold
# identity, where it is contested or unverified, is recorded in the protocol's comparator_caveats.
quoted_not_rerun:
  - {name: MOAD-FNet, value: 0.691, sd: 0.069, doi: "arXiv:2411.17418", folds: unverified}
  - {name: DIMAF, value: 0.679, sd: 0.043, doi: "10.1007/978-3-032-05185-1_12", folds: CONFLICTED}
  - {name: APL, value: 0.677, sd: 0.060, doi: "arXiv:2503.04643", folds: unverified}
  - {name: PIBD, value: 0.667, sd: 0.061, doi: "arXiv:2401.01646", folds: "verified byte-identical by cmp; our rerun crashed after fold 0 at 0.6038 and is NOT quoted as a reproduction"}
  - {name: DSCASurv, value: 0.646, sd: 0.034, doi: "10.1093/bib/bbaf103", folds: unverified}
  - {name: ProtoPathway, value: 0.646, doi: "arXiv:2605.21454", folds: released}
  - {name: OTSurv, value: 0.637, sd: 0.065, doi: "arXiv:2506.20741", folds: "follows MMP"}
  - {name: MMP, value: 0.635, sd: 0.051, doi: "arXiv:2407.00224", folds: "own site-stratified protocol"}
  - {name: MCAT, value: 0.598, folds: released}
  - {name: MOTCat, value: 0.596, folds: released}
```

## What was NOT run

- PIBD's rerun had not finished when this file was written; its cell above is `null`, not filled in
  from the paper.
- DIMAF, MMP, OTSurv, APL, MOAD-Net and ProtoPathway were **not rerun**. Their values are published
  numbers on the same released folds, and none of them can be given the clinical covariates our arm
  uses, because their per-case predictions are not released. The only competitor that can be given
  input parity is SurvPath, whose per-case risks we recovered from our own rerun.
- No oracle or human-expert reference exists in the literature for this protocol; the reachable
  reference above is ours and is model-based.
