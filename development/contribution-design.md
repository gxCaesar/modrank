# Contribution design — S5 slate, TCGA-BLCA WSI survival

2026-08-17. Every candidate here carries the five parts G-M6 requires — mechanism, predicted
slice, falsifier fixed before the run, removal ablation, matched control — and the ones that are
already dead keep theirs, because a kill is only informative if the thing it killed was specified.

`phase_of_origin: exploratory` throughout. **Nothing here is reportable.**

## What counts as one component

One named thing whose **removal can actually be run**. In this campaign that means: a change that
can be expressed as a parameter value the baseline already has (`s = 1`, `n_bag = 1`), or an arm
that can be dropped from the rank average, so the ablation is the same code path rather than a
second implementation that might differ for unrelated reasons.

## The budget, fixed before the first candidate ran

| | |
|---|---|
| total exploratory allocation | one working day of local CPU, plus whatever GPU time PIBD's rerun is not using |
| per-candidate first round | one full five-fold pass on the released folds, with its control, on this session's local machine |
| promotion | beat the candidate's own baseline arm by more than **0.0145** — 2 SD of the paired reseed delta measured on this code |
| abandonment | miss the bar on the first allocation, or need a third engineering rescue |
| rescues allowed | two per candidate |

## The bar, and why it moved

**0.0145.** Measured 2026-08-17 by re-splitting the same 359 cases 24 times under identical code
(`development/baseline-ledger.md`). This project had been carrying **0.0312**, which came from a
*simulation* at 412 patients and 181 events and does not describe this cohort. Both numbers are
kept visible so the change is auditable; every candidate below is judged against 0.0145.

---

## The measurement that pointed the slate — and the correction that moved it

**As first measured (2026-08-17, morning).** An in-sample oracle at 16 principal components reached
0.7666 for the image arm against an out-of-fold 0.6496, with a permuted-outcome control at 0.4980
over 1,000 draws. Read as `estimation_loss = 0.1170` and `recoverable = 0.2686`, this said the
image representation was rich and the estimator was failing to reach it. **Candidates C1, C1' and
C2–C4 were all drafted against that reading.**

**The control was wrong, and C2's falsified mechanism is what exposed it.** That control fitted on
permuted outcomes and scored on *real* ones — which measures whether two different targets are
related, not how much a fit is inflated by being scored on the very outcomes it was fitted to.
Sixteen parameters against ~23 events fit a great deal of *any* target. The corrected control fits
the permutation and scores **that same permutation**:

| arm | real in-sample | **corrected null** | old control | out-of-fold | real − null |
|---|---|---|---|---|---|
| TITAN d=16 | 0.7666 | **0.7474** | 0.4968 | 0.6496 | **+0.0192** |
| TITAN d=8 | 0.7067 | 0.6650 | 0.4959 | 0.6094 | +0.0417 |
| TITAN d=2 | 0.5730 | 0.5793 | 0.5026 | 0.5515 | −0.0063 |
| clinical d=16 | 0.7560 | 0.6990 | 0.4906 | 0.6741 | +0.0570 |
| clinical d=2 | 0.6882 | 0.5756 | 0.5014 | 0.6930 | **+0.1127** |

**A Cox fit to pure noise, scored in-sample at 16 parameters, reaches 0.7474.** So the reported
0.1170 estimation loss was the difference between an inflated number and an honest one, and the
image representation was never shown to hold more than it delivers. The excess *shrinks* from
+0.0417 at d=8 to +0.0192 at d=16, which is what added noise dimensions look like. The clinical
block's out-of-fold value (0.6930) exceeds its own in-sample value (0.6882), because out-of-fold
fits on 287 patients and in-sample on 72 — that arm has no estimation bottleneck at all.

This correction was made because a candidate's **pre-registered mechanism prediction failed**, not
because anyone doubted the arithmetic. That is what pre-registering the prediction bought.

---

## C1 — pan-cancer survival subspace as a prior on the image head · **DEAD**

**Shape.** A subspace constraint on the head's coefficient vector, estimated from out-of-domain
supervision.
**Gain family.** External structure prior.

**Mechanism.** Locating a survival direction in 768 dimensions from 113 events is the difficulty.
If the survival-relevant part of that space is partly shared across cancers, directions estimated
where events are plentiful cut the bladder head's effective dimension. **The sharing was verified
first, as its own falsifier:** a Cox direction fitted on colorectal slides, zero case overlap with
bladder, ranks bladder patients at **0.6492** against a permuted-label control of 0.4727 ± 0.032
(max over 20 draws 0.5235); the mean of four donor directions reaches 0.6211 against 0.4923. The
transfer is real.

**Predicted slice.** The image arm, and not the clinical arm.
**Falsifier (pre-registered).** Beat the plain ridge image arm (0.6596) by more than 0.0145.
**Removal ablation.** `s = 1`, the same code path with the prior off.
**Matched control.** Donor fits on permuted donor labels — same donors, same bootstraps, same k and
s, no survival signal.

**Result.**

| arm | C-index |
|---|---|
| plain ridge baseline | **0.6596** |
| pancancer_prior | 0.6554 (−0.0042) |
| permuted-donor CONTROL | 0.6455 (−0.0141) |
| unsupervised-subspace ablation | 0.6638 (+0.0042) |
| own-cohort-subspace ablation | 0.6548 (−0.0048) |

**Killed on its falsifier.** Every subspace arm sits within ±0.005 of the baseline: subspace
restriction does nothing in either direction, and even label-free dimension reduction beats the
donor-supervised subspace. The control was read first and behaved.

**Why, in a sentence that generalises:** the donor direction and the bladder head's own direction
have cosine **0.1748** — nearly orthogonal, of nearly equal strength. Projecting one onto the
other discards whichever is not kept.

## C1' — the pan-cancer direction as an ADDITIONAL arm · **DEAD**

**Shape.** An extra independent predictor combined at the score level. *Distinct from C1: C1
replaced the bladder head's direction, C1' keeps both.* That sentence is what buys it a budget
allocation under the slate's shape rule.
**Gain family.** External structure prior / second data source.

**Mechanism.** Two nearly orthogonal directions of comparable strength should be combined, not
substituted. The arm needs **no bladder training data at all**, which no other arm on this
benchmark can say.

**Falsifier (pre-registered).** `image + zero-shot` beats `image` by more than 0.0145, and the
permuted-donor control does not move.
**Removal ablation.** Drop the arm from the rank average.
**Matched control.** The same construction with every donor's labels permuted.

**Result.** zero-shot arm alone **0.6259**; permuted control **0.5531**; `image + zero-shot`
**0.6629** against image alone 0.6596 — **gain +0.0033 against a bar of 0.0145. Killed.** Adding it
to the best three-arm combination made that worse (0.7281 → 0.7245).

**One reading worth keeping:** the permuted-donor control reaches 0.5531, not 0.50. A direction
fitted to shuffled labels still lands on the embedding space's leading principal directions, and
those carry some association with survival in bladder. So the honest zero-shot signal is ~0.073,
not ~0.126 — which is exactly why the control exists.

**Family closed.** C1 and C1' share the shape *use non-bladder supervision to improve the bladder
image score*. A third of that shape is not funded.

## C2–C4 — the estimator library · **PARTIALLY ALIVE**

**Shape.** C2 resamples the estimator; C3 changes the penalty geometry; C4 fixes hyperparameter
selection noise. Three different objects, none of them the coefficient-constraint shape C1 died on.
**Gain family.** A different estimator, not a different representation.

**Mechanism (as drafted).** The deficit was measured as estimation variance in a high-dimensional
coefficient, so variance reduction should help.

**Predicted slice, written before the run and this is the matched control:** gains must scale with
block dimension — TITAN (768) > omics (275) > clinical (23). **A uniform gain falsifies the
mechanism even if the numbers rise.**

**Result.**

| block | baseline | C2 bagged | C3 elastic net | C4 dense grid |
|---|---|---|---|---|
| TITAN, 768-d | 0.6596 | **0.6035 (−0.056)** | 0.5480 | 0.6554 (−0.004) |
| omics, 275-d | 0.6510 | **0.6750 (+0.024)** | 0.5821 | 0.6676 (+0.017) |
| clinical, 23-d | 0.6856 | 0.6830 (−0.003) | 0.3417 | 0.6745 (−0.011) |

**The mechanism is falsified.** The gain ordering is the reverse of the prediction: bagging costs
the highest-dimensional block the most and gains only the middle one. This is what sent the oracle
back for re-examination, above.

**What survives is narrower and still useful:** *no single estimator is right for every block*.
Bagging is worth +0.024 on the omics block and −0.056 on the image block, so the estimator has to
be **chosen per block, inside the training fold**. That is the surviving component and it is what
the method arm implements.

**C3 is an implementation failure, not a result.** A 23-dimensional clinical Cox cannot honestly
score 0.3417; the proximal-gradient elastic net did not converge. It is recorded as such and is not
reported as a finding. One rescue remains available; it has not been spent, because sparse Cox is
not where the evidence points.

## C5 — gated fusion, image weight conditioned on clinical risk · **VOID ON ITS CONTROL**

**Shape.** A conditional combination rule.
**Gain family.** Changing the combination, not the inputs.

**Mechanism.** The error atlas measured image gain at **+0.1120** in the middle clinical-risk
tertile and **−0.0452** in the lowest, a spread of 2.9 fold-to-fold standard deviations. No fixed
weight can express that. The fusion oracle put a ceiling on it: gated 0.7492 against linear 0.7268,
so +0.0224, with the oracle's own control difference at −0.0029.

**Falsifier.** Beat the flat rank average by more than 0.0145.
**Matched control.** The identical gating structure driven by a **permuted** stratum label — same
number of strata, same sizes, only membership destroyed.

**Result.** flat 0.7447 · gated **0.7493** · **permuted-stratum control 0.7502**.

**The control scored higher than the real gating.** Per G-M6 the run is `status: void` and nothing
else in it is reported as a gated result. The structure alone buys whatever the gating appeared to
buy, and the atlas's tertile pattern does not survive being made into a component.

---

## What the surviving evidence supports

Three things, and only these:

1. **Per-modality heads beat jointly-trained cross-modal models at this sample size.** Published
   joint models run 0.612–0.679 on these folds; separate ridge Cox heads combined by rank average
   run higher.
2. **The estimator must be chosen per block, in-fold.** Measured directly: +0.024 on one block,
   −0.056 on another, from the same change.
3. **The combination rule should stay equal-weight.** The sweep peaks at w = 0.50, an oracle
   linear fusion is worth +0.008, and gating is void on its control.

## What it does NOT support, stated because it is the objection

The cheap clinical Cox reaches **0.6856**, above every published method on this benchmark. Most of
this campaign's lead over the published numbers is therefore *information* the competitors did not
use, not *method*. The comparison that decides the difference — our arm on **WSI + omics only**,
against DIMAF's 0.679 on the same folds — is declared in `s5_method.py` as comparison C and is
reported whichever way it lands.

---

## Machine-readable record

Read by **G-M6**. The prose above is the reasoning; this is the part a gate can check. Every
`observed_discrepancy` and `effect.observed` below was filled in after the run, and every
`falsifier`, `tolerance` and `effect.predicted` was fixed before it.

`unit: absolute` throughout — these are C-index differences, so no effect is 0, and the tolerance
is the paired split-reseed bar of 0.0145 measured on this code unless a component states otherwise.


### Why four components are `void` rather than `graded`

G-M6 marked C1, C1', C4 and C6 void and it is right to, but the reason differs between the two
pairs and both are worth stating.

**C1 and C4 — the artifact was larger than the signal.** C1's effect moved 0.0042 while its
permuted-donor control moved 0.0141; C4's moved 0.0042 against 0.0111. That is G-M2's arithmetic
one level in: a margin under the task's noise floor carries no information, and an effect under the
apparatus's own null movement carries none either. These two are not "small negative results" —
they are unresolvable on the apparatus that produced them.

**C1' and C6 — the control moved far more than its tolerance, and part of that was expected.**
Their controls add a *permuted* arm to an equal-weight ensemble, and adding a pure-noise arm to an
ensemble SHOULD hurt: dilution is the null's correct behaviour, not a malfunction. So −0.0311 and
−0.0699 are not evidence that the apparatus broke, and the informative quantity is the real arm
against its permuted twin (+0.0033 against −0.0311 for C1'). But G-M6's schema models a control as
an *equivalence pair* whose movement is pure artifact, with no effect at 0, and a dilution control
does not fit that shape. Rather than argue the gate down, the components are marked void: under
either reading neither produced an interpretable positive result, and both failed their own
pre-registered falsifiers regardless (+0.0033 and +0.0105 against a bar of 0.0145).

**The lesson for the next slate**, and it is a design lesson rather than a bookkeeping one: a
matched control should be built as two conditions that CANNOT differ, so that any gap is artifact
by construction. "The same thing with the signal removed" is a null, not an equivalence, and its
expected movement is not zero whenever the component enters an ensemble.

```yaml
components:
  - name: C1_pancancer_subspace_prior
    mechanism: "Locating a survival direction in a 768-dimensional slide embedding from 113 events
      is the difficulty. If the survival-relevant part of that space is partly shared across cancer
      types, directions estimated where events are plentiful cut the bladder head's effective
      dimension from 768 to k, and the estimation loss falls with it. The sharing was verified
      first: a Cox direction fitted on colorectal slides, zero case overlap, ranks bladder patients
      at 0.6492 against a permuted-label control of 0.4727 +/- 0.032."
    predicted_slice: collection_site=XF
    falsifier: "beat the plain ridge image arm (0.6596) by more than 0.0145"
    removal_ablation: "s = 1, the identical code path with the prior switched off"
    status: void
    matched_control:
      name: "the same donor cohorts, bootstraps, k and s, with every donor's survival labels permuted"
      why: "a coefficient vector fitted on ANY 800-patient matrix aligns with the leading principal
        directions of the embedding space, and those are shared across cancers for reasons -- scanner,
        stain, magnification, submitting site -- that have nothing to do with survival"
      unit: absolute
      tolerance: 0.0145
      observed_discrepancy: -0.0141
    effect:
      predicted: 0.03
      observed: -0.0042

  - name: C1p_pancancer_zeroshot_arm
    mechanism: "C1 REPLACED the bladder head's own direction with the donor subspace; the two are
      nearly orthogonal (cosine 0.1748) and of comparable strength (0.6492 against 0.6596), so
      projecting one onto the other discards whichever is not kept. Keeping both and combining at
      the score level should recover what C1 threw away. This arm needs no bladder training data at
      all, which no other arm on this benchmark can say."
    predicted_slice: collection_site=XF
    falsifier: "image + zero-shot beats image alone by more than 0.0145, and the permuted-donor
      control does not move"
    removal_ablation: "drop the arm from the rank average"
    status: void
    matched_control:
      name: "the identical construction with every donor's survival labels permuted"
      why: "same design matrices, same estimator, same aggregation; only the survival signal removed"
      unit: absolute
      tolerance: 0.0145
      observed_discrepancy: -0.0311
    effect:
      predicted: 0.025
      observed: 0.0033

  - name: C2_bagged_ridge
    mechanism: "The image arm's deficit was read as estimation variance in a high-dimensional
      coefficient, and bagging reduces estimation variance directly. The prediction that makes this
      more than a sweep: gains must scale with block dimension -- image 768 > omics 275 > clinical
      23 -- so a UNIFORM gain falsifies the mechanism even if the numbers rise."
    predicted_slice: collection_site=XF
    falsifier: "beat the plain ridge on the image block by more than 0.0145 AND show the
      dimension-ordered gain pattern"
    removal_ablation: "n_bag = 1, which is the plain ridge"
    status: graded
    matched_control:
      name: "the clinical block, 23-d, is the condition that must NOT move"
      why: "whatever fixes high-dimensional estimation variance cannot help a 23-dimensional Cox as
        much as a 768-dimensional one; the control is built into the design rather than bolted on"
      unit: absolute
      tolerance: 0.0145
      observed_discrepancy: -0.0026
    effect:
      predicted: 0.04
      observed: -0.0561

  - name: C4_dense_regularisation_path
    mechanism: "The baseline picks among five alphas spanning 1 to 4096 by one 3-fold inner CV on
      287 patients. A coarse grid plus a noisy criterion is a way to lose performance that looks
      like a modelling limit rather than a selection one."
    predicted_slice: collection_site=XF
    falsifier: "beat the coarse-grid ridge on the image block by more than 0.0145"
    removal_ablation: "revert to the five-alpha grid and 3 inner folds"
    status: void
    matched_control:
      name: "the clinical block, which has almost no alpha sensitivity, must not move"
      why: "a 23-dimensional Cox at 287 patients is insensitive to the penalty over this range"
      unit: absolute
      tolerance: 0.0145
      observed_discrepancy: -0.0111
    effect:
      predicted: 0.02
      observed: -0.0042

  - name: C5_gated_fusion_on_clinical_risk
    mechanism: "The atlas measures image gain at +0.1120 in the middle clinical-risk tertile and
      -0.0452 in the lowest, a spread of 2.9 fold-to-fold standard deviations, and no fixed weight
      can express that. An oracle gated fusion reached 0.7492 against a linear oracle's 0.7268."
    predicted_slice: clinical_risk_tertile=clin_low
    falsifier: "beat the flat rank average by more than 0.0145"
    removal_ablation: "one stratum for everybody, which is the flat rank average"
    status: void
    matched_control:
      name: "the identical gating structure driven by a PERMUTED stratum label"
      why: "same number of strata, same sizes; only which patient is in which is destroyed, so the
        two differ in the clinical information and in nothing else"
      unit: absolute
      tolerance: 0.0145
      observed_discrepancy: 0.0055
    effect:
      predicted: 0.02
      observed: 0.0046

  - name: C6_inner_cv_skill_weighting
    mechanism: "Equal weighting is optimal among forecasters of COMPARABLE accuracy, which is where
      the flat weight sweep's optimum at exactly 0.50 came from -- the image and clinical arms
      differ by 0.026. It says nothing about a pool whose members range from 0.5224 to 0.6856. A
      view at 0.55 given the same weight as one at 0.69 is one seventh of the ensemble spent on
      something barely above chance. Weight each view by max(0, C_inner - 0.5), estimated inside the
      training fold; equal weighting is its special case."
    predicted_slice: collection_site=XF
    falsifier: "beat the uniform rule on the SAME pool by more than 0.0145, on the seven-view pool
      AND without losing on the three-view pool where uniform is at its best"
    removal_ablation: "uniform weights, w = 1 for every view"
    status: void
    matched_control:
      name: "weights computed from inner CV on a PERMUTED training outcome"
      why: "the permuted skills have the same shape and spread of magnitudes; only which view earns
        which weight is destroyed"
      unit: absolute
      tolerance: 0.0145
      observed_discrepancy: -0.0699
    effect:
      predicted: 0.02
      observed: 0.0105

  - name: C7_stratified_partial_likelihood
    mechanism: "Cox's partial likelihood compares each event against everyone still at risk, so a
      molecular head can lower its loss by re-deriving stage -- the easiest thing in the image --
      and part of its fitted capacity goes there. Restricting each risk set to patients in the SAME
      clinical stratum removes that route: within a stratum stage carries no information, so the
      only way to lower the loss is with signal the clinical model does not have."
    predicted_slice: clinical_risk_tertile=clin_low
    falsifier: "beat the marginally-trained arm ON CLINICALLY TIED PAIRS by more than 0.0145. Its
      overall C-index is allowed to fall -- that is what the mechanism predicts, and a rise in both
      would be a different result than the one hypothesised."
    removal_ablation: "one stratum for everybody, which is the ordinary partial likelihood; verified
      by known-answer test to reproduce the unstratified fit to 0.000e+00"
    status: graded
    matched_control:
      name: "strata permuted across patients inside the training fold"
      why: "same number of strata, same sizes, membership destroyed"
      unit: absolute
      tolerance: 0.0145
      observed_discrepancy: -0.0095
    effect:
      predicted: 0.03
      observed: -0.0199
```
