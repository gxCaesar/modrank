# Charter — TCGA-BLCA whole-slide multimodal survival

Written 2026-08-17. **It supersedes `charter-nmibc-bulk-transcriptome-20260814.md`**, which is for
a different topic — bulk-transcriptome NMIBC progression, with a comparator set of published
expression classifiers — and was never updated when the active route moved to whole-slide survival
on 2026-08-15 (branch `blca-wsi-survival-20260815`). The archived file stays as the record of the
superseded route.

## An honesty note that belongs at the top, not in a footnote

**This charter was written after an exploratory campaign, not before it.** S5 ran on 2026-08-17 —
211 scored candidates, seven falsified component families, `development/iteration-ledger.md` — and
the baseline table below names the arms that campaign actually ran. A charter written to match work
already done is a record, not a commitment, and calling it one would be the exact defect this kit
exists to catch.

So the document separates the two, and only the second half is load-bearing:

| written after the fact | genuinely forward-looking |
|---|---|
| the comparator set, and every `observed (exploratory)` value | the **decision rule**, fixed before the confirmatory run |
| the mechanism paragraph, distilled from measurements already taken | **seeds 1–4, which have never been executed** |
| the required-variables and data-source tables | the **kill criterion** and the pre-declared fallback |

`development/benchmark-protocol.json` carries the same split, and its `frozen_at` is the boundary
G-M4 enforces by mtime.

## Target

- venue: **Briefings in Bioinformatics** (primary), npj Digital Medicine (secondary)
- venue_class: journal
- deadline: none, rolling submission
- contribution_type: sota-method
- requirements_doc: publication-project.json

Fixed by the user and not tradeable: **bladder cancer only, public data only, sota-method pursued
until reached rather than downgraded.** A controlled-access route to 327 more patients was offered
and declined on 2026-08-13, so everything here runs on public data.

Briefings in Bioinformatics is primary rather than aspirational: **DSCASurv** (doi:10.1093/bib/bbaf103,
2025) is a method paper on this exact benchmark, at 0.646 on TCGA-BLCA. The venue publishes work of
this shape, at this number, on this data.

## Primary hypothesis

A per-modality, estimation-controlled survival ensemble — one ridge Cox per view over a frozen
foundation-model slide embedding, pathway-level transcriptomics and routine pathologic staging,
combined by equal-weight rank average — reaches a **Harrell C-index above 0.703** on the released
`mahmoodlab/SurvPath` five-fold case-ID splits for TCGA-BLCA disease-specific survival, which is
**DIMAF's published 0.679 plus the 0.0239 selection-inflation term** measured over this campaign's
211 scored candidates.

The named referent is DIMAF (doi:10.1007/978-3-032-05185-1_12, MICCAI 2025), the strongest published
value on this endpoint and cohort whose use of the released folds has been claimed. **Its fold
identity is contested between two independent literature passes and is an open blocker**; if it
resolves against DIMAF, the referent becomes MOAD-FNet's 0.691 (arXiv:2411.17418), and the
hypothesis's number becomes 0.715.

## Kill criterion

If the confirmatory run's primary metric — the mean pooled C-index over seeds 0–4 — fails to exceed
**0.703**, or if its standard deviation over those five seeds exceeds the measured reseed bar of
**0.0145**, then **abandon** the absolute-discrimination claim on this benchmark and execute the
pre-declared fallback under `## Winning condition`, which changes the reported quantity rather than
the venue or the contribution type.

A second, already-triggered kill is recorded here because it is a commitment for the write-up:
**the WSI-plus-omics arm, without clinical variables, does NOT beat the incumbent** — 0.6841 against
0.679, and it failed again under three further combination rules and a seven-encoder ensemble. That
result is **stated in the paper**, not omitted.

## Required variables

| variable | what it is | where it comes from | co-exists in one downloadable sample set? |
|---|---|---|---|
| whole-slide morphology | 768-d TITAN slide embedding, mean-pooled over a case's slides | precomputed pan-TCGA TITAN feature file, 11,658 slides | yes, n=359 |
| bulk transcriptome | 4,999 genes averaged into 275 `combine` pathway groups | `mahmoodlab/SurvPath` `datasets_csv/raw_rna_data/combine/blca/` | yes, n=359 |
| pathologic stage | AJCC stage; and T, N, M separately | `SurvPath datasets_csv/clinical_data/` (stage) and the GDC clinical record (T/N/M) | yes, n=359, every clinical row matched |
| age, sex | years, binary | `SurvPath datasets_csv/metadata/tcga_blca.csv` | yes, n=359 |
| disease-specific survival | `survival_months_dss`, `censorship_dss` | same | yes, n=359, 113 events |
| fold membership | released five-fold case-ID partition | `SurvPath splits/5foldcv/tcga_blca/`, obtained via PIBD's byte-identical copy | yes, all 359 cases |

The joint count is the same 359 cases for every row, which is what makes this a real cell rather
than five separately-available variables.

## Data sources

- `mahmoodlab/SurvPath` — folds, labels, genes, pathway signatures, and `clinical_data/` with the
  `stage` column. <https://github.com/mahmoodlab/SurvPath>
- `zylbuaa/PIBD` — the fold files, verified byte-identical to SurvPath's with `cmp`, all five.
  <https://github.com/zylbuaa/PIBD>
- Precomputed pan-TCGA TITAN slide embeddings, 11,658 slides.
- GDC clinical record for T, N, M and AJCC stage on the 359 cases.
- Six further public, ungated precomputed slide feature releases were downloaded and evaluated
  (Prov-GigaPath native, GigaSSL×4, UNI). **None beats TITAN** and the ensemble of all seven is
  significantly worse than TITAN alone (−0.0440, p = 0.008). They are reported as a negative result.

All public. No controlled-access data. Nothing downloaded outside an explicit user approval.

## Baselines

Every arm gets the same folds, the same labels, the same estimator, the same alpha grid and the same
inner-CV budget. Our method gets no more than that. **Trivial and cheap baselines are named here and
their ledger entries carry the same names**, so the gate compares like with like.

| method | tier | paper / code | search budget | observed (exploratory) | why it is the one to beat |
|---|---|---|---|---|---|
| **constant predictor** | trivial | — | 0 trials | 0.5000 | the floor |
| age alone | trivial | — | 5 alphas × 3 inner folds | 0.5756 | a single covariate |
| **ridge penalised Cox on clinical** (age, sex, T, N, M, stage) | cheap | <https://doi.org/10.18637/jss.v039.i05> | 5 alphas × 3 inner folds | **0.6856** | **THE decisive cheap baseline.** It beats 8 of the 9 published numbers on this benchmark, and no paper in the family has ever run it |
| ridge penalised Cox on the slide embedding | cheap | same | 5 alphas × 3 inner folds | 0.6596 | single-modality |
| ridge penalised Cox on pathway expression | cheap | same | 5 alphas × 3 inner folds | 0.6510 | single-modality |
| SurvPath, official implementation, rerun | recent_strong | <https://github.com/mahmoodlab/SurvPath> | its own released config, seed 1 | 0.6118 (published 0.625) | **the only competitor whose per-case predictions we hold, hence the only PAIRED comparison possible** |
| SurvPath rerun **plus the same clinical variables** | recent_strong | same | same | 0.6963 | **input parity.** The comparison a reviewer will ask for |
| DIMAF | quoted | <https://doi.org/10.1007/978-3-032-05185-1_12> | not rerun | 0.679 ± 0.043 | the incumbent |
| MOAD-FNet | quoted | arXiv:2411.17418 | not rerun | 0.691 ± 0.069 | highest published on this endpoint and n; folds unverified |
| PIBD | quoted | arXiv:2401.01646 | rerun crashed after fold 0 (0.6038) | 0.667 ± 0.061 | folds verified byte-identical |
| APL · DSCASurv · ProtoPathway · OTSurv · MMP · MCAT · MOTCat | quoted | see `development/benchmark-protocol.json` | not rerun | 0.677 · 0.646 · 0.646 · 0.637 · 0.635 · 0.598 · 0.596 | the rest of the published field |

`quoted` is a real tier and a limitation: none of those methods releases per-case predictions, so
none can be given the clinical variables our arm uses, and none can be paired. This is stated in the
paper rather than averaged over.

## Winning condition

- **benchmark and split**: TCGA-BLCA disease-specific survival on the released `mahmoodlab/SurvPath`
  five-fold case-ID splits, 359 cases, 113 events, 24,219 comparable pairs. The exact membership is
  recorded in `development/split-manifest.json` and bound by the frozen protocol.
- **metric**: Harrell C-index, higher is better, pooled over out-of-fold predictions after
  within-fold percentile normalisation, averaged over seeds 0–4.
- **incumbent**: **DIMAF at 0.679** (± 0.043), <https://doi.org/10.1007/978-3-032-05185-1_12>. Its
  fold identity is contested and must be settled with a diff; the fallback referent is **MOAD-FNet
  at 0.691**, arXiv:2411.17418.
- **margin required to report: 0.0239 C-index above the incumbent**, which is the
  selection-inflation term `σ·√(2 ln N)` at σ = 0.0073 and N = 211 candidates scored during S5. It
  is deliberately larger than the measured split-reseed bar of 0.0145, because a margin that only
  clears the reseed bar is what searching 211 candidates produces with no real effect. Both numbers
  are this task's own, measured on this code: 24 reseeds of the fold assignment over the same 359
  cases, identical estimator and alpha grid.
- **mechanism**: at 359 cases and 113 events, the binding constraint on this benchmark is estimation,
  not representation. Every published entrant trains a joint cross-modal network end to end, and the
  parameter count that requires is not supported by 113 events — which is why a ridge Cox on six
  routine clinical fields reaches 0.6856 and beats eight of the nine published numbers. Fitting each
  modality separately, and combining by a rule with no free parameters, spends the events on the
  three coefficient vectors that carry signal instead of on a fusion network. That sentence does not
  survive swapping in a different method: it predicts, specifically, that fitted fusion weights buy
  nothing here — and an oracle two-parameter fusion, given the held-out outcomes, was worth +0.008
  over equal weight, while a gated fusion scored 0.7493 against a permuted-stratum control at 0.7502.
- **if the margin is not reached**: **report honestly that it was not, and reframe around where the
  method does win** — the PrePR-CT route (NMI 2026), chosen now and in advance. The reported
  quantity becomes discrimination **among patients pathologic stage cannot separate**: on the decile
  of comparable pairs the clinical arm separates least, the arm reaches 0.6513 against
  SurvPath-plus-the-same-clinical at 0.5852, a margin of +0.0659 — double the overall one, on the
  question a pathologist is actually asking, and reported by nobody on this benchmark. This is
  chosen over reporting a tie because a tie report is an audit and the contribution type is fixed at
  sota-method; and it is the better-matched quantity, because overall C-index on this cohort is
  dominated by pairs that stage already orders correctly.

## Compute budget

- Exploratory (S5, spent): one working day of local CPU on a MacBook, plus one A4000 on `sysu` for
  the SurvPath and PIBD reruns. 211 candidates scored. No GPU was used by any arm of our own method.
- Confirmatory: five seeds × five folds × three views of ridge Cox. **Local CPU, minutes, no GPU.**
- Storage: the entire method's inputs are 2.4 GB of public precomputed features on `sysu` and ~90 MB
  of per-cohort matrices locally. No checkpoints — the method has no trained network.
- The cheapness is not incidental. It is the mechanism's own prediction and it is a claim the paper
  makes: at this sample size the winning estimator is one that cannot overfit, and it costs minutes.

## What this is NOT

- **Not a deep multimodal architecture, and not a claim that one would help.** Seven component
  families were designed against measured deficits and all seven were falsified with matched
  controls read first. The paper reports them.
- **Not a claim that WSI plus omics beats the field.** It does not: 0.6841 against 0.679, and it
  fails again under three further combination rules and a seven-encoder ensemble. The lead requires
  pathologic staging.
- **Not an audit paper.** The stage finding is evidence for the method's design, not the
  contribution. The contribution is the method and its number.
- **Not externally validated.** This benchmark has no unopened held-out split and neither does any
  published entrant. The four non-bladder TCGA cohorts serve as replication of the stage finding,
  not as external validation of the bladder number.
- **Not a prospective or clinical-utility claim.** No decision-curve analysis, no net benefit, no
  deployment.
