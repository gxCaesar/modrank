# Pathologic stage is in this benchmark's own release, and no paper in it has used one

2026-08-17. `phase_of_origin: exploratory`. **Not reportable** — every number here needs a frozen
confirmatory run before it can carry a claim. Scripts and raw JSON are in the session scratchpad;
the ledger of all 211 scored candidates is `development/iteration-ledger.md`.

## The claim, in one paragraph

Eleven primary papers define the TCGA multimodal survival benchmark that `mahmoodlab/SurvPath`
established — SurvPath, MMP, PIBD, DIMAF, MCAT, MOTCat, OTSurv, APL, MOAD-FNet, DSCASurv,
ProtoPathway. A full-text pass over all eleven on 2026-08-17 found **no clinical-variable baseline
anywhere that uses pathologic stage, overall stage, or T/N/M.** The only clinical covariates
reported are age, sex and cancer grade. Meanwhile SurvPath's own repository ships
`datasets_csv/clinical_data/tcga_<cohort>_clinical.csv` for all five studies, and its columns are
`case_id, stage, grade, subtype`. **Stage was in the release, next to the grade everyone used.**
Swapping that one column, with everything else held fixed, raises the clinical baseline in all five
cohorts, by +0.042 to +0.221 C-index.

## The published clinical baselines, and why they should have prompted a second look

Verbatim from the primary tables (SurvPath Supplementary Table 3; MMP Table 6; DIMAF Table 1):

| cohort | age | sex | grade | age+sex+grade |
|---|---|---|---|---|
| BLCA | 0.578 | 0.489 | 0.515 | 0.570 |
| BRCA | 0.496 | 0.490 | 0.597 | 0.563 |
| COADREAD | **0.357** | 0.542 | not reported | 0.655 |
| HNSC | 0.517 | 0.486 | 0.547 | 0.512 |
| STAD | 0.499 | 0.529 | 0.552 | 0.592 |

A C-index of **0.357** for age in colorectal cancer is not a property of age; it is 0.643 with the
sign the other way. And DIMAF's own BLCA clinical baseline is 0.519 — a coin flip in a disease
whose staging drives every treatment decision made about it.

## The measurement: same pipeline, same cases, one column swapped

Ridge Cox, alphas (1, 8, 64, 512, 4096) by 3-fold inner CV inside each training fold, percentile
normalisation within fold, the released five-fold case-ID splits. Age and sex come from
`datasets_csv/metadata/`; stage and grade come from `datasets_csv/clinical_data/`; every clinical
row matched (359/359, 871/871, 298/298, 394/394, 319/319).

| cohort | n | events | age+sex | **+grade** | **+stage** | **stage − grade** |
|---|---|---|---|---|---|---|
| BLCA | 359 | 113 | 0.5636 | 0.5658 | **0.6218** | **+0.0560** |
| BRCA | 871 | 60 | 0.4828 | 0.4828 | **0.5686** | **+0.0858** |
| COADREAD | 298 | 37 | 0.4690 | 0.4690 | **0.6898** | **+0.2208** |
| HNSC | 394 | 117 | 0.4850 | 0.5027 | **0.5450** | **+0.0423** |
| STAD | 319 | 84 | 0.5150 | 0.5358 | **0.5819** | **+0.0461** |

**Five cohorts, five wins**, against a paired split-reseed bar of 0.0145 measured on this code.
The comparison is internal — our grade arm against our stage arm — which is what makes it a
statement about the covariate rather than about anyone's pipeline.

### The caveat that has to travel with it

**Our rebuilt age+sex+grade baseline does not reproduce the published one**, and the gap is not
small: −0.0042 (BLCA), −0.0802 (BRCA), −0.1860 (COADREAD), −0.0093 (HNSC), −0.0562 (STAD). For
COADREAD the reason is visible in the data — SurvPath's clinical file carries no usable grade for
that cohort, so our `age+sex+grade` arm is numerically identical to `age+sex` (0.4690 both), and
whatever produced the published 0.655 is not the age+sex+grade this file can build. The published
numbers are therefore **context, not the comparator.** The load-bearing comparison is the internal
one, and it is the only one this document rests on.

## What the method does with it

Same estimator, one view per modality, equal-weight rank average, no selection:

| cohort | WSI+omics | **+ stage** | same arm **with grade instead** | best published, any protocol |
|---|---|---|---|---|
| BLCA | 0.6841 | **0.6992** | 0.6861 | MOAD-FNet 0.691 · DIMAF 0.679 |
| BRCA | 0.7101 | 0.6892 | 0.6614 | APL 0.794 |
| COADREAD | 0.6819 | **0.7141** | 0.6231 | DSCASurv 0.832 |
| HNSC | 0.6169 | 0.6179 | 0.6015 | DSCASurv 0.666 |
| STAD | 0.6108 | 0.6311 | 0.6125 | DSCASurv 0.698 |

Stage beats grade inside the full method in all five as well. Against the strongest method whose
use of the released folds is **verified** (SurvPath: 0.625 / 0.655 / 0.673 / 0.600 / 0.592), the
arm is ahead in all five. Against the best number reported on **any** protocol it is ahead on BLCA
only — and APL, DSCASurv and MOAD-FNet all have unverified fold identity, which is exactly the
distinction that has to be made rather than averaged over.

### On bladder specifically, with the richer clinical block

`tnm.json` — T, N, M and AJCC stage fetched from GDC for BLCA in an earlier session — gives a
stronger clinical arm than SurvPath's single `stage` column, and it is what the headline bladder
number uses:

| arm | C-index |
|---|---|
| clinical, SurvPath's `stage` column (age+sex+stage) | 0.6218 |
| clinical, GDC T+N+M+stage (age+sex+T+N+M+stage) | **0.6856** |
| OURS with SurvPath's stage column | 0.6992 |
| **OURS with GDC T/N/M** | **0.7291** |

Both versions are public data; the first uses nothing outside the benchmark's own release, and the
second uses the TCGA clinical record any clinician would have. **Report both.**

## The bladder headline, with its paired uncertainty

Pre-specified arm, released folds, case-level paired bootstrap (6,000 replicates, cases resampled,
not comparable pairs):

| comparison | ours | baseline | Δ | SE | 95% CI | p |
|---|---|---|---|---|---|---|
| vs TITAN alone | 0.7291 | 0.6596 | +0.0696 | 0.0206 | [+0.0286, +0.1103] | **0.0010** |
| vs DIMAF 0.679 (published) | 0.7291 | 0.679 | +0.0501 | — | unpaired | — |
| vs cheap clinical Cox | 0.7291 | 0.6856 | +0.0435 | 0.0231 | [−0.0022, +0.0883] | 0.063 |
| vs SurvPath + same clinical (paired) | 0.7291 | 0.6963 | +0.0331 | 0.0200 | [−0.0061, +0.0721] | 0.097 |
| **without clinical, vs DIMAF** | 0.6841 | 0.679 | +0.0051 | — | **fails the bar** | — |

Fold mean ± SD 0.7285 ± 0.0644; 4 of 5 folds beat SurvPath+clinical.

**Three uncertainties, and they are not interchangeable.** The split-reseed SD is **0.0073** (2 SD
= 0.0145) — what moves if the same 359 patients are re-split, which is what a reader reproducing
the study varies. The selection-inflation term over the campaign's **211 scored candidates** is
**0.0239**. The case-level bootstrap SE is **≈0.020** — what moves under a *new cohort of 359*. The
first two are the right bars for a claim on a fixed benchmark and the margins clear both. The third
is the one a reviewer will ask for, and under it the two most important margins sit at p = 0.06 and
p = 0.10. **This cohort cannot resolve margins of 0.03–0.05 — and that applies to every published
number on it, not only to ours.**

## What was tried and failed

Seven candidate families, each with a pre-registered falsifier, a removal ablation and a matched
control read **before** the result. Full specifications in `development/contribution-design.md`.

| candidate | result | against a bar of 0.0145 |
|---|---|---|
| C1 pan-cancer survival subspace as a prior on the image head | −0.0042 | dead |
| C1' the same direction added as a zero-shot arm | +0.0033 | dead |
| C2 bagged ridge | −0.056 image / +0.024 omics | mechanism falsified |
| C3 elastic-net Cox | implementation failed to converge | not a result |
| C4 dense regularisation path | −0.004 to +0.017 | dead |
| C5 gated fusion on clinical risk tertile | 0.7493 vs **control 0.7502** | **VOID** |
| C6 inner-CV skill weighting | +0.0105 / +0.0037 / +0.0003 | dead |
| C7 stratified partial likelihood | −0.0199 on tied pairs | dead |

One positive result came out of the failures and is worth its own line: a Cox direction fitted on
**colorectal** slides, with zero case overlap, ranks **bladder** patients at 0.6492 against a
permuted-label control of 0.4727 ± 0.032. Cross-cancer survival directions transfer. They just do
not improve a bladder head that already finds a better one — the cosine between them is 0.1748.

## The one slice where the separation is largest

Restrict the comparable pairs to the decile the clinical score separates least — 2,423 of 24,219,
two patients whose stage cannot tell them apart:

| arm | all pairs | clinically tied |
|---|---|---|
| clinical alone | 0.6856 | 0.5167 |
| SurvPath + clinical | 0.6963 | 0.5852 |
| **OURS** | 0.7291 | **0.6513** |
| omics (Xena grouping) alone | 0.6749 | 0.6529 |

+0.0659 over SurvPath+clinical, SE 0.0363, 95% CI [−0.0033, +0.1367], p = 0.065 — double the
overall margin, on the question a pathologist is actually asking, and reported by nobody on this
benchmark. The tie threshold is recomputed inside every bootstrap replicate.

## Open

- **DIMAF's fold identity is contested** between two independent literature passes and must be
  settled with a diff before any claim is written against it. PIBD's was checked here directly
  (`cmp`, all five files byte-identical).
- **APL, DSCASurv and MOAD-FNet report higher numbers on some cohorts with unverified folds.** Their
  status has to be resolved, not averaged over.
- **PIBD's rerun crashed after fold 0** (`FileNotFoundError` on `splits_1.csv`, exit 1, 13:14; the
  file is present and readable from that working directory, so the cause is not configuration).
  Fold 0 reached best val C-index 0.6038 against a published five-fold mean of 0.667. Not restarted:
  a GPU launch needs fresh confirmation.
- **T/N/M for the four non-bladder cohorts** would make the richer clinical arm five-cohort. That is
  a GDC metadata fetch and needs confirmation under `data_download`.

---

# Addendum, 2026-08-17 evening — the finding got harder, and one of our own numbers was wrong

Two things were established after the confirmatory run, by diffing DIMAF's released fold files
against ours. Both belong in the paper.

## 1. DIMAF is scored on the released folds. Settled.

The contested status is resolved. `Trustworthy-AI-UU-NKI/DIMAF` ships
`src/data/data_files/tcga_blca/splits/{0..4}/{train,test}.csv` at **slide** granularity; collapsed
to `case_id` they are the released partition exactly:

| fold | released train/val | DIMAF train/test | identical? |
|---|---|---|---|
| 0 | 289 / 70 | 289 / 70 | **yes** |
| 1 | 287 / 72 | 287 / 72 | **yes** |
| 2 | 284 / 75 | 284 / 75 | **yes** |
| 3 | 288 / 71 | 288 / 71 | **yes** |
| 4 | 288 / 71 | 288 / 71 | **yes** |

359-case universe identical, zero cases in one but not the other, every case in exactly one test
fold. **DIMAF is a valid incumbent and the comparison against 0.679 stands.** The first of the two
literature passes was right; the second was wrong.

## 2. The stage column is in DIMAF's own fold files, beside the grade it used

The columns DIMAF ships in the files it reads to build its own folds:

```
case_id, slide_id, project_id, dss_survival_days, dss_censorship,
tissue_source_site, sex, ajcc_pathologic_tumor_stage, birth_days_to, histological_grade
```

`ajcc_pathologic_tumor_stage` is populated for 2,105 of 2,115 slide rows — Stage I 15, II 720,
III 680, IV 690. It sits directly beside `histological_grade`, which is the covariate DIMAF's
published clinical baseline uses, at 0.519.

So the earlier framing — *"stage is shipped in the benchmark's own release"* — understates it. **It
is in the incumbent's own split files.** There is no version of "the data was not available".

Rebuilding the entire clinical block from those files, so every clinical variable comes from the
incumbent's own repository:

| arm, all covariates from DIMAF's own split files | C-index |
|---|---|
| age + sex + **grade** — the covariate set DIMAF's published baseline uses | 0.5666 |
| age + sex + **stage** — the column beside it | **0.6638** |
| **difference** | **+0.0972** |
| full method with grade | 0.6856 |
| full method with stage | **0.7225** |
| **difference** | **+0.0369** |

## 3. Our own clinical block was defective, and fixing it moved the number UP

The same diff exposed a defect in our arm. The frozen run's clinical block took stage and T/N/M
from a cached GDC query that selected the **first** diagnosis record per case. **217 of 359 cases
carry more than one diagnosis**, and 61 had their T or stage change once the right one was chosen:

```
TCGA-SY-A9G0   Stage IIB / T2 N0 M0 / Adenocarcinoma, NOS           <- what we used
               Stage IV  / T4 N1 M0 / Transitional cell carcinoma   <- the bladder tumour
```

The tell was present and unread: `Stage IIA`, `IIB` and `IIC` appear in our record and **are not
AJCC bladder stages at all**. The errors were systematic — they understated stage (III → I,
IV → I) — because the first-listed diagnosis tends to be the earlier or the other tumour.

Recorded as **amendment A1** in `development/benchmark-protocol.json`, with the two causes
separated:

| clinical block | dim | clinical | full method |
|---|---|---|---|
| A frozen — GDC, first-diagnosis defect | 23 | 0.6863 | 0.7285 |
| B amended — DIMAF's own files, age+sex+stage | 6 | 0.6638 | **0.7225** |
| C corrected — GDC, right diagnosis, with T/N/M | 17 | 0.6992 | **0.7335** |
| D corrected block, grade instead of stage | 4 | 0.5666 | 0.6856 |

- **the defect alone** (A → C, same variables, right diagnosis): clinical **+0.0129**, full **+0.0050**
- **T/N/M alone** (B → C, same diagnosis, adds T/N/M): clinical **+0.0354**, full **+0.0110**
- frozen versus corrected, paired: **+0.0050, SE 0.0085, 95% CI [−0.0117, +0.0219], p = 0.56**

**Fixing the defect improved the number**, so the frozen 0.7260 was not inflated by it, and the
change to the headline is not significant.

**Arm B is primary going forward at 0.7212** (mean over seeds 0–4, SD 0.0048), not the higher C.
B has no provenance question at all — every column is in the file the incumbent itself reads —
while C uses a GDC query under a selection rule we designed. Both clear the pre-registered target
of 0.7029. C is reported beside it as the with-full-TNM arm.
