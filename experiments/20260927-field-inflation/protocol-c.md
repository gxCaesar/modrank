# 20260927-field-inflation, analysis C: the reference gap in external bladder cohorts

Written after the systematic search and before any download or outcome analysis. Post-freeze and
exploratory, like A, B and D. Nothing here changes ModRank or any committed value.

## Why C changed

C was specified to replicate the inflation D. Analyses A, B and D then showed that D is positive for a
score without information under the rank combination (README, finding 1), so a replication of D
against zero would carry little information. C therefore tests the part of the paper's first finding
that no combination rule produces: whether a clinical reference built from pathologic stage separates
patients better than one built from grade, and whether that gap tracks how little grade varies, as the
paper's mechanism says. This needs grade, stage and survival for the same patients, not expression.

## Cohorts

From the search of 2026-09-27 (GEO, ArrayExpress/BioStudies, cBioPortal, source papers; every field
marked verified or unverified in the search record). A cohort is admitted when, after download and
before any outcome is analysed, a row-level count shows at least 50 patients with grade (two or more
levels), stage (two or more levels), a time and an event indicator, and at least 20 events.

| cohort | setting | endpoint in the record | reported size | status |
|---|---|---|---|---|
| GSE5479 | NMIBC | progression or bladder-cancer death | 300 with complete covariates, 83 events | candidate |
| GSE1827 | mixed, 50 cystectomies | overall survival | 80, 44 deaths | candidate |
| GSE19915 (GPL5186 subset) | mixed | overall survival | 87 with OS, 26 deaths | candidate |
| E-MTAB-1803 | MIBC | overall survival | 73, 43 deaths; grade G2 4, G3 69 | candidate |
| GSE13507 | mixed | overall and cancer-specific survival | 165 primary tumours | PI decision: a retained bundle under a project closure |
| E-MTAB-4321 (UROMOL 2016) | NMIBC | progression-free survival | 454 evaluable | PI decision: a retained bundle under a project closure |

Excluded before download: GSE32548 (overlaps GSE32894, already analysed), GSE83586 and GSE57933
(grade not established per patient), GSE120736 and GSE128959 (repeat sampling, survival not linked),
E-MTAB-1940 (no survival time found), E-TABM-147 (fewer than 50 with survival), UROMOL 2021 (overlaps
E-MTAB-4321), treatment cohorts (GSE48276, GSE69795, GSE87304, IMvigor210, UC-GENOME, which stays
unopened), and the MSK/TCGA and PDX studies (overlap with TCGA, or fewer than 50). A lookup that failed
is not evidence that a field is absent, and a cohort excluded on an unverified field is re-admitted if
its source file turns out to carry it.

## Analysis (analysis/s36_reference_gap_external.py, written after the downloads are approved)

Per cohort, on the patients the row-level count admits:

1. Two clinical models by the paper's recipe: ridge Cox with the penalty chosen by inner three-fold
   cross-validation, five repeats of stratified five-fold cross-validation (the GEO protocol's
   construction), out-of-fold scores averaged over repeats. The grade model has age, sex and one-hot
   grade; the stage model has age, sex and one-hot stage (T category or AJCC group, whichever the
   record carries, stated per cohort). Where age or sex is missing for the whole cohort both models
   drop it.
2. The gap C(stage model) − C(grade model), Harrell's concordance, with a paired case bootstrap of
   6,000 resamples, seed 20260911.
3. The normalised Shannon entropy of grade and of stage, from the admitted patients.

Across cohorts: the gaps of the new cohorts together with the five TCGA studies (committed, released
folds) and the two GEO cohorts (committed), plotted against grade's normalised entropy, with a
Spearman correlation reported as a description over at most eleven points, not as a test.

Endpoints differ between cohorts (disease-specific, overall, progression). Each gap is reported with
its own endpoint and no pooled estimate is formed across endpoints.

## Predictions, written before any download

1. The gap is positive in the cohorts where grade barely varies (E-MTAB-1803; any cystectomy cohort
   whose grade is above 85% one level).
2. The gap is smaller, and may be negative, in the non-muscle-invasive cohort (GSE5479), where grade
   separates low from high risk and stage takes only Ta and T1.
3. Across all cohorts the gap is larger where grade's normalised entropy is lower (Spearman below
   zero).

If prediction 1 fails, the paper's mechanism sentence is narrowed to the TCGA studies. If 2 or 3 fails,
the paper says so. Every admitted cohort is reported.

## Downloads, for approval (none made)

Clinical tables only, to sysu `/data1/guanxing/bladder_cancer/data/geo-c/`. Sizes are from indexed
listings and are unverified until the transfer.

| cohort | file | approximate size |
|---|---|---|
| GSE5479 | `GSE5479_clinical_information.txt` (GEO supplementary) | under 1 MB |
| GSE1827 | series matrix (characteristics, carries the expression table too) | about 12 MB |
| GSE19915 | GPL5186 series matrix | a few MB |
| E-MTAB-1803 | SDRF and the clinical sample table (BioStudies) | under 1 MB |

GSE13507 and E-MTAB-4321 are not on this list; they need the PI's decision first.

## Amendment, 2026-09-27, before any download

The PI approved the four downloads and enabled GSE13507 and E-MTAB-4321 for C as a new use (decision
H4-2026-09-27-c-downloads-and-retained-bundles). Both join the candidate list under the same admission
rule and analysis. GSE13507 is mixed-stage with cancer-specific survival; E-MTAB-4321 is
non-muscle-invasive with progression-free survival, so prediction 2 applies to it as to GSE5479. Their
clinical fields are read from the copies already on sysu where those carry them, and only the
clinical table is used. The predictions above are unchanged.

## Gate 0, 2026-09-27, before any outcome was analysed

Row-level counts (`analysis/s36_reference_gap_external.py gate0`, output
`analysis-results/reference-gap-gate0.json`): GSE19915 GPL5186 85 complete patients, 25
disease-specific deaths, no age or sex field (its models omit both); E-MTAB-1803 146, 86 deaths
(overall survival; one row per individual after collapsing the SDRF's repeated assays); GSE13507 165
primary tumours, 32 cancer-specific deaths; E-MTAB-4321 462, 31 progressions to T2+. All four meet
the admission rule. GSE5479 and GSE1827 carry no time or event field in their deposited files and are
not admitted (a source outside those files was not sought, as it would need another download).

One coding rule fixed here: GSE19915 records four tumours as "G3^Nested", grade 3 with the nested
variant appended, and these are coded G3. Grade is otherwise used as recorded (E-MTAB-4321 keeps
PUNLMP as its own level).
