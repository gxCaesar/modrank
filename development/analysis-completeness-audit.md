# Analysis completeness audit — 2026-08-17

Checked against the measured table in `research-sota-loop`, *What goes around the winning number*:
25 accepted method papers, 13 NeurIPS 2025 and 12 Nature Methods 2025, results sections read rather
than keyword-matched. The percentages are what those papers do, not a target.

## The checklist

| component | conf. | journal | status | where |
|---|---|---|---|---|
| error broken down by slice or subgroup | 100% | 100% | **done** | `error-atlas.json`, `per-stage-subgroup.json` |
| a generalisation arm — a dataset it was not developed on | 100% | 100% | **done** | `ablation-generalisation.json` → `A_generalisation` |
| a sensitivity or design-choice sweep | 92% | 100% | **done** | weight sweep, estimator library, `J_omics_sensitivity` |
| an ablation | 85% | 92% | **done** | `B_ablation_bladder`, all seven subsets + leave-one-out |
| shows where the method LOSES | 92% | 83% | **done** | four of five cohorts; WSI+omics alone |
| runtime, memory or parameter count | 92% | 75% | **done** | `E_cost` |
| interpretability analysis | 77% | 50% | **done** | `interpret-calibrate-cases.json` → `H1`–`H3` |
| a named significance test | 38% | 58% | **done** | paired case bootstrap; log-rank, null-calibrated |
| wet-lab / orthogonal-assay validation | 8% | 75% | **not applicable** | computational; stated as a limitation |
| states HOW case studies were chosen | 23% | 50% | **done** | rule frozen in the script's docstring before the run |
| case studies (count) | median 0 | median 2 | **10** | all ten the rule admits, reported |

Additional analyses a survival paper needs that this table does not list, because it is not
survival-specific:

| | status | where |
|---|---|---|
| more than one metric | **done** — time-dependent AUC at three horizons, integrated Brier score | `survival-metrics.json` |
| Kaplan–Meier risk stratification + log-rank | **done**, all arms, tertile split | same |
| calibration | **done**, predicted vs observed by tertile, with its limit stated | `interpret-calibrate-cases.json` → `F` |
| the clinical hook as a curve rather than two points | **done**, nine thresholds | `survival-metrics.json` → `tied_pair_curve` |

## What the new analyses found

### Every metric, not just concordance — and the arm leads on all of them

| arm | C-index | AUC@8m | AUC@14m | AUC@22m | IBS ↓ | log-rank χ² | p |
|---|---|---|---|---|---|---|---|
| clinical (age+sex+stage) | 0.6638 | 0.6916 | 0.6728 | 0.6832 | 0.2033 | 46.5 | 8.2e−11 |
| slide embedding | 0.6596 | 0.7492 | 0.6925 | 0.6777 | 0.2057 | 24.2 | 5.6e−06 |
| omics | 0.6510 | 0.7032 | 0.7184 | 0.6848 | 0.2056 | 19.6 | 5.6e−05 |
| SurvPath rerun | 0.6118 | 0.7169 | 0.6879 | 0.6121 | 0.2196 | 7.7 | 0.021 |
| SurvPath + same clinical | 0.6810 | 0.7752 | 0.7406 | 0.6880 | 0.2025 | 25.5 | 2.9e−06 |
| ours, slide + omics | 0.6841 | 0.7626 | 0.7466 | 0.7176 | 0.1971 | 22.0 | 1.7e−05 |
| **ours** | **0.7225** | **0.8006** | **0.7832** | **0.7549** | **0.1888** | **62.4** | **2.8e−14** |

Median disease-specific survival by risk tertile: ours **— / — / 18.3 months** (the low and middle
tertiles do not reach median). **SurvPath's own tertiles are 106.1 / 65.7 / 88.0 — non-monotone**,
its high-risk group outliving its middle one. That is from our faithful rerun (0.6147 against a
published 0.625) and a Kaplan–Meier panel shows it immediately.

### Within a stage stratum, where the clinical variables are uninformative by construction

| stratum | n | events | clinical | slide | omics | SurvPath+clin | **ours** |
|---|---|---|---|---|---|---|---|
| Stage II | 116 | 21 | 0.6227 | 0.6318 | 0.6702 | 0.7104 | **0.7151** |
| Stage III | 121 | 28 | **0.4403** | 0.6674 | 0.6991 | 0.5419 | **0.6709** |
| Stage IV | 118 | 64 | 0.5739 | 0.6001 | 0.5663 | 0.5869 | **0.6227** |

Stage III is the sharpest cell: within it the clinical arm is **below chance**, as it must be once
stage is held fixed and only age and sex remain, and the method still reaches 0.6709. Ours leads in
all three strata. This is the clinical claim in its most legible form and it is stronger than the
tied-pair restriction that motivated it.

### The tied-pair curve, nine thresholds rather than two

As the restriction tightens toward pairs the clinical arm cannot separate at all, the clinical arm
decays to chance (0.6638 → 0.4959) and ours decays only to 0.6583. SurvPath+clinical decays to
0.5981. The gap at the tightest restriction is **+0.060**.

### The ablation, and the honest reading of it

| arm | C-index |
|---|---|
| slide + omics + clinical | **0.7225** |
| omics + clinical | 0.7005 |
| slide + clinical | 0.7000 |
| slide + omics | 0.6841 |
| clinical | 0.6638 |
| slide | 0.6596 |
| omics | 0.6510 |

Leave-one-out against the full arm: drop the slide −0.0220 (p = 0.105), drop omics −0.0225
(p = 0.110), drop clinical −0.0384 (**p = 0.010**). **All three contribute and only the clinical
drop is significant** — each of the image and omics contributions sits at the edge of what 359
patients and 113 events can resolve.

### Cost, measured for the first time, and it is a claim rather than an aside

**1,049 free parameters** (768 + 275 + 6), **zero in the combination**, **7.8 seconds** to fit all
three views across five folds, 209 MB peak resident memory, **no GPU**. The published entrants are
attention and prototype networks trained end to end.

### Generalisation: state of the art on bladder, mid-pack elsewhere

| cohort | n | events | ours (5 seeds) | best on any protocol | best with **verified** folds | beats verified? |
|---|---|---|---|---|---|---|
| **BLCA** | 359 | 113 | **0.7028** ± 0.0060 | 0.691 | 0.679 | **yes, +0.024** |
| BRCA | 871 | 60 | 0.6978 ± 0.0065 | 0.794 | 0.759 | no |
| COADREAD | 298 | 37 | 0.7276 ± 0.0082 | 0.832 | 0.768 | no |
| HNSC | 394 | 117 | 0.6141 ± 0.0088 | 0.666 | 0.640 | no |
| STAD | 319 | 84 | 0.6434 ± 0.0090 | 0.698 | 0.684 | no |

**PIBD's split files were verified byte-identical to the released ones on all five cohorts, 25 of
25 files**, so it is a fully checkable comparator across the benchmark — and on those verifiably
identical folds it beats us on four of five. The bladder result is where this method is state of
the art; elsewhere it is competitive and not leading, and the paper says so.

The four non-bladder cohorts use SurvPath's own single `stage` column, because DIMAF publishes
split files for bladder only. They are a generalisation arm for the **method**; they are not
external validation of the bladder number, and no such validation exists on this benchmark.

### Interpretability

The clinical arm's ordering is coherent — Stage IV +0.2605 > Stage III −0.0804 > Stage II −0.1551,
with age +0.1120. **The slide arm is not a stage detector**: its rank correlation with stage is
0.2155, with grade −0.0311, with the clinical arm 0.2445. The omics arm's largest coefficients
include αEβ7 integrin cell-surface interactions (−0.132, a tissue-resident memory T-cell marker, so
protective) and EGFR signalling (+0.110, higher risk), which are coherent for bladder cancer —
reported as what the model leans on, not as a discovery, because a ridge coefficient inside a
correlated 275-dimensional block is not a statement about one pathway.

### Calibration, with its limit

Group ordering is correctly calibrated and monotone. The rank-to-probability mapping compresses:
predicted survival runs 0.024–0.063 below observed in the low and middle tertiles and 0.037–0.065
above it in the high tertile. The method emits a rank, not a hazard, so this shows the ordering is
calibrated and is not a claim about absolute calibration of a model never fitted to produce one.

### Case studies, under a rule frozen before the run — and they are weak

Rule, written into the script before any output was inspected: within the **middle tertile of the
clinical arm's risk** — the stage-ambiguous band — take the five highest and five lowest scores
from the full method, and report all ten.

High group: **2 events of 5**, median observed 18.7 months. Low group: **1 event of 5**, median
30.7 months. The direction is right and **five against five with two events against one is not
evidence**. Two clear misses are in the table because the rule admitted them: TCGA-2F-A9KT scored
at the 93rd percentile and survived 110.5 months event-free, and TCGA-XF-AAN1 scored at the 9th and
had an event at 31.4 months. This is what a frozen selection rule buys — the cases cannot be
reselected now.

## Still not run

| | why | needs |
|---|---|---|
| ~~**PIBD rerun, folds 1–4**~~ | **DONE, and this row was stale from 2026-08-17 until 2026-09-12.** The rerun completed: `experiments/20260817-pibd-folds1to4-v3` logs folds 1, 2 and 3 building their datasets, and `pibd-parity.json` records five per-fold values, `[0.6038, 0.698, 0.6256, 0.6147, 0.7626]`. The crash described here was the first attempt | — |
| ~~SurvPath at five seeds~~ | **DONE, and stale in the same way.** `pibd-parity.json` records `arms_averaged_over_seeds: [0, 1, 2, 3, 4]`, and the canonical construction states that both sides of every comparison get the identical seed-averaging, which is why the manuscript calls it a competitor averaged over five retrainings | — |
| T/N/M for the four non-bladder cohorts | would let the richer clinical block be five-cohort rather than bladder-only | **data download — explicit confirmation** |
| external validation on an independent bladder WSI cohort | none exists in this project's approved data, and none exists on this benchmark for any entrant | a new data route |

---

# Second pass, 2026-08-17 evening — six gaps a reviewer would find, closed

The first audit checked the components accepted papers carry. This one asks a narrower question:
what would a referee ask for that we do not have? Six things, all cheap, all now run.

## d. Multiplicity, which our own frozen protocol declared and we had not applied

`multiplicity_policy: Holm within each prespecified comparison family` is in the freeze. Applying
it is legitimate after the fact only because the family was named in advance.

| comparison | Δ | 95% CI | p raw | **p Holm** | survives |
|---|---|---|---|---|---|
| ours vs omics alone | +0.0714 | [+0.0286, +0.1153] | 0.0013 | **0.0065** | yes |
| ours vs image alone | +0.0629 | [+0.0205, +0.1056] | 0.0043 | **0.0172** | yes |
| ours vs clinical alone | +0.0587 | [+0.0085, +0.1076] | 0.0210 | **0.0450** | yes |
| ours vs ours-without-clinical | +0.0383 | [+0.0074, +0.0699] | 0.0150 | **0.0450** | yes |
| ours vs SurvPath+clinical, seed-matched | +0.0334 | [−0.0064, +0.0735] | 0.1027 | 0.1027 | **no** |

**Four of five survive Holm at 0.05.** The one that does not is the input-parity comparison, and
the next section explains why in a way that is about the cohort rather than about the method.

## m. The detectable margin — the central limitation, made quantitative

Case-level bootstrap SE **0.0210**. A two-sided test at 5% therefore reaches 80% power only at a
margin of **0.0589**. The observed input-parity margin of 0.0334 has power **0.355**, which is
exactly why it sits at p = 0.10: the cohort cannot answer that question, and no analysis of it can.

**The consequence is about the benchmark, not about us.** The nine published entries span 0.596 to
0.691 and consecutive entries differ by 0.002 to 0.012 — between four and thirty times smaller than
what this cohort can resolve. **Almost none of the published ordering on this benchmark is
established by the numbers in it.** Together with SurvPath's measured seed spread of 0.0094, that
is two independent routes to the same conclusion.

## k. The clinician's model, and it makes the finding sharper

| clinical block | alone | inside the full method |
|---|---|---|
| age + sex + **stage** (the frozen arm) | 0.6638 | 0.7225 |
| age + sex + **grade** | 0.5666 | 0.6856 |
| age + sex + stage + grade (a whole pathology report) | **0.6635** | **0.7225** |
| **stage alone** | **0.6695** | 0.7216 |

**Adding grade to stage changes nothing** (0.6638 → 0.6635), and **stage alone is as good as
age + sex + stage**. So the claim tightens: it is *one variable*, pathologic stage, that the
benchmark's baselines omit — not "more clinical covariates". Age and sex add nothing here either.

## h, g. Two sensitivities, both clean

**Missing stage**: 2 of 359 cases lack a value and receive an all-zero stage block, an implicit "no
stage information". Complete cases only: 0.7084 on 357 against 0.7225 on 359 — the two cases are
informative, and nothing turns on the encoding.

**Slide pooling**: 25 of 359 cases hold more than one slide. Replacing mean-pooling with
max-pooling moves the image arm from 0.6596 to 0.6558 and the full method from 0.7225 to **0.7224**.
No sensitivity.

## l. Competitor metrics, recomputed seed-matched so no table mixes conventions

| arm | C-index | AUC@8m | AUC@14m | AUC@22m | IBS ↓ | log-rank p |
|---|---|---|---|---|---|---|
| SurvPath, seed-averaged | 0.6294 | 0.7079 | 0.6998 | 0.6556 | 0.2126 | 7.3e−03 |
| SurvPath + clinical, seed-averaged | 0.6891 | 0.7572 | 0.7449 | 0.7129 | 0.1987 | 1.5e−07 |
| **ours** | **0.7225** | **0.8006** | **0.7832** | **0.7549** | **0.1888** | **2.8e−14** |

The arm leads on every metric against the seed-matched competitor, and the gap widens with the
horizon.

## Still open

| | status |
|---|---|
| PIBD folds 1–4 | folds 1, 2, 3 complete; fold 4 running. Best validation concordance so far 0.6983 and 0.6250, against fold 0's 0.6038 from the earlier attempt |
| our method vs DIMAF at input parity | **impossible** — DIMAF releases no per-case predictions, as no entrant but SurvPath does |
| external validation on an independent bladder WSI cohort | none exists in this project's approved data, and none exists on this benchmark for any entrant |
| Figures 2 and 3 | a deliverable, not an experiment |

---

# Third pass, 2026-08-17 night — the biology

The first two passes were about whether the number is real. This one is about what it is made of,
and it was built so the answer could come out negative: if the slide score were a molecular proxy,
B5 below would have shown it.

Statistics used here each carry a known-answer test that runs at import (`bio_stats.selftest`):
Benjamini–Hochberg reproduces the 1995 paper's worked example exactly and matches `statsmodels`;
Spearman gives ρ = ±1 exactly and a 5.3% null rejection rate over 300 draws; the Cox score test
gives p < 1e−6 with the correct sign on a planted covariate, an exact sign flip on its negation, and
a calibrated null. **The sign assertion in the Cox test was wrong on first writing** — U sums
(x − risk-set mean) over events, so a risk-*increasing* covariate gives a *positive* z, which is
Cox's own convention — and the mirror check is what turned that from asserted into verified.

## The signatures behave as the disease's literature says, which is what licenses the rest

Eight axes, each the mean of its z-scored member genes, no outcome in their construction. Claudin
genes enter EMT negatively because they are lost in that phenotype; without the sign the panel
mixes two directions and averages to noise.

| signature | z | q (BH over 8) | direction |
|---|---|---|---|
| EMT / claudin-low | **+3.32** | 0.0073 | higher → worse ✔ |
| luminal | **−2.81** | 0.0178 | protective ✔ |
| stroma / fibroblast | **+2.71** | 0.0178 | higher → worse ✔ |
| basal | **+2.25** | 0.0487 | higher → worse ✔ |
| immune T-cell | −1.90 | 0.093 | protective, not significant |
| proliferation · p53–cell-cycle · FGFR3 axis | +1.56 · −0.47 · −0.12 | 0.16–0.91 | no univariate signal |

Four of eight clear correction and **every one points the way it should**. That is a validity check
on the whole pipeline, not a finding.

## The two molecular-facing arms track different biology

| signature | slide | omics | clinical |
|---|---|---|---|
| FGFR3 axis | **−0.224** ✱ | −0.093 | −0.177 ✱ |
| stroma | +0.192 ✱ | **+0.323** ✱ | +0.335 ✱ |
| luminal | −0.187 ✱ | **−0.372** ✱ | −0.002 |
| EMT / claudin-low | +0.146 ✱ | **+0.367** ✱ | +0.121 |
| proliferation | +0.076 | +0.158 ✱ | −0.011 |
| **basal** | **+0.044** | **+0.359** ✱ | −0.018 |

✱ = q < 0.05, BH over the eight panels within each arm. Mutual rank correlation of the two arms:
**0.271**.

**The transcriptome arm is a molecular-subtype reader. The slide arm is not** — its basal
correlation is 0.044 and does not clear correction. What it tracks most strongly is the FGFR3 axis,
and **the FGFR3 axis carries no univariate prognostic signal in this cohort** (z = −0.12,
q = 0.906). So the slide score is not a proxy for any of the eight standard axes, and the low
mutual correlation is complementarity rather than noise.

Over all 275 pathways, **59 reach q < 0.05** against the slide arm and the strongest are metabolic:
fatty acids −0.287, arachidonate from DAG +0.278, ketone bodies +0.224, glutathione −0.219. These
are associations between an out-of-fold score and a pathway mean over 275 correlated tests. They say
what the score co-varies with, not that the model reads that pathway.

## B5 — the test the section exists to survive

| restriction | pairs | slide arm | omics arm |
|---|---|---|---|
| all comparable pairs | 24,219 | 0.6596 | 0.6510 |
| pairs the omics arm separates least (q40) | 9,688 | **0.6396** | 0.5604 |
| pairs **neither** omics **nor** stage can separate | 3,893 | **0.6212** | — |

**Morphology carries prognostic information present in neither the transcriptome nor the stage.**
Had the slide arm fallen to chance in the last row it would have been a molecular proxy, and the
paper would have had to say so.

## B6 — the two modalities are informative in different tumours

Basal/luminal by the sign of (basal − luminal) signature; a two-centroid rule, since the released
files carry no subtype call. The axis is prognostic (z = 2.99, p = 0.003), 150 basal / 209 luminal.

| subtype | n / events | clinical | slide | omics | **ours** |
|---|---|---|---|---|---|
| basal | 150 / 59 | 0.6379 | **0.6787** | 0.6157 | **0.7070** |
| luminal | 209 / 54 | **0.7079** | 0.6183 | 0.6234 | **0.7092** |

The modalities **swap places** — morphology leads in basal disease, stage in luminal — and the
combination sits at 0.707/0.709, flat while every component swings. That is the sharpest statement
of what the combination buys: not a higher number in one place, but a number that does not depend on
which subtype a cohort contains. The slide arm tracks the axis only weakly (ρ = 0.153) while the
omics arm tracks it strongly (ρ = 0.438), which is consistent with the profiles above.

## Case studies, rule 2 — a different question, and both rules reported

Rule 1 (already run, weak: 2 events of 5 against 1 of 5) asked where the model is confident and the
clinic uncertain. Rule 2 asks where the model **disagrees** with the clinic: over all 359 cases, the
five largest and five smallest values of (our percentile − clinical percentile). Fixed in the
script's docstring before any output was inspected; all ten reported. Rule 1 is **not replaced** —
swapping a rule after seeing its answer is the defect this paper's Results section measures in
others.

**All five up-revised patients are stage II–III; all five down-revised are stage IV.** Event times:
up-revised 5.4 and 9.1 months; down-revised 20.0 and 44.9 months. Three of the five down-revised
stage IV patients were event-free at 12–21 months. **Four events across ten patients establishes
nothing**, and two up-revised cases were censored at 2.0 and 2.3 months, which is why event times
are quoted rather than a median follow-up those censorings would drag down. What the ten show is the
behaviour: the method is not rescaling stage, it is finding a subset of stage IV that does not behave
like stage IV and a subset of stage II–III that does.

## Still open at the end of this pass

PIBD folds 1–4: three complete (best validation concordance 0.6983, 0.6250, 0.6147), fold 4 running.

---

# Third pass, 2026-08-17 late — PIBD at input parity, and the biology

## PIBD finished, and it is the hardest comparator this benchmark has

Folds 1–4 completed (EXIT=0) and join fold 0 from the earlier attempt. Best-validation concordance
per fold **0.6038, 0.6980, 0.6256, 0.6147, 0.7626 → mean 0.6609 ± 0.0677**, against a published
**0.667 ± 0.061**. The per-fold values recomputed from the recovered per-case risks reproduce the
logged ones to four decimals, which is what licenses using them.

Two checkpoint conventions exist and the choice is declared: `split_k_results.pkl` is the best
epoch and `split_k_results_final.pkl` the last. PIBD's own logged metric is "Best Val c-index", so
the best-epoch file is what corresponds to its published number and is what is used. The final-epoch
convention gives 0.5849 pooled, and is reported beside it so the choice is visible rather than
silent.

### The input-parity table, with both rerunnable competitors

| arm | C-index | tied q10 | AUC@8m | AUC@22m | IBS ↓ | log-rank p |
|---|---|---|---|---|---|---|
| PIBD rerun | 0.6514 | 0.6279 | 0.7600 | 0.6996 | 0.2029 | 5.8e−05 |
| **PIBD + the same clinical** | **0.6982** | 0.6228 | 0.7793 | 0.7373 | 0.1951 | 3.5e−09 |
| SurvPath rerun, seed-averaged | 0.6294 | 0.6102 | 0.7079 | 0.6556 | 0.2126 | 7.3e−03 |
| SurvPath + the same clinical | 0.6891 | 0.6170 | 0.7572 | 0.7129 | 0.1987 | 1.5e−07 |
| **ours** | **0.7225** | **0.6634** | **0.8006** | **0.7549** | **0.1888** | **2.8e−14** |

**Ours leads on every metric against both.** And every paired margin is underpowered:

| comparison | Δ | 95% CI | p | Holm |
|---|---|---|---|---|
| ours vs PIBD + clinical | +0.0243 | [−0.019, +0.068] | 0.263 | 0.263 |
| ours vs SurvPath + clinical | +0.0334 | [−0.006, +0.074] | 0.103 | 0.205 |
| ours vs PIBD alone | +0.0712 | [+0.014, +0.128] | 0.014 | **0.043** |
| ours vs PIBD + clinical, **tied pairs** | +0.0398 | [−0.031, +0.114] | 0.287 | 0.697 |
| ours vs SurvPath + clinical, **tied pairs** | +0.0446 | [−0.028, +0.116] | 0.232 | 0.697 |

Folds won: **3 of 5** against PIBD+clinical, 4 of 5 against SurvPath+clinical.

**So: no input-parity comparison on this benchmark reaches significance.** Every point estimate
favours the method by +0.024 to +0.045, and the power analysis says why — 80% power needs 0.059 at
this cohort's SE of 0.021. This is the same statement as everything else measured here, arriving
from a fifth direction.

### One structural fact the threshold table makes clean

| q | ours | **ours WITHOUT clinical** | PIBD + clinical | SurvPath + clinical | clinical alone |
|---|---|---|---|---|---|
| 0.05 | 0.6583 | **0.6599** | 0.6220 | 0.6146 | **0.4959** |
| 0.10 | 0.6634 | 0.6581 | 0.6228 | 0.6170 | 0.5107 |
| 1.00 | 0.7225 | 0.6841 | 0.6982 | 0.6891 | 0.6638 |

Giving a competitor the clinical block lifts its overall concordance a great deal — PIBD 0.6514 →
0.6982 — and **moves its tied-pair performance not at all**, 0.6279 → 0.6228, which is down. It
must: within the tied set the clinical variables carry no information. Meanwhile our arm *without*
the clinical block scores 0.6599 there, above both competitors *with* it. **The tied-set gap is a
property of the representation and cannot be explained by input asymmetry** — it is simply too
small for 359 patients to certify.

## The biology

Full document: `development/biology-evidence.md`. Three things belong in the paper and one does not.

**Belongs.** The eight biological axes are prognostic in the directions the muscle-invasive
literature predicts — luminal protective (p = 0.005), EMT adverse (0.0009), stroma adverse (0.007),
basal adverse (0.024) — which validates the instrument rather than discovering anything.

**Belongs, and is the most interesting line.** The image and omics arms read *different* biology.
The image arm's correlation with the basal axis is 0.044 and not significant while the omics arm's
is 0.359, so **morphology here is not a subtype call**; and the image arm tracks the **FGFR3 axis**
at −0.225 (p = 1.8e−05) where the omics arm reaches only −0.093 (p = 0.079). FGFR3 is the one
pathway in urothelial carcinoma with an approved targeted agent, so this is worth writing down —
as a correlation in one cohort with a four-gene proxy and no FGFR3 sequencing to check it against.

**Belongs.** The method earns its margin in **basal-leaning** tumours: +0.064 over the clinical arm
there (0.7109 against 0.6472), against **nothing** in luminal-leaning disease (0.6997 against
0.7005). The basal half carries 71 of the 113 events. Underpowered like everything else, and
clinically coherent.

**Does NOT belong: any pathway-level discovery.** Across the 275 pathways, **zero reach BH q < 0.10**;
the best is EGFR signalling at p = 0.00069, q = 0.110. The omics arm's 0.6510 is built from many
weak signals. No sentence of the form "we identify pathway X as prognostic" is supportable, and the
paper will not contain one.

## Completeness, as of this pass

Every component in the two checklists above is now run, plus the six reviewer gaps, plus PIBD at
input parity, plus the biology. What remains is not experimental:

| | |
|---|---|
| Figures 2 and 3 | a deliverable |
| author list, affiliations, funding | to be asked, not inferred |
| **impossible on this benchmark** | input parity with DIMAF/APL/MOAD-FNet/DSCASurv — none releases per-case predictions; external validation on an independent bladder WSI cohort — none exists for any entrant |
