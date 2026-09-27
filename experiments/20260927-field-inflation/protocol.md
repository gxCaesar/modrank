# 20260927-field-inflation: protocol, written and committed before any of these runs

Corrected before any run (second commit): MCAT's upstream script is `mcat_pathways.sh` (lr 0.001, two
epochs), not the MOTCat script whose values the first commit gave it.

Three analyses requested on 2026-09-27 to test how far the paper's two findings extend. All three are
specified after the 2026-08-17 freeze of the bladder protocol and will be reported as such. Nothing
here changes ModRank, its primary value (0.7212) or any existing comparison.

## A. The inflation in five further multimodal architectures (TCGA-BLCA)

**Question.** The paper measures the inflation D for three constructions (ModRank, and SurvPath and
PIBD given the clinical block). Is it a property of how the field's models are compared, rather than of
three models?

**Models.** The multimodal baselines implemented in the SurvPath code base (the same tarball fetched on
2026-08-15 that the SurvPath rerun used): MCAT (`--modality coattn`), MOTCat (`coattn_motcat`),
ABMIL with pathways (`abmil_wsi_pathways`, `--fusion concat`), TransMIL with pathways
(`transmil_wsi_pathways`, `--fusion concat`) and DeepMISL with pathways (`deepmisl_wsi_pathways`,
`--fusion concat`).

**Inputs and training, fixed.** The released five folds; `survival_months_dss`; the same 768-d CHIEF
(CTransPath) tile features the SurvPath rerun read (`--encoding_dim 768`); the `combine` pathways.
Hyperparameters are the upstream scripts' own, not tuned here: MOTCat `--lr 0.0005 --max_epochs 5`
(`scripts/mcat_motcat.sh`); MCAT, ABMIL, TransMIL and DeepMISL `--lr 0.001 --max_epochs 2`
(`scripts/mcat_pathways.sh`, `abmil.sh`, `transmil.sh`, `deepmisl.sh`), MCAT and MOTCat with
`--fusion concat` as upstream, ABMIL with
`--encoding_layer_1_dim 8 --encoding_layer_2_dim 16 --encoder_dropout 0.25`; for all, `--batch_size 1
--opt radam --reg 0.0001 --alpha_surv 0.5 --weighted_sample --bag_loss nll_surv --n_classes 4
--num_patches 4096 --wsi_projection_dim 256`. Five retrainings per model, seeds 1 to 5, as for
SurvPath. The per-case risk each run writes (`split_k_results.pkl`) is what is analysed. An
engineering failure (memory, a missing import) may be fixed and is recorded; a model that cannot be
run is reported as not run and is not replaced.

**Analysis (`analysis/s33_field_inflation.py`).** For each model, its five runs are percentile-ranked
within fold and averaged, and the average is ranked again (M). M is combined with each clinical
reference by the paper's rank rule: M + grade = fold percentile of (M + the canonical grade arm), and
M + stage likewise with the canonical stage arm, both arms read from the committed per-case vectors
(`experiments/20260911-blca-posthoc/results/percase-canonical-vectors.jsonl`). Added value over each
reference is C(M + R) − C(R), and D is added value over grade minus added value over stage. All four
concordances are recomputed on each of 6,000 case resamples drawn with seed 20260911, the draws the
paper's existing D used. Before reporting anything the script must reproduce SurvPath's committed
values (added value +0.0670 over grade, +0.0253 over stage, D +0.0417) from the committed SurvPath
vectors by the same code path. Holm correction is applied across the five new D tests.

**Predictions, written before the runs.** (1) D is positive in all five. (2) D's interval excludes
zero in at least three of the five. (3) The added value over the stage-based reference has an
interval covering zero in at least three of the five. Each model's own concordance is reported beside
its published value on these folds where one exists.

**What each outcome means.** If (1) and (2) hold, the inflation holds for every multimodal architecture
that could be rerun on the benchmark, eight in all, and the paper says so. If D is positive but mostly
unresolved, the paper reports the direction and the count. If D is negative for any model, the paper
reports it and the claim is scoped to the architectures where it holds. Every result is reported.

## B. The inflation in the other TCGA studies with a usable grade

**Question.** Does D shrink where grade carries more information, as the paper's mechanism predicts?

**Analysis (`analysis/s34_five_study_extensions.py`).** Under the five-study protocol (released folds,
the benchmark's clinical file, seed 0), for bladder, head and neck and stomach, the studies whose files
carry a usable grade: D = [C(ModRank with grade) − C(age, sex, grade)] − [C(ModRank) − C(age, sex,
stage)], all four concordances recomputed on each of 6,000 case resamples (seed 0) drawn once per
study. The four concordances of each D are already committed; B adds D's interval.

**Prediction.** D is smaller in head and neck and stomach, where grade's normalised entropy is 0.643
and 0.623, than in bladder (0.311), and their intervals may cover zero.

## D. Further fitted fusions on identical inputs (all five studies)

**Question.** The paper tests two fitted fusions (concatenation, stacking) against the rank average.
Does any other common way of fitting the combination separate from it?

**Fusions (same script),** each built from the same three arms, folds and seed as the five-study fusion
analysis, and each fitted only inside the training fold on inner out-of-fold percentiles:

1. simplex weights: three non-negative weights summing to one on the arms' percentiles, chosen from a
   grid of step 0.1 (66 points) by inner three-fold concordance;
2. interaction stacking: the stacked Cox with the three pairwise products of arm percentiles added
   (ridge 1.0, as the existing stacking);
3. risk-gated weights: separate simplex weights inside each tertile of the clinical arm's inner
   out-of-fold score, chosen as in 1, applied by the validation patient's clinical tertile.

Each is compared with the rank average by the paired case bootstrap (6,000, seed 0) in each study, 15
comparisons in all, and the rank average must first reproduce its committed value in every study.

**Prediction.** None of the 15 separates from the rank average in either direction. If any fitted
fusion is ahead with an interval excluding zero, the paper says so and narrows its fusion statement to
where it holds.

## C. External replication of D in further public cohorts

Specified separately, after a systematic search and before any download, in `protocol-c.md`.
