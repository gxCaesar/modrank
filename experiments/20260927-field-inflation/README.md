# 20260927-field-inflation: analyses A, B and D of protocol.md, and two findings they exposed

Post-freeze, exploratory. `protocol.md` was committed before any run (7bbeb71, corrected pre-run in
fb5f296). Nothing here changes ModRank or any committed number.

## Observed

### A. Four further architectures on TCGA-BLCA (sysu training, Mac analysis)

Training: `bash run_a.sh <model> <gpu>` on sysu, 2026-09-27, seeds 1-5, SurvPath code base with the
upstream scripts' settings. Engineering fix recorded: the four baselines hard-code a 1,024-d WSI input
and the CHIEF features are 768-d, so `core_utils.py` and the four model files were patched to pass
`--encoding_dim` (`survpath-baselines.patch`, originals in `*.before`). The first attempt's logs and
result directories carry the suffix `attempt1-encoding-dim`. MOTCat (`coattn_motcat`) is a stub in the
released code base (`model_motcat.py` is not shipped) and was not run.

Analysis: `analysis/s33_field_inflation.py`, run on the Mac on the 100 pulled `split_k_results.pkl`
(log `logs/analysis-a-mac.log`):

```
known answer reproduced: SurvPath over grade +0.0670, over stage +0.0253, D +0.0417
coattn                   alone 0.5620 | over grade +0.0195 [-0.0243, 0.0634] | over stage -0.0122 [-0.0573, 0.0318] | D +0.0317 [0.0104, 0.054]
abmil_wsi_pathways       alone 0.5930 | over grade +0.0457 [-0.0026, 0.0938] | over stage +0.0155 [-0.0292, 0.0597] | D +0.0302 [0.0072, 0.0551]
transmil_wsi_pathways    alone 0.5964 | over grade +0.0420 [-0.0048, 0.0878] | over stage +0.0136 [-0.0305, 0.055] | D +0.0284 [0.0047, 0.0542]
deepmisl_wsi_pathways    alone 0.5074 | over grade +0.0007 [-0.0465, 0.0479] | over stage -0.0333 [-0.0766, 0.0091] | D +0.0340 [0.0111, 0.0574]
```

The known-answer check reproduces every field of SurvPath's committed block, including all three
bootstrap intervals, and was shown to reject: 500 draws instead of 6,000, and one swapped pair of stage
percentiles in fold 0, each exit 2 with nothing written.

### B and D. Five-study extensions (Mac)

`analysis/s34_five_study_extensions.py`, log `logs/analysis-bd-mac.log`, output
`analysis-results/five-study-extensions-bd.json`. Every committed known answer reproduced (the five
ModRank values, the four concordances behind each B study, and the committed stacked value rebuilt
from the shared inner out-of-fold path).

```
B blca      D +0.0429 [0.0099, 0.0767] | over grade +0.1203 | over stage +0.0774
B hnsc      D +0.0259 [-0.0155, 0.0667] | over grade +0.0988 | over stage +0.0729
B stad      D +0.0275 [-0.0082, 0.0628] | over grade +0.0767 | over stage +0.0491
```

D: none of the 15 comparisons (simplex weights, interaction stacking, risk-gated weights, in five
studies) has an interval excluding zero in either direction. Largest gaps: stomach interaction
stacking behind by 0.0539 ([-0.0099, +0.1211]); breast simplex ahead by 0.0290 ([-0.0815, +0.0239]);
bladder risk-gated behind by 0.0275 ([-0.0002, +0.0567]).

## Finding 1: D has a positive null under the rank combination

A model with no signal (random within-fold ranks), combined with each reference by the paper's rank
rule, has D > 0, because an equal-weight combination dilutes the stronger reference more
(`diagnostics/nullD.txt`: added value over grade -0.0218, over stage -0.0505, D +0.0287). The first
diagnostic (`diagnostics/excessD.py`, one random draw per bootstrap resample) put no construction's
excess over this benchmark beyond its interval. That estimator carries the random arm's own draw
variance into the interval and is conservative. The committed analysis, `analysis/s35_inflation_null.py`
(20 fixed random draws averaged on each resample, 6,000 resamples, seed 20260911, known answers
first), gives:

```
ModRank                D +0.0610 | no-signal D +0.0417 [0.0289, 0.0538] | excess +0.0193 [0.0022, 0.0359] p 0.0263
SurvPath               D +0.0417 | no-signal D +0.0288 [0.016, 0.0413] | excess +0.0129 [-0.0035, 0.027] p 0.1373
PIBD                   D +0.0496 | no-signal D +0.0287 [0.0154, 0.0406] | excess +0.0210 [0.0066, 0.0393] p 0.0033
coattn                 D +0.0317 | no-signal D +0.0287 [0.0158, 0.0409] | excess +0.0031 [-0.0122, 0.0149] p 0.7670
abmil_wsi_pathways     D +0.0302 | no-signal D +0.0288 [0.0165, 0.0412] | excess +0.0014 [-0.0139, 0.0196] p 0.7310
transmil_wsi_pathways  D +0.0284 | no-signal D +0.0286 [0.016, 0.041] | excess -0.0001 [-0.0152, 0.0195] p 0.7873
deepmisl_wsi_pathways  D +0.0340 | no-signal D +0.0287 [0.0157, 0.0412] | excess +0.0053 [-0.009, 0.0249] p 0.3500
GSE31684  n 93 ev 38 | D +0.0253 | no-signal D +0.0150 [-0.0148, 0.0436] | excess +0.0103 [-0.0325, 0.0335] p 0.9813
GSE32894  n 224 ev 25 | D +0.0581 | no-signal D +0.0183 [-0.0069, 0.042] | excess +0.0397 [0.004, 0.0775] p 0.0317
```

(BLCA on the Mac, `logs/s35-blca-mac.log`; GEO on sysu with the interpreter that produced the
committed GEO values, `logs/s35-geo-sysu.log`. The approximate GEO null in `diagnostics/excessD.txt`
used one cohort-wide percentile of repeat-averaged vectors and is superseded.) Reading: D splits into
a part the combination produces for any score and an excess that appears only for models with signal.
The four weak architectures of A (alone 0.51-0.60) have no excess. ModRank's and PIBD's excess
intervals exclude zero, SurvPath's does not, and in GEO the Lund cohort's does and the cystectomy
cohort's does not. The no-signal comparison is post hoc.

## Finding 2: the tie order of the rank function depends on the NumPy build

`pct()` ranks by `np.argsort(np.argsort(v))`, whose tie order is unspecified. Rank sums carry many
ties (fold 0: 58 distinct sums among 70 cases), and so do clinical arms with discrete covariates. On
the same input (`diagnostics/tieprobe-*.txt`) NumPy 1.24.4 and 2.0.2 on the Mac (arm64) and NumPy
1.24.2 on sysu (x86-64) return one order, and NumPy 1.26.4 on sysu, whose argsort dispatches to
AVX-512 code, returns another. Under 1.26.4, s33 and s34 refused to write
(`diagnostics/s33-sysu-numpy126.txt`, `diagnostics/s34-sysu-numpy126.txt`, largest drift 0.0035 on
the head and neck age + sex + grade arm); under 1.24.2 on the same machine s33 reproduced its known
answer and MCAT's values exactly (`diagnostics/s33-sysu-numpy124.txt`). Over 1,000 random tie
orders the concordance of seed-0 ModRank ranges 0.7203-0.7229 (sd 0.00043, committed 0.7225,
mid-rank 0.7215, `diagnostics/tiespread.txt`). Tie order is unrelated to outcome, so this is noise,
not bias, but the fourth decimal reproduces only under the builds that produced the values.

## Not run

- MOTCat: no implementation in the released code base. Fetching the official MOTCat code needs a
  download approval; it was not requested.
- C (external cohorts): the search is running; no protocol written and nothing downloaded.
- No manuscript, claim-registry or release change has been made on the basis of anything here.
