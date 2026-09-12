# 20260911-geo-external: the combination rule and D in public GEO bladder cohorts

Pre-registered: `protocol.json`, frozen at f0fc6bf before any GEO outcome was analysed. Amendment
01 (`amendment-01.md`, 0d89861) adds the evidence-gated variant for GSE32894 and GSE48075 only.
Code: `analysis/s24_geo_external.py`, `analysis/s24b_geo_gated.py`. Run on sysu with the system
python3; the logs are `run.log` and `run_s24b.log` on the host.

## Observed (s24)

Two-arm ModRank uses the transcriptome pathways and the clinical stage block. Values are
5×5-CV means. The bar is 2 × the SD of C(ModRank) − C(stage) over 24 resplits.

| cohort | n / DSS deaths | ModRank | stage block | stacked | concatenated | transcriptome | bar | clears all? |
|---|---|---|---|---|---|---|---|---|
| GSE31684 (cystectomy, primary) | 93 / 38 | 0.5898 | 0.6590 | 0.6317 | 0.4895 | 0.4742 | 0.0423 | no |
| GSE32894 (Lund, mixed) | 224 / 25 | 0.8779 | 0.8672 | 0.8747 | 0.8216 | 0.8163 | 0.0217 | no |
| GSE48075 (MDA) | 73 / 35 | 0.6384 | 0.6322 | 0.6201 | 0.5864 | 0.5683 | 0.0540 | no |

Paired comparisons use repeat-averaged vectors, 6,000 case resamples, and Holm within each cohort.

- **GSE31684.** ModRank vs stage −0.0537 [−0.1164, 0.0131]. vs transcriptome +0.1401, Holm 0.0028.
  vs concatenated +0.1154, Holm 0.012.
- **GSE32894.** ModRank vs stage +0.0185 [−0.0245, 0.0622]. vs stacked +0.0034 [−0.0044, 0.0129].
  vs transcriptome +0.0641, Holm 0.004. vs concatenated +0.0543, Holm 0.004.
- **GSE48075.** No comparison significant, with n = 73.

**D, added-value inflation** from a grade-based reference:

| cohort | D | 95% interval | p |
|---|---|---|---|
| GSE32894 | +0.0581 | [0.0175, 0.0971] | 0.0073 |
| GSE31684 | +0.0253 | [−0.018, 0.068] | 0.258 |
| GSE31684 without the 3 pre-cystectomy chemotherapy cases | +0.1094 | — | — |

The protocol predicted a *smaller* D in GSE32894 than in GSE31684. It was larger. The prediction was
wrong and is reported as written.

## Reading, stated once

- The inflation that a grade-based clinical reference produces **replicates in an independent
  cohort** (GSE32894, p = 0.007) and is positive in the other.
- The two-arm combination rule **does not replicate beyond the floor** in any of the three cohorts.
  - Where the microarray transcriptome arm is informative (GSE32894, GSE48075), ModRank has the best
    point estimate but sits within the floor of the clinical stage block.
  - Where the arm is uninformative (GSE31684, 0.4742), the equal-weight average is diluted below
    the stage block. That is the dilution an equal-weight rule has by construction when one arm
    carries nothing.

## Observed (s24b, amendment 01: the evidence-gated variant)

| cohort | gated | ModRank | stage | margin over stage | margin over stacked | bar | succeeds |
|---|---|---|---|---|---|---|---|
| GSE32894 | 0.8779 | 0.8779 | 0.8672 | +0.0107 | +0.0032 | 0.0217 | no |
| GSE48075 | 0.6326 | 0.6384 | 0.6322 | +0.0004 | +0.0125 | 0.0540 | no |

- In GSE32894 both arms pass the gate in every fold, so the variant equals ModRank.
- In GSE48075 the gate sometimes drops the transcriptome arm, and the variant moves toward the stage
  block.
- The gate prevents dilution. It cannot make the combination exceed its best informative arm, and on
  these cohorts nothing exceeds the clinical stage block beyond the floor.

## Amendment 01 in the public archive

`amendment-01.md` is withheld from the public archive while a publication embargo on another cohort
holds, because its motivating paragraph refers to a result from that cohort. Its operative text,
reproduced here unchanged in substance:

- **When.** Recorded 2026-09-11 (commit 0d89861), after the GSE31684 results of `s24` were read and
  before any GSE32894 or GSE48075 result was inspected.
- **What.** The evidence-gated variant, with its definition copied unchanged from an external
  protocol frozen at 9e60ada, before any GEO data existed. Inside each outer training fold an arm
  enters the rank average only if its inner out-of-fold concordance exceeds the 95th percentile of
  1,000 permutations of the training outcomes against that fixed score
  (`default_rng(10000 + 100*repeat + fold)`). If no arm passes, all arms enter.
- **Where.** GSE32894 and GSE48075 only. GSE31684 is excluded because its arms' results were already
  seen.
- **Prediction.** If a transcriptome arm is uninformative the gate drops it and the variant approaches
  the clinical arm. If both arms are informative the variant equals two-arm ModRank.
- **Decision rule.** The variant succeeds in a cohort if its mean concordance exceeds every comparator
  (clinical stage, transcriptome, concatenated, stacked and two-arm ModRank) by more than that
  cohort's bar, 2 x the split-reseed SD from `s24`.
- **Status.** A hypothesis tested here, not the frozen method. On TCGA-BLCA it is bit-identical to the
  frozen method (s19).
