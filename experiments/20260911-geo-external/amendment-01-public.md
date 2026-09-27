<!-- Public copy of amendment-01.md, the file geo-gated.json cites. Two references to an external
cohort whose results are embargoed are withheld, as marked; every other word is the record as
written on 2026-09-11. -->

# GEO protocol amendment 01: the evidence-gated variant, tested only on cohorts not yet seen

**Recorded** 2026-09-11, after the GSE31684 results of `s24` were read and **before** any GSE32894 or
GSE48075 result was inspected. At that moment s24 was still running on those two cohorts; their log
lines and results file had not been opened.

**What was seen that motivates it.** In GSE31684 the transcriptome arm scored 0.4742, below chance,
and the equal-weight two-arm average (0.5898) fell below the clinical stage block alone (0.6590).
[One sentence withheld: it compares this with a result from an external cohort whose terms of use
embargo its results.]

**What is added.** The evidence-gated variant, with its definition copied **unchanged** from an
earlier protocol for that external cohort (frozen at 9e60ada, before any GEO data existed). Inside each outer training fold, an arm enters the rank average only if its inner
out-of-fold concordance exceeds the 95th percentile of 1,000 permutations of the training outcomes
against that fixed score. The permutation seed is `default_rng(10000 + 100*repeat + fold)`. If no
arm passes, all arms enter. It uses the same splits, estimator and comparators as `s24`.

**Where it is tested.** GSE32894 and GSE48075 only. **GSE31684 is excluded from this test** because
its arms' results were already seen. On GSE31684 the gate's behaviour can be predicted from those
results, so running it there would not be a test.

**Prediction, written now.** If a cohort's transcriptome arm is uninformative, the gate drops it and
the variant approaches the clinical arm. If both arms are informative, the variant equals the
two-arm ModRank.

**Decision rule, the same form as the protocol's.** The variant succeeds in a cohort if its mean
concordance exceeds every comparator's (clinical stage, transcriptome, concatenated, stacked, and
the equal-weight two-arm ModRank) by more than that cohort's bar, which is 2 × the split-reseed SD
from s24.

**Status of the variant.** It was defined after an observation in that embargoed cohort, so it is a
hypothesis being tested here, not the frozen method. On TCGA-BLCA it is bit-identical to the frozen
method (s19).
