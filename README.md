# ModRank

Code, prespecified protocols and results for *Parameter-free multimodal fusion with ModRank and a
corrected clinical reference for bladder cancer survival*.

Two findings, of different kinds.

**A measurement.** On the TCGA bladder whole-slide survival benchmark, every method paper we could
read that reports a clinical-variable baseline builds it from tumour grade, none uses pathologic
stage, and half the papers whose text could be read report no clinical baseline at all. Stage is in
the released files. Rebuilding the clinical reference with it raises its concordance from 0.5664 to
0.6637 on identical cases, above five of the nine published entries, and more than halves the
apparent added value of three multimodal constructions. Muscle-invasive urothelial carcinoma is
94.4% one grade level in this cohort, so a grade-based reference separates few patients. Part of the
difference in added value is produced by the rank combination itself, since a score carrying no
information also gains less over stage than over grade under the same rule, so every model is
measured against that benchmark as well (`analysis/s35_inflation_null.py`).

**A method.** At 359 patients and 113 events an in-sample fit at sixteen parameters reaches 0.7474
on pure noise, so ModRank spends its events on one penalised Cox model per modality and none on the
combination, which is an equal-weight rank average. It reaches 0.7214 with 1,049 coefficients
against two competitors' 25 million, the highest value reported on these folds by point estimate,
and does not separate from either competitor. Across the benchmark's five studies, none of six
fitted fusions of the same inputs significantly improves on it, and it exceeds the corrected
clinical reference in three of the five.

**And one thing this repository reports against itself.** The study's pre-registered decision rule
used the standard deviation of a paired between-arm difference where a maximum-of-N correction needs
that of an absolute score. Two hundred outcome-permuted runs put the corrected bar at 0.7708, and
0.7214 does not clear it. That measurement is here, in
`experiments/20260818-selection-null/`, and it is reported in the paper rather than left out.

## Start here

```bash
python3 verify_release.py
```

It needs Python 3.9 and numpy, nothing else, and it recomputes 45 quantities from the result
files in this repository rather than from anything written down. `FOR_REVIEWERS.md` says what each
one establishes, and the last two check that every artefact the paper's Availability statement
promises is actually here.

## What is where

| | |
|---|---|
| `analysis/` | the scripts that produced every reported number, one per experiment |
| `development/benchmark-protocol.json` | the protocol frozen before the confirmatory run, with amendment A1 recorded in it |
| `experiments/20260817-blca-confirm/` | the confirmatory run: command, environment, results |
| `experiments/20260818-reporting-dump/` | the per-case and per-pathway tables the figures are built from |
| `experiments/20260818-selection-null/` | the permuted-outcome selection null and the encoder swap |
| `paper/figures/*.py` | one builder per figure, each reading committed result files |
| `paper/check_numbers.py` | the gate that fails if the manuscript and the results disagree |
| `docs/DATA_CARD.md` | the cohort, the three inputs, and what one cohort at 113 events cannot support |
| `docs/MODEL_CARD.md` | what the model is, what it must not be used for, and the corrected selection bar it fails |
| `MANIFEST.md` | every file in this repository, checked complete in both directions |

## Inputs are public and are not redistributed here

Every input is a public release and this repository points at them rather than shipping them: the
TITAN slide embeddings, SurvPath's gene matrix and pathway signatures, SurvPath's released per-fold
predictions, and DIMAF's released split files. Put them in one directory and the confirmatory run
executes from this tree:

```bash
BLCA_DATA=/path/to/inputs bash experiments/20260817-blca-confirm/run.sh
```

That script is the frozen confirmatory run and writes `results/confirmatory.json`, whose primary is
**0.7260**. Amendment A1 rebuilt the clinical block from DIMAF's released split files (0.7212), and
amendment A2 corrected the Cox risk sets and tied ranks and regenerated every later result
(**0.7214**, the value the manuscript reports). `experiments/20260818-reporting-dump/run.sh`
reproduces 0.7214, and with `BLCA_A2=0` it reproduces 0.7212 exactly.

## Amendment A2 and reproducing to the fourth decimal

The Cox fit gives every event the full risk set of its tied time (Breslow), and tied scores share
their average rank, so no reported value depends on the order in which rows or ties are sorted. Both
were corrected after review; `experiments/20260928-amendment-a2/` holds the protocol written before
anything was regenerated, the reproduction of every committed pre-A2 file under the original
functions (`BLCA_A2=0`), the driver that ran both passes (`run_jobs.py`) and the headline values
before and after (`results/before-after.json`). The corrected functions have their own tests:

```bash
python3 tests/test_a2_survival_primitives.py
```

Before A2, within-fold ranks broke exact ties in NumPy's unspecified sort order, and
`experiments/20260927-field-inflation/diagnostics/` records how much that mattered (over 1,000 random
tie orders the pre-A2 seed-0 ModRank concordance ranged from 0.7203 to 0.7229).

## Citing

See `CITATION.cff`. The results published here are in whole or part based upon data generated by
the TCGA Research Network: <https://www.cancer.gov/tcga>.

## Licence

MIT, see `LICENSE`. The licence covers this code. It does not extend to the public datasets, which
carry their own terms.
