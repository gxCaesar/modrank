# 20260927-geo-percase: the external-cohort analysis, re-run with its per-case scores

`analysis/s24_geo_external.py` originally wrote only the summary of each external cohort
(`experiments/20260911-geo-external/results/geo-external.json`). This run repeats the same analysis
with `--percase-out`, which also writes the repeat-averaged score of every arm for every patient, so
that the added-value inflation D can be recomputed from the rows rather than read from the summary.
The option only writes output; no computation reads it.

Inputs: the three GEO series matrices and platform annotations (GSE31684, GSE32894, GSE48075), the
SurvPath combine pathway signatures, and the protocol fixed before any GEO outcome was analysed
(`experiments/20260911-geo-external/protocol.json`). Command: `run.sh`. Environment: `env.txt`
(Python 3.8.10, NumPy 1.24.2, SciPy 1.8.1, CPU only).

## Result

The last lines of `logs/run.log`:

```
  GSE32894: n 224 events 25 | modrank 0.8779 stage 0.8672 omics 0.8163 concat 0.8216 stacked 0.8747 | D +0.0581
  GSE32894 rule block: n 224 events 25 | modrank 0.8651 stage 0.8458 omics 0.8163 concat 0.8235 stacked 0.8598
wrote 480 per-case rows to experiments/20260927-geo-percase/results/geo-percase.jsonl
wrote experiments/20260927-geo-percase/results/geo-external-rerun.json (323 s)
exit 0
```

`results/geo-external-rerun.json` equals the original summary in every field except
`runtime_seconds` (369.2 then, 323.5 now). `results/geo-percase.jsonl` has 480 rows, one per patient
in each of the four analyses that carry sample identifiers (93 + 90 + 224 + 73; the GSE32894
rule-block analysis is not exported). Recomputed from these rows with `analysis/blca_common.cindex`,
D is +0.0581 in GSE32894 and +0.0253 in GSE31684, the values the manuscript reports.
`paper/emit_claim.py geo_inflation_replication --rows experiments/20260927-geo-percase/results/geo-percase.jsonl`
does this and refuses to print if the two disagree.

Not re-run: the permutation-screened variant (`analysis/s24b_geo_gated.py`), whose results are
unchanged.
