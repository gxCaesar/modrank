# 20260911-five-study-repro: the five-study numbers reproduced from restored inputs

The inputs to `analysis/s5_stage_five_cohorts.py` for the four non-bladder studies had lived only
in a session scratchpad that was later deleted, so the supplement's five-study numbers could no
longer be recomputed. On 2026-09-11 they were restored from the SurvPath clone on sysu into
`data/frozen-inputs/pan/`, which is gitignored and not redistributed:

- `rna_<study>.csv`, SurvPath's `raw_rna_data/combine/<study>/rna_clean.csv`. sha256 was checked
  against sysu, and all five match;
- `splits/<study>/`, SurvPath's five-fold case-ID splits.

**Observed:** `run.sh` rewrote `results/stage5.json`, and all 80 numeric fields equal
`development/s5-results/stage5.json`, with 0 differences.

**Not done:** `why-grade-fails.json` and `decomp.json` still have no committed producer. Their
values match stage5.json and amendment A1 respectively, but no script in the repository writes
them.
