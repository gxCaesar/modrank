# 20260817-blca-confirm — the frozen confirmatory run

One run, against `development/benchmark-protocol.json` (`frozen_at: 2026-08-17`,
`split_hash 0d0b8f82e7f5aa6e…`). No iteration after the freeze. `phase_of_origin: frozen`.

## Observed

```
$ bash experiments/20260817-blca-confirm/run.sh
seed 0  OURS=0.7291  titan=0.6596 omics=0.6510 clinical=0.6856
seed 1  OURS=0.7264  titan=0.6515 omics=0.6733 clinical=0.6751
seed 2  OURS=0.7186  titan=0.6550 omics=0.6470 clinical=0.6871
seed 3  OURS=0.7290  titan=0.6553 omics=0.6664 clinical=0.6803
seed 4  OURS=0.7268  titan=0.6472 omics=0.6681 clinical=0.6748
{
 "rule": "primary > 0.679 + 0.0239 = 0.7029 AND sd over seeds <= 0.0145",
 "target": 0.7029,
 "primary": 0.726,
 "primary_clears_target": true,
 "seed_sd": 0.0043,
 "seed_sd_within_bar": true,
 "VERDICT": "PASS"
}
```

**Seeds 1–4 had never been executed before this run.** They are the only genuinely unseen axis this
benchmark offers, and they were put in the protocol to catch an arm that only works under one
inner-CV partition. The spread came back at **0.0043**, a third of the reseed bar.

| | value |
|---|---|
| **primary** — mean pooled C-index over seeds 0–4 | **0.7260** |
| SD over seeds | 0.0043 (bar 0.0145) |
| target — DIMAF 0.679 + selection inflation 0.0239 | 0.7029 |
| **clears target by** | **+0.0231** |

### Comparisons, case-level paired bootstrap, 6,000 replicates

| comparison | ours | baseline | Δ | SE | 95% CI | p |
|---|---|---|---|---|---|---|
| vs TITAN alone | 0.7291 | 0.6596 | +0.0696 | 0.0206 | [+0.0286, +0.1103] | **0.0010** |
| vs DIMAF 0.679, published | 0.7260 | 0.679 | +0.0470 | — | unpaired | — |
| vs cheap clinical Cox | 0.7291 | 0.6856 | +0.0435 | 0.0231 | [−0.0022, +0.0883] | 0.063 |
| vs SurvPath + same clinical, **input parity** | 0.7291 | 0.6963 | +0.0331 | 0.0200 | [−0.0061, +0.0721] | 0.097 |
| **without clinical, vs DIMAF** | 0.6841 | 0.679 | +0.0051 | — | **declared in advance as not claimed** | — |

Per fold: 0.6923 · 0.7880 · 0.7815 · 0.6355 · 0.7453. **4 of 5 folds** beat SurvPath+clinical.

On the decile of comparable pairs the clinical arm separates least: **0.6513** against
SurvPath-plus-the-same-clinical at **0.5852**.

## Interpretation

The pre-registered decision rule is met on both of its conditions. The margin over the incumbent
clears the selection-inflation term over 269 scored candidates, not merely the split-reseed bar, and
the arm is stable across the one axis that had never been run.

Two things the number does **not** say, both declared before the run:

1. **The lead requires the clinical staging variables.** WSI + omics alone is 0.6841 against 0.679 —
   inside the bar. This is in the protocol as `declared_in_advance_as_NOT_claimed` and belongs in
   the paper, not in a footnote.
2. **Under case-level resampling the two most important margins sit at p = 0.06 and p = 0.10.**
   At 359 patients and 113 events the bootstrap SE of a C-index difference is ≈0.020, so this cohort
   cannot resolve margins of 0.03–0.05 — which is true of every published number on it, not only of
   ours. The reseed and selection-inflation bars are the right instruments for a claim on a *fixed*
   benchmark; the bootstrap is the right instrument for "would this hold in a new cohort", and it
   says: not established.

## Not run

- **No held-out test.** This benchmark has none, and neither does any published entrant; G14 is
  declared inapplicable in `research_gates.yaml` with the reason, rather than satisfied by pointing
  at a file that does not exist.
- **DIMAF, MMP, OTSurv, APL, MOAD-FNet, DSCASurv, ProtoPathway, MCAT, MOTCat were not rerun.** Their
  values are published. None releases per-case predictions, so none can be given the clinical
  variables our arm uses and none can be paired. SurvPath is the only competitor at input parity.
- **DIMAF's fold identity is contested** between two independent literature passes and is unresolved.
  PIBD's was verified here directly (`cmp`, all five files byte-identical).
- **PIBD's rerun crashed after fold 0** (`FileNotFoundError` on `splits_1.csv`, exit 1); fold 0
  reached 0.6038 against a published five-fold mean of 0.667. Not restarted — a GPU launch needs
  fresh confirmation. Its published value, not our partial rerun, is what the tables quote.

## Provenance

- `run.sh` — the exact command, written before the run
- `env.txt` — host, UTC time, git HEAD, Python 3.8.20 / numpy 1.24.4 / scipy 1.10.1, the frozen
  protocol date, the split hash, and the leakage-falsifier verdict
- `code/` — a copy of the four modules the run imported, so the run is readable without the
  scratchpad
- `results/confirmatory.json` — the full record
- `logs/confirm.log` — stdout and stderr

Gates at the time of the run: G01, G-M1, G-M2, G-M3, G-M4, G-M5, G-M6, G-M11, G-M12 **PASS**;
G14 **FAIL, declared inapplicable**. G-M4 verified by mtime that the freeze (16:03:32) predates the
results (16:20:01) rather than taking it on assertion.
