# For reviewers

Start by running the check. It needs nothing installed beyond Python 3.9 and numpy:

```bash
python3 verify_release.py
```

It recomputes, from the result files in this bundle and not from anything written down: the gaps
between consecutive published entries, the capacity sweep's excess over its matched null, the seven
components' arithmetic against their controls, and the six-comparison family's deltas against the
canonical arm table, the selection bar under an outcome-permuted null, and the rank of the frozen
arm re-executed from the ledger. Every one of those is a number the manuscript quotes, and the bar
is one the reported result does not clear.

## What is here

| | |
|---|---|
| `analysis/` | 63 scripts, the code that produced every reported number. `analysis/README.md` maps each to what it establishes |
| `development/benchmark-protocol.json` | the protocol, frozen before the confirmatory run, with amendment A1 recorded in it |
| `development/split-manifest.json` | the split manifest whose SHA-256 the protocol binds |
| `development/iteration-ledger.md` | all 269 scored candidates, which is what the selection-inflation term is computed over |
| `development/contribution-design.md` | the seven components as pre-registered, each with its falsifier and matched control |
| `experiments/20260817-blca-confirm/` | the confirmatory run: command, environment, and 20 result files |
| `experiments/20260818-selection-null/` | the 200 outcome-permuted runs behind the selection bar, and the rank of the frozen arm |
| `paper/figures/` | the figure builders and their source data |

## What is deliberately NOT here

No slides, no expression matrices, no split files: none of the inputs is ours to redistribute. They
are the released SurvPath split files, the public precomputed TITAN embeddings, and DIMAF's released
split files, each cited in the manuscript's Data availability.

Also absent by construction: this project's own control plane. The bundle is built from an explicit
allowlist rather than by exporting the working tree, so no approval receipt, pipeline state or
review trace can reach it.

## The three things worth checking first

**The pre-registered selection correction was mis-specified, and the paper reports that rather than
the version that passes.** The rule used the standard deviation of a paired between-arm difference
where the maximum-of-N formula needs the standard deviation of a candidate's own absolute score. An
outcome-permuted null measures the second quantity at 0.0397 against the 0.0073 assumed, and under
either correct construction the reported 0.7212 does not clear the bar. `verify_release.py`
recomputes both bars and asserts the failure.

**The protocol's freeze predates the confirmatory run**, and that is the claim everything else rests
on. `development/benchmark-protocol.json` carries its own timestamps and amendment A1 is recorded in
it rather than applied silently.

**Two of the three self-reported corrections moved the headline number down.** They are in the
manuscript's supplementary and the code that produced both the before and the after is in
`analysis/`.
