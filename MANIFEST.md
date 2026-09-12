# Manifest

199 files. Every one of them is listed here, and every file listed here is in the tree; the builder asserts both directions and refuses to write this file otherwise.

SHA-256 digests for every file except the digest list itself are in `results_manifest.sha256`.

| directory | files | what it holds |
|---|---|---|
| `analysis/` | 55 | the code that produced every reported number, one script per question |
| `development/` | 33 | the frozen protocol, the split manifest it binds, the ledger of all 269 scored candidates, the error atlas, and the component pre-registrations |
| `experiments/20260817-blca-confirm/` | 29 | the confirmatory run: command, environment, log, code and result files |
| `experiments/20260818-selection-null/` | 4 | the 200 outcome-permuted runs behind the corrected selection bar, and the encoder swap |
| `experiments/20260818-reporting-dump/` | 2 | the reporting dump the manuscript's tables read from |
| `paper/` | 39 | the manuscript sources, the figure builders and their source data, and the number-consistency checker |
| `docs/` | 2 | data card and model card |

## At the top level

| file | what it is |
|---|---|
| `.gitignore` | nothing generated belongs in the repository |
| `.zenodo.json` |  |
| `CITATION.cff` | citation metadata, machine readable |
| `FOR_REVIEWERS.md` | the first file to open, opening with a runnable command |
| `LICENSE` | MIT |
| `MANIFEST.md` | this file |
| `README.md` | what this is, and how to run it |
| `charter.md` | what the study committed to before it ran, including its kill criterion |
| `experiments/20260820-encoder-parity/README.md` |  |
| `experiments/20260820-encoder-parity/results/encoder-parity.json` |  |
| `experiments/20260820-encoder-parity/run.sh` |  |
| `experiments/20260911-blca-posthoc/README.md` |  |
| `experiments/20260911-blca-posthoc/env.txt` |  |
| `experiments/20260911-blca-posthoc/logs/run.log` |  |
| `experiments/20260911-blca-posthoc/logs/s22.log` |  |
| `experiments/20260911-blca-posthoc/logs/s22b.log` |  |
| `experiments/20260911-blca-posthoc/logs/s23.log` |  |
| `experiments/20260911-blca-posthoc/results/calibration-and-decision-curve.json` |  |
| `experiments/20260911-blca-posthoc/results/double-tied-probe.json` |  |
| `experiments/20260911-blca-posthoc/results/modality-atlas-amended.json` |  |
| `experiments/20260911-blca-posthoc/results/percase-canonical-vectors.json` |  |
| `experiments/20260911-blca-posthoc/results/resplit-and-site-cv.json` |  |
| `experiments/20260911-blca-posthoc/results/unified-fusion-and-added-value.json` |  |
| `experiments/20260911-blca-posthoc/results/utility-intervals.json` |  |
| `experiments/20260911-blca-posthoc/run.sh` |  |
| `experiments/20260911-five-study-repro/README.md` |  |
| `experiments/20260911-five-study-repro/results/stage5.json` |  |
| `experiments/20260911-five-study-repro/run.sh` |  |
| `experiments/20260911-geo-external/README.md` |  |
| `experiments/20260911-geo-external/protocol.json` |  |
| `experiments/20260911-geo-external/results/geo-external.json` |  |
| `experiments/20260911-geo-external/results/geo-gated.json` |  |
| `experiments/20260911-geo-gate0/results/gate0.json` |  |
| `results_manifest.sha256` | one SHA-256 per file, C-sorted so two platforms agree |
| `verify_release.py` | recomputes the manuscript's load-bearing numbers from these files |

