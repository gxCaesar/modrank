# Manifest

149 files. Every one of them is listed here, and every file listed here is in the tree.

| directory | files | what it holds |
|---|---|---|
| `analysis/` | 44 | the code that produced every reported number, one script per question |
| `development/` | 33 | the frozen protocol, the split manifest it binds, the ledger of all 269 scored candidates, the error atlas, and the component pre-registrations |
| `experiments/20260817-blca-confirm/` | 29 | the confirmatory run: command, environment, log, code and result files |
| `experiments/20260818-selection-null/` | 4 | the 200 outcome-permuted runs behind the corrected selection bar, and the encoder swap |
| `experiments/20260818-reporting-dump/` | 2 | the reporting dump the manuscript's tables read from |
| `paper/` | 26 | the manuscript sources, the figure builders and their source data, and the number-consistency checker |
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
| `verify_release.py` | recomputes the manuscript's load-bearing numbers from these files |
