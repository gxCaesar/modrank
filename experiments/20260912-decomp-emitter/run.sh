#!/usr/bin/env bash
# The clinical-block decomposition, rebuilt from committed code. `development/s5-results/decomp.json`
# was the one file whose numbers reached a display item with no committed producer; this replaces it
# for the rows that can be regenerated and reports the ones that cannot.
#
# B and D are the known answers, because their pairs also appear in the amendment record, which has
# its own producer. A is rebuilt and compared: it does NOT match the old file, and the rebuild is the
# value that agrees with the rest of the project. C is not attempted at all -- the rule that chose
# one diagnosis per case was never recorded, so rebuilding it would mean inventing a rule and seeing
# whether the invention matched.
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${PYTHON:=/usr/bin/python3}"   # this run used Python 3.9.6, NumPy 2.0.2, SciPy 1.13.1
D=data/frozen-inputs
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 "$PYTHON" -W ignore \
  analysis/s31_decomp_emitter.py \
  --protocol development/benchmark-protocol.json \
  --root "$D" --titan "$D/TCGA_TITAN_features.pkl" --sp-dir "$D/survpath_preds" \
  --omics-dir "$D/omics" --dimaf-dir "$D/dimaf_splits" \
  --frozen development/s5-results/decomp.json \
  --out experiments/20260912-decomp-emitter/results/clinical-block-decomposition.json
