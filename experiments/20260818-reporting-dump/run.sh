#!/usr/bin/env bash
# Reporting dump against the SAME frozen protocol and the SAME inputs as the confirmatory run.
# Not a new experiment: it re-executes the frozen arm construction and writes out the per-case and
# per-pathway tables that run summarised away. It ASSERTS the recomputed primary equals 0.7212,
# which is the amended value, and refuses to write anything if it does not.
#
# Inputs are public releases and are not redistributed. Set BLCA_DATA to a directory holding them.
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${BLCA_DATA:?set BLCA_DATA to the directory holding the public inputs}"
: "${PYTHON:=python3}"
"$PYTHON" analysis/s12_reporting_dump.py \
  --protocol development/benchmark-protocol.json \
  --root "$BLCA_DATA" \
  --titan "$BLCA_DATA/TCGA_TITAN_features.pkl" \
  --sp-dir "$BLCA_DATA/survpath_preds" \
  --omics-dir "$BLCA_DATA/omics" \
  --dimaf-dir "$BLCA_DATA/dimaf_splits" \
  --out experiments/20260818-reporting-dump/results/reporting-dump.json
