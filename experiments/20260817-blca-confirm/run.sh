#!/usr/bin/env bash
# Confirmatory run against the frozen protocol. Executed once. No iteration.
#
# INPUTS ARE NOT REDISTRIBUTED HERE. They are public releases and this bundle points at them
# rather than shipping them: the TITAN slide embeddings, SurvPath's rna_combine.csv and
# combine_signatures.csv, SurvPath's released per-fold predictions, and DIMAF's released split
# files. Set BLCA_DATA to a directory holding them and this script runs from the bundle:
#
#   BLCA_DATA=/path/to/inputs bash experiments/20260817-blca-confirm/run.sh
#
# NOTE ON THE NUMBER THIS PRODUCES. It writes results/confirmatory.json, whose primary is 0.7260.
# The manuscript reports 0.7212, which is the value AFTER amendment A1 rebuilt the clinical block
# from DIMAF's released split files. A1 is a separate, also-shipped step; see
# analysis/s6_amend_clinical.py and results/amendment-A1-clinical-provenance.json. The amended
# value is what experiments/20260818-reporting-dump/run.sh reproduces, and that script asserts it
# equals 0.7212 before writing anything.
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${BLCA_DATA:?set BLCA_DATA to the directory holding the public inputs}"
: "${PYTHON:=python3}"
"$PYTHON" experiments/20260817-blca-confirm/code/s6_confirmatory.py \
  --protocol development/benchmark-protocol.json \
  --root "$BLCA_DATA" \
  --titan "$BLCA_DATA/TCGA_TITAN_features.pkl" \
  --sp-dir "$BLCA_DATA/survpath_preds" \
  --omics-dir "$BLCA_DATA/omics" \
  --out-dir experiments/20260817-blca-confirm/results
