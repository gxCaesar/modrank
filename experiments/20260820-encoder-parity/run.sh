#!/usr/bin/env bash
# ModRank on the encoder the comparators themselves consumed. Answers cold panel round 2's G8.
#
# The SurvPath rerun in this study read CHIEF tile features. This pools those same features to the
# case and puts them where ModRank's slide arm goes, so the comparison happens on the comparator's
# own ground rather than on a 2025 slide encoder the comparators never saw.
#
# Post-freeze and descriptive. It changes no reported number and the script refuses to report its
# own result unless it first reproduces the shipped sweep's TITAN values through the identical code
# path, which is the only thing that makes a new code path's number worth reading.
#
# CPU only, about 15 seconds, no download and no GPU. All inputs are already local.
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${PYTHON:=/usr/bin/python3}"
D=data/frozen-inputs

"$PYTHON" -W ignore analysis/s16_encoder_parity_chief.py \
  --protocol development/benchmark-protocol.json \
  --root "$D" --titan "$D/TCGA_TITAN_features.pkl" \
  --sp-dir "$D/survpath_preds" --omics-dir "$D/omics" --dimaf-dir "$D/dimaf_splits" \
  --chief "$D/chief_pooled.npz" \
  --sweep experiments/20260818-selection-null/results/selection-null.json \
  --out experiments/20260820-encoder-parity/results/encoder-parity.json
