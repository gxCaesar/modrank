#!/usr/bin/env bash
# The two experiments cold panel round 2 named: an outcome-permuted selection null, and the
# encoder swap that separates "method" from "2025 encoder". Post-freeze, descriptive, and neither
# changes a reported number.
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${PERMS:=200}"
: "${PYTHON:=python3}"
D=data/frozen-inputs
"$PYTHON" analysis/s13_selection_null_and_encoder.py \
  --protocol development/benchmark-protocol.json \
  --root "$D" --titan "$D/TCGA_TITAN_features.pkl" \
  --sp-dir "$D/survpath_preds" --omics-dir "$D/omics" --dimaf-dir "$D/dimaf_splits" \
  --fm-dir "$D/fmfeat" --perms "$PERMS" \
  --out experiments/20260818-selection-null/results/selection-null.json

# Where the frozen arm sat among the campaign's comparable candidates, recomputed from the ledger
# rather than remembered. No data, no model: it reads one markdown file.
"$PYTHON" analysis/s15_selection_rank.py \
  --ledger development/iteration-ledger.md \
  --out experiments/20260818-selection-null/results/selection-rank.json
