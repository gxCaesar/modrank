#!/usr/bin/env bash
# Reproducibility repair: the five-study stage-vs-grade comparison and the method on five studies,
# re-run from inputs restored on 2026-09-11 (they had lived only in a deleted scratchpad). The run
# must reproduce development/s5-results/stage5.json exactly; check.py compares them field by field.
# Inputs (public, not redistributed): TITAN pan-TCGA slide embeddings; SurvPath's metadata, clinical
# and combine RNA files; SurvPath's five-fold splits (byte-identical to PIBD's copies).
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${PYTHON:=python3}"   # the run used Python 3.8.20 with NumPy 1.24.4 and SciPy 1.10.1; env.txt records it
D=data/frozen-inputs
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 "$PYTHON" -W ignore analysis/s5_stage_five_cohorts.py \
  --titan "$D/TCGA_TITAN_features.pkl" --meta-dir "$D/meta" --clin-dir "$D/clin" \
  --pan-dir "$D/pan" --sig "$D/omics/combine_signatures.csv" \
  --out experiments/20260911-five-study-repro/results/stage5.json
