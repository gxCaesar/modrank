#!/usr/bin/env bash
# Do the fusion alternatives beat the equal-weight rank average on the benchmark's other four
# studies? The development cohort had this comparison and the other four never did, so a claim about
# fitting the fusion rested on one study of five.
#
# Concatenation is one ridge Cox on all three blocks with the penalty grid extended by 32,768 in its
# favour. Stacking learns three weights per fold from the arms' INNER out-of-fold percentiles, which
# is the part that keeps the weights from seeing the scores they weight; both constructions are
# copied from s19 rather than re-derived.
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${PYTHON:=/usr/bin/python3}"   # this run used Python 3.9.6, NumPy 2.0.2, SciPy 1.13.1
D=data/frozen-inputs
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 "$PYTHON" -W ignore \
  analysis/s32_five_cohort_fusion.py \
  --titan "$D/TCGA_TITAN_features.pkl" --meta-dir "$D/meta" --clin-dir "$D/clin" \
  --pan-dir "$D/pan" --sig "$D/omics/combine_signatures.csv" \
  --frozen experiments/20260911-five-study-repro/results/stage5.json \
  --out experiments/20260913-five-cohort-fusion/results/five-cohort-fusion.json
