#!/usr/bin/env bash
# The five-study arms with the intervals they had never carried. Every point estimate here already
# existed in experiments/20260911-five-study-repro/results/stage5.json; this run adds the paired
# case-level bootstrap the manuscript uses for every other quantity it states, and refuses to write
# its output unless it first reproduces those point estimates for all five cohorts.
#
# Inputs are public and are not redistributed: the pan-TCGA TITAN slide embeddings, SurvPath's
# metadata, clinical and combine RNA files, and SurvPath's five-fold case-ID splits. They live under
# data/frozen-inputs/, which is gitignored.
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${PYTHON:=/usr/bin/python3}"   # this run used Python 3.9.6, NumPy 2.0.2, SciPy 1.13.1; env.txt records it
D=data/frozen-inputs
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 "$PYTHON" -W ignore \
  analysis/s27_five_cohort_intervals.py \
  --titan "$D/TCGA_TITAN_features.pkl" --meta-dir "$D/meta" --clin-dir "$D/clin" \
  --pan-dir "$D/pan" --sig "$D/omics/combine_signatures.csv" \
  --frozen experiments/20260911-five-study-repro/results/stage5.json \
  --out experiments/20260912-five-cohort-intervals/results/five-cohort-intervals.json
