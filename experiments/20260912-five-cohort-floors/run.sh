#!/usr/bin/env bash
# A noise floor for each of the benchmark's five studies, built the way the bladder floor was built:
# 24 re-partitions of the same patients, every arm refitted inside each one, and the floor is the
# standard deviation of the PAIRED difference C(ModRank) - C(clinical arm) across them.
#
# The bladder row here is NOT the main text's 0.0073. That number uses the incumbent's clinical
# block and five seeds; this run uses the five-study protocol, the benchmark's released clinical
# file and one seed. Having bladder measured both ways is what makes the other four comparable to
# something already in the paper.
#
# Partition 0 is the released folds and must reproduce the committed point estimates before any
# re-partition is drawn. About fifteen minutes on one CPU host.
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${PYTHON:=/usr/bin/python3}"   # this run used Python 3.9.6, NumPy 2.0.2, SciPy 1.13.1
D=data/frozen-inputs
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 "$PYTHON" -W ignore \
  analysis/s29_five_cohort_floors.py \
  --titan "$D/TCGA_TITAN_features.pkl" --meta-dir "$D/meta" --clin-dir "$D/clin" \
  --pan-dir "$D/pan" --sig "$D/omics/combine_signatures.csv" \
  --frozen experiments/20260911-five-study-repro/results/stage5.json \
  --out experiments/20260912-five-cohort-floors/results/five-cohort-floors.json
