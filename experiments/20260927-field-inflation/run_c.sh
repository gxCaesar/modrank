#!/bin/bash
# Analysis C of protocol-c.md, on sysu with the system python3 (NumPy 1.24.2), which reproduces the
# committed GEO values. gate0 counts before any outcome; analyze reproduces the two earlier GEO
# cohorts' clinical blocks first; summary runs anywhere from the committed JSONs.
set -eu
ROOT=/data1/guanxing/bladder_cancer
RUN=$ROOT/experiments/20260927-field-inflation
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=/usr/bin/python3
cd $ROOT
$PY -W ignore $RUN/code/s36_reference_gap_external.py gate0 --geo-c data/geo-c \
    --retained data/public-nmibc-progression --out $RUN/analysis-results/reference-gap-gate0.json
$PY -W ignore $RUN/code/s36_reference_gap_external.py analyze --geo-c data/geo-c \
    --retained data/public-nmibc-progression --gate0 $RUN/analysis-results/reference-gap-gate0.json \
    --geo-dir data/geo-gate0 --committed-geo $RUN/inputs/geo-external.json \
    --out $RUN/analysis-results/reference-gap-external.json
# then, from the project root on any machine:
#   python3 analysis/s36_reference_gap_external.py summary \
#     --gap experiments/20260927-field-inflation/analysis-results/reference-gap-external.json \
#     --stage5 experiments/20260911-five-study-repro/results/stage5.json \
#     --intervals experiments/20260912-five-cohort-intervals/results/five-cohort-intervals.json \
#     --why experiments/20260817-blca-confirm/results/why-grade-fails.json \
#     --out experiments/20260927-field-inflation/analysis-results/reference-gap-summary.json
