#!/usr/bin/env bash
# Subgroup performance on the axes the released data carries, answering the two TRIPOD+AI rows this
# work had marked "Not addressed": 3c (health inequalities) and 14 (fairness).
#
# It refits nothing. Every score it reads is a committed seed-averaged out-of-fold value from the
# confirmatory run, so the strata are partitions of the same numbers the manuscript reports, and the
# run refuses to write unless the pooled values still equal 0.7237 and 0.6644.
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${PYTHON:=/usr/bin/python3}"   # this run used Python 3.9.6, NumPy 2.0.2, SciPy 1.13.1
# Race and ethnicity are not in the benchmark's released files. They are fetched from the Genomic
# Data Commons, frozen into the gitignored inputs tree, and summarised as counts in a committed
# file. Re-running this line is what makes the race strata below reproducible.
"$PYTHON" -W ignore analysis/s30_fetch_gdc_race.py \
  --meta data/frozen-inputs/meta/tcga_blca.csv \
  --out-map data/frozen-inputs/gdc_demographics.json \
  --out-counts experiments/20260912-subgroup-fairness/results/gdc-demographics-counts.json

OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 "$PYTHON" -W ignore \
  analysis/s28_subgroup_fairness.py \
  --vectors experiments/20260911-blca-posthoc/results/percase-canonical-vectors.json \
  --meta data/frozen-inputs/meta/tcga_blca.csv \
  --race data/frozen-inputs/gdc_demographics.json \
  --out experiments/20260912-subgroup-fairness/results/subgroup-fairness.json
