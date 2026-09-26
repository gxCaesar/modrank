#!/usr/bin/env bash
# The five-study validation, rerun with a different pretrained slide encoder and nothing else
# changed. It answers the question the validation otherwise invites: whether three of five studies
# separate because of the method or because of TITAN.
#
# GigaSSL-GigaPath is the alternative because it is the one encoder on disk that covers all five
# studies (blca 359, brca 863, coadread 298, hnsc 394, stad 318). It is a different architecture
# trained by a different group under a different objective, so agreement between the two is not
# agreement between two versions of one thing.
#
# Known answers are deliberately skipped here: they were measured on TITAN, and this run replaces
# the arm they describe.
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${PYTHON:=/usr/bin/python3}"   # this run used Python 3.9.6, NumPy 2.0.2, SciPy 1.13.1
D=data/frozen-inputs
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 "$PYTHON" -W ignore \
  analysis/s27_five_cohort_intervals.py \
  --titan "$D/TCGA_TITAN_features.pkl" \
  --slide-npz "$D/fmfeat/gigassl_gigapath.npz" \
  --meta-dir "$D/meta" --clin-dir "$D/clin" --pan-dir "$D/pan" \
  --sig "$D/omics/combine_signatures.csv" \
  --frozen experiments/20260911-five-study-repro/results/stage5.json \
  --out experiments/20260913-encoder-sensitivity/results/five-cohort-gigassl.json
