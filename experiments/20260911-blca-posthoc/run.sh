#!/usr/bin/env bash
# Post-freeze analyses on TCGA-BLCA: fusion on identical inputs (concatenated and stacked Cox),
# the evidence-gated variant's known-answer check, the added-value inflation D for three multimodal
# constructions, and the modality paragraph recomputed with the amended clinical block.
#
# CPU only, local, on the stack the frozen run recorded (Python 3.8.20, NumPy 1.24.4, SciPy 1.10.1).
# Refuses to report unless eight values already in the manuscript reproduce first.
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${PYTHON:=python3}"   # the run used Python 3.8.20 with NumPy 1.24.4 and SciPy 1.10.1; env.txt records it
D=data/frozen-inputs
E=experiments
SP=tcga_blca__nll_surv_a0.5_lr5e-04_l2Weight_0.0001_5foldcv_b1_survival_months_dss_dim1_768_patches_4096_wsiDim_256_epochs_5_fusion_None_modality_survpath_pathT_combine
PB=tcga_blca_b32_survival_months_dss_wsiDim_256_epochs_30_omics_pathways_pathT_combine_s1
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 "$PYTHON" -W ignore \
  analysis/s19_unified_fusion_and_added_value.py \
  --protocol development/benchmark-protocol.json \
  --root "$D" --titan "$D/TCGA_TITAN_features.pkl" --sp-dir "$D/survpath_preds" \
  --omics-dir "$D/omics" --dimaf-dir "$D/dimaf_splits" \
  --sp-seed-dirs "$E/20260817-survpath-seeds23-v2/results_seed2/$SP" \
                 "$E/20260817-survpath-seeds23-v2/results_seed3/$SP" \
                 "$E/20260817-survpath-seeds45-v2/results_seed4/$SP" \
                 "$E/20260817-survpath-seeds45-v2/results_seed5/$SP" \
  --pibd-dirs "$E/20260817-pibd-blca-repro/results/$PB" "$E/20260817-pibd-folds1to4-v3/results/$PB" \
  --out-dir experiments/20260911-blca-posthoc/results

# s19c: the double-tied probe (slide arm on pairs neither the transcriptome arm nor the clinical
# model separates). Gives the August audit note's four values a committed emitter; refuses unless
# they and five known answers reproduce. Added 2026-09-11.
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 "$PYTHON" -W ignore \
  analysis/s19c_double_tied_probe.py \
  --vectors experiments/20260911-blca-posthoc/results/percase-canonical-vectors.json \
  --root "$D" --titan "$D/TCGA_TITAN_features.pkl" --sp-dir "$D/survpath_preds" \
  --dimaf-dir "$D/dimaf_splits" --tnm "$D/tnm.json" \
  --out experiments/20260911-blca-posthoc/results/double-tied-probe.json

# s23b: intervals for the accuracy and net-benefit differences s23 reports as point estimates.
# Reproduces every value s23 committed before resampling. Added 2026-09-11.
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 "$PYTHON" -W ignore \
  analysis/s23b_utility_intervals.py --root "$D" --titan "$D/TCGA_TITAN_features.pkl" \
  --sp-dir "$D/survpath_preds" --omics-dir "$D/omics" --dimaf-dir "$D/dimaf_splits" \
  --committed experiments/20260911-blca-posthoc/results/calibration-and-decision-curve.json \
  --out experiments/20260911-blca-posthoc/results/utility-intervals.json
