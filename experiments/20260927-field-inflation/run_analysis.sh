#!/bin/bash
# Analyses A (s33) and B + D (s34) of protocol.md, on sysu CPU, from copies of the committed code and
# inputs placed in code/ and inputs/ of this run directory. Each script refuses to write unless it first
# reproduces its committed known answers.
#   bash run_analysis.sh a      # after run_a.sh has finished every model
#   bash run_analysis.sh bd
set -eu
RUN=/data1/guanxing/bladder_cancer/experiments/20260927-field-inflation
D=/data1/guanxing/bladder_cancer/data/frozen-inputs
PY=~/miniconda3/envs/hacl/bin/python
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONWARNINGS=ignore
mkdir -p $RUN/analysis-results
case $1 in
  a)  $PY $RUN/code/s33_field_inflation.py \
        --percase $RUN/inputs/percase-canonical-vectors.jsonl \
        --committed $RUN/inputs/unified-fusion-and-added-value.json \
        --runs $RUN/results \
        --models coattn coattn_motcat abmil_wsi_pathways transmil_wsi_pathways deepmisl_wsi_pathways \
        --out $RUN/analysis-results/field-inflation-a.json ;;
  bd) $PY $RUN/code/s34_five_study_extensions.py \
        --titan $D/TCGA_TITAN_features.pkl --meta-dir $D/meta --clin-dir $D/clin \
        --pan-dir $D/pan --sig $D/omics/combine_signatures.csv \
        --frozen $RUN/inputs/stage5.json --frozen-fusion $RUN/inputs/five-cohort-fusion.json \
        --out $RUN/analysis-results/five-study-extensions-bd.json ;;
  *) echo "usage: run_analysis.sh a|bd"; exit 2 ;;
esac
