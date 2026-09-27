#!/bin/bash
# Analysis A of protocol.md: one SurvPath-code-base baseline, five retrainings (seeds 1-5), on one card.
#   bash run_a.sh <model> <nvidia-smi index>
# The SurvPath rerun's environment, features, folds and flags; only --modality, --lr, --max_epochs,
# --fusion and the ABMIL encoder settings follow the model's own upstream script (protocol.md, A).
set -u
MODEL=$1
GPU=$2
ROOT=/data1/guanxing/bladder_cancer
RUN=$ROOT/experiments/20260927-field-inflation
mkdir -p $RUN/results/$MODEL $RUN/logs
ulimit -n 65536
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES=$GPU
export PYTHONPATH=$ROOT/envs/survpath-deps
export PYTHONWARNINGS=ignore
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
PY=~/miniconda3/envs/hacl/bin/python

case $MODEL in
  coattn)                EXTRA="--lr 0.001 --max_epochs 2 --fusion concat" ;;
  coattn_motcat)         EXTRA="--lr 0.0005 --max_epochs 5 --fusion concat" ;;
  abmil_wsi_pathways)    EXTRA="--lr 0.001 --max_epochs 2 --fusion concat --encoding_layer_1_dim 8 --encoding_layer_2_dim 16 --encoder_dropout 0.25" ;;
  transmil_wsi_pathways) EXTRA="--lr 0.001 --max_epochs 2 --fusion concat" ;;
  deepmisl_wsi_pathways) EXTRA="--lr 0.001 --max_epochs 2 --fusion concat" ;;
  *) echo "unknown model $MODEL"; exit 2 ;;
esac

# refuse to train unless torch holds the intended card and it is free (the SurvPath rerun's lesson:
# this box mixes 4x A4000 with 1x A6000 and CUDA's default ordering is not nvidia-smi's)
$PY - <<'PY' || exit 9
import torch
assert torch.cuda.is_available(), "no CUDA visible to torch"
free = torch.cuda.mem_get_info(0)[0] / 1e9
print("PLACEMENT %s, %.1f GB free" % (torch.cuda.get_device_name(0), free))
if free < 10:
    raise SystemExit("REFUSING: the intended card has under 10 GB free")
PY

cd $ROOT/SurvPath
for SEED in 1 2 3 4 5; do
  echo "=== $MODEL seed $SEED START $(date -Is) gpu $GPU"
  $PY main.py \
    --study tcga_blca --task survival --split_dir splits --which_splits 5foldcv \
    --type_of_path combine --modality $MODEL \
    --data_root_dir $ROOT/data/tcga-blca-wsi-pilot/chief \
    --label_file datasets_csv/metadata/tcga_blca.csv \
    --omics_dir datasets_csv/raw_rna_data/combine/blca \
    --results_dir $RUN/results/$MODEL/seed$SEED \
    --batch_size 1 --opt radam --reg 0.0001 \
    --alpha_surv 0.5 --weighted_sample --encoding_dim 768 \
    --label_col survival_months_dss --k 5 --bag_loss nll_surv --n_classes 4 \
    --num_patches 4096 --wsi_projection_dim 256 --seed $SEED $EXTRA
  echo "=== $MODEL seed $SEED EXIT=$? $(date -Is)"
done
