#!/usr/bin/env bash
# Re-run of s24 with a per-case export, on the Linux CPU host and with the system python3 that ran
# the 2026-09-11 analysis. The project root below is that host's.
set -euo pipefail
cd /data1/guanxing/bladder_cancer
OUT=experiments/20260927-geo-percase
mkdir -p $OUT/results $OUT/logs
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 CUDA_VISIBLE_DEVICES=""
{ hostname; date -u +%FT%TZ; /usr/bin/python3 -c "import sys,numpy,scipy; print(sys.version.split()[0], numpy.__version__, scipy.__version__)"; } > $OUT/env.txt
/usr/bin/python3 analysis/s24_geo_external.py \
  --protocol experiments/20260911-geo-external/protocol.json \
  --geo-dir data/geo-gate0 \
  --signatures SurvPath/datasets_csv/metadata/combine_signatures.csv \
  --out $OUT/results/geo-external-rerun.json \
  --percase-out $OUT/results/geo-percase.jsonl \
  > $OUT/logs/run.log 2>&1
echo "exit $?" >> $OUT/logs/run.log
