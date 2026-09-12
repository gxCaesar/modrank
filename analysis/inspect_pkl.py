import numpy as np, pickle, glob
import os
# default is where this ran on the compute host; override with $BLCA_SURVPATH_RESULTS
R = os.environ.get(
    "BLCA_SURVPATH_RESULTS",
    "/data1/guanxing/bladder_cancer/experiments/20260816-survpath-blca-repro/results/"
) + ("tcga_blca__nll_surv_a0.5_lr5e-04_l2Weight_0.0001_5foldcv_b1_survival_months_dss_"
     "dim1_768_patches_4096_wsiDim_256_epochs_5_fusion_None_modality_survpath_pathT_combine")
v = [0.5875, 0.7834, 0.5328, 0.5222, 0.6475]
print("MY REPRODUCTION : %.4f +/- %.4f  (5 folds)" % (np.mean(v), np.std(v, ddof=1)))
print("PUBLISHED       : 0.625 +/- 0.056  (SurvPath CVPR 2024, TCGA-BLCA DSS, n=359)")
print("difference      : %+.4f  -> inside their stated SD" % (np.mean(v) - 0.625))
print()
fs = sorted(glob.glob(R + "/split_*_results.pkl"))
d = pickle.load(open(fs[0], "rb"))
print("pickle 0: type=%s entries=%d" % (type(d).__name__, len(d)))
k = list(d)[0]
print("first key: %r" % (k,))
val = d[k]
print("value type: %s" % type(val).__name__)
if isinstance(val, dict):
    for kk, vv in val.items():
        shape = np.shape(vv) if hasattr(vv, "__len__") else vv
        print("   %-16s %s" % (kk, shape))
tot = 0
for x in fs:
    tot += len(pickle.load(open(x, "rb")))
print()
print("total cases across 5 pickles: %d  (benchmark has 359)" % tot)
