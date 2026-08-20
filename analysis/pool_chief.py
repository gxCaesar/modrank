#!/usr/bin/env python3
"""Per-slide distributional summaries of the CHIEF/CTransPath patch features.

Every statistic is computed from ONE slide's own patches and nothing else, so there is no fold,
no fitting and no leakage question to answer -- which is why this is the first patch-level thing
to run rather than a clustering that would need a fold-aware basis.

Mean pooling is what every slide-level arm in this campaign has used so far, and it is the one
summary that cannot see heterogeneity. std, and the gap between the 90th and 10th percentile, can.
"""
import glob, os, sys
import numpy as np
import torch

src, out = sys.argv[1], sys.argv[2]
files = sorted(glob.glob(os.path.join(src, "*.pt")))
names, M, S, Q10, Q90, NP = [], [], [], [], [], []
for i, f in enumerate(files):
    x = torch.load(f, map_location="cpu")
    a = x.numpy().astype(np.float32) if hasattr(x, "numpy") else np.asarray(x, dtype=np.float32)
    if a.ndim != 2:
        print("skip %s shape %s" % (f, a.shape)); continue
    names.append(os.path.basename(f)[:-3])
    M.append(a.mean(0)); S.append(a.std(0))
    q = np.quantile(a, [0.1, 0.9], axis=0)
    Q10.append(q[0]); Q90.append(q[1]); NP.append(a.shape[0])
    if (i + 1) % 50 == 0:
        print("%d/%d" % (i + 1, len(files)), flush=True)
np.savez_compressed(out, names=np.array(names), mean=np.vstack(M), std=np.vstack(S),
                    q10=np.vstack(Q10), q90=np.vstack(Q90), n_patches=np.array(NP))
print("wrote", out, len(names), "slides")
