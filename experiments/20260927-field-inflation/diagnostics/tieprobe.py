import json, sys, platform
import numpy as np
rows=[json.loads(x) for x in open(sys.argv[1])]
fold=np.array([r["fold"] for r in rows])
def exact(v):
    out=np.full(len(v),np.nan)
    for k in np.unique(fold):
        m=fold==k; out[m]=np.round(v[m]*(m.sum()-1))/(m.sum()-1.0)
    return out
sp=exact(np.array([r["survpath"] for r in rows])); zc=exact(np.array([r["clinical_stage"] for r in rows]))
s=sp+zc
m=fold==0; v=s[m]
print(platform.machine(), np.__version__)
print("fold0 n", m.sum(), "distinct sums", len(np.unique(v)))
print("argsort default :", np.argsort(v)[:25].tolist())
print("argsort stable  :", np.argsort(v, kind="stable")[:25].tolist())
try:
    np.show_runtime()
except Exception as e: print(e)
