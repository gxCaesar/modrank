import json, os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
import numpy as np
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from blca_common import cidx, cpairs, pct
from scipy.stats import rankdata
rows=[json.loads(x) for x in open(os.path.join(ROOT, "experiments", "20260911-blca-posthoc", "results", "percase-canonical-vectors.jsonl"))]
fold=np.array([r["fold"] for r in rows]); t=np.array([r["months"] for r in rows]); e=np.array([r["event"] for r in rows],float)
ii,jj=cpairs(t,e)
def exact(name):
    v=np.array([r[name] for r in rows]); out=np.full(len(v),np.nan)
    for k in np.unique(fold):
        m=fold==k; out[m]=np.round(v[m]*(m.sum()-1))/(m.sum()-1.0)
    return out
def fp(v, mode, rng=None):
    out=np.full(len(v),np.nan)
    for k in np.unique(fold):
        m=fold==k; x=v[m]; n=m.sum()
        if mode=="default": r=np.argsort(np.argsort(x))
        elif mode=="stable": r=np.argsort(np.argsort(x,kind="stable"),kind="stable")
        elif mode=="random": o=np.lexsort((rng.random(n),x)); r=np.empty(n); r[o]=np.arange(n)
        elif mode=="mid": r=rankdata(x)-1
        out[m]=r/(n-1.0)
    return out
C=lambda v: cidx(v,ii,jj)
rng=np.random.default_rng(1)
cases={
 "seed0 ModRank (slide+omics+stage clin)": (exact("seed0_slide")+exact("seed0_omics")+exact("seed0_clinical_stage"), "seed0_ours"),
 "SurvPath + stage": (exact("survpath")+exact("clinical_stage"), "survpath_plus_stage"),
 "SurvPath + grade": (exact("survpath")+exact("clinical_grade"), "survpath_plus_grade"),
 "PIBD(best val) + stage": (exact("pibd_best_val")+exact("clinical_stage"), "pibd_best_val_plus_stage"),
}
for name,(s,col) in cases.items():
    com=C(exact(col))
    d=fp(s,"default"); st=fp(s,"stable"); mid=fp(s,"mid")
    rr=np.array([C(fp(s,"random",rng)) for _ in range(1000)])
    nt=sum(len(x)-len(np.unique(x)) for x in (s[fold==k] for k in np.unique(fold)))
    print("%-40s committed %.4f | mac-default %.4f | stable %.4f | midrank %.4f | random: min %.4f max %.4f sd %.5f | tied cases %d"
          %(name,com,C(d),C(st),C(mid),rr.min(),rr.max(),rr.std(),nt))
