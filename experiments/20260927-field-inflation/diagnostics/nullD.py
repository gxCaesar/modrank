import json, os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
import numpy as np
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from blca_common import cidx, cpairs, pct
rows=[json.loads(x) for x in open(os.path.join(ROOT, "experiments", "20260911-blca-posthoc", "results", "percase-canonical-vectors.jsonl"))]
fold=np.array([r["fold"] for r in rows]); t=np.array([r["months"] for r in rows]); e=np.array([r["event"] for r in rows],float)
ii,jj=cpairs(t,e)
def exact(name):
    v=np.array([r[name] for r in rows]); out=np.full(len(v),np.nan)
    for k in np.unique(fold):
        m=fold==k; out[m]=np.round(v[m]*(m.sum()-1))/(m.sum()-1.0)
    return out
def fp(v):
    out=np.full(len(v),np.nan)
    for k in np.unique(fold):
        m=fold==k; out[m]=pct(v[m])
    return out
zg,zc=exact("clinical_grade"),exact("clinical_stage")
Cg,Cs=cidx(zg,ii,jj),cidx(zc,ii,jj)
rng=np.random.default_rng(7)
res=[]
for _ in range(2000):
    M=fp(rng.random(len(rows)))
    aw=cidx(fp(M+zg),ii,jj)-Cg; as_=cidx(fp(M+zc),ii,jj)-Cs
    res.append((cidx(M,ii,jj),aw,as_,aw-as_))
r=np.array(res)
print("C grade %.4f  C stage %.4f"%(Cg,Cs))
print("random M: C alone mean %.4f | added over grade mean %+.4f | over stage mean %+.4f | D mean %+.4f sd %.4f  q2.5 %+.4f q97.5 %+.4f"
      %(r[:,0].mean(),r[:,1].mean(),r[:,2].mean(),r[:,3].mean(),r[:,3].std(),np.quantile(r[:,3],.025),np.quantile(r[:,3],.975)))
# D as a function of how informative M is: blend the true outcome-free signal? use ModRank-without-clinical scaled by noise
ours_nc=exact("ours_without_clinical")
for w in (0.0,0.25,0.5,1.0,2.0,1e9):
    ds=[]; cs=[]
    for _ in range(300):
        M=fp(ours_nc + (rng.random(len(rows))/w if w>0 else rng.random(len(rows))*1e9) ) if w<1e8 else fp(ours_nc)
        aw=cidx(fp(M+zg),ii,jj)-Cg; as_=cidx(fp(M+zc),ii,jj)-Cs; ds.append(aw-as_); cs.append(cidx(M,ii,jj))
        if w>=1e8: break
    print("signal weight %-6s C alone %.4f  D %+.4f"%(w,np.mean(cs),np.mean(ds)))
