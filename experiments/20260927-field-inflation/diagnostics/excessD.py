import json, os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
import numpy as np
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from blca_common import cidx, cpairs, pct
P=os.path.join(ROOT, "experiments") + "/"
rows=[json.loads(x) for x in open(P+"20260911-blca-posthoc/results/percase-canonical-vectors.jsonl")]
fold=np.array([r["fold"] for r in rows]); t=np.array([r["months"] for r in rows]); e=np.array([r["event"] for r in rows],float)
n=len(rows); ii,jj=cpairs(t,e)
def exact(name):
    v=np.array([r[name] for r in rows]); out=np.full(n,np.nan)
    for k in np.unique(fold):
        m=fold==k; out[m]=np.round(v[m]*(m.sum()-1))/(m.sum()-1.0)
    return out
def fp(v, f=fold):
    out=np.full(len(v),np.nan)
    for k in np.unique(f):
        m=f==k; out[m]=pct(v[m])
    return out
zg,zc=exact("clinical_grade"),exact("clinical_stage")
cons={"ModRank":(exact("ours_grade"),exact("ours"),3),"SurvPath":(exact("survpath_plus_grade"),exact("survpath_plus_stage"),2),
      "PIBD":(exact("pibd_best_val_plus_grade"),exact("pibd_best_val_plus_stage"),2)}
rng=np.random.default_rng(11)
def noiseD(k, idx=None):
    z=[rng.random(n) for _ in range(k-1)]
    Ms=fp(sum(fp(x) for x in z)) if k==3 else fp(z[0])
    if k==3:
        mw=fp(fp(z[0])+fp(z[1])+zg); ms=fp(fp(z[0])+fp(z[1])+zc)
    else:
        mw=fp(Ms+zg); ms=fp(Ms+zc)
    if idx is None: idx=np.arange(n)
    i2,j2=cpairs(t[idx],e[idx])
    return (cidx(mw[idx],i2,j2)-cidx(zg[idx],i2,j2))-(cidx(ms[idx],i2,j2)-cidx(zc[idx],i2,j2))
for k in (2,3):
    nd=[noiseD(k) for _ in range(1000)]
    print("BLCA noise-null D, %d-way sum: mean %+.4f sd %.4f q97.5 %+.4f"%(k,np.mean(nd),np.std(nd),np.quantile(nd,.975)))
draws=[np.random.default_rng(20260911).choice(n,size=n,replace=True) for _ in range(1)]
rngb=np.random.default_rng(20260911); draws=[rngb.choice(n,size=n,replace=True) for _ in range(2000)]
for name,(w,s,k) in cons.items():
    D=(cidx(w,ii,jj)-cidx(zg,ii,jj))-(cidx(s,ii,jj)-cidx(zc,ii,jj))
    ex=[]
    for b in draws:
        i2,j2=cpairs(t[b],e[b])
        Db=(cidx(w[b],i2,j2)-cidx(zg[b],i2,j2))-(cidx(s[b],i2,j2)-cidx(zc[b],i2,j2))
        ex.append(Db-noiseD(k,b))
    ex=np.array(ex)
    print("%-9s D %+.4f | D minus no-signal D: mean %+.4f 95%% [%+.4f, %+.4f] p %.4f"%(name,D,ex.mean(),np.quantile(ex,.025),np.quantile(ex,.975),2*min((ex<=0).mean(),(ex>=0).mean())))
# GEO, approximate: cohort-wide percentiles on repeat-averaged vectors
g=[json.loads(x) for x in open(P+"20260927-geo-percase/results/geo-percase.jsonl")]
for coh in ("GSE31684","GSE32894"):
    R=[r for r in g if r["analysis"]==coh]; m=len(R)
    tt=np.array([r["months"] for r in R]); ee=np.array([r["event"] for r in R],float); one=np.zeros(m,int)
    a={k:np.array([r[k] for r in R]) for k in ("clinical_stage","clinical_weak","modrank","modrank_with_weak","transcriptome")}
    i2,j2=cpairs(tt,ee)
    D=(cidx(a["modrank_with_weak"],i2,j2)-cidx(a["clinical_weak"],i2,j2))-(cidx(a["modrank"],i2,j2)-cidx(a["clinical_stage"],i2,j2))
    nd=[]
    for _ in range(1000):
        M=pct(rng.random(m)); 
        nd.append((cidx(pct(M+pct(a["clinical_weak"])),i2,j2)-cidx(a["clinical_weak"],i2,j2))-(cidx(pct(M+pct(a["clinical_stage"])),i2,j2)-cidx(a["clinical_stage"],i2,j2)))
    # same approximation applied to the real transcriptome arm, to calibrate the approximation
    Dapprox=(cidx(pct(pct(a["transcriptome"])+pct(a["clinical_weak"])),i2,j2)-cidx(a["clinical_weak"],i2,j2))-(cidx(pct(pct(a["transcriptome"])+pct(a["clinical_stage"])),i2,j2)-cidx(a["clinical_stage"],i2,j2))
    print("%s n %d ev %d | D committed-vectors %+.4f | same-approximation D %+.4f | approx no-signal D mean %+.4f sd %.4f q97.5 %+.4f"
          %(coh,m,int(ee.sum()),D,Dapprox,np.mean(nd),np.std(nd),np.quantile(nd,.975)))
