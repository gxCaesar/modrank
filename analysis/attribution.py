#!/usr/bin/env python3
"""Is the advantage the clinical variables, or the representation and the fusion?

Pre-H1. NOT reportable. This is the strongest attack available on our own result, and it needs no
GPU because SurvPath's per-case scores are already on disk.

THE OBJECTION. Our arm reaches 0.7150 against SurvPath's 0.6118 on the released folds, but the
inputs differ: ours is a TITAN slide embedding plus clinical staging, theirs is a slide embedding
plus transcriptomics. A reviewer's first question is whether we win because the method is better or
because we brought T, N, M and stage to a fight where the opponent had none. That question has a
direct answer rather than a rhetorical one.

THE TEST. Give SurvPath the same clinical variables, in the same way we use them. Its out-of-fold
risk score is rank-averaged, within fold, with the SAME clinical Cox model our arm uses, fitted
out-of-fold on the same folds. Then:

  survpath + clinical  ~=  ours     -> the advantage is the clinical data, not the method. The
                                       contribution shrinks to "we added stage", and it must be
                                       written that way.
  survpath + clinical   <  ours     -> the TITAN representation is carrying something SurvPath's
                                       WSI-plus-omics score does not, and the fusion is not merely
                                       a clinical wrapper.

Both outcomes are reported. The point of running it is that the first one is entirely possible.

Two further decompositions are included because they bound the answer from the other side:
clinical alone, which is the floor any multimodal claim must clear, and TITAN alone, which says how
much the image contributes without any clinical help.
"""

from __future__ import annotations

import collections
import csv
import glob
import json
import pickle
import sys

import numpy as np


def cpairs(time, event):
    ii, jj = [], []
    for i in np.flatnonzero(event == 1):
        later = np.flatnonzero(time > time[i])
        if later.size:
            ii.append(np.full(later.size, i))
            jj.append(later)
    if not ii:
        return np.empty(0, int), np.empty(0, int)
    return np.concatenate(ii), np.concatenate(jj)


def cidx(risk, ii, jj):
    if ii.size == 0:
        return float("nan")
    a, b = risk[ii], risk[jj]
    return float((np.sum(a > b) + 0.5 * np.sum(a == b)) / ii.size)


def cindex(risk, t, e):
    ii, jj = cpairs(t, e)
    return cidx(risk, ii, jj)


def cox_fit(X, t, e, alpha):
    from scipy.optimize import minimize
    o = np.argsort(t)
    Xo, eo = X[o], e[o]
    ev = np.flatnonzero(eo == 1)
    if ev.size < 2:
        return np.zeros(X.shape[1])

    def f(b):
        eta = Xo @ b
        m = float(eta.max())
        ex = np.exp(eta - m)
        rs = np.cumsum(ex[::-1])[::-1]
        rx = np.cumsum((Xo * ex[:, None])[::-1], axis=0)[::-1]
        ll = float(np.sum(eta[ev] - (np.log(rs[ev]) + m)))
        g = (Xo[ev] - rx[ev] / rs[ev, None]).sum(0)
        return -ll + 0.5 * alpha * float(b @ b), -g + alpha * b

    return minimize(f, np.zeros(X.shape[1]), jac=True, method="L-BFGS-B",
                    options={"maxiter": 400, "gtol": 1e-8}).x


ALPHAS = [1.0, 8.0, 64.0, 512.0, 4096.0]


def fitapply(Xtr, ttr, etr, Xte, seed=0):
    best, ba = -1.0, ALPHAS[0]
    inner = np.array_split(np.random.default_rng(seed).permutation(len(ttr)), 3)
    for al in ALPHAS:
        ip = np.zeros(len(ttr))
        for ite in inner:
            itr = np.setdiff1d(np.arange(len(ttr)), ite)
            mu, sd = Xtr[itr].mean(0), Xtr[itr].std(0) + 1e-9
            ip[ite] = ((Xtr[ite] - mu) / sd) @ cox_fit((Xtr[itr] - mu) / sd, ttr[itr], etr[itr], al)
        c = cindex(ip, ttr, etr)
        if c == c and c > best:
            best, ba = c, al
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
    return ((Xte - mu) / sd) @ cox_fit((Xtr - mu) / sd, ttr, etr, ba)


def pct(v):
    return np.argsort(np.argsort(v)).astype(float) / max(1.0, len(v) - 1.0)


def main():
    root, titan_pkl, sp_dir = sys.argv[1], sys.argv[2], sys.argv[3]

    sp = {}
    for k, f in enumerate(sorted(glob.glob(sp_dir + "/split_*_results.pkl"))):
        d = pickle.load(open(f, "rb"))
        for c in d:
            sp[c] = {"risk": float(d[c]["risk"]), "time": float(d[c]["time"]),
                     "event": 1.0 if float(d[c]["censorship"]) == 0.0 else 0.0, "fold": k}

    meta = list(csv.DictReader(open(root + "/tcga_blca_meta.csv")))
    stem2case = {r["slide_id"].replace(".svs", ""): r["case_id"] for r in meta}
    lab = {}
    for r in meta:
        lab.setdefault(r["case_id"], r)
    folds = []
    for k in range(5):
        tr, va = set(), set()
        for r in csv.DictReader(open(root + "/splits/splits_%d.csv" % k)):
            if r.get("train", "").strip():
                tr.add(r["train"].strip())
            if r.get("val", "").strip():
                va.add(r["val"].strip())
        folds.append((tr, va))

    def norm(x):
        x = str(x)
        for s in (".svs", ".h5", ".pt"):
            if x.endswith(s):
                x = x[:-len(s)]
        return x

    d = pickle.load(open(titan_pkl, "rb"))
    E, F = d["embeddings"], [norm(x) for x in d["filenames"]]
    per = collections.defaultdict(list)
    for i, f in enumerate(F):
        c = stem2case.get(f)
        if c:
            per[c].append(i)
    tnm = json.load(open(root + "/tnm.json"))

    keep = sorted(set(sp) & set(per))
    idx = {c: i for i, c in enumerate(keep)}
    T = np.vstack([E[per[c]].mean(0) for c in keep]).astype(float)
    t = np.array([sp[c]["time"] for c in keep])
    e = np.array([sp[c]["event"] for c in keep])
    spr = np.array([sp[c]["risk"] for c in keep])
    spfold = np.array([sp[c]["fold"] for c in keep])

    age = np.array([float(lab[c]["age"]) if lab[c]["age"].strip() else 65.0 for c in keep])[:, None]
    fem = np.array([float(lab[c]["is_female"]) if lab[c]["is_female"].strip() else 0.0
                    for c in keep])[:, None]

    def oh(v):
        lv = sorted({x for x in v if x})
        return (np.array([[1.0 if x == l else 0.0 for l in lv] for x in v])
                if lv else np.zeros((len(v), 0)))

    def g(c, k):
        x = (tnm.get(c) or {}).get(k)
        return "" if x in (None, "NX", "MX", "TX") else str(x)

    CLIN = np.hstack([age, fem,
                      oh([g(c, "ajcc_pathologic_t")[:2] for c in keep]),
                      oh([g(c, "ajcc_pathologic_n") for c in keep]),
                      oh([g(c, "ajcc_pathologic_m") for c in keep]),
                      oh([str((tnm.get(c) or {}).get("ajcc_pathologic_stage") or "") for c in keep])])

    # out-of-fold component scores on the identical folds
    p_titan = np.full(len(keep), np.nan)
    p_clin = np.full(len(keep), np.nan)
    for k, (tr, va) in enumerate(folds):
        tri = np.array([idx[c] for c in tr if c in idx])
        vai = np.array([idx[c] for c in va if c in idx])
        if len(tri) < 30 or not len(vai):
            continue
        p_titan[vai] = fitapply(T[tri], t[tri], e[tri], T[vai], 0)
        p_clin[vai] = fitapply(CLIN[tri], t[tri], e[tri], CLIN[vai], 0)

    def fold_pct(v):
        out = np.full(len(v), np.nan)
        for k in range(5):
            m = spfold == k
            out[m] = pct(v[m])
        return out

    z_titan, z_clin, z_sp = fold_pct(p_titan), fold_pct(p_clin), fold_pct(spr)

    ii, jj = cpairs(t, e)
    arms = {
        "clinical alone (the floor any multimodal claim must clear)": cidx(z_clin, ii, jj),
        "TITAN alone": cidx(z_titan, ii, jj),
        "SurvPath alone (WSI + transcriptomics)": cidx(z_sp, ii, jj),
        "OURS = TITAN + clinical, rank-averaged": cidx(fold_pct(z_titan + z_clin), ii, jj),
        "SURVPATH + clinical, same rank-average": cidx(fold_pct(z_sp + z_clin), ii, jj),
        "SurvPath + TITAN (no clinical)": cidx(fold_pct(z_sp + z_titan), ii, jj),
        "all three": cidx(fold_pct(z_sp + z_titan + z_clin), ii, jj),
    }

    rep = {"artifact_type": "attribution_clinical_vs_representation", "reportable": False,
           "phase_of_origin": "pre_h1", "n": len(keep), "events": int(e.sum()),
           "arms": {k: round(v, 4) for k, v in arms.items()}}

    ours = fold_pct(z_titan + z_clin)
    spc = fold_pct(z_sp + z_clin)
    rng = np.random.default_rng(0)
    boots = []
    n = len(keep)
    for _ in range(4000):
        bs = rng.choice(n, size=n, replace=True)
        i2, j2 = cpairs(t[bs], e[bs])
        if i2.size == 0:
            continue
        boots.append(cidx(ours[bs], i2, j2) - cidx(spc[bs], i2, j2))
    b = np.asarray(boots)
    rep["ours_minus_survpath_plus_clinical"] = {
        "mean": round(float(b.mean()), 4), "se": round(float(b.std(ddof=1)), 4),
        "ci95": [round(float(np.quantile(b, 0.025)), 4), round(float(np.quantile(b, 0.975)), 4)],
        "reading": ("if this interval contains zero, the advantage over SurvPath is explained by the "
                    "clinical variables rather than by the TITAN representation, and the claim must "
                    "be rewritten accordingly"),
    }
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
