#!/usr/bin/env python3
"""The controlled head-to-head an audit said was missing.

Pre-H1. NOT reportable.

WHAT WAS WRONG BEFORE. Our number could only be set against SurvPath's PUBLISHED MEAN of 0.625,
because the repository ships the folds but no per-fold or per-case outputs. An adversarial audit
put it plainly: "Aggregate mean comparisons cannot support paired superiority, and cross-paper
numbers should not be described as a controlled head-to-head result."

WHAT CHANGED. SurvPath was rerun on the released `splits/5foldcv/tcga_blca` case-ID folds with the
768-d CTransPath-family features its own paper specifies, reproducing 0.6147 +/- 0.1067 against the
published 0.625 +/- 0.056. Its per-case risk scores were then recovered from split_{k}_results.pkl,
and the per-fold C-index recomputed from those scores matches the run's own log to within 0.0005 on
every fold, which is what licenses using them.

So both models can now be scored on THE SAME 359 patients, in THE SAME held-out folds, and the
difference has a paired uncertainty instead of a difference of two published averages.

DISCIPLINE.
  * Our arm uses in-fold rule selection: the choice among {titan, clinical, concat, rank-average}
    is made by inner cross-validation inside each training fold, never after seeing the outer fold.
  * Predictions are percentile-normalised WITHIN each fold before pooling, for both arms equally.
    Without this, a Cox log-hazard on a fold-specific scale is pooled against a rank on [0,1] and
    the comparison measures scale, not skill.
  * The bootstrap resamples CASES, not comparable pairs: pairs from one patient are not independent
    and resampling them would understate the variance.
  * Fold membership is asserted identical between the two arms before anything is computed.
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

    # ---- SurvPath per-case risks and their fold assignment
    sp = {}
    for k, f in enumerate(sorted(glob.glob(sp_dir + "/split_*_results.pkl"))):
        d = pickle.load(open(f, "rb"))
        for c in d:
            sp[c] = {"risk": float(d[c]["risk"]), "time": float(d[c]["time"]),
                     "event": 1.0 if float(d[c]["censorship"]) == 0.0 else 0.0, "fold": k}

    # ---- our side, on the same released folds
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

    # fold membership must be identical, or this is not a paired comparison
    myfold = np.full(len(keep), -1)
    for k, (_, va) in enumerate(folds):
        for c in va:
            if c in idx:
                myfold[idx[c]] = k
    assert (myfold == spfold).all(), "fold assignment differs between the two arms"

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

    CAND = ("titan", "clinical", "concat", "rank_avg")

    def score(kind, tri, vai, seed):
        if kind == "titan":
            return fitapply(T[tri], t[tri], e[tri], T[vai], seed)
        if kind == "clinical":
            return fitapply(CLIN[tri], t[tri], e[tri], CLIN[vai], seed)
        if kind == "concat":
            Z = np.hstack([T, CLIN])
            return fitapply(Z[tri], t[tri], e[tri], Z[vai], seed)
        a = fitapply(T[tri], t[tri], e[tri], T[vai], seed)
        b = fitapply(CLIN[tri], t[tri], e[tri], CLIN[vai], seed)
        r = lambda v: np.argsort(np.argsort(v)).astype(float)
        return r(a) + r(b)

    ours = np.full(len(keep), np.nan)
    picks = collections.Counter()
    SEED = 0
    for k, (tr, va) in enumerate(folds):
        tri = np.array([idx[c] for c in tr if c in idx])
        vai = np.array([idx[c] for c in va if c in idx])
        if len(tri) < 30 or not len(vai):
            continue
        best, bk = -1.0, CAND[0]
        inner = np.array_split(np.random.default_rng(SEED + 11).permutation(len(tri)), 3)
        for kind in CAND:
            ip = np.zeros(len(tri))
            for ite in inner:
                itr = np.setdiff1d(np.arange(len(tri)), ite)
                ip[ite] = score(kind, tri[itr], tri[ite], SEED)
            c = cindex(ip, t[tri], e[tri])
            if c == c and c > best:
                best, bk = c, kind
        picks[bk] += 1
        ours[vai] = pct(score(bk, tri, vai, SEED))

    # theirs percentile-normalised within the same folds, so both arms are treated identically
    theirs = np.full(len(keep), np.nan)
    for k in range(5):
        m = spfold == k
        theirs[m] = pct(spr[m])

    rep = {"artifact_type": "paired_head_to_head", "reportable": False, "phase_of_origin": "pre_h1",
           "n_cases": len(keep), "n_events": int(e.sum()),
           "fold_assignment_identical": True,
           "rules_chosen_in_fold": dict(picks)}

    ii, jj = cpairs(t, e)
    c_ours, c_theirs = cidx(ours, ii, jj), cidx(theirs, ii, jj)
    rep["pooled"] = {"ours": round(c_ours, 4), "survpath_rerun": round(c_theirs, 4),
                     "difference": round(c_ours - c_theirs, 4)}

    perfold = {}
    for k in range(5):
        m = spfold == k
        i2, j2 = cpairs(t[m], e[m])
        perfold["fold_%d" % k] = {
            "n": int(m.sum()), "events": int(e[m].sum()),
            "ours": round(cidx(ours[m], i2, j2), 4),
            "survpath": round(cidx(theirs[m], i2, j2), 4),
            "diff": round(cidx(ours[m], i2, j2) - cidx(theirs[m], i2, j2), 4)}
    rep["per_fold"] = perfold
    rep["folds_we_win"] = sum(1 for v in perfold.values() if v["diff"] > 0)

    rng = np.random.default_rng(0)
    boots = []
    n = len(keep)
    for _ in range(4000):
        bs = rng.choice(n, size=n, replace=True)
        i3, j3 = cpairs(t[bs], e[bs])
        if i3.size == 0:
            continue
        boots.append(cidx(ours[bs], i3, j3) - cidx(theirs[bs], i3, j3))
    b = np.asarray(boots)
    rep["paired_bootstrap"] = {
        "mean_difference": round(float(b.mean()), 4),
        "se": round(float(b.std(ddof=1)), 4),
        "ci95": [round(float(np.quantile(b, 0.025)), 4), round(float(np.quantile(b, 0.975)), 4)],
        "p_two_sided_bootstrap": round(float(2 * min((b <= 0).mean(), (b >= 0).mean())), 4),
        "resampled": "cases, not comparable pairs",
    }
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
