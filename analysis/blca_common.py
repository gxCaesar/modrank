#!/usr/bin/env python3
"""Shared loading and survival primitives for the TCGA-BLCA S5 loop.

Pre-freeze, exploratory. Nothing produced through this module is reportable.

Everything here was already written and exercised inside `paired_head_to_head.py` and
`attribution.py`; it is factored out so the headroom probe, the error atlas and every candidate
run share ONE definition of the folds, the labels, the clinical design matrix and the Cox fit.
Two copies of a Cox solver is how a candidate wins by using a different alpha grid than the
baseline it is compared against.
"""

from __future__ import annotations

import collections
import csv
import glob
import pickle

import numpy as np

ALPHAS = (1.0, 8.0, 64.0, 512.0, 4096.0)


# --------------------------------------------------------------------------- survival primitives

def cpairs(time, event):
    """Index pairs (i, j) with i's event strictly before j's last known time."""
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
    return cidx(risk, *cpairs(t, e))


def cox_fit(X, t, e, alpha):
    """Breslow-tie ridge Cox by L-BFGS on the exact partial likelihood and its gradient."""
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


def fitapply(Xtr, ttr, etr, Xte, seed=0, alphas=ALPHAS):
    """Ridge Cox with alpha chosen by 3-fold inner CV inside the training set only."""
    best, ba = -1.0, alphas[0]
    inner = np.array_split(np.random.default_rng(seed).permutation(len(ttr)), 3)
    for al in alphas:
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
    """Rank to [0, 1]. Applied within a fold before pooling, for every arm equally."""
    return np.argsort(np.argsort(v)).astype(float) / max(1.0, len(v) - 1.0)


# --------------------------------------------------------------------------- the cohort

class Cohort:
    """The 359 cases the released SurvPath folds define, with every arm's inputs aligned to them.

    `keep` is the case order every array in this object is indexed by. SurvPath's rerun supplies
    the labels, because those are the labels its own evaluation used; taking them from anywhere
    else would make the head-to-head compare two different endpoints.
    """

    def __init__(self, root, titan_pkl, sp_dir):
        import json

        sp = {}
        for k, f in enumerate(sorted(glob.glob(sp_dir + "/split_*_results.pkl"))):
            for c, v in pickle.load(open(f, "rb")).items():
                sp[c] = {"risk": float(v["risk"]), "time": float(v["time"]),
                         "event": 1.0 if float(v["censorship"]) == 0.0 else 0.0, "fold": k}

        meta = list(csv.DictReader(open(root + "/tcga_blca_meta.csv")))
        stem2case = {r["slide_id"].replace(".svs", ""): r["case_id"] for r in meta}
        lab = {}
        for r in meta:
            lab.setdefault(r["case_id"], r)

        self.folds = []
        for k in range(5):
            tr, va = set(), set()
            for r in csv.DictReader(open(root + "/splits/splits_%d.csv" % k)):
                if r.get("train", "").strip():
                    tr.add(r["train"].strip())
                if r.get("val", "").strip():
                    va.add(r["val"].strip())
            self.folds.append((tr, va))

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

        self.keep = keep
        self.idx = {c: i for i, c in enumerate(keep)}
        self.n_slides = np.array([len(per[c]) for c in keep])
        self.T = np.vstack([E[per[c]].mean(0) for c in keep]).astype(float)
        self.t = np.array([sp[c]["time"] for c in keep])
        self.e = np.array([sp[c]["event"] for c in keep])
        self.spr = np.array([sp[c]["risk"] for c in keep])
        self.fold = np.array([sp[c]["fold"] for c in keep])
        self.tnm, self.lab = tnm, lab

        # fold membership must be identical to SurvPath's, or nothing here is paired
        mine = np.full(len(keep), -1)
        for k, (_, va) in enumerate(self.folds):
            for c in va:
                if c in self.idx:
                    mine[self.idx[c]] = k
        assert (mine == self.fold).all(), "fold assignment differs from SurvPath's"

        self.age = np.array([float(lab[c]["age"]) if lab[c]["age"].strip() else 65.0 for c in keep])
        self.fem = np.array([float(lab[c]["is_female"]) if lab[c]["is_female"].strip() else 0.0
                             for c in keep])
        self.site = np.array([c.split("-")[1] if len(c.split("-")) > 1 else "??" for c in keep])
        self.CLIN = self._clinical()

    def _g(self, c, k):
        x = (self.tnm.get(c) or {}).get(k)
        return "" if x in (None, "NX", "MX", "TX") else str(x)

    def _clinical(self):
        def oh(v):
            lv = sorted({x for x in v if x})
            return (np.array([[1.0 if x == l else 0.0 for l in lv] for x in v])
                    if lv else np.zeros((len(v), 0)))
        return np.hstack([
            self.age[:, None], self.fem[:, None],
            oh([self._g(c, "ajcc_pathologic_t")[:2] for c in self.keep]),
            oh([self._g(c, "ajcc_pathologic_n") for c in self.keep]),
            oh([self._g(c, "ajcc_pathologic_m") for c in self.keep]),
            oh([str((self.tnm.get(c) or {}).get("ajcc_pathologic_stage") or "")
                for c in self.keep]),
        ])

    def fold_indices(self, folds=None):
        """(train_idx, val_idx) per fold, restricted to cases this cohort actually holds."""
        out = []
        for tr, va in (folds if folds is not None else self.folds):
            out.append((np.array(sorted(self.idx[c] for c in tr if c in self.idx)),
                        np.array(sorted(self.idx[c] for c in va if c in self.idx))))
        return out

    def fold_pct(self, v, assign=None):
        """Percentile-normalise within each fold. Both arms get this or the comparison is scale."""
        a = self.fold if assign is None else assign
        out = np.full(len(v), np.nan)
        for k in np.unique(a):
            m = a == k
            out[m] = pct(v[m])
        return out
