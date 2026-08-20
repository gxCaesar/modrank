#!/usr/bin/env python3
"""C1 -- a pan-cancer survival subspace as a prior on the bladder image head.

Pre-freeze, exploratory. NOT reportable.

THE DEFICIT THIS AIMS AT, measured today and not assumed. At 16 principal components the TITAN
representation reaches 0.7666 when the Cox head is fitted IN-SAMPLE on the held-out fold, 0.6496
when it is fitted out-of-fold, and 0.4980 when the same in-sample fit is given permuted outcomes
(1,000 draws). So the representation holds 0.2686 of recoverable structure above its own control
and gives up 0.1170 of it to ESTIMATION at 287 training patients and 113 events. The clinical
block, at a fraction of the dimension, gives up 0.0819. Estimation is what caps the image arm.

MECHANISM. Locating a survival direction inside a 768-dimensional embedding space from 113 events
is the whole difficulty. If the survival-relevant part of that space is partly shared across
cancer types, then directions estimated where events are plentiful can cut the bladder head's
effective dimension from 768 to a handful, and the estimation loss falls with it. Today's
falsifier says the sharing is real: a Cox direction fitted on COLORECTAL slides, with zero case
overlap, ranks bladder patients at 0.6492 against a permuted-label control of 0.4727 +/- 0.032
(max over 20 draws 0.5235), and the mean of four donor directions reaches 0.6211 against 0.4923.

THE COMPONENT. Estimate a k-dimensional subspace from bootstrap Cox fits on the four non-bladder
cohorts SurvPath ships (BRCA, COADREAD, HNSC, STAD -- 1,882 cases, 298 events, all with TITAN
embeddings already on disk), then fit the bladder head under an anisotropic ridge that penalises
off-subspace directions more than in-subspace ones. The prior strength s and the dimension k are
chosen by inner cross-validation INSIDE each training fold; s = 1 reproduces the plain ridge
exactly, so the component contains its own null.

FALSIFIER, FIXED BEFORE THE RUN. The prior arm must beat the plain ridge image arm (0.6596 on the
released folds) by more than the paired noise bar measured on this code today -- 2 SD of the
reseeded paired delta, 0.0145. Anything smaller is what a reseed would have produced anyway.

REMOVAL ABLATION. `s = 1`, which is the same code path with the prior switched off.

THREE MORE ARMS, because the obvious ways for this to be right for the wrong reason each have one:

  unsupervised_subspace   the top-k principal components of the bladder TRAINING embeddings, same
                          k, same anisotropic ridge. If this matches, the gain is dimension
                          reduction and nothing to do with donor survival labels.
  own_cohort_subspace     the same bootstrap-Cox subspace construction fitted on the bladder
                          TRAINING fold instead of the donors. If this matches, supervision helps
                          but it does not have to come from other cancers.
  permuted_donor_subspace THE MATCHED CONTROL. Identical donors, identical bootstraps, identical
                          k and s, with each donor's survival labels permuted before fitting. Same
                          dimensionality, same estimation noise, no survival signal. A gain here
                          voids the run.

AND A SECOND CONTROL FOR A DIFFERENT SHORTCUT. A permuted-label control cannot rule out a donor
signal that is mediated by a technical factor -- scanner, stain, magnification, submitting site --
which also happens to track survival in bladder. So the transfer is additionally scored WITHIN
bladder collection site: only pairs of bladder patients from the same submitting site contribute.
If the transfer is a batch artifact it collapses there, and if it survives it is not one.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cindex, cpairs, cox_fit  # noqa: E402

DONORS = ("brca", "coadread", "hnsc", "stad")
ALPHAS = (1.0, 8.0, 64.0, 512.0, 4096.0)


def load_donor(path, stem2emb):
    per, lab = collections.defaultdict(list), {}
    for r in csv.DictReader(open(path)):
        stem = r["slide_id"].replace(".svs", "")
        if stem in stem2emb:
            per[r["case_id"]].append(stem2emb[stem])
        if r.get("survival_months_dss", "").strip() and r.get("censorship_dss", "").strip():
            lab[r["case_id"]] = (float(r["survival_months_dss"]),
                                 1.0 - float(r["censorship_dss"]))
    cases = sorted(c for c in per if c in lab)
    if not cases:
        return None
    X = np.vstack([np.mean(per[c], axis=0) for c in cases]).astype(float)
    return (X, np.array([lab[c][0] for c in cases]),
            np.array([lab[c][1] for c in cases]), cases)


def bootstrap_directions(X, t, e, n_boot, rng, alpha=64.0, permute=False):
    """Unit ridge-Cox coefficient vectors over bootstrap resamples of one cohort."""
    mu, sd = X.mean(0), X.std(0) + 1e-9
    Z = (X - mu) / sd
    out = []
    for _ in range(n_boot):
        bs = rng.choice(len(t), size=len(t), replace=True)
        tt, ee = t[bs], e[bs]
        if permute:
            q = rng.permutation(len(tt))
            tt, ee = tt[q], ee[q]
        if ee.sum() < 5:
            continue
        b = cox_fit(Z[bs], tt, ee, alpha)
        n = float(np.linalg.norm(b))
        if n > 0:
            out.append(b / n)
    return np.vstack(out) if out else np.zeros((0, X.shape[1]))


def subspace(dirs, k):
    """Top-k right singular directions of a stack of coefficient vectors."""
    if dirs.shape[0] == 0:
        return np.zeros((0, dirs.shape[1]))
    _, _, Vt = np.linalg.svd(dirs, full_matrices=False)
    return Vt[:min(k, Vt.shape[0])]


def aniso_fitapply(Xtr, ttr, etr, Xte, V, s, seed=0, alpha=None):
    """Ridge Cox with off-subspace directions penalised harder than in-subspace ones.

    Implemented by amplifying the in-subspace coordinates before an ordinary ridge:
    W = I + (s-1) V^T V, then plain ridge on X W. At s = 1 the transform is the identity and this
    is exactly the baseline estimator, which is why the removal ablation is a parameter value
    rather than a second code path.

    `alpha=None` selects the ridge strength by inner CV, which is what the FINAL refit does.
    Passing an explicit alpha skips that, and the (k, s) search below uses it: a full alpha search
    nested inside a (k, s) search inside a fold is a triple nesting that costs 16x for a choice the
    baseline has already made, and the arms are compared to each other at matched alpha rather than
    each being allowed its own -- which is the fairer comparison anyway.
    """
    d = Xtr.shape[1]
    if V.shape[0] == 0 or s == 1.0:
        Wtr, Wte = Xtr, Xte
    else:
        P = V.T @ V
        W = np.eye(d) + (s - 1.0) * P
        Wtr, Wte = Xtr @ W, Xte @ W
    grid = ALPHAS if alpha is None else (alpha,)
    if len(grid) == 1:
        ba = grid[0]
    else:
        best, ba = -1.0, grid[0]
        inner = np.array_split(np.random.default_rng(seed).permutation(len(ttr)), 3)
        for al in grid:
            ip = np.zeros(len(ttr))
            for ite in inner:
                itr = np.setdiff1d(np.arange(len(ttr)), ite)
                mu, sd = Wtr[itr].mean(0), Wtr[itr].std(0) + 1e-9
                ip[ite] = ((Wtr[ite] - mu) / sd) @ cox_fit((Wtr[itr] - mu) / sd,
                                                           ttr[itr], etr[itr], al)
            c = cindex(ip, ttr, etr)
            if c == c and c > best:
                best, ba = c, al
    mu, sd = Wtr.mean(0), Wtr.std(0) + 1e-9
    return ((Wte - mu) / sd) @ cox_fit((Wtr - mu) / sd, ttr, etr, ba), ba


def pick_alpha(Xtr, ttr, etr, seed=0):
    """The plain-ridge alpha for this training fold; every (k, s) candidate is searched at it."""
    best, ba = -1.0, ALPHAS[0]
    inner = np.array_split(np.random.default_rng(seed).permutation(len(ttr)), 3)
    for al in ALPHAS:
        ip = np.zeros(len(ttr))
        for ite in inner:
            itr = np.setdiff1d(np.arange(len(ttr)), ite)
            mu, sd = Xtr[itr].mean(0), Xtr[itr].std(0) + 1e-9
            ip[ite] = ((Xtr[ite] - mu) / sd) @ cox_fit((Xtr[itr] - mu) / sd,
                                                       ttr[itr], etr[itr], al)
        c = cindex(ip, ttr, etr)
        if c == c and c > best:
            best, ba = c, al
    return ba


def inner_select_and_fit(co, tri, vai, basis_fn, ks, ss, seed=0):
    """Choose (k, s) by inner CV inside the TRAINING fold, then refit on the whole training fold.

    `basis_fn(train_index_subset)` returns the stack of direction vectors the subspace is built
    from. For the donor arms it ignores its argument -- no bladder patient takes part -- and for
    the bladder-own arm it must see only the inner-training patients, which is why it is a callback
    rather than a precomputed matrix.
    """
    inner = list(np.array_split(np.random.default_rng(seed + 3).permutation(len(tri)), 3))
    al = pick_alpha(co.T[tri], co.t[tri], co.e[tri], seed)
    bases = [basis_fn(tri[np.setdiff1d(np.arange(len(tri)), ite)]) for ite in inner]
    best, bk, bs_ = -1.0, ks[0], ss[0]
    for k in ks:
        for s in ss:
            ip = np.zeros(len(tri))
            for ite, B in zip(inner, bases):
                itr = tri[np.setdiff1d(np.arange(len(tri)), ite)]
                ip[ite] = aniso_fitapply(co.T[itr], co.t[itr], co.e[itr], co.T[tri[ite]],
                                         subspace(B, k), s, seed, alpha=al)[0]
            c = cindex(ip, co.t[tri], co.e[tri])
            if c == c and c > best:
                best, bk, bs_ = c, k, s
    V = subspace(basis_fn(tri), bk)
    sc, _ = aniso_fitapply(co.T[tri], co.t[tri], co.e[tri], co.T[vai], V, bs_, seed)
    return sc, bk, bs_


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--pan-dir", required=True)
    ap.add_argument("--boots", type=int, default=20)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    d = pickle.load(open(a.titan, "rb"))
    E = np.asarray(d["embeddings"])
    stem2emb = {}
    for i, f in enumerate(d["filenames"]):
        s = str(f)
        for suf in (".svs", ".h5", ".pt"):
            if s.endswith(suf):
                s = s[:-len(suf)]
        stem2emb[s] = E[i]

    rng = np.random.default_rng(0)
    donor_real, donor_perm = [], []
    for nm in DONORS:
        got = load_donor(os.path.join(a.pan_dir, "tcga_%s.csv" % nm), stem2emb)
        if got is None:
            continue
        X, t, e, _ = got
        donor_real.append(bootstrap_directions(X, t, e, a.boots, np.random.default_rng(11), 64.0))
        donor_perm.append(bootstrap_directions(X, t, e, a.boots, np.random.default_rng(12), 64.0,
                                               permute=True))
        print("donor %s: %d real dirs, %d permuted" % (nm, len(donor_real[-1]),
                                                       len(donor_perm[-1])),
              file=sys.stderr, flush=True)
    D_real = np.vstack(donor_real)
    D_perm = np.vstack(donor_perm)

    KS, SS = (2, 4, 8, 16, 32), (1.0, 3.0, 10.0, 30.0)
    arms = {
        "plain_ridge_baseline": (lambda idx: np.zeros((0, co.T.shape[1])), (0,), (1.0,)),
        "pancancer_prior": (lambda idx: D_real, KS, SS),
        "permuted_donor_subspace_CONTROL": (lambda idx: D_perm, KS, SS),
        "unsupervised_subspace_ABLATION": (
            lambda idx: np.linalg.svd((co.T[idx] - co.T[idx].mean(0))
                                      / (co.T[idx].std(0) + 1e-9),
                                      full_matrices=False)[2][:32], KS, SS),
        "own_cohort_subspace_ABLATION": (
            lambda idx: bootstrap_directions(co.T[idx], co.t[idx], co.e[idx], 20,
                                             np.random.default_rng(13), 64.0), KS, SS),
    }

    fi = co.fold_indices()
    ii, jj = cpairs(co.t, co.e)
    same_site = co.site[ii] == co.site[jj]

    rep = {"artifact_type": "s5_c1_pancancer_prior", "reportable": False,
           "phase_of_origin": "exploratory",
           "n": len(co.keep), "events": int(co.e.sum()),
           "donor_directions": {"real": int(D_real.shape[0]), "permuted": int(D_perm.shape[0])},
           "grid": {"k": list(KS), "s": list(SS)},
           "falsifier": "the pancancer_prior arm must beat plain_ridge_baseline by more than "
                        "0.0145, the 2-SD paired reseed bar measured on this code today",
           "within_site_pairs": int(same_site.sum()), "all_pairs": int(ii.size)}

    out = {}
    for nm, (bf, ks, ss) in arms.items():
        s = np.full(len(co.keep), np.nan)
        picks = []
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            if nm == "plain_ridge_baseline":
                s[vai] = aniso_fitapply(co.T[tri], co.t[tri], co.e[tri], co.T[vai],
                                        np.zeros((0, co.T.shape[1])), 1.0, 0)[0]
                picks.append((0, 1.0))
            else:
                s[vai], k, sv = inner_select_and_fit(co, tri, vai, bf, ks, ss, 0)
                picks.append((k, sv))
        z = co.fold_pct(s)
        out[nm] = {"cindex": round(cidx(z, ii, jj), 4),
                   "cindex_within_site_pairs": round(cidx(z, ii[same_site], jj[same_site]), 4),
                   "selected_k_s_per_fold": [[int(k), float(v)] for k, v in picks]}
        print("%-34s C=%.4f  within-site=%.4f  picks=%s"
              % (nm, out[nm]["cindex"], out[nm]["cindex_within_site_pairs"], picks),
              file=sys.stderr, flush=True)
    rep["arms"] = out

    base = out["plain_ridge_baseline"]["cindex"]
    rep["verdict"] = {
        "prior_minus_baseline": round(out["pancancer_prior"]["cindex"] - base, 4),
        "control_minus_baseline": round(out["permuted_donor_subspace_CONTROL"]["cindex"] - base, 4),
        "unsupervised_minus_baseline": round(
            out["unsupervised_subspace_ABLATION"]["cindex"] - base, 4),
        "own_cohort_minus_baseline": round(
            out["own_cohort_subspace_ABLATION"]["cindex"] - base, 4),
        "bar": 0.0145,
        "read_the_control_first": "if permuted_donor_subspace_CONTROL moved as much as "
                                  "pancancer_prior did, this run is VOID and nothing else in it "
                                  "is reportable",
    }

    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
