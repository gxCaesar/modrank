#!/usr/bin/env python3
"""C2-C5 -- attack the ESTIMATION loss directly, with a mechanism prediction that can be wrong.

Pre-freeze, exploratory. NOT reportable.

WHERE THIS COMES FROM. The oracle probe put the image arm's deficit at 0.1170 of C-index lost
between an in-sample fit at matched capacity (0.7666, permuted-outcome control 0.4980 over 1,000
draws) and its out-of-fold twin (0.6496). Two attempts to import that missing strength from other
cancers have since been killed on their own falsifiers -- a donor subspace constraint at -0.0042
against a bar of +0.0145, and the donor direction added as a separate arm at +0.0033 -- with the
matched controls behaving correctly in both. Under the slate's shape rule those two share one
shape, `use non-bladder supervision to improve the bladder image score`, and a third of that shape
is not funded.

So the remaining move is the estimator itself, and these four candidates are four different
objects:

  C2  bagged ridge Cox        resample the ESTIMATOR       average B ridge fits over bootstrap
                              resamples of the training fold. Estimation variance is exactly what
                              bagging reduces, and the deficit was measured as variance.
  C3  elastic net Cox         change the PENALTY GEOMETRY  an L1 component selects a sparse set of
                              embedding coordinates instead of shrinking all 768 equally.
  C4  dense regularisation    fix the SELECTION NOISE      the baseline picks among five alphas
      path                    spanning 1 to 4096 by one 3-fold inner CV on 287 patients. That
                              choice is itself noisy, and a coarse grid plus a noisy criterion is
                              a way to lose performance that looks like a modelling limit.
  C5  gated fusion            change the COMBINATION        the image weight is allowed to differ
                              across clinical risk strata, selected in-fold. The oracle for this
                              was +0.0224 with its control at -0.0029.

THE MECHANISM PREDICTION, WRITTEN BEFORE THE RUN, AND IT IS WHAT MAKES THIS MORE THAN A SWEEP.
If the deficit really is estimation variance in a high-dimensional coefficient, then C2-C4 must
help the arms IN PROPORTION TO THEIR DIMENSION:

    TITAN     768 features, measured estimation loss 0.1170   -> largest gain
    omics     275 pathway features                            -> intermediate gain
    clinical   ~15 features, measured estimation loss 0.0819  -> little or no gain

A gain that is UNIFORM across the three arms falsifies the mechanism even if the numbers go up,
because something that helps a 15-dimensional Cox as much as a 768-dimensional one is not fixing
high-dimensional estimation variance. This is the matched control for C2-C4 and it is built into
the design rather than bolted on: the clinical arm is the condition that must NOT move.

C5 carries its own control instead, since it is a different object: the same gating structure
driven by a PERMUTED stratum label, which preserves the number of strata and their sizes and
destroys only which patient is in which.

FALSIFIER for every candidate: beat its own baseline arm by more than 0.0145, the 2-SD paired
reseed bar measured on this code today.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cindex, cpairs, cox_fit, fitapply, pct  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402

COARSE = (1.0, 8.0, 64.0, 512.0, 4096.0)
DENSE = tuple(float(x) for x in np.round(np.geomspace(0.25, 16384.0, 25), 4))


def _std(Xtr, Xte):
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
    return (Xtr - mu) / sd, (Xte - mu) / sd


def ridge_arm(Xtr, ttr, etr, Xte, alphas, folds=3, seed=0):
    """Ridge Cox, alpha by inner CV. `folds` and `alphas` are the two knobs C4 turns."""
    best, ba = -1.0, alphas[0]
    inner = np.array_split(np.random.default_rng(seed).permutation(len(ttr)), folds)
    for al in alphas:
        ip = np.zeros(len(ttr))
        for ite in inner:
            itr = np.setdiff1d(np.arange(len(ttr)), ite)
            a, b = _std(Xtr[itr], Xtr[ite])
            ip[ite] = b @ cox_fit(a, ttr[itr], etr[itr], al)
        c = cindex(ip, ttr, etr)
        if c == c and c > best:
            best, ba = c, al
    a, b = _std(Xtr, Xte)
    return b @ cox_fit(a, ttr, etr, ba)


def bagged_arm(Xtr, ttr, etr, Xte, alphas, n_bag=40, seed=0):
    """C2. Average the RANKS of B ridge fits on bootstrap resamples of the training fold.

    Ranks rather than raw log-hazards, because the resampled fits sit on different scales and
    averaging them directly lets whichever bootstrap happened to fit the largest coefficients
    dominate the ensemble.

    THE RIDGE STRENGTH IS CHOSEN ONCE, on the full training fold, and reused inside every
    bootstrap. A first version re-ran the whole alpha search inside each of the 40 resamples, which
    is 16x the work for a hyperparameter the fold has already decided, and it made the in-fold
    estimator-selection loop cost about four hours. It is also the statistically cleaner choice:
    bagging is variance reduction over the COEFFICIENT estimate, and letting alpha jump between
    resamples averages over a different estimator each time.
    """
    rng = np.random.default_rng(seed + 77)
    al = alphas[0] if len(alphas) == 1 else _pick_alpha(Xtr, ttr, etr, alphas, seed)
    acc = np.zeros(len(Xte))
    used = 0
    for _ in range(n_bag):
        bs = rng.choice(len(ttr), size=len(ttr), replace=True)
        if etr[bs].sum() < 8:
            continue
        a, b = _std(Xtr[bs], Xte)
        acc += np.argsort(np.argsort(b @ cox_fit(a, ttr[bs], etr[bs], al))).astype(float)
        used += 1
    return acc / max(used, 1)


def _pick_alpha(Xtr, ttr, etr, alphas, seed=0, folds=3):
    best, ba = -1.0, alphas[0]
    inner = np.array_split(np.random.default_rng(seed).permutation(len(ttr)), folds)
    for al in alphas:
        ip = np.zeros(len(ttr))
        for ite in inner:
            itr = np.setdiff1d(np.arange(len(ttr)), ite)
            a, b = _std(Xtr[itr], Xtr[ite])
            ip[ite] = b @ cox_fit(a, ttr[itr], etr[itr], al)
        c = cindex(ip, ttr, etr)
        if c == c and c > best:
            best, ba = c, al
    return ba


def enet_arm(Xtr, ttr, etr, Xte, seed=0):
    """C3. Elastic-net Cox by proximal gradient, l1_ratio and lambda by inner CV."""
    def fit(X, t, e, lam, l1r, iters=300):
        o = np.argsort(t)
        Xo, eo = X[o], e[o]
        ev = np.flatnonzero(eo == 1)
        if ev.size < 2:
            return np.zeros(X.shape[1])
        b = np.zeros(X.shape[1])
        L = float(np.linalg.norm(Xo, 2) ** 2) / max(len(t), 1) + 1e-9
        step = 1.0 / (L + lam * (1 - l1r))
        for _ in range(iters):
            eta = Xo @ b
            m = float(eta.max())
            ex = np.exp(eta - m)
            rs = np.cumsum(ex[::-1])[::-1]
            rx = np.cumsum((Xo * ex[:, None])[::-1], axis=0)[::-1]
            g = -(Xo[ev] - rx[ev] / rs[ev, None]).sum(0) + lam * (1 - l1r) * b
            z = b - step * g
            thr = step * lam * l1r
            b = np.sign(z) * np.maximum(np.abs(z) - thr, 0.0)
        return b

    grid = [(lam, l1r) for lam in (0.5, 2.0, 8.0, 32.0, 128.0) for l1r in (0.1, 0.5, 0.9)]
    best, bp = -1.0, grid[0]
    inner = np.array_split(np.random.default_rng(seed).permutation(len(ttr)), 3)
    for lam, l1r in grid:
        ip = np.zeros(len(ttr))
        for ite in inner:
            itr = np.setdiff1d(np.arange(len(ttr)), ite)
            a, b_ = _std(Xtr[itr], Xtr[ite])
            ip[ite] = b_ @ fit(a, ttr[itr], etr[itr], lam, l1r)
        c = cindex(ip, ttr, etr)
        if c == c and c > best:
            best, bp = c, (lam, l1r)
    a, b_ = _std(Xtr, Xte)
    return b_ @ fit(a, ttr, etr, bp[0], bp[1])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--rna", required=True)
    ap.add_argument("--sig", required=True)
    ap.add_argument("--bag", type=int, default=40)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    G, genes, _ = load_rna(a.rna, co.keep)
    P, pnames = pathway_matrix(G, genes, a.sig)
    blocks = {"titan_768d": co.T, "omics_%dd" % P.shape[1]: P,
              "clinical_%dd" % co.CLIN.shape[1]: co.CLIN}
    ii, jj = cpairs(co.t, co.e)

    estimators = {
        "baseline_ridge_coarse5_inner3": lambda Xtr, t, e, Xte: ridge_arm(Xtr, t, e, Xte, COARSE,
                                                                          3, 0),
        "C2_bagged_ridge": lambda Xtr, t, e, Xte: bagged_arm(Xtr, t, e, Xte, COARSE, a.bag, 0),
        "C3_elastic_net": lambda Xtr, t, e, Xte: enet_arm(Xtr, t, e, Xte, 0),
        "C4_dense25_inner5": lambda Xtr, t, e, Xte: ridge_arm(Xtr, t, e, Xte, DENSE, 5, 0),
    }

    rep = {"artifact_type": "s5_c2_c5_estimator_slate", "reportable": False,
           "phase_of_origin": "exploratory", "n": len(co.keep), "events": int(co.e.sum()),
           "bar": 0.0145, "bag_size": a.bag,
           "mechanism_prediction": "gains must scale with block dimension: titan(768) > "
                                   "omics(%d) > clinical(%d). A uniform gain falsifies the "
                                   "estimation-variance mechanism." % (P.shape[1],
                                                                       co.CLIN.shape[1])}

    scores, table = {}, {}
    for bname, X in blocks.items():
        table[bname] = {}
        for ename, fn in estimators.items():
            s = np.full(len(co.keep), np.nan)
            for tri, vai in fi:
                if len(tri) < 30 or not len(vai):
                    continue
                s[vai] = fn(X[tri], co.t[tri], co.e[tri], X[vai])
            z = co.fold_pct(s)
            scores[(bname, ename)] = z
            table[bname][ename] = round(cidx(z, ii, jj), 4)
            print("%-16s %-32s %.4f" % (bname, ename, table[bname][ename]),
                  file=sys.stderr, flush=True)
        b0 = table[bname]["baseline_ridge_coarse5_inner3"]
        table[bname]["_gains_vs_baseline"] = {
            k: round(v - b0, 4) for k, v in table[bname].items()
            if k != "baseline_ridge_coarse5_inner3" and not k.startswith("_")}
    rep["per_block"] = table

    # ---- C5: gated fusion, strata and gating chosen INSIDE the training fold
    def best_estimator_for(bname):
        d = {k: v for k, v in table[bname].items() if not k.startswith("_")}
        return max(d, key=lambda k: d[k])

    picks = {b: best_estimator_for(b) for b in blocks}
    rep["estimator_picked_per_block_on_pooled_outcome"] = picks
    rep["pick_caveat"] = ("four estimators were scored per block on the pooled outcome; taking the "
                          "best carries a selection term of about sd*sqrt(2*ln 4). The gated arm "
                          "below inherits it and is a reference, not a reportable number.")

    zt = scores[("titan_768d", picks["titan_768d"])]
    zo = scores[[k for k in blocks if k.startswith("omics")][0],
                picks[[k for k in blocks if k.startswith("omics")][0]]]
    zc = scores[[k for k in blocks if k.startswith("clinical")][0],
                picks[[k for k in blocks if k.startswith("clinical")][0]]]

    def gated(zimg, zom, zcl, strata):
        out = np.zeros(len(zcl))
        for g in np.unique(strata):
            m = strata == g
            if m.sum() < 10:
                out[m] = zimg[m] + zom[m] + zcl[m]
                continue
            wi = cidx(zimg[m], *cpairs(co.t[m], co.e[m]))
            wo = cidx(zom[m], *cpairs(co.t[m], co.e[m]))
            wi = max(wi - 0.5, 0.0)
            wo = max(wo - 0.5, 0.0)
            tot = wi + wo + 0.25
            out[m] = (wi * zimg[m] + wo * zom[m] + 0.25 * zcl[m]) / max(tot, 1e-9)
        return out

    q = np.quantile(zc, [1 / 3, 2 / 3])
    tert = np.digitize(zc, q)
    rngp = np.random.default_rng(5)
    tert_perm = tert[rngp.permutation(len(tert))]
    flat = co.fold_pct(zt + zo + zc)
    rep["C5_gated_fusion"] = {
        "flat_rank_average": round(cidx(flat, ii, jj), 4),
        "gated_on_clinical_tertile": round(cidx(co.fold_pct(gated(zt, zo, zc, tert)), ii, jj), 4),
        "gated_on_PERMUTED_tertile_CONTROL": round(
            cidx(co.fold_pct(gated(zt, zo, zc, tert_perm)), ii, jj), 4),
        "note": "the gating weights here are read off the SAME pooled outcome they are scored on, "
                "so this is an oracle-flavoured reference for whether gating is worth building "
                "in-fold, not an arm. Its control says how much the structure alone buys.",
    }

    rep["best_flat_combination"] = round(cidx(flat, ii, jj), 4)
    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
