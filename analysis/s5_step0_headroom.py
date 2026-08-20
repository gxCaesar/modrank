#!/usr/bin/env python3
"""S5 Step 0. Is there room above the cheap baselines, and where does the room live?

Pre-freeze, exploratory. NOT reportable.

Three questions, and the third is the one that decides what to build.

1. CHEAP BASELINES. Constant, age alone, stage alone, the full clinical Cox, TITAN alone. The
   loop's own rule: a method that does not beat these has not been shown to do anything. Clinical
   alone already sits at 0.6754 in an earlier run, which is above every published WSI method on
   this benchmark, so this is not a formality here -- it is the bar.

2. THE NOISE FLOOR, MEASURED ON THIS CODE. Not the simulation from 2026-08-14 and not a number
   quoted from a paper: the same 359 patients, the same arms, the same estimator, re-split under
   new fold assignments. What comes back is the standard deviation of the number a paper would
   report. A margin under twice that does not survive somebody else's reseed.

3. THE ORACLE, AND ITS CONTROL. Fit the SAME estimator in-sample on the validation fold at matched
   capacity (d principal components, so the fit is not degenerate), and score it there. That is an
   input the deployed system can never have, and it says how much of the outcome this
   representation can express once the estimation problem is removed.

   An oracle overshoots by overfitting, so it is paired with a structural control: the identical
   in-sample fit with the OUTCOME PERMUTED inside the fold. Whatever the permuted arm reaches is
   pure overfitting at that capacity. The excess of the real oracle over it is the recoverable
   structure, and the difference between the two readings is exactly the difference between "this
   representation is exhausted" and "the estimator cannot reach what is there".

   The two have opposite correct responses, which is why this runs before any component is
   designed:
     representation exhausted -> patch level, another modality, a different target
     estimator short          -> regularisation, capacity, a better head on the same features
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cindex, cpairs, cidx, fitapply, pct  # noqa: E402


def reseeded_folds(n, rng, k=5):
    parts = np.array_split(rng.permutation(n), k)
    return [(np.setdiff1d(np.arange(n), p), p) for p in parts]


def run_arms(co, fold_idx, seed=0):
    """Out-of-fold scores for the three arms, percentile-normalised inside each fold."""
    assign = np.full(len(co.keep), -1)
    ti = np.full(len(co.keep), np.nan)
    cl = np.full(len(co.keep), np.nan)
    for k, (tri, vai) in enumerate(fold_idx):
        if len(tri) < 30 or not len(vai):
            continue
        assign[vai] = k
        ti[vai] = fitapply(co.T[tri], co.t[tri], co.e[tri], co.T[vai], seed)
        cl[vai] = fitapply(co.CLIN[tri], co.t[tri], co.e[tri], co.CLIN[vai], seed)
    z_ti, z_cl = co.fold_pct(ti, assign), co.fold_pct(cl, assign)
    return {"titan": z_ti, "clinical": z_cl,
            "rank_avg": co.fold_pct(z_ti + z_cl, assign)}, assign


def oracle_block(co, fold_idx, dims, rng, n_perm=20):
    """In-sample fit on the held-out fold at matched capacity, against a permuted-outcome control.

    The PCA basis is fitted on the TRAINING fold, so the representation itself is not chosen with
    the validation patients; only the Cox coefficients are, which is the leak the oracle is
    deliberately allowed. The control repeats that leak with the outcome shuffled inside the same
    fold, holding n, the event count and the censoring pattern fixed.

    THE CONTROL NEEDS MANY DRAWS, and a first version of this took ONE per fold. At d=2 the fitted
    direction lives in a two-dimensional span of a genuinely prognostic block, so a permuted fit
    still points somewhere in that span: it scores far from 0.5 in each fold and its SIGN is what
    is random. Five draws then returned 0.6467 for a control that must average 0.5, which would
    have made a real oracle look like an artifact. `n_perm` draws per fold fixes it, and the
    control's own spread is reported so the reading can be checked rather than trusted.

    The three representations are kept SEPARATE and there is deliberately no concatenated arm:
    standardising 768 TITAN dimensions beside ~15 clinical ones and then taking principal
    components yields a basis the clinical block cannot influence, so a "titan+clinical" oracle
    would be a TITAN oracle wearing a label.
    """
    from blca_common import cox_fit
    out = {}
    for name, Xall in (("titan", co.T), ("clinical", co.CLIN)):
        # the PCA basis depends only on the fold, not on d -- hoisted out of the d loop
        bases = []
        for tri, vai in fold_idx:
            if len(tri) < 30 or len(vai) < 20:
                bases.append(None)
                continue
            mu, sd = Xall[tri].mean(0), Xall[tri].std(0) + 1e-9
            Z = (Xall - mu) / sd
            _, _, Vt = np.linalg.svd(Z[tri], full_matrices=False)
            bases.append((Z, Vt))
        for d in dims:
            real, perm, oof = [], [], []
            for (tri, vai), bs in zip(fold_idx, bases):
                if bs is None:
                    continue
                Z, Vt = bs
                dd = min(d, Z.shape[1], len(tri) - 1)
                P = Z @ Vt[:dd].T
                tv, ev = co.t[vai], co.e[vai]
                ii, jj = cpairs(tv, ev)
                if ii.size == 0:
                    continue
                real.append(cidx(P[vai] @ cox_fit(P[vai], tv, ev, 1.0), ii, jj))
                for _ in range(n_perm):
                    q = rng.permutation(len(vai))
                    perm.append(cidx(P[vai] @ cox_fit(P[vai], tv[q], ev[q], 1.0), ii, jj))
                oof.append(cidx(P[vai] @ cox_fit(P[tri], co.t[tri], co.e[tri], 1.0), ii, jj))
            if real:
                out["%s_d%d" % (name, d)] = {
                    "oracle_in_sample": round(float(np.mean(real)), 4),
                    "control_permuted_outcome": round(float(np.mean(perm)), 4),
                    "control_sd": round(float(np.std(perm, ddof=1)), 4),
                    "control_draws": len(perm),
                    "out_of_fold": round(float(np.mean(oof)), 4),
                    "recoverable_above_control": round(float(np.mean(real) - np.mean(perm)), 4),
                    "estimation_loss": round(float(np.mean(real) - np.mean(oof)), 4),
                }
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--reseeds", type=int, default=20)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    rng = np.random.default_rng(0)
    ii, jj = cpairs(co.t, co.e)
    rep = {"artifact_type": "s5_step0_headroom", "reportable": False,
           "phase_of_origin": "exploratory",
           "n": len(co.keep), "events": int(co.e.sum()),
           "comparable_pairs": int(ii.size)}

    # ---- 1. cheap baselines, on the released folds
    fi = co.fold_indices()
    cheap = {"constant (every patient the same risk)": cidx(np.zeros(len(co.keep)), ii, jj)}
    for nm, X in (("age alone", co.age[:, None]),
                  ("sex alone", co.fem[:, None]),
                  ("stage block alone (T,N,M,stage, no age/sex)", co.CLIN[:, 2:]),
                  ("clinical Cox (age+sex+T+N+M+stage)", co.CLIN),
                  ("TITAN slide embedding Cox", co.T)):
        s = np.full(len(co.keep), np.nan)
        assign = np.full(len(co.keep), -1)
        for k, (tri, vai) in enumerate(fi):
            if len(tri) < 30 or not len(vai):
                continue
            assign[vai] = k
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], 0)
        cheap[nm] = cidx(co.fold_pct(s, assign), ii, jj)
    rep["cheap_baselines"] = {k: round(v, 4) for k, v in cheap.items()}

    # ---- 2. the noise floor, measured by reseeding THIS code on THESE patients
    base, _ = run_arms(co, fi, 0)
    pooled = {k: [] for k in base}
    deltas = []
    for r in range(a.reseeds):
        f2 = reseeded_folds(len(co.keep), np.random.default_rng(1000 + r))
        arms, assign = run_arms(co, f2, 0)
        for k, v in arms.items():
            pooled[k].append(cidx(v, ii, jj))
        deltas.append(cidx(arms["rank_avg"], ii, jj) - cidx(arms["clinical"], ii, jj))
        print("reseed %2d/%d  rank_avg=%.4f clinical=%.4f titan=%.4f"
              % (r + 1, a.reseeds, pooled["rank_avg"][-1], pooled["clinical"][-1],
                 pooled["titan"][-1]), file=sys.stderr, flush=True)
    rep["noise_floor_measured_on_this_code"] = {
        "reseeds": a.reseeds,
        "what_was_varied": "the 5-fold assignment of the SAME 359 patients; arms, estimator, "
                           "alpha grid and inner-CV seed all held fixed",
        "per_arm_pooled_cindex": {k: {"mean": round(float(np.mean(v)), 4),
                                      "sd": round(float(np.std(v, ddof=1)), 4)}
                                  for k, v in pooled.items()},
        "delta_rank_avg_minus_clinical": {
            "mean": round(float(np.mean(deltas)), 4),
            "sd": round(float(np.std(deltas, ddof=1)), 4),
            "bar_2sd": round(float(2 * np.std(deltas, ddof=1)), 4)},
        "released_fold_values": {k: round(cidx(v, ii, jj), 4) for k, v in base.items()},
    }

    # ---- 3. the oracle and its structural control
    rep["oracle"] = oracle_block(co, fi, (2, 4, 8, 16), rng)
    rep["oracle_reading"] = (
        "recoverable_above_control is how much of the outcome this representation expresses once "
        "pure overfitting at the same capacity is subtracted. estimation_loss is what the "
        "estimator gives up going from in-sample to out-of-fold. If the first is small the "
        "representation is exhausted and the next component must change the inputs; if the second "
        "is large the information is present and unreached.")

    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
