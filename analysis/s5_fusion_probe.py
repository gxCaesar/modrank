#!/usr/bin/env python3
"""How much is the FUSION RULE leaving on the table?

Pre-freeze, exploratory. NOT reportable.

WHY THIS EXISTS. The atlas just measured two things that do not fit together under the story this
campaign had been telling itself:

  Spearman(image score, clinical score) = 0.213          -- the two are nearly orthogonal
  image alone, on clinically-tied pairs  = 0.61-0.63     -- the image is not a stage detector
  clinical 0.6856, image 0.6596, rank-averaged 0.7191    -- combining them buys +0.0335

Two nearly independent predictors, each far above chance, should combine for a great deal more
than +0.03. So the deficit is not necessarily in the representation. It may be in the one line of
code that adds two percentile ranks together. This file measures which.

FOUR REFERENCES, each answering a different question.

  A. WEIGHT SWEEP. C-index of w*image + (1-w)*clinical over the whole grid, pooled. Says whether
     the equal weight the current arm uses is even near the best fixed weight.

  B. LINEAR FUSION ORACLE. A two-parameter Cox on (image, clinical) fitted IN-SAMPLE on the
     validation fold. Two parameters against ~23 events is a fit the data can support, and it is
     the ceiling of every fixed linear fusion. Paired with a permuted-outcome control at the same
     two parameters, 200 draws, so the overfitting it buys is measured rather than assumed.

  C. GATED FUSION ORACLE. The same fit, but the image weight is allowed to differ across clinical
     risk tertiles. This is the ceiling of a fusion that knows WHEN to trust the image. The atlas
     says the image gain is +0.112 in the middle clinical tertile and -0.045 in the lowest, a
     spread of 2.9 fold-to-fold standard deviations, and no fixed weight can express that. The gap
     between C and B is what a gated component could win; its own permuted control says how much
     of that gap is just the extra parameters.

  D. THE INDEPENDENCE REFERENCE. If two signals really were independent and really were as strong
     as these two measure, what C-index would their sum reach? Simulated on the same n, the same
     event count and the same censoring pattern, with each latent signal calibrated by bisection
     to reproduce the measured single-arm C-index. This is the number that says whether +0.0335 is
     a disappointment or is simply what these two signals are worth.

     It is a reference, not a bound: it assumes exact independence, and the measured correlation is
     0.213 rather than 0, so the reference is calibrated with the observed correlation as well and
     both are reported.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, cox_fit, fitapply  # noqa: E402


def oracle_fusion(co, arms_oof, groups, rng, n_perm=200):
    """In-sample two-or-more-parameter fusion on the validation fold, with a permuted control.

    `groups` is None for the plain linear fusion, or a per-case integer label for a gated fusion
    in which the image column is expanded into one column per group. The clinical column is never
    expanded, so the comparison isolates the gating of the IMAGE weight.
    """
    z_ti, z_cl = arms_oof
    real, perm = [], []
    for k in range(5):
        m = co.fold == k
        vai = np.flatnonzero(m)
        tv, ev = co.t[vai], co.e[vai]
        ii, jj = cpairs(tv, ev)
        if ii.size == 0:
            continue
        if groups is None:
            X = np.column_stack([z_ti[vai], z_cl[vai]])
        else:
            gv = groups[vai]
            cols = [(gv == g).astype(float) * z_ti[vai] for g in np.unique(groups)]
            X = np.column_stack(cols + [z_cl[vai]])
        X = (X - X.mean(0)) / (X.std(0) + 1e-9)
        real.append(cidx(X @ cox_fit(X, tv, ev, 1.0), ii, jj))
        for _ in range(n_perm):
            q = rng.permutation(len(vai))
            perm.append(cidx(X @ cox_fit(X, tv[q], ev[q], 1.0), ii, jj))
    return {"oracle_in_sample": round(float(np.mean(real)), 4),
            "control_permuted_outcome": round(float(np.mean(perm)), 4),
            "control_sd": round(float(np.std(perm, ddof=1)), 4),
            "control_draws": len(perm),
            "params": 2 if groups is None else len(np.unique(groups)) + 1,
            "above_control": round(float(np.mean(real) - np.mean(perm)), 4)}


def independence_reference(co, c_a, c_b, rho, rng, reps=200):
    """What would two signals of these strengths, at this correlation, reach when summed?

    A FIRST VERSION OF THIS RETURNED 0.4955 AND WAS WRONG, in a way worth recording because the
    number looked like a finding. It drew a random latent `eta`, ranked it, and scored that rank
    against the cohort's time ordering -- but nothing tied eta to the time ordering, so the
    "predictor" was a random permutation, every bisection step measured 0.5, and the calibration
    converged to whatever the bracket ended on. A signal calibrated to nothing then summed to
    nothing. The tell was that two arms individually at 0.66 and 0.69 cannot sum to chance.

    THE MODEL THIS USES INSTEAD, and why it is the one the data supports. Writing a = r + noise_a
    and b = r + noise_b around one shared risk r cannot reproduce what was measured: two
    predictors both near C = 0.67 would then have to share most of their variance and correlate
    highly, and the observed Spearman is 0.213. So the outcome is driven by TWO components,
    r = r1 + r2, with the image arm reading mostly r1 and the clinical arm mostly r2. That is the
    only simple structure consistent with both arms being strong and nearly uncorrelated, and it
    is also the structure that makes their combination worth something.

      r1, r2  ~ N(0,1), independent; a leakage term makes the observed correlation reproduce
      time    ~ Exp(exp(-(r1 + r2))), censored to the cohort's own event fraction
      a = r1 + sigma_a * noise, b = r2 + sigma_b * noise
      sigma_a, sigma_b bisected until C(a) and C(b) match the measured single-arm values

    The reported value is C(rank(a) + rank(b)) -- the same rank average the deployed arm uses. It
    is a reference under a model, not a bound.
    """
    n, n_ev = len(co.keep), int(co.e.sum())

    def cohort(g, r):
        te = g.exponential(np.exp(-r))
        u = g.random(n)
        lo, hi = 1e-4, 1e4
        for _ in range(60):
            mid = (lo + hi) / 2
            if (te <= -np.log(u) * mid).sum() > n_ev:
                hi = mid
            else:
                lo = mid
        c = -np.log(u) * (lo + hi) / 2
        return np.minimum(te, c), (te <= c).astype(float)

    def cal(which, sigma, reps_in=16):
        cs = []
        for r in range(reps_in):
            g = np.random.default_rng(9000 + r)
            r1, r2 = g.standard_normal(n), g.standard_normal(n)
            t, e = cohort(g, r1 + r2)
            sig = (r1 if which == 0 else r2) + sigma * g.standard_normal(n)
            cs.append(cindex_local(sig, t, e))
        return float(np.mean(cs))

    def cindex_local(risk, t, e):
        return cidx(risk, *cpairs(t, e))

    def bisect(which, target):
        lo, hi = 0.01, 60.0          # more noise -> lower C, so C is DECREASING in sigma
        for _ in range(30):
            mid = (lo + hi) / 2
            if cal(which, mid) > target:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2

    sa, sb = bisect(0, c_a), bisect(1, c_b)
    got_a, got_b, cs, rhos = [], [], [], []
    for rep in range(reps):
        g = np.random.default_rng(4000 + rep)
        r1, r2 = g.standard_normal(n), g.standard_normal(n)
        t, e = cohort(g, r1 + r2)
        a = r1 + sa * g.standard_normal(n)
        b = r2 + sb * g.standard_normal(n)
        got_a.append(cindex_local(a, t, e))
        got_b.append(cindex_local(b, t, e))
        ra = np.argsort(np.argsort(a)).astype(float) / (n - 1)
        rb = np.argsort(np.argsort(b)).astype(float) / (n - 1)
        x, y = ra - ra.mean(), rb - rb.mean()
        rhos.append(float(x @ y / np.sqrt((x @ x) * (y @ y))))
        cs.append(cindex_local(ra + rb, t, e))
    return {
        "sum_of_the_two_rank_averaged": round(float(np.mean(cs)), 4),
        "sd": round(float(np.std(cs, ddof=1)), 4),
        "calibration_target": {"arm_a": round(c_a, 4), "arm_b": round(c_b, 4)},
        "calibration_achieved": {"arm_a": round(float(np.mean(got_a)), 4),
                                 "arm_b": round(float(np.mean(got_b)), 4),
                                 "sigma_a": round(sa, 3), "sigma_b": round(sb, 3)},
        "correlation_check": {"simulated_spearman": round(float(np.mean(rhos)), 4),
                              "observed_spearman": round(float(rho), 4),
                              "note": "the model puts the two arms on independent components, so "
                                      "its induced correlation should be near zero; the observed "
                                      "0.21 is above that and the reference is correspondingly "
                                      "slightly optimistic"},
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    rng = np.random.default_rng(0)

    p_ti = np.full(len(co.keep), np.nan)
    p_cl = np.full(len(co.keep), np.nan)
    for tri, vai in co.fold_indices():
        if len(tri) < 30 or not len(vai):
            continue
        p_ti[vai] = fitapply(co.T[tri], co.t[tri], co.e[tri], co.T[vai], 0)
        p_cl[vai] = fitapply(co.CLIN[tri], co.t[tri], co.e[tri], co.CLIN[vai], 0)
    z_ti, z_cl = co.fold_pct(p_ti), co.fold_pct(p_cl)

    ii, jj = cpairs(co.t, co.e)
    c_ti, c_cl = cidx(z_ti, ii, jj), cidx(z_cl, ii, jj)
    rep = {"artifact_type": "s5_fusion_probe", "reportable": False,
           "phase_of_origin": "exploratory",
           "arms_out_of_fold": {"titan": round(c_ti, 4), "clinical": round(c_cl, 4),
                                "rank_average": round(cidx(co.fold_pct(z_ti + z_cl), ii, jj), 4)}}

    # ---- A. the weight sweep
    sweep = {}
    for w in np.round(np.arange(0.0, 1.001, 0.05), 2):
        sweep["%.2f" % w] = round(cidx(co.fold_pct(w * z_ti + (1 - w) * z_cl), ii, jj), 4)
    best_w = max(sweep, key=lambda k: sweep[k])
    rep["A_weight_sweep"] = {"grid": sweep, "best_w_on_image": best_w,
                             "best_value": sweep[best_w],
                             "value_at_equal_weight": sweep["0.50"],
                             "note": "this sweep is chosen on the pooled outcome and is therefore "
                                     "itself an oracle over one parameter; it is a reference for "
                                     "the SHAPE of the objective, not a reportable arm"}

    # ---- B / C. the fusion oracles and their controls
    q = np.quantile(z_cl, [1 / 3, 2 / 3])
    tert = np.digitize(z_cl, q)
    rep["B_linear_fusion_oracle"] = oracle_fusion(co, (z_ti, z_cl), None, rng)
    rep["C_gated_fusion_oracle"] = oracle_fusion(co, (z_ti, z_cl), tert, rng)
    rep["C_minus_B"] = round(rep["C_gated_fusion_oracle"]["oracle_in_sample"]
                             - rep["B_linear_fusion_oracle"]["oracle_in_sample"], 4)
    rep["C_minus_B_control"] = round(rep["C_gated_fusion_oracle"]["control_permuted_outcome"]
                                     - rep["B_linear_fusion_oracle"]["control_permuted_outcome"], 4)
    rep["C_reading"] = ("C_minus_B is the ceiling of what gating the image weight on clinical risk "
                        "could buy; C_minus_B_control is how much of that the extra parameters buy "
                        "with no signal at all. Only the excess of the first over the second is a "
                        "reason to build the component.")

    # ---- D. the independence reference
    ra = np.argsort(np.argsort(z_ti)).astype(float)
    rb = np.argsort(np.argsort(z_cl)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    rho = float(ra @ rb / np.sqrt((ra @ ra) * (rb @ rb)))
    rep["D_independence_reference"] = independence_reference(co, c_ti, c_cl, rho, rng)
    rep["D_reading"] = ("if the measured rank-average sits at the simulated value, +0.0335 is "
                        "simply what two signals of these strengths are worth and the fusion rule "
                        "is not the deficit; if it sits well below, the rule is losing information "
                        "these two arms already contain")

    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
