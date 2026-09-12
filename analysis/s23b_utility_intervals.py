#!/usr/bin/env python3
"""Intervals for the differences s23 reports as point estimates: accuracy and net benefit.

WHY. The calibration paragraph says ModRank's index of prediction accuracy (IPA) is higher than the
clinical model's at one and two years and that its net benefit is higher at two years across a
range of thresholds. s23 gives each model's value, not the uncertainty of the difference, and a
difference between two point estimates on 359 patients is not a finding until it has an interval.

HOW. The per-fold construction is s23's, line for line, and this file refuses to report unless it
reproduces every value s23 committed (the two known concordances, the calibration slopes, and the
Brier score and IPA of both models at all three horizons). The predictions are then held fixed and
the evaluation is resampled: 6,000 draws of patients with replacement (rng 20260911), and in each
draw the censoring distribution, both Brier scores, the null Brier score, both IPAs and both net
benefits are recomputed on the drawn patients. The interval is the 2.5th and 97.5th percentile of
the difference, and p is twice the smaller tail. Resampling the evaluation only is the standard
conditional interval: it does not include the variation from refitting the models.

Post-freeze and descriptive. Nothing here changes a reported primary.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cox_fit, cpairs, fitapply, pct   # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix                     # noqa: E402
from s6_amend_clinical import dimaf_clinical                          # noqa: E402
from s7_survival_metrics import km_censoring                          # noqa: E402
from s23_calibration_dca import HORIZONS, H_at, breslow, brier_graf, km_at   # noqa: E402

REPS, SEED = 6000, 20260911
NB_THRESHOLDS = (0.2, 0.3, 0.4, 0.5, 0.6)


def net_benefit(pr, t, e, th, h=24.0):
    n = len(t)
    hi = pr >= th
    n_hi = int(hi.sum())
    s_hi = km_at(t[hi], e[hi], h) if n_hi else 1.0
    return (n_hi * (1.0 - s_hi) - n_hi * s_hi * th / (1.0 - th)) / n


def summary(v):
    v = np.asarray(v)
    return {"mean": round(float(v.mean()), 4),
            "ci95": [round(float(np.quantile(v, 0.025)), 4), round(float(np.quantile(v, 0.975)), 4)],
            "p_two_sided": round(float(2 * min((v <= 0).mean(), (v >= 0).mean())), 4)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("root", "titan", "sp-dir", "omics-dir", "dimaf-dir", "committed", "out"):
        ap.add_argument("--" + f, required=True)
    a = ap.parse_args()
    co = Cohort(a.root, a.titan, a.sp_dir)
    n, t, e = len(co.keep), co.t, co.e
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    B = {"slide": co.T, "omics": P, "clinical": np.hstack([dc["age"], dc["fem"], dc["stage"]])}

    # ---- s23's construction, unchanged
    score_val = {"ours": np.full(n, np.nan), "clinical": np.full(n, np.nan)}
    surv = {k: {h: np.full(n, np.nan) for h in HORIZONS} for k in score_val}
    surv_null = {h: np.full(n, np.nan) for h in HORIZONS}
    lp = {k: np.full(n, np.nan) for k in score_val}
    for tri, vai in co.fold_indices():
        zv = {b: pct(fitapply(X[tri], t[tri], e[tri], X[vai], 0)) for b, X in B.items()}
        v_scores = {"ours": pct(zv["slide"] + zv["omics"] + zv["clinical"]), "clinical": zv["clinical"]}
        m = len(tri)
        inner = np.array_split(np.random.default_rng(0).permutation(m), 3)
        zi = {b: np.full(m, np.nan) for b in B}
        i_ours = np.full(m, np.nan)
        for ite in inner:
            itr = np.setdiff1d(np.arange(m), ite)
            for b, X in B.items():
                Xt = X[tri]
                zi[b][ite] = pct(fitapply(Xt[itr], t[tri][itr], e[tri][itr], Xt[ite], 0))
            i_ours[ite] = pct(zi["slide"][ite] + zi["omics"][ite] + zi["clinical"][ite])
        i_scores = {"ours": i_ours, "clinical": zi["clinical"]}
        for k in score_val:
            beta = float(cox_fit(i_scores[k][:, None], t[tri], e[tri], 1e-4)[0])
            taus, H = breslow(t[tri], e[tri], beta * i_scores[k])
            score_val[k][vai] = v_scores[k]
            lp[k][vai] = beta * v_scores[k]
            for h in HORIZONS:
                surv[k][h][vai] = np.exp(-H_at(taus, H, h) * np.exp(beta * v_scores[k]))
        for h in HORIZONS:
            surv_null[h][vai] = km_at(t[tri], e[tri], h)

    # ---- known answers: every value s23 committed
    com = json.load(open(a.committed))
    ii, jj = cpairs(t, e)
    uniq, g = km_censoring(t, e)
    bad = []
    for k in score_val:
        if round(cidx(score_val[k], ii, jj), 4) != com["known_answer"]["%s_seed0" % k]:
            bad.append((k, "concordance"))
        if round(float(cox_fit(lp[k][:, None], t, e, 1e-6)[0]), 4) != com[k]["calibration_slope"]:
            bad.append((k, "slope"))
        for h in HORIZONS:
            bs = brier_graf(surv[k][h], t, e, h, uniq, g)
            bs0 = brier_graf(surv_null[h], t, e, h, uniq, g)
            c = com[k]["by_horizon"][str(h)]
            if round(bs, 4) != c["brier"] or round(1 - bs / bs0, 4) != c["ipa"]:
                bad.append((k, h))
    if bad:
        print(json.dumps({"status": "error", "error_code": "s23_not_reproduced", "which": str(bad)}),
              file=sys.stderr)
        return 2
    print("s23 reproduced: concordances, slopes, Brier and IPA at 12, 24, 36 months", file=sys.stderr)

    # ---- the evaluation, resampled
    pr = {k: 1.0 - surv[k][24.0] for k in score_val}
    rng = np.random.default_rng(SEED)
    d_ipa = {h: [] for h in HORIZONS}
    d_nb = {th: [] for th in NB_THRESHOLDS}
    for _ in range(REPS):
        bs_ = rng.choice(n, size=n, replace=True)
        tb, eb = t[bs_], e[bs_]
        ub, gb = km_censoring(tb, eb)
        for h in HORIZONS:
            b0 = brier_graf(surv_null[h][bs_], tb, eb, h, ub, gb)
            ipa = {k: 1.0 - brier_graf(surv[k][h][bs_], tb, eb, h, ub, gb) / b0 for k in score_val}
            d_ipa[h].append(ipa["ours"] - ipa["clinical"])
        for th in NB_THRESHOLDS:
            d_nb[th].append(net_benefit(pr["ours"][bs_], tb, eb, th) - net_benefit(pr["clinical"][bs_], tb, eb, th))
    out = {"artifact_type": "s23b_utility_intervals", "phase_of_origin": "post_freeze_2026-09-11",
           "reportable": True, "reproduces": os.path.basename(a.committed),
           "bootstrap": {"replicates": REPS, "unit": "patient", "rng_seed": SEED,
                         "what_is_resampled": "the evaluation; predictions are held fixed"},
           "ipa_modrank_minus_clinical": {str(h): summary(d_ipa[h]) for h in HORIZONS},
           "net_benefit_24m_modrank_minus_clinical": {str(th): summary(d_nb[th]) for th in NB_THRESHOLDS}}
    json.dump(out, open(a.out, "w"), indent=1)
    for h in HORIZONS:
        s = out["ipa_modrank_minus_clinical"][str(h)]
        print("IPA difference at %2d months: %+.4f [%+.4f, %+.4f] p=%.4f"
              % (h, s["mean"], s["ci95"][0], s["ci95"][1], s["p_two_sided"]), file=sys.stderr)
    for th in NB_THRESHOLDS:
        s = out["net_benefit_24m_modrank_minus_clinical"][str(th)]
        print("net benefit difference at %.0f%%: %+.4f [%+.4f, %+.4f] p=%.4f"
              % (100 * th, s["mean"], s["ci95"][0], s["ci95"][1], s["p_two_sided"]), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
