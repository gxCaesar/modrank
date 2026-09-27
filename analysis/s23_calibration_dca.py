#!/usr/bin/env python3
"""Calibration and clinical utility, with every quantity fitted inside the training fold.

POST-FREEZE AND DESCRIPTIVE. WHY IT REPLACES WHAT WAS THERE. The manuscript's calibration paragraph
came from a heuristic transform of the percentile score (s7_interpret_calibrate_cases.py) with both
Kaplan-Meier curves fitted on all 359 patients. That is not a calibration of anything that could be
deployed. TRIPOD-AI, which Nature Communications requires for prognostic models, asks for
calibration of the model as it would be used.

HOW, per outer fold of the released partition, seed 0 (the protocol's reference seed):
  * the validation cases get the frozen construction's score: each arm fitted on the whole training
    fold, validation percentiles, their sum re-ranked (this reproduces 0.7225 at seed 0, checked);
  * the TRAINING cases get the same construction's score from a 3-fold inner split, so the
    recalibration never sees a training case scored by a model that was fitted on it;
  * one Cox coefficient on that inner score, and a Breslow baseline hazard, both from training cases
    only, turn a validation score into a survival probability at 12, 24 and 36 months.
The clinical arm (age, sex, stage) gets the identical treatment, so the two are compared on equal
terms. Nothing here uses a validation outcome.

REPORTED, pooled over the five validation folds: the calibration slope (Cox on the predicted log
relative hazard; 1 is ideal), calibration-in-the-large and grouped calibration at 24 months, the
IPCW Brier score (Graf et al.: divided by n, unlike the helper in s7, which divided by the uncensored
count and fed no reported number), the index of prediction accuracy against a training-fold
Kaplan-Meier null, and net benefit across thresholds at 24 months against treat-all and treat-none.
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
from s7_survival_metrics import g_at, km_censoring                    # noqa: E402

HORIZONS = (12.0, 24.0, 36.0)
THRESHOLDS = tuple(round(x, 2) for x in np.arange(0.05, 0.61, 0.05))


def km_at(t, e, h):
    """Kaplan-Meier survival at h."""
    s = 1.0
    for tau in np.unique(t[(e == 1) & (t <= h)]):
        at_risk = float((t >= tau).sum())
        d = float(((t == tau) & (e == 1)).sum())
        s *= 1.0 - d / at_risk
    return s


def breslow(t, e, lp):
    taus = np.unique(t[e == 1])
    w = np.exp(lp)
    H = np.cumsum([((t == tau) & (e == 1)).sum() / w[t >= tau].sum() for tau in taus])
    return taus, H


def H_at(taus, H, h):
    k = np.searchsorted(taus, h, side="right") - 1
    return float(H[k]) if k >= 0 else 0.0


def brier_graf(surv, t, e, h, uniq, g):
    died = (t <= h) & (e == 1)
    alive = t > h
    w = np.zeros(len(t))
    w[died] = 1.0 / g_at(uniq, g, t[died])
    w[alive] = 1.0 / float(g_at(uniq, g, np.array([h]))[0])
    return float(np.sum(w * (alive.astype(float) - surv) ** 2) / len(t))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("root", "titan", "sp-dir", "omics-dir", "dimaf-dir", "out"):
        ap.add_argument("--" + f, required=True)
    a = ap.parse_args()
    co = Cohort(a.root, a.titan, a.sp_dir)
    n, t, e = len(co.keep), co.t, co.e
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    B = {"slide": co.T, "omics": P, "clinical": np.hstack([dc["age"], dc["fem"], dc["stage"]])}

    score_val = {"ours": np.full(n, np.nan), "clinical": np.full(n, np.nan)}
    surv = {k: {h: np.full(n, np.nan) for h in HORIZONS} for k in score_val}
    surv_null = {h: np.full(n, np.nan) for h in HORIZONS}
    lp = {k: np.full(n, np.nan) for k in score_val}
    betas = {k: [] for k in score_val}
    for tri, vai in co.fold_indices():
        zv = {b: pct(fitapply(X[tri], t[tri], e[tri], X[vai], 0)) for b, X in B.items()}
        v_scores = {"ours": pct(zv["slide"] + zv["omics"] + zv["clinical"]),
                    "clinical": zv["clinical"]}
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
            betas[k].append(round(beta, 4))
            taus, H = breslow(t[tri], e[tri], beta * i_scores[k])
            score_val[k][vai] = v_scores[k]
            lp[k][vai] = beta * v_scores[k]
            for h in HORIZONS:
                surv[k][h][vai] = np.exp(-H_at(taus, H, h) * np.exp(beta * v_scores[k]))
        for h in HORIZONS:
            surv_null[h][vai] = km_at(t[tri], e[tri], h)

    ii, jj = cpairs(t, e)
    known = {"ours_seed0": round(cidx(score_val["ours"], ii, jj), 4),
             "clinical_seed0": round(cidx(score_val["clinical"], ii, jj), 4)}
    import blca_common
    _tol = blca_common.A2_SANITY if blca_common.A2 else 5e-5   # amendment A2: exact check under BLCA_A2=0
    if abs(known["ours_seed0"] - 0.7225) > _tol or abs(known["clinical_seed0"] - 0.6638) > _tol:
        print(json.dumps({"status": "error", "error_code": "known_answer_failed", "got": known}),
              file=sys.stderr)
        return 2

    uniq, g = km_censoring(t, e)
    out = {"artifact_type": "s23_calibration_and_decision_curve",
           "phase_of_origin": "post_freeze_2026-09-11", "reportable": True,
           "known_answer": known, "horizons_months": list(HORIZONS),
           "at_risk": {str(h): int((t > h).sum()) for h in HORIZONS},
           "events_by": {str(h): int(((t <= h) & (e == 1)).sum()) for h in HORIZONS},
           "recalibration_slopes_per_fold": betas}
    for k in score_val:
        slope = float(cox_fit(lp[k][:, None], t, e, 1e-6)[0])
        rec = {"calibration_slope": round(slope, 4), "by_horizon": {}}
        for h in HORIZONS:
            pred_risk = 1.0 - surv[k][h]
            bs = brier_graf(surv[k][h], t, e, h, uniq, g)
            bs0 = brier_graf(surv_null[h], t, e, h, uniq, g)
            rec["by_horizon"][str(h)] = {
                "mean_predicted_risk": round(float(pred_risk.mean()), 4),
                "observed_km_risk": round(1.0 - km_at(t, e, h), 4),
                "brier": round(bs, 4), "brier_null_km": round(bs0, 4),
                "ipa": round(1.0 - bs / bs0, 4)}
        # grouped calibration and net benefit at 24 months
        h = 24.0
        pr = 1.0 - surv[k][h]
        q = np.quantile(pr, [0.2, 0.4, 0.6, 0.8])
        grp = np.digitize(pr, q)
        rec["grouped_24m"] = [{"quintile": int(j + 1), "n": int((grp == j).sum()),
                               "events": int(e[grp == j].sum()),
                               "predicted": round(float(pr[grp == j].mean()), 4),
                               "observed_km": round(1.0 - km_at(t[grp == j], e[grp == j], h), 4)}
                              for j in range(5)]
        nb = []
        for th in THRESHOLDS:
            hi = pr >= th
            n_hi = int(hi.sum())
            s_hi = km_at(t[hi], e[hi], h) if n_hi else 1.0
            tp, fp = n_hi * (1.0 - s_hi), n_hi * s_hi
            nb.append(round((tp - fp * th / (1.0 - th)) / n, 4))
        rec["net_benefit_24m"] = dict(zip([str(x) for x in THRESHOLDS], nb))
        out[k] = rec
    s_all = km_at(t, e, 24.0)
    out["treat_all_net_benefit_24m"] = {str(th): round((1 - s_all) - s_all * th / (1 - th), 4)
                                        for th in THRESHOLDS}
    json.dump(out, open(a.out, "w"), indent=1)
    for k in score_val:
        r = out[k]
        print("%-9s slope %.3f | 24m pred %.3f obs %.3f Brier %.4f IPA %.4f"
              % (k, r["calibration_slope"], r["by_horizon"]["24.0"]["mean_predicted_risk"],
                 r["by_horizon"]["24.0"]["observed_km_risk"], r["by_horizon"]["24.0"]["brier"],
                 r["by_horizon"]["24.0"]["ipa"]), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
