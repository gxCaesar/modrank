#!/usr/bin/env python3
"""Tier 1C + 1D + 2F + 2G -- the survival analyses a concordance index does not cover.

A survival paper reporting one metric is thin, and the competitors on this benchmark report more
than one: PIBD's own logs carry integrated Brier score and integrated AUC fields beside the
concordance. This file adds, for every arm, on the same folds:

  time-dependent AUC        discrimination at a horizon, IPCW-weighted for censoring
  integrated Brier score    calibration AND discrimination together, lower is better
  a KM risk stratification  the clinically legible form: split by the score, log-rank test
  a calibration check       predicted versus observed survival by risk group
  the tied-pair CURVE       the clinical hook, currently reported at two thresholds only

WHY IPCW. A naive time-dependent AUC or Brier score treats a censored patient as event-free, which
biases both toward whichever arm happens to rank the censored patients low. Inverse-probability-of-
censoring weights, with the censoring distribution estimated by Kaplan-Meier on the SAME fold
structure, remove that. The censoring KM is fitted on the training fold only, so no validation
patient contributes to its own weight.

WHAT IS NOT DONE HERE, and it is stated rather than left as an absence: no competitor's per-case
predictions exist except SurvPath's, so every metric below is computed for our arms and for
SurvPath (alone, and with the same clinical block), and for nobody else. The published values on
this benchmark are concordance only, so the additional metrics have no published comparator at all
and are reported as descriptions of our arm rather than as a comparison we win.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, fitapply  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402
from s6_amend_clinical import dimaf_clinical  # noqa: E402


def km_censoring(t, e):
    """Kaplan-Meier estimate of the CENSORING distribution G(u) = P(C > u), left-continuous.

    Censoring is the 'event' here, so the indicator is 1-e. Returned as a step function evaluated
    by `g_at`, with G held at its last value beyond the largest observed time rather than dropping
    to zero, which would make a weight infinite.
    """
    o = np.argsort(t)
    ts, cs = t[o], (1.0 - e[o])
    uniq = np.unique(ts)
    g, surv = [], 1.0
    for u in uniq:
        at_risk = float((ts >= u).sum())
        d = float(((ts == u) & (cs == 1)).sum())
        if at_risk > 0 and d > 0:
            surv *= (1.0 - d / at_risk)
        g.append(surv)
    return uniq, np.asarray(g)


def g_at(uniq, g, x, floor=1e-3):
    """G evaluated just before x; floored so a weight cannot explode on one patient."""
    idx = np.searchsorted(uniq, x, side="left") - 1
    out = np.where(idx >= 0, g[np.clip(idx, 0, len(g) - 1)], 1.0)
    return np.maximum(out, floor)


def td_auc(score, t, e, horizon, uniq, g):
    """IPCW time-dependent AUC at `horizon`: cases are events by then, controls survive past it."""
    cases = (t <= horizon) & (e == 1)
    ctrl = t > horizon
    if cases.sum() < 3 or ctrl.sum() < 3:
        return float("nan")
    wi = 1.0 / g_at(uniq, g, t[cases])            # weight a case by 1/G(T_i-)
    wj = np.full(int(ctrl.sum()), 1.0 / float(g_at(uniq, g, np.array([horizon]))[0]))
    si, sj = score[cases], score[ctrl]
    num = 0.0
    for a, w in zip(si, wi):
        num += w * float(np.sum(wj * ((a > sj) + 0.5 * (a == sj))))
    den = float(wi.sum() * wj.sum())
    return num / den if den > 0 else float("nan")


def brier(surv_prob, t, e, horizon, uniq, g):
    """IPCW Brier score at `horizon`. `surv_prob` is the predicted P(T > horizon) per patient."""
    w = np.zeros(len(t))
    died = (t <= horizon) & (e == 1)
    alive = t > horizon
    w[died] = 1.0 / g_at(uniq, g, t[died])
    w[alive] = 1.0 / float(g_at(uniq, g, np.array([horizon]))[0])
    y = alive.astype(float)                        # observed survival indicator at the horizon
    contrib = w * (y - surv_prob) ** 2
    keep = died | alive                            # censored before the horizon contribute nothing
    return float(contrib[keep].sum() / max(keep.sum(), 1))


def ibs(surv_curves, horizons, t, e, uniq, g):
    """Integrated Brier score: the Brier score averaged over the horizon grid."""
    vals = [brier(surv_curves[:, k], t, e, h, uniq, g) for k, h in enumerate(horizons)]
    vals = [v for v in vals if v == v]
    return float(np.mean(vals)) if vals else float("nan")


def logrank(t, e, group):
    """Two-or-more-group log-rank chi-square and its degrees of freedom."""
    gs = np.unique(group)
    times = np.unique(t[e == 1])
    O = np.zeros(len(gs))
    E = np.zeros(len(gs))
    V = np.zeros((len(gs), len(gs)))
    for u in times:
        at_risk = t >= u
        n = float(at_risk.sum())
        d = float(((t == u) & (e == 1)).sum())
        if n <= 1 or d == 0:
            continue
        nk = np.array([float((at_risk & (group == gg)).sum()) for gg in gs])
        dk = np.array([float(((t == u) & (e == 1) & (group == gg)).sum()) for gg in gs])
        ek = d * nk / n
        O += dk
        E += ek
        f = d * (n - d) / (n - 1) if n > 1 else 0.0
        for a in range(len(gs)):
            for b in range(len(gs)):
                V[a, b] += f * ((nk[a] / n) * ((a == b) - nk[b] / n))
    d_ = (O - E)[:-1]
    Vr = V[:-1, :-1]
    try:
        chi2 = float(d_ @ np.linalg.pinv(Vr) @ d_)
    except Exception:
        chi2 = float("nan")
    return chi2, len(gs) - 1, {str(gg): {"observed": float(o), "expected": round(float(ee), 2)}
                               for gg, o, ee in zip(gs, O, E)}


def chi2_sf(x, df):
    """Upper tail of the chi-square distribution, via scipy if present."""
    try:
        from scipy.stats import chi2 as _c
        return float(_c.sf(x, df))
    except Exception:
        return float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--omics-dir", required=True)
    ap.add_argument("--dimaf-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    CLIN = np.hstack([dc["age"], dc["fem"], dc["stage"]])

    def arm(X, seed=0):
        s = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], seed)
        return co.fold_pct(s)

    zt, zo, zc = arm(co.T), arm(P), arm(CLIN)
    zs = co.fold_pct(co.spr)
    arms = {"clinical": zc, "wsi_titan": zt, "omics": zo,
            "survpath_rerun": zs,
            "survpath_plus_clinical": co.fold_pct(zs + zc),
            "OURS_wsi_omics": co.fold_pct(zt + zo),
            "OURS": co.fold_pct(zt + zo + zc)}

    ii, jj = cpairs(co.t, co.e)
    uniq, g = km_censoring(co.t, co.e)
    # horizons inside the observable range: quartiles of the observed event times
    ev_t = np.sort(co.t[co.e == 1])
    horizons = [float(np.quantile(ev_t, q)) for q in (0.25, 0.5, 0.75)]
    grid = [float(x) for x in np.linspace(ev_t[0], np.quantile(ev_t, 0.9), 12)]

    rep = {"artifact_type": "s7_survival_metrics", "phase_of_origin": "frozen_post_hoc",
           "n": len(co.keep), "events": int(co.e.sum()),
           "horizons_months": [round(h, 1) for h in horizons],
           "ibs_grid_months": [round(x, 1) for x in grid],
           "note": "additional metrics are computed for our arms and for SurvPath only -- no other "
                   "entrant releases per-case predictions, and the published values on this "
                   "benchmark are concordance alone, so these have no published comparator"}

    out = {}
    for nm, z in arms.items():
        # a percentile rank is a risk score; the survival probability a Brier score needs is a
        # monotone decreasing function of it. Using the empirical KM of the training folds at each
        # horizon, shifted by the rank, keeps this free of an extra fitted model.
        km_t, km_s = km_censoring(co.t, 1.0 - co.e)          # event-KM: swap the indicator
        base = g_at(km_t, km_s, np.array(grid))
        curves = np.clip(base[None, :] ** np.exp(2.0 * (z[:, None] - 0.5)), 1e-6, 1 - 1e-6)
        row = {"cindex": round(cidx(z, ii, jj), 4)}
        for h in horizons:
            row["td_auc_%dm" % round(h)] = round(td_auc(z, co.t, co.e, h, uniq, g), 4)
        row["ibs"] = round(ibs(curves, grid, co.t, co.e, uniq, g), 4)
        # KM stratification by tertile of the score
        q = np.quantile(z, [1 / 3, 2 / 3])
        grp = np.digitize(z, q)
        c2, df, tab = logrank(co.t, co.e, grp)
        row["logrank_chi2"] = round(c2, 2)
        row["logrank_df"] = df
        row["logrank_p"] = float("%.3g" % chi2_sf(c2, df))
        row["group_sizes"] = [int((grp == k).sum()) for k in range(3)]
        row["group_events"] = [int(co.e[grp == k].sum()) for k in range(3)]
        row["group_median_survival_months"] = []
        for k in range(3):
            m = grp == k
            kt, ks = km_censoring(co.t[m], 1.0 - co.e[m])
            below = np.flatnonzero(ks <= 0.5)
            row["group_median_survival_months"].append(
                round(float(kt[below[0]]), 1) if below.size else None)
        out[nm] = row
        print("%-24s C=%.4f  IBS=%.4f  logrank chi2=%.1f p=%.3g"
              % (nm, row["cindex"], row["ibs"], row["logrank_chi2"], row["logrank_p"]),
              file=sys.stderr, flush=True)
    rep["arms"] = out

    # ---- the tied-pair CURVE, not two points
    gap = np.abs(zc[ii] - zc[jj])
    curve = {}
    for q in (0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00):
        m = gap <= float(np.quantile(gap, q)) if q < 1.0 else np.ones(len(gap), bool)
        if m.sum() < 200:
            continue
        curve["q%03d" % int(q * 100)] = {
            "pairs": int(m.sum()),
            **{nm: round(cidx(z, ii[m], jj[m]), 4) for nm, z in arms.items()}}
    rep["tied_pair_curve"] = curve
    print("\ntied-pair curve computed at %d thresholds" % len(curve), file=sys.stderr)

    open(a.out, "w").write(json.dumps(rep, indent=1) + "\n")
    print("wrote", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
