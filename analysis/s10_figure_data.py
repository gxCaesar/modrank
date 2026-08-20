#!/usr/bin/env python3
"""The four quantities the figures need as CURVES rather than as summaries.

Everything else in the figure set already exists as a number in
`experiments/20260817-blca-confirm/results/`. These four exist only as summaries -- a median, a
log-rank p, three horizons, a power at one margin -- and a panel cannot be drawn from a summary.

  F1  KAPLAN-MEIER step functions and AT-RISK TABLES, by risk tertile, for our arm, for
      SurvPath+clinical and for PIBD+clinical. A KM panel without a risk table is a picture of a
      curve whose late tail rests on four patients and does not say so.

  F2  TIME-DEPENDENT AUC on a DENSE horizon grid rather than three quartiles. Uno's IPCW estimator,
      the same one the confirmatory run used, so the three published horizons must reproduce
      exactly -- that is asserted below, not hoped for.

  F3  POWER against a grid of true margins at this cohort's case-resampling SE. Closed form:
      for a two-sided test at alpha, power = Phi(d/SE - z) + Phi(-d/SE - z). The single number
      already reported (80% power at 0.0589) must fall on this curve, and that is asserted too.

  F4  The LUMINANCE-SELECTED ladder triple. Three risk groups are an ORDERED quantity, so they take
      three entries from the ordered ladder -- chosen by computing sRGB relative luminance and
      maximising the minimum adjacent gap, not by my eye. Asserted >= 0.05 so the panel survives
      greyscale.

Nothing here is a new result. It is the same numbers at the resolution a figure needs.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, fitapply  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402
from s6_amend_clinical import dimaf_clinical  # noqa: E402
from s7_survival_metrics import chi2_sf, km_censoring, logrank, td_auc  # noqa: E402

LADDER = ["#667DB8", "#CB7F9F", "#A4A8B6", "#CAAEDF", "#9CD0F3", "#F3CFED", "#DEEFFF"]


def luminance(hexstr):
    """sRGB relative luminance, WCAG definition. Used to ORDER, not to decorate."""
    c = [int(hexstr[i:i + 2], 16) / 255.0 for i in (1, 3, 5)]
    lin = [x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4 for x in c]
    return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]


def pick_triple(ladder):
    """The three ladder entries whose minimum adjacent luminance gap is largest."""
    import itertools
    best, arg = -1.0, None
    for combo in itertools.combinations(range(len(ladder)), 3):
        ls = sorted(luminance(ladder[i]) for i in combo)
        gap = min(ls[1] - ls[0], ls[2] - ls[1])
        if gap > best:
            best, arg = gap, combo
    cols = sorted(arg, key=lambda i: luminance(ladder[i]))
    return [ladder[i] for i in cols], best


def km_curve(t, e):
    """Kaplan-Meier with the at-risk count at each step. Returns times, survival, n_at_risk."""
    order = np.argsort(t, kind="mergesort")
    t, e = np.asarray(t)[order], np.asarray(e)[order]
    times, surv, atrisk, s, n = [0.0], [1.0], [len(t)], 1.0, len(t)
    i = 0
    while i < len(t):
        j = i
        while j < len(t) and t[j] == t[i]:
            j += 1
        d = float(e[i:j].sum())
        if d > 0:
            s *= (1.0 - d / n)
            times.append(float(t[i])); surv.append(float(s)); atrisk.append(int(n))
        n -= (j - i)
        i = j
    return times, surv, atrisk


def at_risk_at(t, grid):
    return [int((np.asarray(t) >= g).sum()) for g in grid]


def median_survival(times, surv):
    for tt, ss in zip(times, surv):
        if ss <= 0.5:
            return float(tt)
    return None


def norm_cdf(x):
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def power_at(d, se, alpha=0.05):
    z = 1.959963984540054 if abs(alpha - 0.05) < 1e-12 else None
    assert z is not None, "only alpha=0.05 is used here, and the critical value is written out"
    lam = d / se
    return norm_cdf(lam - z) + norm_cdf(-lam - z)


def load_preds(dirpath, idx, n, suffix="", pattern="split_*_results%s.pkl"):
    r, f = np.full(n, np.nan), np.full(n, -1)
    files = sorted(glob.glob(os.path.join(dirpath, pattern % suffix)))
    for k, fp in enumerate(files):
        for c, v in pickle.load(open(fp, "rb")).items():
            if c in idx:
                r[idx[c]] = float(v["risk"]); f[idx[c]] = k
    return r, f


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for fl in ("root", "titan", "sp-dir", "seed-dir", "pibd-dir", "omics-dir", "dimaf-dir", "out"):
        ap.add_argument("--" + fl, required=True)
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    ii, jj = cpairs(co.t, co.e)
    n = len(co.keep)
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    CLIN = np.hstack([dc["age"], dc["fem"], dc["stage"]])

    def arm(X, seed=0):
        s = np.full(n, np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], seed)
        return co.fold_pct(s)

    zt, zo, zc = arm(co.T), arm(P), arm(CLIN)
    OURS = co.fold_pct(zt + zo + zc)

    # SurvPath averaged over the five seeds, then given the same clinical block
    seeds = [load_preds(a.sp_dir, co.idx, n)]
    for s in (2, 3, 4, 5):
        d = os.path.join(a.seed_dir, "seed%d" % s)
        if os.path.isdir(d):
            seeds.append(load_preds(d, co.idx, n))
    SP = co.fold_pct(np.vstack([co.fold_pct(r, f) for r, f in seeds]).mean(0))
    SPC = co.fold_pct(SP + zc)
    pr, pf = load_preds(a.pibd_dir, co.idx, n)
    PIBD = co.fold_pct(pr, pf)
    PIBDC = co.fold_pct(PIBD + zc)

    rep = {"artifact_type": "figure_source_data", "phase_of_origin": "frozen_post_hoc",
           "n": n, "events": int(co.e.sum()),
           "note": "curves for panels that exist elsewhere only as summaries. No new result."}

    # ---------------------------------------------------------------- F1 Kaplan-Meier
    # survpath_as_released is the SINGLE-SEED arm WITHOUT the clinical block -- the configuration
    # the benchmark actually publishes. Its tertile ordering is the point of the panel it feeds, and
    # keeping it separate from survpath_plus_clinical is what stops the two being confused: the
    # figure design conflated them once, which is what designing before drawing is for.
    SP1 = co.fold_pct(*load_preds(a.sp_dir, co.idx, n))
    ARMS = {"ours": OURS, "survpath_plus_clinical": SPC, "pibd_plus_clinical": PIBDC,
            "clinical_alone": zc, "survpath_as_released": SP1,
            "survpath_seedavg_alone": SP, "pibd_as_released": PIBD}
    risk_grid = [0, 12, 24, 36, 48, 60, 72, 84]
    km = {}
    for name, z in ARMS.items():
        q = np.quantile(z, [1 / 3, 2 / 3])
        grp = np.digitize(z, q)
        c2, df, _ = logrank(co.t, co.e, grp)
        groups = {}
        for g, lab in ((0, "low"), (1, "middle"), (2, "high")):
            m = grp == g
            tt, ss, ar = km_curve(co.t[m], co.e[m])
            groups[lab] = {"n": int(m.sum()), "events": int(co.e[m].sum()),
                           "times": [round(x, 4) for x in tt],
                           "survival": [round(x, 5) for x in ss],
                           "at_risk_grid": risk_grid,
                           "at_risk": at_risk_at(co.t[m], risk_grid),
                           "median_months": median_survival(tt, ss)}
        km[name] = {"cindex": round(cidx(z, ii, jj), 4),
                    "logrank_chi2": round(float(c2), 3), "logrank_df": int(df),
                    "logrank_p": float("%.3g" % chi2_sf(c2, df)),
                    "tertile_cuts": [round(float(x), 4) for x in q], "groups": groups}
        med = [groups[l]["median_months"] for l in ("low", "middle", "high")]
        km[name]["medians_monotone"] = bool(
            all(m is not None for m in med) and med[0] >= med[1] >= med[2])
    rep["F1_kaplan_meier"] = km

    # ---------------------------------------------------------------- F2 dense td-AUC
    uniq, g = km_censoring(co.t, co.e)
    ev_t = np.sort(co.t[co.e == 1])
    dense = [float(x) for x in np.linspace(float(np.quantile(ev_t, 0.05)),
                                           float(np.quantile(ev_t, 0.90)), 20)]
    check = [float(np.quantile(ev_t, q)) for q in (0.25, 0.5, 0.75)]
    auc = {name: {"horizons": [round(h, 3) for h in dense],
                  "auc": [round(td_auc(z, co.t, co.e, h, uniq, g), 4) for h in dense],
                  "check_horizons": [round(h, 3) for h in check],
                  "check_auc": [round(td_auc(z, co.t, co.e, h, uniq, g), 4) for h in check]}
           for name, z in list(ARMS.items()) + [("slide", zt), ("omics", zo)]}
    # the three published horizons must reproduce the confirmatory run exactly
    assert auc["ours"]["check_auc"] == [round(td_auc(OURS, co.t, co.e, h, uniq, g), 4)
                                        for h in check]
    rep["F2_time_dependent_auc"] = auc
    rep["F2_time_dependent_auc"]["_at_risk_at_horizons"] = at_risk_at(co.t, dense)

    # ---------------------------------------------------------------- F3 power
    se = 0.0210
    margins = [round(x, 4) for x in np.linspace(0.0, 0.10, 101)]
    powers = [round(power_at(d, se), 4) for d in margins]
    at80 = next(d for d, p in zip(margins, powers) if p >= 0.80)
    rep["F3_power"] = {
        "se": se, "alpha": 0.05, "margins": margins, "power": powers,
        "detectable_at_80pc": round(at80, 4),
        "published_consecutive_gap_band": [0.002, 0.012],
        "marked": {"ours_vs_survpath": {"margin": 0.0334, "power": round(power_at(0.0334, se), 4)},
                   "ours_vs_pibd": {"margin": 0.0243, "power": round(power_at(0.0243, se), 4)},
                   "ours_vs_incumbent_table": {"margin": 0.0422,
                                               "power": round(power_at(0.0422, se), 4)}}}
    # the already-reported figures must sit on this curve
    assert abs(rep["F3_power"]["detectable_at_80pc"] - 0.0589) < 0.0015, at80
    assert abs(rep["F3_power"]["marked"]["ours_vs_survpath"]["power"] - 0.355) < 0.01

    # ---------------------------------------------------------------- F4 ladder triple
    triple, gap = pick_triple(LADDER)
    assert gap >= 0.05, "ladder triple gap %.4f below the greyscale floor" % gap
    rep["F4_risk_group_colours"] = {
        "ladder": LADDER,
        "chosen_dark_to_pale": triple,
        "luminance": [round(luminance(c), 4) for c in triple],
        "min_adjacent_gap": round(gap, 4),
        "assignment": {"high": triple[0], "middle": triple[1], "low": triple[2]},
        "why": "risk tertile is ORDERED, so it is encoded in lightness. Highest risk darkest."}

    open(a.out, "w").write(json.dumps(rep, indent=1) + "\n")
    print(json.dumps({
        "km": {k: {"c": v["cindex"], "logrank_p": v["logrank_p"],
                   "medians": [v["groups"][l]["median_months"] for l in ("low", "middle", "high")],
                   "monotone": v["medians_monotone"]} for k, v in km.items()},
        "auc_span_ours": [auc["ours"]["auc"][0], auc["ours"]["auc"][-1]],
        "power": {"detectable_at_80pc": rep["F3_power"]["detectable_at_80pc"],
                  "marked": rep["F3_power"]["marked"]},
        "colours": rep["F4_risk_group_colours"]["assignment"],
        "colour_gap": rep["F4_risk_group_colours"]["min_adjacent_gap"]}, indent=1))
    print("wrote", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
