#!/usr/bin/env python3
"""Tier 2F + 3H + 3I -- calibration, interpretability, and case studies under a frozen rule.

=========================== THE CASE-STUDY SELECTION RULE, FIXED BEFORE ANY OUTPUT IS INSPECTED

Written here, in the file, before the run. Every case the rule admits is reported, whatever it
shows -- a narrative fitted to whichever examples happened to look good is the easiest defect for a
reviewer to find and the hardest to answer.

  POPULATION   patients in the MIDDLE tertile of the clinical arm's out-of-fold risk. That is the
               stage-ambiguous band: the clinical model has placed them neither clearly high nor
               clearly low, so it is where an image-and-molecule model is being asked to add
               something.
  SELECTION    within that band, the 5 patients with the HIGHEST and the 5 with the LOWEST
               out-of-fold score from the full method. Ten cases, no discretion.
  REPORTED     for each: the three component scores, the clinical stage, the observed time and
               whether it was an event. Nothing is excluded for being inconvenient.
  READ AS      if the rule works, the five high-score patients should have shorter observed times
               and more events than the five low-score ones. If they do not, that is the result.

The rule uses no outcome, so it can be applied before the outcomes are looked at, which is what
makes it a rule rather than a description of what was found.

=========================== CALIBRATION

Predicted versus observed survival by risk group. The method emits a RANK, not a hazard, so a
calibration curve needs the rank turned into a survival probability; that is done with the
training-fold Kaplan-Meier baseline raised to a power monotone in the rank, and the resulting
prediction is compared with the group's own Kaplan-Meier estimate at three horizons. What this can
show is miscalibration in the ordering of groups; what it cannot show is absolute calibration of a
model that was never fitted to produce one, and that limit is stated rather than glossed.

=========================== INTERPRETABILITY

Three questions a reader will ask, answered without a saliency map:
  1. which stage levels carry the clinical arm's risk, as fitted coefficients;
  2. which pathways carry the omics arm's, as the largest coefficients averaged over folds --
     reported with the caveat that a ridge coefficient in a correlated 275-dimensional block is
     not a claim about a pathway in isolation;
  3. what the IMAGE arm is, in terms a clinician can check: its rank correlation with stage,
     grade, age and sex. A slide score that is merely a stage detector would show it here.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, cox_fit, fitapply  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402
from s6_amend_clinical import dimaf_clinical  # noqa: E402
from s7_survival_metrics import km_censoring, g_at  # noqa: E402


def spearman(a, b):
    from blca_common import ranks            # average ranks for ties under amendment A2
    ra, rb = ranks(a), ranks(b)
    ra -= ra.mean(); rb -= rb.mean()
    d = float(np.sqrt((ra @ ra) * (rb @ rb)))
    return float(ra @ rb / d) if d else float("nan")


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
    ii, jj = cpairs(co.t, co.e)
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, pnames = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    CLIN = np.hstack([dc["age"], dc["fem"], dc["stage"]])
    clin_names = ["age", "female"] + ["stage=" + s for s in dc["stage_levels"]]

    def arm(X, seed=0):
        s = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], seed)
        return co.fold_pct(s)

    zt, zo, zc = arm(co.T), arm(P), arm(CLIN)
    OURS = co.fold_pct(zt + zo + zc)
    rep = {"artifact_type": "s7_interpret_calibrate_cases", "phase_of_origin": "frozen_post_hoc",
           "n": len(co.keep), "events": int(co.e.sum())}

    # ============================================================ H1  clinical coefficients
    def mean_coef(X, names):
        acc = []
        for tri, _ in fi:
            if len(tri) < 30:
                continue
            mu, sd = X[tri].mean(0), X[tri].std(0) + 1e-9
            acc.append(cox_fit((X[tri] - mu) / sd, co.t[tri], co.e[tri], 64.0))
        b = np.mean(acc, axis=0)
        s = np.std(acc, axis=0, ddof=1)
        return [{"name": n, "coef": round(float(x), 4), "sd_over_folds": round(float(y), 4)}
                for n, x, y in zip(names, b, s)]

    rep["H1_clinical_coefficients"] = {
        "note": "standardised ridge Cox coefficients at alpha=64, averaged over the five training "
                "folds, with the spread across folds. Positive is higher risk.",
        "coefficients": sorted(mean_coef(CLIN, clin_names), key=lambda r: -abs(r["coef"]))}

    # ============================================================ H2  omics coefficients
    oc = mean_coef(P, pnames)
    oc.sort(key=lambda r: -abs(r["coef"]))
    rep["H2_omics_top_pathways"] = {
        "note": "a ridge coefficient inside a correlated 275-dimensional block is not a claim "
                "about a pathway in isolation; neighbouring pathways share genes and the penalty "
                "spreads weight across them. Reported as what the model leans on, not as a "
                "discovery.",
        "n_pathways": len(pnames),
        "top_15_by_absolute_coefficient": oc[:15]}

    # ============================================================ H3  what the image arm IS
    stage_ord = {s: k for k, s in enumerate(dc["stage_levels"])}
    stage_num = np.array([np.argmax(r) if r.sum() else -1 for r in dc["stage"]], dtype=float)
    grade_num = np.array([np.argmax(r) if r.sum() else -1 for r in dc["grade"]], dtype=float)
    have_s, have_g = stage_num >= 0, grade_num >= 0
    rep["H3_what_the_image_arm_is"] = {
        "spearman_with_stage": round(spearman(zt[have_s], stage_num[have_s]), 4),
        "spearman_with_grade": round(spearman(zt[have_g], grade_num[have_g]), 4),
        "spearman_with_age": round(spearman(zt, dc["age"].ravel()), 4),
        "spearman_with_sex": round(spearman(zt, dc["fem"].ravel()), 4),
        "spearman_with_the_clinical_arm": round(spearman(zt, zc), 4),
        "spearman_with_the_omics_arm": round(spearman(zt, zo), 4),
        "reading": "a slide score that were merely a stage detector would correlate strongly with "
                   "stage and with the clinical arm. Its concordance restricted to pairs the "
                   "clinical arm cannot separate is the direct version of the same question."}

    # ============================================================ F  calibration
    ev_t = np.sort(co.t[co.e == 1])
    horizons = [float(np.quantile(ev_t, q)) for q in (0.25, 0.5, 0.75)]
    q = np.quantile(OURS, [1 / 3, 2 / 3])
    grp = np.digitize(OURS, q)
    kt, ks = km_censoring(co.t, 1.0 - co.e)          # cohort event-free survival
    cal = []
    for k in range(3):
        m = grp == k
        gt, gs = km_censoring(co.t[m], 1.0 - co.e[m])
        row = {"group": ["low", "mid", "high"][k], "n": int(m.sum()),
               "events": int(co.e[m].sum())}
        for h in horizons:
            obs = float(g_at(gt, gs, np.array([h]))[0])
            base = float(g_at(kt, ks, np.array([h]))[0])
            pred = float(np.mean(np.clip(base ** np.exp(2.0 * (OURS[m] - 0.5)), 1e-6, 1 - 1e-6)))
            row["h%dm" % round(h)] = {"observed_KM": round(obs, 3), "predicted": round(pred, 3),
                                      "difference": round(pred - obs, 3)}
        cal.append(row)
    rep["F_calibration"] = {
        "horizons_months": [round(h, 1) for h in horizons],
        "by_risk_tertile": cal,
        "limit": "the method emits a RANK, not a hazard. The predicted column turns that rank into "
                 "a survival probability through the cohort baseline, so this shows whether the "
                 "GROUP ORDERING is calibrated; it is not a claim about absolute calibration of a "
                 "model never fitted to produce one."}

    # ============================================================ I  case studies, frozen rule
    qc = np.quantile(zc, [1 / 3, 2 / 3])
    band = (zc > qc[0]) & (zc <= qc[1])              # middle clinical tertile
    idx = np.flatnonzero(band)
    order = idx[np.argsort(OURS[idx])]
    picked = list(order[-5:][::-1]) + list(order[:5])
    stage_lab = {k: s for s, k in stage_ord.items()}
    cases = []
    for i in picked:
        sj = int(stage_num[i]) if stage_num[i] >= 0 else None
        cases.append({
            "case_id": co.keep[i],
            "selected_as": "high" if i in order[-5:] else "low",
            "our_score_percentile": round(float(OURS[i]), 3),
            "image": round(float(zt[i]), 3), "omics": round(float(zo[i]), 3),
            "clinical": round(float(zc[i]), 3),
            "stage": stage_lab.get(sj), "age": round(float(dc["age"][i, 0]), 1),
            "observed_months": round(float(co.t[i]), 1),
            "event": bool(co.e[i] == 1)})
    hi = [c for c in cases if c["selected_as"] == "high"]
    lo = [c for c in cases if c["selected_as"] == "low"]
    rep["I_case_studies"] = {
        "selection_rule": "FIXED BEFORE THE RUN, in the file's docstring: within the MIDDLE tertile "
                          "of the clinical arm's out-of-fold risk, the 5 highest and 5 lowest "
                          "scores from the full method. All ten are reported.",
        "band_size": int(band.sum()),
        "cases": cases,
        "summary": {
            "high_group": {"events": sum(c["event"] for c in hi),
                           "median_observed_months": round(float(np.median([c["observed_months"] for c in hi])), 1)},
            "low_group": {"events": sum(c["event"] for c in lo),
                          "median_observed_months": round(float(np.median([c["observed_months"] for c in lo])), 1)}},
        "read_as": "if the rule works the high group should carry more events and shorter observed "
                   "times. Reported whichever way it lands."}

    open(a.out, "w").write(json.dumps(rep, indent=1) + "\n")
    print(json.dumps({"H3": rep["H3_what_the_image_arm_is"],
                      "cases_summary": rep["I_case_studies"]["summary"],
                      "calibration": rep["F_calibration"]["by_risk_tertile"]}, indent=1))
    print("wrote", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
