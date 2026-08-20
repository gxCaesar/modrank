#!/usr/bin/env python3
"""The pre-specified arm, its three comparisons, and every paired uncertainty.

Pre-freeze, exploratory. NOT reportable -- this is the arm intended for the freeze.

THE ARM, and why this specification is not a selection. One view per modality, and each view is the
one the benchmark's own reference implementation uses:

    WSI       TITAN slide embedding, mean-pooled over a case's slides
    omics     SurvPath's 4,999 genes averaged into ITS OWN 275 `combine` pathway groups -- the
              representation `mahmoodlab/SurvPath` ships and its omics branch consumes
    clinical  age, sex, T, N, M, AJCC stage

Estimator: ridge Cox, alphas (1, 8, 64, 512, 4096) by 3-fold inner CV inside each training fold,
identical for every view. Combination: percentile-normalise within fold, equal weight.

A better-scoring omics view exists -- the Xena grouping reaches 0.6749 against `combine`'s 0.6510 --
and it is NOT used here, because it was found by scoring three groupings against the pooled outcome
and taking the best. That is a selection worth about sigma*sqrt(2 ln 3) = 0.011 at this task's
sigma, and it is reported below as a sensitivity rather than promoted into the arm.

FIVE COMPARISONS, each with a case-level paired bootstrap where a pairing is possible:

  1  vs the cheap clinical Cox                  -- the baseline S5 exists to make you beat
  2  vs the TITAN image arm alone
  3  vs SurvPath rerun + THE SAME clinical      -- input parity, paired, the strongest available
  4  vs DIMAF 0.679 on the same released folds  -- unpaired; the verified incumbent
  5  ours on WSI+omics only vs DIMAF            -- unpaired; whether the win survives without the
                                                   clinical variables

Comparison 5 has already failed under three different combination rules (0.6770 flat, 0.6767
skill-weighted, 0.6776 thresholded, against 0.679). It is recomputed here so the arm's own number
sits next to it rather than in a different file.
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


def boot(a, b, t, e, reps=6000, seed=0):
    rng = np.random.default_rng(seed)
    n = len(t)
    out = []
    for _ in range(reps):
        bs = rng.choice(n, size=n, replace=True)
        ii, jj = cpairs(t[bs], e[bs])
        if ii.size:
            out.append(cidx(a[bs], ii, jj) - cidx(b[bs], ii, jj))
    v = np.asarray(out)
    return {"mean": round(float(v.mean()), 4), "se": round(float(v.std(ddof=1)), 4),
            "ci95": [round(float(np.quantile(v, .025)), 4), round(float(np.quantile(v, .975)), 4)],
            "p_two_sided": round(float(2 * min((v <= 0).mean(), (v >= 0).mean())), 4),
            "reps": int(v.size), "resampled": "cases, not comparable pairs"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--omics-dir", required=True)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, pn = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    Gx, gx, _ = load_rna(os.path.join(a.omics_dir, "rna_xena.csv"), co.keep)
    Px, pnx = pathway_matrix(Gx, gx, os.path.join(a.omics_dir, "xena_signatures.csv"))

    def arm(X):
        s = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], 0)
        return co.fold_pct(s)

    z_ti, z_om, z_cl = arm(co.T), arm(P), arm(co.CLIN)
    z_omx = arm(Px)
    z_sp = co.fold_pct(co.spr)

    OURS = co.fold_pct(z_ti + z_om + z_cl)
    OURS_NOCLIN = co.fold_pct(z_ti + z_om)
    SPC = co.fold_pct(z_sp + z_cl)
    SENS = co.fold_pct(z_ti + z_omx + z_cl)

    ii, jj = cpairs(co.t, co.e)
    gap = np.abs(z_cl[ii] - z_cl[jj])
    tied = gap <= float(np.quantile(gap, 0.10))
    ss = co.site[ii] == co.site[jj]

    def row(z):
        return {"cindex": round(cidx(z, ii, jj), 4),
                "clinically_tied_q10": round(cidx(z, ii[tied], jj[tied]), 4),
                "within_site": round(cidx(z, ii[ss], jj[ss]), 4)}

    rep = {"artifact_type": "s5_headline_prespecified_arm", "reportable": False,
           "phase_of_origin": "exploratory",
           "n": len(co.keep), "events": int(co.e.sum()), "comparable_pairs": int(ii.size),
           "bar_paired_2sd": 0.0145,
           "arm_specification": {
               "views": ["TITAN slide embedding (768-d)",
                         "SurvPath `combine` pathway means (%d-d)" % P.shape[1],
                         "clinical age/sex/T/N/M/stage (%d-d)" % co.CLIN.shape[1]],
               "estimator": "ridge Cox, alphas (1,8,64,512,4096), 3-fold inner CV in-fold",
               "combination": "percentile within fold, equal weight",
               "selection_performed": "none -- one view per modality, the omics grouping taken "
                                      "from the benchmark's own reference implementation"},
           "components": {"titan": row(z_ti), "omics_combine": row(z_om),
                          "clinical": row(z_cl), "survpath_rerun": row(z_sp),
                          "omics_xena_NOT_IN_ARM": row(z_omx)},
           "arms": {"OURS": row(OURS), "OURS_no_clinical": row(OURS_NOCLIN),
                    "survpath_plus_clinical": row(SPC),
                    "SENSITIVITY_xena_instead_of_combine": row(SENS)}}

    rep["comparisons"] = {
        "1_vs_cheap_clinical_cox": {"ours": row(OURS)["cindex"], "baseline": row(z_cl)["cindex"],
                                    "paired": boot(OURS, z_cl, co.t, co.e)},
        "2_vs_titan_alone": {"ours": row(OURS)["cindex"], "baseline": row(z_ti)["cindex"],
                             "paired": boot(OURS, z_ti, co.t, co.e)},
        "3_vs_survpath_plus_same_clinical_INPUT_PARITY": {
            "ours": row(OURS)["cindex"], "baseline": row(SPC)["cindex"],
            "paired": boot(OURS, SPC, co.t, co.e)},
        "4_vs_DIMAF_published_same_folds": {
            "ours": row(OURS)["cindex"], "incumbent": 0.679, "incumbent_sd": 0.043,
            "gap": round(row(OURS)["cindex"] - 0.679, 4),
            "paired": "NOT POSSIBLE -- DIMAF's per-case predictions are not released"},
        "5_ours_without_clinical_vs_DIMAF": {
            "ours_wsi_plus_omics": row(OURS_NOCLIN)["cindex"], "incumbent": 0.679,
            "gap": round(row(OURS_NOCLIN)["cindex"] - 0.679, 4),
            "passes": bool(row(OURS_NOCLIN)["cindex"] - 0.679 > 0.0145),
            "already_failed_under": {"flat_equal_weight_7_views": 0.6770,
                                     "inner_skill_weighted_5_views": 0.6767,
                                     "threshold_0.55_5_views": 0.6776}},
    }

    pf = {}
    for k in range(5):
        m = co.fold == k
        i2, j2 = cpairs(co.t[m], co.e[m])
        pf["fold_%d" % k] = {"n": int(m.sum()), "events": int(co.e[m].sum()),
                             "ours": round(cidx(OURS[m], i2, j2), 4),
                             "survpath_plus_clinical": round(cidx(SPC[m], i2, j2), 4),
                             "clinical": round(cidx(z_cl[m], i2, j2), 4)}
    rep["per_fold"] = pf
    v = [x["ours"] for x in pf.values()]
    rep["ours_fold_mean_sd"] = [round(float(np.mean(v)), 4), round(float(np.std(v, ddof=1)), 4)]
    rep["folds_won_vs_survpath_plus_clinical"] = sum(
        1 for x in pf.values() if x["ours"] > x["survpath_plus_clinical"])

    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
