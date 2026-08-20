#!/usr/bin/env python3
"""AMENDMENT A1 -- the clinical block is rebuilt from the incumbent's own released split files.

WHY THIS EXISTS, and it is a defect in our own arm rather than a refinement of it. The frozen run's
clinical block took T, N, M and AJCC stage from a GDC query cached as `tnm.json`. Diffing DIMAF's
released BLCA folds against ours -- to settle whether DIMAF is a valid comparator, which it is,
5 of 5 validation folds identical at case level -- exposed that DIMAF ships `ajcc_pathologic_tumor_stage`
in those same files, and that it disagrees with ours on 26 of 357 cases with a further 35 missing:
**61 of 359, 17%.**

The cause was found by querying GDC directly for four of the disagreeing cases. They carry MULTIPLE
diagnosis records, and `tnm.json` took the first one:

    TCGA-FD-A6TG   Stage IIB / T2 N0 MX / Adenocarcinoma, NOS          <- what we used
                   Stage IV  / T3a N2 MX / Transitional cell carcinoma <- the bladder tumour
    TCGA-SY-A9G0   Stage IIB / T2 N0 M0  / Adenocarcinoma, NOS         <- what we used
                   Stage IV  / T4 N1 M0  / Transitional cell carcinoma <- the bladder tumour

The tell was there before the GDC query and was not read: `Stage IIA`, `IIB` and `IIC` appear in our
record, and those are not AJCC bladder stages at all. The errors are also SYSTEMATIC rather than
random -- they understate the stage (III -> I, IV -> I), because the first-listed diagnosis tends to
be the earlier or the other tumour.

THE FIX, and it is stronger than the thing it replaces. Every clinical variable is taken from
DIMAF's own released split files: `birth_days_to` for age, `sex`, `ajcc_pathologic_tumor_stage`,
and `histological_grade` for the grade-versus-stage contrast. So the clinical block now comes from
the INCUMBENT'S OWN REPOSITORY. There is no longer any version of the objection "you brought data
the competitors did not have" -- the column is in the file DIMAF reads to build the folds it is
scored on.

WHAT IS REPORTED. Both numbers, the frozen one and the amended one, whichever way the amendment
moves it. An amendment that is only reported when it helps is not an amendment.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, fitapply  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402


def dimaf_clinical(dimaf_dir, keep):
    """age, sex, stage and grade for each case, from DIMAF's own released split files."""
    rec = {}
    for k in range(5):
        for f in ("train", "test"):
            p = os.path.join(dimaf_dir, "%d_%s.csv" % (k, f))
            for r in csv.DictReader(open(p)):
                rec.setdefault(r["case_id"], r)
    have = sum(1 for c in keep if c in rec)

    def num(c, key, default):
        v = (rec.get(c) or {}).get(key, "")
        try:
            return float(v)
        except (TypeError, ValueError):
            return default

    age = np.array([-num(c, "birth_days_to", -65 * 365.25) / 365.25 for c in keep])[:, None]
    fem = np.array([1.0 if (rec.get(c) or {}).get("sex", "") == "F" else 0.0 for c in keep])[:, None]

    def oh(vals, drop=("[Not Available]", "[Unknown]", "", "[Not Evaluated]")):
        lv = sorted({v for v in vals if v not in drop})
        return (np.array([[1.0 if v == l else 0.0 for l in lv] for v in vals])
                if lv else np.zeros((len(vals), 0))), lv

    stage, slv = oh([(rec.get(c) or {}).get("ajcc_pathologic_tumor_stage", "") for c in keep])
    grade, glv = oh([(rec.get(c) or {}).get("histological_grade", "") for c in keep])
    return {"age": age, "fem": fem, "stage": stage, "grade": grade,
            "stage_levels": slv, "grade_levels": glv, "cases_matched": have}


def boot(x, y, t, e, reps=6000, seed=0):
    rng = np.random.default_rng(seed)
    n = len(t)
    v = []
    for _ in range(reps):
        bs = rng.choice(n, size=n, replace=True)
        ii, jj = cpairs(t[bs], e[bs])
        if ii.size:
            v.append(cidx(x[bs], ii, jj) - cidx(y[bs], ii, jj))
    v = np.asarray(v)
    return {"mean": round(float(v.mean()), 4), "se": round(float(v.std(ddof=1)), 4),
            "ci95": [round(float(np.quantile(v, .025)), 4),
                     round(float(np.quantile(v, .975)), 4)],
            "p_two_sided": round(float(2 * min((v <= 0).mean(), (v >= 0).mean())), 4)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--protocol", required=True)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--omics-dir", required=True)
    ap.add_argument("--dimaf-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    proto = json.load(open(a.protocol))
    seeds = proto["method"]["seeds"]["values"]
    alphas = tuple(proto["method"]["estimator"]["alpha_grid"])
    target = 0.679 + 0.0239

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    CLIN_NEW = np.hstack([dc["age"], dc["fem"], dc["stage"]])
    CLIN_GRADE = np.hstack([dc["age"], dc["fem"], dc["grade"]])
    ii, jj = cpairs(co.t, co.e)

    def arm(X, seed):
        s = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], seed, alphas)
        return co.fold_pct(s)

    rep = {"artifact_type": "s6_amendment_A1_clinical_provenance",
           "amends": "development/benchmark-protocol.json",
           "executed_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "phase_of_origin": "frozen_amended",
           "defect": {
             "what": "the frozen run's clinical block took stage and T/N/M from a cached GDC query "
                     "that selected the FIRST diagnosis record; cases with several diagnoses got "
                     "the wrong tumour",
             "cases_disagreeing_with_the_incumbents_own_file": 26,
             "cases_missing_a_stage_in_our_record": 35,
             "of_total": 359,
             "direction": "systematic understatement of stage (III->I, IV->I)",
             "tell_that_was_present_and_unread": "Stage IIA/IIB/IIC appear in our record and are "
                                                 "not AJCC bladder stages"},
           "fix": "every clinical variable now comes from DIMAF's own released split files: "
                  "birth_days_to, sex, ajcc_pathologic_tumor_stage, histological_grade",
           "cases_matched_in_dimaf_files": dc["cases_matched"],
           "stage_levels": dc["stage_levels"], "grade_levels": dc["grade_levels"]}

    per_seed = []
    for sd in seeds:
        zt, zo, zc = arm(co.T, sd), arm(P, sd), arm(CLIN_NEW, sd)
        ours = co.fold_pct(zt + zo + zc)
        per_seed.append({"seed": sd, "OURS": round(cidx(ours, ii, jj), 4),
                         "wsi_titan": round(cidx(zt, ii, jj), 4),
                         "omics_combine": round(cidx(zo, ii, jj), 4),
                         "clinical": round(cidx(zc, ii, jj), 4)})
        print("seed %d  OURS=%.4f  clinical=%.4f" % (sd, per_seed[-1]["OURS"],
                                                     per_seed[-1]["clinical"]),
              file=sys.stderr, flush=True)
    vals = [r["OURS"] for r in per_seed]
    primary, primary_sd = float(np.mean(vals)), float(np.std(vals, ddof=1))

    zt0, zo0, zc0 = arm(co.T, 0), arm(P, 0), arm(CLIN_NEW, 0)
    zg0 = arm(CLIN_GRADE, 0)
    ours0 = co.fold_pct(zt0 + zo0 + zc0)
    spc0 = co.fold_pct(co.fold_pct(co.spr) + zc0)
    ours_grade = co.fold_pct(zt0 + zo0 + zg0)

    rep["per_seed"] = per_seed
    rep["primary"] = {"value": round(primary, 4), "sd_over_seeds": round(primary_sd, 4),
                      "seeds": seeds}
    rep["frozen_run_for_comparison"] = {"primary": 0.7260, "sd_over_seeds": 0.0043,
                                        "clinical_arm": 0.6856}
    rep["arms_seed0"] = {
        "clinical_age_sex_stage_from_DIMAF_file": round(cidx(zc0, ii, jj), 4),
        "clinical_age_sex_GRADE_from_DIMAF_file": round(cidx(zg0, ii, jj), 4),
        "stage_minus_grade": round(cidx(zc0, ii, jj) - cidx(zg0, ii, jj), 4),
        "wsi_titan": round(cidx(zt0, ii, jj), 4),
        "omics_combine": round(cidx(zo0, ii, jj), 4),
        "OURS": round(cidx(ours0, ii, jj), 4),
        "OURS_with_grade_instead_of_stage": round(cidx(ours_grade, ii, jj), 4),
        "survpath_plus_same_clinical": round(cidx(spc0, ii, jj), 4)}
    rep["comparisons"] = {
        "vs_clinical_alone": boot(ours0, zc0, co.t, co.e),
        "vs_titan_alone": boot(ours0, zt0, co.t, co.e),
        "vs_survpath_plus_same_clinical": boot(ours0, spc0, co.t, co.e),
        "vs_DIMAF_published": {"ours_primary": round(primary, 4), "incumbent": 0.679,
                               "gap": round(primary - 0.679, 4),
                               "fold_identity": "SETTLED -- 5 of 5 validation folds identical at "
                                                "case level, 359-case universe identical"}}
    rep["decision"] = {
        "rule": "primary > %.4f AND sd over seeds <= 0.0145" % target,
        "target": round(target, 4), "primary": round(primary, 4),
        "primary_clears_target": bool(primary > target),
        "seed_sd": round(primary_sd, 4), "seed_sd_within_bar": bool(primary_sd <= 0.0145),
        "VERDICT": "PASS" if (primary > target and primary_sd <= 0.0145) else "MISS",
        "moved_from_frozen_run": round(primary - 0.7260, 4)}

    os.makedirs(a.out_dir, exist_ok=True)
    out = os.path.join(a.out_dir, "amendment-A1-clinical-provenance.json")
    open(out, "w").write(json.dumps(rep, indent=1) + "\n")
    print(json.dumps({"arms": rep["arms_seed0"], "decision": rep["decision"]}, indent=1))
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
