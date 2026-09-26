#!/usr/bin/env python3
"""The five-study numbers, with the intervals they have never carried.

WHY THIS EXISTS. `s5_stage_five_cohorts.py` produced, for each of the benchmark's five TCGA studies,
the three-modality arm and the same arm with grade substituted for stage. Every one of those numbers
is a point estimate. The manuscript reports intervals for every other quantity it states, and an
external review on 2026-09-12 recommended desk rejection partly because the method has no
independent validation with all three modalities. Promoting these five cohorts from a replication of
the stage finding to a validation of the method is a change to the pre-registered boundary
(charter.md: "The four non-bladder TCGA cohorts serve as replication of the stage finding, not as
external validation of the bladder number"), and it may not be made on point estimates.

Two of these cohorts are small in the way that matters: BRCA carries 60 events and COADREAD 37. A
+0.0910 difference on 37 events is not evidence until its interval is known, and if the interval
covers zero then the honest reading is a consistency check rather than a validation. This script
settles that either way, before any manuscript text is written.

WHAT IT MEASURES, per cohort, each by the same paired case-level bootstrap the confirmatory run
uses (6,000 resamples of CASES, comparable pairs rebuilt inside every resample, seed 0):

  ours_minus_grade_control   the added-value inflation, replicated: the three-modality arm against
                             the identical arm with grade in place of stage
  ours_minus_clinical_stage  what the two molecular modalities add over the corrected clinical
                             reference. This is the quantity a validation claim rests on
  ours_minus_wsi_omics       what the clinical block adds to the two molecular arms
  stage_minus_grade          the one-column swap in the clinical arm alone, which the manuscript
                             already reports as a point estimate in all five studies

WHAT IT DOES NOT DO. It changes no arm, no fold, no penalty grid and no input. Every score vector is
built by the same functions `s5_stage_five_cohorts` used, imported rather than copied, and the run
refuses to write its output unless it first reproduces that script's committed point estimates for
all five cohorts. A bootstrap around numbers that have drifted would be worse than no bootstrap.

THE COMPARISON THIS DOES NOT LICENSE. Beating `published_best_verified_released_folds` (SurvPath's
own released-fold value per cohort) is not the same as beating the best published value under any
protocol: on BRCA that is APL at 0.794 and on COADREAD DSCASurv at 0.832, both far above these arms.
The output carries both reference values per cohort so the narrower claim cannot be widened by
accident.

Post-freeze, exploratory. Reportable only if the charter amendment is recorded.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from blca_common import cidx, cpairs                                          # noqa: E402
from s5_stage_five_cohorts import (COHORTS, PUBLISHED_BEST,                   # noqa: E402
                                   PUBLISHED_RELEASED_FOLD_BEST, arm,
                                   fold_pct, load_cohort)

REPS = 6000


def boot(x, y, t, e, reps=REPS, seed=0):
    """Copied verbatim from s6_confirmatory.py, which is what every interval in the paper uses.

    Cases are resampled, never comparable pairs, and the pairs are rebuilt from the resampled
    survival times inside each replicate. Reusing the original pairs would make every interval too
    narrow, which is the one way this function can be wrong without looking wrong.
    """
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
            "p_two_sided": round(float(2 * min((v <= 0).mean(), (v >= 0).mean())), 4),
            "reps": int(v.size), "resampled": "cases, not comparable pairs"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("titan", "meta-dir", "clin-dir", "pan-dir", "sig", "frozen", "out"):
        ap.add_argument("--" + f, required=True)
    # A second encoder answers the question this run otherwise invites: whether the validation is a
    # property of the method or of TITAN. With --slide-npz the slide arm is rebuilt from a different
    # pretrained encoder and everything else is held fixed, so a conclusion that survives the swap is
    # about the recipe. Known answers are skipped in that mode, since they were measured on TITAN.
    ap.add_argument("--slide-npz", default="", help="alternative encoder, npz with names and X")
    a = ap.parse_args()
    t0 = time.time()

    frozen = json.load(open(a.frozen))["cohorts"]

    if a.slide_npz:
        z = np.load(a.slide_npz, allow_pickle=True)
        names, X = [str(x) for x in z["names"]], np.asarray(z["X"])
    else:
        d0 = pickle.load(open(a.titan, "rb"))
        names, X = [str(x) for x in d0["filenames"]], np.asarray(d0["embeddings"])
    stem2emb = {}
    for i, f in enumerate(names):
        s = f
        for suf in (".svs", ".h5", ".pt", ".ndpi"):
            if s.endswith(suf):
                s = s[:-len(suf)]
        stem2emb[s] = X[i]

    out = {"artifact_type": "s27_five_cohort_intervals",
           "phase_of_origin": "post_freeze_2026-09-12",
           "reportable": "only with the charter amendment recorded",
           "amends": "charter.md, which scoped the four non-bladder cohorts to replication of the "
                     "stage finding rather than validation of the method",
           "bootstrap": {"replicates": REPS, "unit": "cases", "seed": 0,
                         "source": "identical to s6_confirmatory.boot"},
           "slide_encoder": os.path.basename(a.slide_npz) if a.slide_npz else os.path.basename(a.titan),
           "slide_encoder_is_the_primary": not a.slide_npz,
           "slides_available": len(stem2emb),
           "cohorts": {}}

    drift = {}
    for c in COHORTS:
        paths = (os.path.join(a.meta_dir, "tcga_%s.csv" % c),
                 os.path.join(a.clin_dir, "tcga_%s_clinical.csv" % c),
                 os.path.join(a.pan_dir, "splits", c),
                 os.path.join(a.pan_dir, "rna_%s.csv" % c))
        if not all(os.path.exists(p) for p in paths):
            out["cohorts"][c] = {"status": "inputs missing"}
            continue
        co = load_cohort(paths[0], paths[1], paths[2], paths[3], a.sig, stem2emb)
        ii, jj = cpairs(co["t"], co["e"])

        AS = np.hstack([co["age"], co["fem"]])
        ASG = np.hstack([AS, co["grade"]])
        ASS = np.hstack([AS, co["stage"]])
        z_ass, z_asg = arm(ASS, co), arm(ASG, co)
        z_tit, z_omi = arm(co["T"], co), arm(co["P"], co)
        wsi_om = fold_pct(z_tit + z_omi, co["assign"])
        ours = fold_pct(z_tit + z_omi + z_ass, co["assign"]) if z_ass is not None else None
        ctrl = fold_pct(z_tit + z_omi + z_asg, co["assign"]) if z_asg is not None else None

        def sc(v):
            return None if v is None else round(cidx(v, ii, jj), 4)

        row = {"n_cases": len(co["keep"]), "events": int(co["e"].sum()),
               "comparable_pairs": int(ii.size),
               "points": {"OURS_wsi_omics_age_sex_stage": sc(ours),
                          "control_same_arm_with_GRADE_instead": sc(ctrl),
                          "age_sex_stage": sc(z_ass), "age_sex_grade": sc(z_asg),
                          "WSI_plus_OMICS": sc(wsi_om)},
               "published_best_verified_released_folds": PUBLISHED_RELEASED_FOLD_BEST.get(c),
               "published_best_any_protocol": PUBLISHED_BEST.get(c)}

        # known answers first: a bootstrap around a drifted point estimate is worse than none. They
        # do not apply when the slide arm has been deliberately replaced.
        if not a.slide_npz:
            for key in ("OURS_wsi_omics_age_sex_stage", "control_same_arm_with_GRADE_instead"):
                want = frozen.get(c, {}).get(key)
                got = row["points"][key]
                if want is not None and got is not None and abs(got - want) > 5e-5:
                    drift.setdefault(c, {})[key] = {"got": got, "want": want}

        row["differences"] = {}
        if ours is not None and ctrl is not None:
            row["differences"]["ours_minus_grade_control"] = boot(ours, ctrl, co["t"], co["e"])
        if ours is not None and z_ass is not None:
            row["differences"]["ours_minus_clinical_stage"] = boot(ours, z_ass, co["t"], co["e"])
        if ours is not None:
            row["differences"]["ours_minus_wsi_omics"] = boot(ours, wsi_om, co["t"], co["e"])
        if z_ass is not None and z_asg is not None:
            row["differences"]["stage_minus_grade"] = boot(z_ass, z_asg, co["t"], co["e"])

        out["cohorts"][c] = row
        d = row["differences"]
        print("%-9s n=%4d ev=%3d | OURS %.4f vs grade-ctrl %.4f  d=%+.4f %s p=%.3f | "
              "vs stage-arm d=%+.4f %s p=%.3f  [%.0f s]"
              % (c, row["n_cases"], row["events"], row["points"]["OURS_wsi_omics_age_sex_stage"],
                 row["points"]["control_same_arm_with_GRADE_instead"],
                 d["ours_minus_grade_control"]["mean"], d["ours_minus_grade_control"]["ci95"],
                 d["ours_minus_grade_control"]["p_two_sided"],
                 d["ours_minus_clinical_stage"]["mean"], d["ours_minus_clinical_stage"]["ci95"],
                 d["ours_minus_clinical_stage"]["p_two_sided"], time.time() - t0),
              file=sys.stderr, flush=True)

    scored = [v for v in out["cohorts"].values() if v.get("differences")]
    out["summary"] = {
        "cohorts_scored": len(scored),
        "total_cases": sum(v["n_cases"] for v in scored),
        "total_events": sum(v["events"] for v in scored),
        "cohorts_where_ours_beats_the_grade_control_with_an_interval_excluding_zero": [
            k for k, v in out["cohorts"].items()
            if v.get("differences", {}).get("ours_minus_grade_control", {}).get("ci95", [0, 0])[0] > 0],
        "cohorts_where_ours_beats_the_corrected_clinical_arm_with_an_interval_excluding_zero": [
            k for k, v in out["cohorts"].items()
            if v.get("differences", {}).get("ours_minus_clinical_stage", {}).get("ci95", [0, 0])[0] > 0],
        "reading": "the second list is what a validation claim rests on. A cohort absent from it "
                   "shows the method not separable from the corrected clinical reference there, "
                   "which is a consistency check and not a validation, and must be reported as one",
    }
    out["runtime_seconds"] = round(time.time() - t0, 1)
    if drift:
        out["known_answer_drift"] = drift
        print("KNOWN ANSWER DRIFT, refusing to write: %s" % drift, file=sys.stderr)
        return 2
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1)
    print("wrote %s" % a.out, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
