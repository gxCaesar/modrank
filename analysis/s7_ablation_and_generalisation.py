#!/usr/bin/env python3
"""Tier 1A + 1B + 2E + 2J -- the generalisation arm, the ablation table, cost, and sensitivity.

FOUR THINGS, in one pass because they share the same fitting machinery.

A. GENERALISATION. The frozen arm, unchanged, on all five TCGA studies the benchmark covers, with
   seeds 0-4, against each study's published values. A generalisation arm -- a dataset the method
   was not developed on -- appears in 100% of the accepted papers this lab surveyed, and ours was
   until now only run for the STAGE finding, not for the method.

   The clinical block differs by necessity and is labelled as such: bladder uses the stage column
   from DIMAF's released split files, and the other four use SurvPath's own clinical_data stage,
   because DIMAF publishes split files for bladder only.

B. ABLATION. Every subset of the three modalities, seven arms, under the identical estimator and
   combination. Leave-one-out is the row that attributes the result; the singletons say what each
   modality is worth alone.

E. COST. Parameter count, wall-clock fit time and peak resident memory for the whole method. The
   published entrants are attention and prototype networks; ours is three penalised Cox models, and
   the difference is a claim the paper makes rather than an aside.

J. SENSITIVITY to the omics representation: SurvPath's `combine` grouping (the one the frozen arm
   uses, chosen because it is the benchmark's own), against its `hallmarks` and `xena` groupings.
   Xena scores higher and is deliberately NOT the frozen choice -- it was found by scoring three
   groupings against the pooled outcome, which is a selection worth about 0.011 at this task's
   sigma, and promoting it after the fact is the exact move this campaign has been measuring in
   other people's work.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import resource
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, fitapply, pct  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402
from s6_amend_clinical import dimaf_clinical, boot  # noqa: E402
from s5_pancohort import load_cohort, fold_pct  # noqa: E402

PUBLISHED = {
    "blca": {"MOAD-FNet": 0.691, "DIMAF": 0.679, "APL": 0.677, "PIBD": 0.667,
             "DSCASurv": 0.646, "OTSurv": 0.637, "MMP": 0.635, "SurvPath": 0.625},
    "brca": {"APL": 0.794, "DSCASurv": 0.765, "DIMAF": 0.759, "MMP": 0.738,
             "PIBD": 0.736, "SurvPath": 0.655},
    "coadread": {"DSCASurv": 0.832, "APL": 0.812, "PIBD": 0.768, "SurvPath": 0.673,
                 "OTSurv": 0.667, "MMP": 0.630},
    "hnsc": {"DSCASurv": 0.666, "APL": 0.653, "PIBD": 0.640, "SurvPath": 0.600},
    "stad": {"DSCASurv": 0.698, "APL": 0.686, "PIBD": 0.684, "MMP": 0.598,
             "SurvPath": 0.592},
}
VERIFIED_FOLDS = {"SurvPath", "PIBD", "DIMAF"}      # checked here; the rest are quoted


def onehot(vals, drop=("", None, "N/A", "[Not Available]", "[Unknown]", "nan", "None")):
    lv = sorted({v for v in vals if v not in drop})
    return (np.array([[1.0 if v == l else 0.0 for l in lv] for v in vals])
            if lv else np.zeros((len(vals), 0)))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--omics-dir", required=True)
    ap.add_argument("--dimaf-dir", required=True)
    ap.add_argument("--meta-dir", required=True)
    ap.add_argument("--pan-dir", required=True)
    ap.add_argument("--clin-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rep = {"artifact_type": "s7_ablation_and_generalisation",
           "phase_of_origin": "frozen_post_hoc",
           "note": "the frozen arm, unchanged; nothing here re-tunes anything"}

    # ================================================================ B, E, J on bladder
    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    ii, jj = cpairs(co.t, co.e)
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    CLIN = np.hstack([dc["age"], dc["fem"], dc["stage"]])
    omics = {}
    for tag, sig in (("combine", "combine_signatures.csv"),
                     ("hallmarks", "hallmarks_signatures.csv"),
                     ("xena", "xena_signatures.csv")):
        rna = os.path.join(a.omics_dir, "rna_%s.csv" % tag)
        sp = os.path.join(a.omics_dir, sig)
        if os.path.exists(rna) and os.path.exists(sp):
            G, genes, _ = load_rna(rna, co.keep)
            P, _ = pathway_matrix(G, genes, sp)
            omics[tag] = P

    def arm(X, seed=0):
        s = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], seed)
        return co.fold_pct(s)

    t0 = time.time()
    m0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    zt, zo, zc = arm(co.T), arm(omics["combine"]), arm(CLIN)
    fit_s = time.time() - t0
    m1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    parts = {"wsi": zt, "omics": zo, "clin": zc}

    # ---- B: every subset
    import itertools
    abl = {}
    for r in range(1, 4):
        for cmb in itertools.combinations(("wsi", "omics", "clin"), r):
            z = co.fold_pct(sum(parts[k] for k in cmb))
            abl[" + ".join(cmb)] = round(cidx(z, ii, jj), 4)
    full = co.fold_pct(zt + zo + zc)
    rep["B_ablation_bladder"] = {
        "arms": abl,
        "leave_one_out": {
            "drop wsi   (omics + clin)": round(abl["omics + clin"] - abl["wsi + omics + clin"], 4),
            "drop omics (wsi + clin)": round(abl["wsi + clin"] - abl["wsi + omics + clin"], 4),
            "drop clin  (wsi + omics)": round(abl["wsi + omics"] - abl["wsi + omics + clin"], 4)},
        "paired_vs_full": {
            k: boot(co.fold_pct(sum(parts[x] for x in k.split(" + "))), full, co.t, co.e, 3000)
            for k in ("omics + clin", "wsi + clin", "wsi + omics")}}

    # ---- E: cost
    n_par = int(co.T.shape[1] + omics["combine"].shape[1] + CLIN.shape[1])
    rep["E_cost"] = {
        "free_parameters_total": n_par,
        "per_view": {"wsi_titan": int(co.T.shape[1]), "omics_combine": int(omics["combine"].shape[1]),
                     "clinical": int(CLIN.shape[1])},
        "parameters_in_the_combination": 0,
        "fit_seconds_all_three_views_five_folds": round(fit_s, 1),
        "peak_rss_mb": round(max(m0, m1) / (1024.0 * 1024.0 if sys.platform == "darwin" else 1024.0), 1),
        "hardware": "one CPU core-set on a laptop; no GPU is used by any arm of this method",
        "contrast": "the published entrants are attention and prototype networks trained end to "
                    "end; this is three penalised Cox models and the combination has no parameters"}

    # ---- J: omics representation sensitivity
    sens = {}
    for tag, P in omics.items():
        zoo = arm(P)
        sens[tag] = {"omics_alone": round(cidx(zoo, ii, jj), 4),
                     "full_arm_with_this_grouping": round(
                         cidx(co.fold_pct(zt + zoo + zc), ii, jj), 4),
                     "n_pathways": int(P.shape[1])}
    rep["J_omics_sensitivity"] = {
        "groupings": sens,
        "frozen_choice": "combine",
        "why_not_the_best": "xena scores higher and is NOT used. It was found by scoring three "
                            "groupings against the pooled outcome, a selection worth about "
                            "sigma*sqrt(2 ln 3) = 0.011 at this task's sigma, and promoting it "
                            "afterwards would be the move this campaign measures in others' work. "
                            "`combine` is the benchmark's own default."}
    print("bladder: ablation and sensitivity done (%.0fs fit)" % fit_s, file=sys.stderr, flush=True)

    # ================================================================ A: five cohorts
    d0 = __import__("pickle").load(open(a.titan, "rb"))
    E = np.asarray(d0["embeddings"])
    stem2emb = {}
    for i, f in enumerate(d0["filenames"]):
        s = str(f)
        for suf in (".svs", ".h5", ".pt"):
            if s.endswith(suf):
                s = s[:-len(suf)]
        stem2emb[s] = E[i]

    gen = {}
    for c in ("blca", "brca", "coadread", "hnsc", "stad"):
        meta = os.path.join(a.meta_dir, "tcga_%s.csv" % c)
        rna = os.path.join(a.pan_dir, "rna_%s.csv" % c)
        spl = os.path.join(a.pan_dir, "splits", c)
        clin_csv = os.path.join(a.clin_dir, "tcga_%s_clinical.csv" % c)
        if not all(os.path.exists(p) for p in (meta, rna, clin_csv)) or not os.path.isdir(spl):
            gen[c] = {"status": "inputs missing"}
            continue
        cc = load_cohort(meta, spl, rna, os.path.join(a.omics_dir, "combine_signatures.csv"),
                         stem2emb)
        stage_of = {r["case_id"]: str(r.get("stage", "")).strip()
                    for r in csv.DictReader(open(clin_csv))}
        CL = np.hstack([cc["AS"], onehot([stage_of.get(k, "") for k in cc["keep"]])])
        i2, j2 = cpairs(cc["t"], cc["e"])

        def arm2(X, seed):
            s = np.full(len(cc["keep"]), np.nan)
            for tri, vai in cc["fi"]:
                if len(tri) < 30 or not len(vai):
                    continue
                s[vai] = fitapply(X[tri], cc["t"][tri], cc["e"][tri], X[vai], seed)
            return fold_pct(s, cc["assign"])

        per_seed = []
        for sd in range(5):
            z = fold_pct(arm2(cc["T"], sd) + arm2(cc["P"], sd) + arm2(CL, sd), cc["assign"])
            per_seed.append(round(cidx(z, i2, j2), 4))
        pub = PUBLISHED.get(c, {})
        best_any = max(pub.values()) if pub else None
        best_ver = max((v for k, v in pub.items() if k in VERIFIED_FOLDS), default=None)
        gen[c] = {
            "n_cases": len(cc["keep"]), "events": int(cc["e"].sum()),
            "clinical_stage_source": ("SurvPath datasets_csv/clinical_data" if c != "blca"
                                      else "SurvPath datasets_csv/clinical_data (the DIMAF-file "
                                           "variant used for the headline is bladder-only)"),
            "per_seed": per_seed,
            "primary_mean_over_seeds": round(float(np.mean(per_seed)), 4),
            "sd_over_seeds": round(float(np.std(per_seed, ddof=1)), 4),
            "published": pub,
            "best_published_any_protocol": best_any,
            "best_published_verified_folds": best_ver,
            "beats_best_any": bool(best_any is not None and np.mean(per_seed) > best_any),
            "beats_best_verified": bool(best_ver is not None and np.mean(per_seed) > best_ver)}
        print("%-9s n=%4d ev=%3d  ours=%.4f+/-%.4f  best-any=%.3f  best-verified=%.3f"
              % (c, gen[c]["n_cases"], gen[c]["events"], gen[c]["primary_mean_over_seeds"],
                 gen[c]["sd_over_seeds"], best_any or float("nan"), best_ver or float("nan")),
              file=sys.stderr, flush=True)
    rep["A_generalisation"] = {
        "cohorts": gen,
        "beats_best_verified_in": [c for c, v in gen.items() if v.get("beats_best_verified")],
        "beats_best_any_in": [c for c, v in gen.items() if v.get("beats_best_any")],
        "caveat": "the four non-bladder cohorts use SurvPath's own single `stage` column, because "
                  "DIMAF publishes split files for bladder only. They are a generalisation arm for "
                  "the METHOD, not external validation of the bladder number -- there is no "
                  "held-out bladder cohort anywhere in this benchmark."}

    open(a.out, "w").write(json.dumps(rep, indent=1) + "\n")
    print("wrote", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
