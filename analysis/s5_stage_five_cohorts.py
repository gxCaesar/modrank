#!/usr/bin/env python3
"""The stage variable was in the benchmark's own release all along. What does using it do?

Pre-freeze, exploratory. NOT reportable.

THE FACT THIS TESTS, and it is a fact about the field rather than about us. A full-text pass over
the eleven primary papers of this benchmark family -- SurvPath, MMP, PIBD, DIMAF, MCAT, MOTCat,
OTSurv, APL, MOAD-FNet, DSCASurv, ProtoPathway -- found **no clinical baseline anywhere that uses
pathologic stage, overall stage, or T/N/M**. The only clinical covariates reported are age, sex and
cancer grade, and the resulting baselines are weak enough to be suspicious on their face:

    BLCA      age 0.578 · sex 0.489 · grade 0.515 · age+sex+grade 0.570
    COADREAD  age 0.357 · sex 0.542 · age+sex+grade 0.655
    HNSC      age 0.517 · sex 0.486 · grade 0.547 · age+sex+grade 0.512
    STAD      age 0.499 · sex 0.529 · grade 0.552 · age+sex+grade 0.592
    BRCA      age 0.496 · sex 0.490 · grade 0.597 · age+sex+grade 0.563
                            (SurvPath Supplementary Table 3; MMP Table 6; DIMAF Table 1)

A C-index of 0.357 for age in colorectal cancer is not a property of age.

AND THE DATA WAS NEVER MISSING. `mahmoodlab/SurvPath` ships `datasets_csv/clinical_data/
tcga_<cohort>_clinical.csv` for all five studies, and its columns are `case_id, stage, grade,
subtype`. Stage is in the benchmark's own release, next to the grade every paper used instead.

So this run is not "we brought extra data". It is the same repository, the same folds, the same
cases, and one different column. Five arms per cohort, and the pair that matters is the third and
the fourth:

    age+sex                     what our cross-cohort arm had before today
    age+sex+grade               THE PUBLISHED BASELINE, rebuilt from the same file
    age+sex+stage               the column nobody used
    WSI+omics                   our method without any clinical variable
    WSI+omics+age+sex+stage     our method

Rebuilding the published baseline from the same file is what makes the comparison honest: if our
age+sex+grade lands near the published age+sex+grade, the pipeline agrees with theirs and the
stage column is the only thing that differs. If it does not, that is the finding instead, and it is
reported.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import cidx, cpairs, fitapply, pct  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402

COHORTS = ("blca", "brca", "coadread", "hnsc", "stad")
PUBLISHED_AGE_SEX_GRADE = {"blca": 0.570, "brca": 0.563, "coadread": 0.655,
                           "hnsc": 0.512, "stad": 0.592}
PUBLISHED_BEST = {"blca": ("DIMAF", 0.679), "brca": ("APL", 0.794),
                  "coadread": ("DSCASurv", 0.832), "hnsc": ("DSCASurv", 0.666),
                  "stad": ("DSCASurv", 0.698)}
PUBLISHED_RELEASED_FOLD_BEST = {"blca": ("SurvPath", 0.625), "brca": ("SurvPath", 0.655),
                                "coadread": ("SurvPath", 0.673), "hnsc": ("SurvPath", 0.600),
                                "stad": ("SurvPath", 0.592)}


def onehot(vals):
    lv = sorted({v for v in vals if v not in ("", "N/A", "nan", "None")})
    if not lv:
        return np.zeros((len(vals), 0)), []
    return (np.array([[1.0 if v == l else 0.0 for l in lv] for v in vals]), lv)


def load_cohort(meta_csv, clin_csv, splits_dir, rna_csv, sig_csv, stem2emb):
    per, lab = collections.defaultdict(list), {}
    for r in csv.DictReader(open(meta_csv)):
        stem = r["slide_id"].replace(".svs", "")
        if stem in stem2emb:
            per[r["case_id"]].append(stem2emb[stem])
        if r.get("survival_months_dss", "").strip() and r.get("censorship_dss", "").strip():
            lab[r["case_id"]] = (float(r["survival_months_dss"]),
                                 1.0 - float(r["censorship_dss"]),
                                 float(r["age"]) if r.get("age", "").strip() else 65.0,
                                 float(r["is_female"]) if r.get("is_female", "").strip() else 0.0)
    clin = {}
    for r in csv.DictReader(open(clin_csv)):
        cid = r.get("case_id", "").strip()
        if cid:
            clin[cid] = (str(r.get("stage", "")).strip(), str(r.get("grade", "")).strip())

    folds = []
    for k in range(5):
        tr, va = set(), set()
        for r in csv.DictReader(open(os.path.join(splits_dir, "splits_%d.csv" % k))):
            if r.get("train", "").strip():
                tr.add(r["train"].strip())
            if r.get("val", "").strip():
                va.add(r["val"].strip())
        folds.append((tr, va))

    keep = sorted(c for c in per if c in lab)
    idx = {c: i for i, c in enumerate(keep)}
    T = np.vstack([np.mean(per[c], axis=0) for c in keep]).astype(float)
    t = np.array([lab[c][0] for c in keep])
    e = np.array([lab[c][1] for c in keep])
    age = np.array([lab[c][2] for c in keep])[:, None]
    fem = np.array([lab[c][3] for c in keep])[:, None]
    stg, stg_lv = onehot([clin.get(c, ("", ""))[0] for c in keep])
    grd, grd_lv = onehot([clin.get(c, ("", ""))[1] for c in keep])
    G, genes, _ = load_rna(rna_csv, keep)
    P, _ = pathway_matrix(G, genes, sig_csv)

    fi, assign = [], np.full(len(keep), -1)
    for k, (tr, va) in enumerate(folds):
        tri = np.array(sorted(idx[c] for c in tr if c in idx))
        vai = np.array(sorted(idx[c] for c in va if c in idx))
        fi.append((tri, vai))
        for i in vai:
            assign[i] = k
    have_clin = sum(1 for c in keep if c in clin)
    return {"keep": keep, "T": T, "P": P, "age": age, "fem": fem, "stage": stg, "grade": grd,
            "stage_levels": stg_lv, "grade_levels": grd_lv, "t": t, "e": e, "fi": fi,
            "assign": assign, "cases_with_clinical_row": have_clin}


def fold_pct(v, assign):
    out = np.full(len(v), np.nan)
    for k in np.unique(assign):
        if k < 0:
            continue
        m = assign == k
        out[m] = pct(v[m])
    return out


def arm(X, d):
    if X.shape[1] == 0:
        return None
    s = np.full(len(d["keep"]), np.nan)
    for tri, vai in d["fi"]:
        if len(tri) < 30 or not len(vai):
            continue
        s[vai] = fitapply(X[tri], d["t"][tri], d["e"][tri], X[vai], 0)
    return fold_pct(s, d["assign"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--meta-dir", required=True)
    ap.add_argument("--clin-dir", required=True)
    ap.add_argument("--pan-dir", required=True)
    ap.add_argument("--sig", required=True)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    d0 = pickle.load(open(a.titan, "rb"))
    E = np.asarray(d0["embeddings"])
    stem2emb = {}
    for i, f in enumerate(d0["filenames"]):
        s = str(f)
        for suf in (".svs", ".h5", ".pt"):
            if s.endswith(suf):
                s = s[:-len(suf)]
        stem2emb[s] = E[i]

    rep = {"artifact_type": "s5_stage_five_cohorts", "reportable": False,
           "phase_of_origin": "exploratory",
           "data_provenance": "every column comes from mahmoodlab/SurvPath's own release: "
                              "datasets_csv/metadata/tcga_<c>.csv for labels, age and sex; "
                              "datasets_csv/clinical_data/tcga_<c>_clinical.csv for stage and "
                              "grade; datasets_csv/raw_rna_data/combine/<c>/rna_clean.csv for "
                              "genes; PIBD's copy of the released 5-fold case-ID splits",
           "cohorts": {}}

    for c in COHORTS:
        paths = (os.path.join(a.meta_dir, "tcga_%s.csv" % c),
                 os.path.join(a.clin_dir, "tcga_%s_clinical.csv" % c),
                 os.path.join(a.pan_dir, "splits", c),
                 os.path.join(a.pan_dir, "rna_%s.csv" % c))
        if not all(os.path.exists(p) for p in paths):
            rep["cohorts"][c] = {"status": "inputs missing"}
            continue
        co = load_cohort(paths[0], paths[1], paths[2], paths[3], a.sig, stem2emb)
        ii, jj = cpairs(co["t"], co["e"])

        AS = np.hstack([co["age"], co["fem"]])
        ASG = np.hstack([AS, co["grade"]])
        ASS = np.hstack([AS, co["stage"]])
        z = {"age_sex": arm(AS, co), "age_sex_grade": arm(ASG, co), "age_sex_stage": arm(ASS, co),
             "titan": arm(co["T"], co), "omics": arm(co["P"], co)}
        wsi_om = fold_pct(z["titan"] + z["omics"], co["assign"])
        full = (fold_pct(z["titan"] + z["omics"] + z["age_sex_stage"], co["assign"])
                if z["age_sex_stage"] is not None else None)
        full_grade = (fold_pct(z["titan"] + z["omics"] + z["age_sex_grade"], co["assign"])
                      if z["age_sex_grade"] is not None else None)

        def sc(v):
            return None if v is None else round(cidx(v, ii, jj), 4)

        row = {"n_cases": len(co["keep"]), "events": int(co["e"].sum()),
               "cases_with_clinical_row": co["cases_with_clinical_row"],
               "stage_levels": co["stage_levels"], "grade_levels": co["grade_levels"],
               "arms": {k: sc(v) for k, v in z.items()},
               "WSI_plus_OMICS": sc(wsi_om),
               "OURS_wsi_omics_age_sex_stage": sc(full),
               "control_same_arm_with_GRADE_instead": sc(full_grade),
               "published_age_sex_grade": PUBLISHED_AGE_SEX_GRADE.get(c),
               "published_best_any_protocol": PUBLISHED_BEST.get(c),
               "published_best_verified_released_folds": PUBLISHED_RELEASED_FOLD_BEST.get(c)}
        row["our_age_sex_grade_minus_published"] = (
            round(row["arms"]["age_sex_grade"] - PUBLISHED_AGE_SEX_GRADE[c], 4)
            if row["arms"]["age_sex_grade"] is not None and c in PUBLISHED_AGE_SEX_GRADE else None)
        row["stage_minus_grade_same_pipeline"] = (
            round(row["arms"]["age_sex_stage"] - row["arms"]["age_sex_grade"], 4)
            if None not in (row["arms"]["age_sex_stage"], row["arms"]["age_sex_grade"]) else None)
        rep["cohorts"][c] = row
        print("%-9s n=%4d ev=%3d | age+sex %.4f | +grade %.4f (pub %.3f) | +stage %.4f | "
              "WSI+OM %.4f | OURS %.4f | best-any %.3f"
              % (c, row["n_cases"], row["events"], row["arms"]["age_sex"],
                 row["arms"]["age_sex_grade"] or float("nan"),
                 PUBLISHED_AGE_SEX_GRADE.get(c, float("nan")),
                 row["arms"]["age_sex_stage"] or float("nan"), row["WSI_plus_OMICS"],
                 row["OURS_wsi_omics_age_sex_stage"] or float("nan"),
                 PUBLISHED_BEST[c][1]), file=sys.stderr, flush=True)

    ok = [v for v in rep["cohorts"].values() if v.get("OURS_wsi_omics_age_sex_stage")]
    if ok:
        rep["summary"] = {
            "cohorts_scored": len(ok),
            "mean_OURS": round(float(np.mean([v["OURS_wsi_omics_age_sex_stage"] for v in ok])), 4),
            "mean_WSI_plus_OMICS": round(float(np.mean([v["WSI_plus_OMICS"] for v in ok])), 4),
            "cohorts_where_OURS_beats_best_published_any_protocol": [
                k for k, v in rep["cohorts"].items()
                if v.get("OURS_wsi_omics_age_sex_stage")
                and v["OURS_wsi_omics_age_sex_stage"] > v["published_best_any_protocol"][1]],
            "cohorts_where_OURS_beats_best_verified_released_fold_method": [
                k for k, v in rep["cohorts"].items()
                if v.get("OURS_wsi_omics_age_sex_stage")
                and v["OURS_wsi_omics_age_sex_stage"]
                > v["published_best_verified_released_folds"][1]],
            "cohorts_where_stage_beats_grade": [
                k for k, v in rep["cohorts"].items()
                if v.get("stage_minus_grade_same_pipeline") is not None
                and v["stage_minus_grade_same_pipeline"] > 0]}

    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
