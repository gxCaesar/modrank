#!/usr/bin/env python3
"""Does the method win on the OTHER FOUR cohorts of this benchmark, or is bladder a special case?

Pre-freeze, exploratory. NOT reportable.

WHY THIS RUN EXISTS. Seven candidate families have now been killed on their own pre-registered
falsifiers, and the pattern in how they died is the finding: at 359 bladder cases and 113 events the
case-level bootstrap SE of a C-index difference is about 0.020-0.023, so margins of 0.03-0.05 -- the
size this whole literature publishes -- land at p = 0.06 to 0.10 whichever method produces them.
The cohort cannot resolve what the field claims on it. That is not a reason to stop; it is a reason
to stop asking one cohort to carry the whole claim.

Every method on this benchmark reports all five TCGA studies SurvPath released -- BLCA, BRCA,
COADREAD, HNSC, STAD -- so evaluating on five is the benchmark's own convention, not a widening of
scope. Bladder stays the lead cohort. And the data is already on disk: the TITAN feature file
covers 11,658 TCGA slides rather than the 457 bladder ones, SurvPath ships each cohort's genes and
labels, and PIBD ships each cohort's released fold files.

THE QUESTION IT SETTLES. On bladder, our WSI+omics arm reaches 0.6841 against DIMAF's 0.679 -- a
+0.005 gap that fails the bar under three different combination rules. If the same arm clears the
published per-cohort numbers on the other four, the bladder result is a cohort quirk and the method
is real. If it fails there too, the honest reading is that the method's advantage needs the
clinical variables, and the paper has to say so.

WHAT IS AND IS NOT MATCHED HERE. The arm is identical to the bladder one: ridge Cox per view,
alphas by 3-fold inner CV inside the training fold, percentile-normalise within fold, equal weight.
The released case-ID folds are used for every cohort. The clinical block is NOT matched -- T, N, M
and AJCC stage were fetched from GDC for bladder only, so the cross-cohort arms carry age and sex
alone and are reported as such rather than quietly labelled `clinical`.
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


def load_cohort(meta_csv, splits_dir, rna_csv, sig_csv, stem2emb):
    rows = list(csv.DictReader(open(meta_csv)))
    per, lab = collections.defaultdict(list), {}
    for r in rows:
        stem = r["slide_id"].replace(".svs", "")
        if stem in stem2emb:
            per[r["case_id"]].append(stem2emb[stem])
        if r.get("survival_months_dss", "").strip() and r.get("censorship_dss", "").strip():
            lab[r["case_id"]] = (float(r["survival_months_dss"]),
                                 1.0 - float(r["censorship_dss"]),
                                 float(r["age"]) if r.get("age", "").strip() else 65.0,
                                 float(r["is_female"]) if r.get("is_female", "").strip() else 0.0)
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
    AS = np.column_stack([[lab[c][2] for c in keep], [lab[c][3] for c in keep]]).astype(float)

    G, genes, _ = load_rna(rna_csv, keep)
    P, _pn = pathway_matrix(G, genes, sig_csv)

    fi, assign = [], np.full(len(keep), -1)
    for k, (tr, va) in enumerate(folds):
        tri = np.array(sorted(idx[c] for c in tr if c in idx))
        vai = np.array(sorted(idx[c] for c in va if c in idx))
        fi.append((tri, vai))
        for i in vai:
            assign[i] = k
    return {"keep": keep, "T": T, "P": P, "AS": AS, "t": t, "e": e, "fi": fi, "assign": assign}


def fold_pct(v, assign):
    out = np.full(len(v), np.nan)
    for k in np.unique(assign):
        if k < 0:
            continue
        m = assign == k
        out[m] = pct(v[m])
    return out


def arm(X, d):
    s = np.full(len(d["keep"]), np.nan)
    for tri, vai in d["fi"]:
        if len(tri) < 30 or not len(vai):
            continue
        s[vai] = fitapply(X[tri], d["t"][tri], d["e"][tri], X[vai], 0)
    return fold_pct(s, d["assign"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--meta-dir", required=True, help="dir holding tcga_<cohort>.csv")
    ap.add_argument("--pan-dir", required=True, help="dir holding rna_<cohort>.csv and splits/")
    ap.add_argument("--sig", required=True)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    d = pickle.load(open(a.titan, "rb"))
    E = np.asarray(d["embeddings"])
    stem2emb = {}
    for i, f in enumerate(d["filenames"]):
        s = str(f)
        for suf in (".svs", ".h5", ".pt"):
            if s.endswith(suf):
                s = s[:-len(suf)]
        stem2emb[s] = E[i]

    rep = {"artifact_type": "s5_pancohort", "reportable": False,
           "phase_of_origin": "exploratory",
           "arm": "identical to the bladder arm: ridge Cox per view, alphas by 3-fold inner CV "
                  "in-fold, percentile within fold, equal weight",
           "clinical_block": "age + sex ONLY -- T/N/M/stage were fetched from GDC for bladder only, "
                             "so the cross-cohort clinical arm is NOT the bladder one",
           "cohorts": {}}

    for c in COHORTS:
        meta = os.path.join(a.meta_dir, "tcga_%s.csv" % c)
        rna = os.path.join(a.pan_dir, "rna_%s.csv" % c)
        spl = os.path.join(a.pan_dir, "splits", c)
        if not (os.path.exists(meta) and os.path.exists(rna) and os.path.isdir(spl)):
            rep["cohorts"][c] = {"status": "inputs missing"}
            continue
        co = load_cohort(meta, spl, rna, a.sig, stem2emb)
        n, ev = len(co["keep"]), int(co["e"].sum())
        ii, jj = cpairs(co["t"], co["e"])
        z_ti, z_om, z_as = arm(co["T"], co), arm(co["P"], co), arm(co["AS"], co)
        wsi_om = fold_pct(z_ti + z_om, co["assign"])
        allthree = fold_pct(z_ti + z_om + z_as, co["assign"])
        pf = []
        for k in range(5):
            m = co["assign"] == k
            i2, j2 = cpairs(co["t"][m], co["e"][m])
            if i2.size:
                pf.append(cidx(wsi_om[m], i2, j2))
        rep["cohorts"][c] = {
            "n_cases": n, "events": ev, "comparable_pairs": int(ii.size),
            "omics_pathways_dim": int(co["P"].shape[1]),
            "titan_only": round(cidx(z_ti, ii, jj), 4),
            "omics_only": round(cidx(z_om, ii, jj), 4),
            "age_sex_only": round(cidx(z_as, ii, jj), 4),
            "WSI_plus_OMICS": round(cidx(wsi_om, ii, jj), 4),
            "WSI_plus_OMICS_plus_age_sex": round(cidx(allthree, ii, jj), 4),
            "WSI_plus_OMICS_fold_mean": round(float(np.mean(pf)), 4),
            "WSI_plus_OMICS_fold_sd": round(float(np.std(pf, ddof=1)), 4)}
        print("%-10s n=%4d ev=%3d  titan=%.4f omics=%.4f age+sex=%.4f  WSI+OM=%.4f (folds %.4f+/-%.4f)"
              % (c, n, ev, rep["cohorts"][c]["titan_only"], rep["cohorts"][c]["omics_only"],
                 rep["cohorts"][c]["age_sex_only"], rep["cohorts"][c]["WSI_plus_OMICS"],
                 rep["cohorts"][c]["WSI_plus_OMICS_fold_mean"],
                 rep["cohorts"][c]["WSI_plus_OMICS_fold_sd"]), file=sys.stderr, flush=True)

    ok = [v for v in rep["cohorts"].values() if "WSI_plus_OMICS" in v]
    if ok:
        rep["mean_across_cohorts_WSI_plus_OMICS"] = round(
            float(np.mean([v["WSI_plus_OMICS"] for v in ok])), 4)
    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
