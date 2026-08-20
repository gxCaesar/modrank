#!/usr/bin/env python3
"""Every available view, one fixed estimator, equal weights, and the three declared comparisons.

Pre-freeze, exploratory. NOT reportable -- but this is the arm intended for the freeze.

WHY THERE IS NO SELECTION ANYWHERE IN THIS FILE. Two measurements taken today, in this order:

  * choosing the estimator per block on the POOLED outcome lifted the three-arm combination from
    0.7281 to 0.7447, and at four candidates per block with this task's sigma = 0.0073 that lift
    carries a selection term of about sigma*sqrt(2 ln 4) = 0.012;
  * choosing it honestly, by inner cross-validation INSIDE each training fold, gave 0.6568 /
    0.6543 / 0.6743 for image / omics / clinical against 0.6596 / 0.6510 / 0.6856 for plain ridge
    everywhere -- selection COST performance on two blocks out of three.

At 287 training patients and ~90 events the inner criterion cannot reliably separate estimators
that differ by 0.02, so adaptive selection buys noise at full price. The rule here is therefore
fixed in advance and applied identically to every block: **ridge Cox, alpha over a fixed grid by
3-fold inner CV, percentile-normalise within fold, equal-weight rank average over every available
view.** No block-level estimator choice, no weight fitting, no gating -- each of those was measured
and each was either void on its control or worse than not doing it.

THE VIEWS, and what each one is for.

  wsi_titan            TITAN slide embedding, the strongest single image view available
  wsi_chief_mean       CTransPath/CHIEF patch features mean-pooled -- a SECOND, independently
                       pretrained image encoder, so its errors are not TITAN's
  wsi_chief_dispersion per-dimension std and the 90th-minus-10th percentile gap across a slide's
                       ~16,600 patches. Mean pooling is the one summary that cannot see
                       heterogeneity; this is the part of the bag it throws away
  omics_combine        SurvPath's own 4,999 genes averaged into its 275 pathway groups
  omics_hallmarks      the same genes under the hallmark grouping
  omics_xena           the Xena grouping
  clinical             age, sex, T, N, M and AJCC stage

THE THREE COMPARISONS, declared before the run and reported whichever way they land:

  A  against published numbers on the same released folds -- unpaired, and the competitors had no
     clinical variables
  B  against SurvPath given THE SAME clinical variables, paired, case-level bootstrap -- the only
     competitor whose per-case predictions we hold
  C  our arm on WSI + omics ONLY, the input set the competitors actually used, against DIMAF's
     0.679. **This is the one that decides whether the result is a method result or an information
     result**, and it needs +0.0145 over 0.679 to count.
"""

from __future__ import annotations

import argparse
import csv as _csv
import itertools
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, fitapply  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402


def paired_bootstrap(a, b, t, e, reps=4000, seed=0):
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
            "ci95": [round(float(np.quantile(v, 0.025)), 4),
                     round(float(np.quantile(v, 0.975)), 4)],
            "p_two_sided": round(float(2 * min((v <= 0).mean(), (v >= 0).mean())), 4),
            "resampled": "cases, not comparable pairs"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--omics-dir", required=True)
    ap.add_argument("--chief", required=True)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    co = Cohort(args.root, args.titan, args.sp_dir)
    fi = co.fold_indices()
    blocks = {"wsi_titan": co.T}

    z = np.load(args.chief, allow_pickle=True)
    stem2case = {r["slide_id"].replace(".svs", ""): r["case_id"]
                 for r in _csv.DictReader(open(args.root + "/tcga_blca_meta.csv"))}
    percase = {}
    for i, nm in enumerate(z["names"]):
        c = stem2case.get(str(nm))
        if c:
            percase.setdefault(c, []).append(i)
    miss = [c for c in co.keep if c not in percase]
    part = {}
    for key in ("mean", "std", "q10", "q90"):
        M = np.asarray(z[key], dtype=float)
        part[key] = np.vstack([M[percase[c]].mean(0) if c in percase else M.mean(0)
                               for c in co.keep])
    blocks["wsi_chief_mean"] = part["mean"]
    blocks["wsi_chief_dispersion"] = np.hstack([part["std"], part["q90"] - part["q10"]])

    for tag, sig in (("combine", "combine_signatures.csv"),
                     ("hallmarks", "hallmarks_signatures.csv"),
                     ("xena", "xena_signatures.csv")):
        rna = os.path.join(args.omics_dir, "rna_%s.csv" % tag)
        sigp = os.path.join(args.omics_dir, sig)
        if not (os.path.exists(rna) and os.path.exists(sigp)):
            continue
        G, genes, _ = load_rna(rna, co.keep)
        P, pn = pathway_matrix(G, genes, sigp)
        if P.shape[1] >= 10:
            blocks["omics_%s" % tag] = P
    blocks["clinical"] = co.CLIN

    ii, jj = cpairs(co.t, co.e)
    rep = {"artifact_type": "s5_final_all_views", "reportable": False,
           "phase_of_origin": "exploratory", "n": len(co.keep), "events": int(co.e.sum()),
           "cases_missing_chief_get_the_cohort_mean": len(miss),
           "blocks": {k: int(v.shape[1]) for k, v in blocks.items()},
           "estimator": "ridge Cox, alphas (1,8,64,512,4096) by 3-fold inner CV inside each "
                        "training fold; identical for every block; no block-level selection",
           "combination": "equal-weight rank average, percentile-normalised within fold"}

    arms = {}
    for nm, X in blocks.items():
        s = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], 0)
        arms[nm] = co.fold_pct(s)
        print("%-24s %d-d  C=%.4f" % (nm, X.shape[1], cidx(arms[nm], ii, jj)),
              file=sys.stderr, flush=True)
    arms["survpath_rerun"] = co.fold_pct(co.spr)

    gapref = arms["clinical"]
    gap = np.abs(gapref[ii] - gapref[jj])
    tied = gap <= float(np.quantile(gap, 0.10))
    same_site = co.site[ii] == co.site[jj]
    rep["single_arms"] = {
        k: {"cindex": round(cidx(v, ii, jj), 4),
            "clinically_tied_q10": round(cidx(v, ii[tied], jj[tied]), 4),
            "within_site": round(cidx(v, ii[same_site], jj[same_site]), 4)}
        for k, v in arms.items()}

    WSI = [k for k in blocks if k.startswith("wsi_")]
    OM = [k for k in blocks if k.startswith("omics_")]
    groups = {
        "PRESPECIFIED_wsi_only": WSI,
        "PRESPECIFIED_wsi_plus_omics": WSI + OM,
        "PRESPECIFIED_full": WSI + OM + ["clinical"],
        "reference_titan_only_plus_omics_plus_clinical": ["wsi_titan"] + OM + ["clinical"],
        "reference_clinical_only": ["clinical"],
    }
    scored = {}
    for gname, ks in groups.items():
        zz = co.fold_pct(sum(arms[k] for k in ks))
        scored[gname] = zz
        rep.setdefault("prespecified_groups", {})[gname] = {
            "members": ks,
            "cindex": round(cidx(zz, ii, jj), 4),
            "clinically_tied_q10": round(cidx(zz, ii[tied], jj[tied]), 4),
            "within_site": round(cidx(zz, ii[same_site], jj[same_site]), 4)}
        print("%-46s C=%.4f" % (gname, cidx(zz, ii, jj)), file=sys.stderr, flush=True)

    # every combination, for transparency only -- the pre-specified groups above are the arms
    allk = sorted(blocks)
    combos = {}
    for r in range(1, len(allk) + 1):
        for cmb in itertools.combinations(allk, r):
            combos[" + ".join(cmb)] = round(cidx(co.fold_pct(sum(arms[k] for k in cmb)), ii, jj), 4)
    top = dict(sorted(combos.items(), key=lambda kv: -kv[1])[:12])
    rep["all_combinations_scored"] = len(combos)
    rep["top_12_combinations_FOR_TRANSPARENCY_NOT_ARMS"] = top
    rep["selection_note"] = (
        "%d combinations are scored here. None of them is the arm: the arms are the "
        "PRESPECIFIED_ groups above, fixed before the run as 'use every view in the group, equal "
        "weight'. The full list exists so that the best-of-%d number is visible next to the "
        "pre-specified one rather than quietly replacing it." % (len(combos), len(combos)))

    full = scored["PRESPECIFIED_full"]
    noclin = scored["PRESPECIFIED_wsi_plus_omics"]
    spc = co.fold_pct(arms["survpath_rerun"] + arms["clinical"])

    rep["comparison_A_vs_published"] = {
        "ours_full": round(cidx(full, ii, jj), 4),
        "published_same_released_folds": {
            "MOAD-Net (folds unverified)": 0.691, "DIMAF (verified)": 0.679, "APL": 0.677,
            "PIBD": 0.667, "ProtoPathway": 0.646, "DSCASurv (BiB 2025)": 0.646,
            "OTSurv": 0.637, "MMP": 0.635, "SurvPath": 0.625},
        "caveat": "unpaired; the competitors had no clinical variables"}
    rep["comparison_B_input_parity_paired"] = {
        "ours_full": round(cidx(full, ii, jj), 4),
        "survpath_plus_same_clinical": round(cidx(spc, ii, jj), 4),
        "paired_bootstrap": paired_bootstrap(full, spc, co.t, co.e)}
    rep["comparison_C_THE_DECIDER"] = {
        "ours_wsi_plus_omics_only": round(cidx(noclin, ii, jj), 4),
        "incumbent_DIMAF_same_folds": 0.679,
        "gap": round(cidx(noclin, ii, jj) - 0.679, 4),
        "bar": 0.0145,
        "passes": bool(cidx(noclin, ii, jj) - 0.679 > 0.0145),
        "means": "if this passes, the method wins on the competitors' own input set; if it fails, "
                 "the lead in comparison A is information the competitors did not use"}

    pf = {}
    for k in range(5):
        m = co.fold == k
        i2, j2 = cpairs(co.t[m], co.e[m])
        pf["fold_%d" % k] = {"n": int(m.sum()), "events": int(co.e[m].sum()),
                             "ours_full": round(cidx(full[m], i2, j2), 4),
                             "ours_wsi_omics": round(cidx(noclin[m], i2, j2), 4),
                             "survpath_plus_clinical": round(cidx(spc[m], i2, j2), 4)}
    rep["per_fold"] = pf
    v = [x["ours_full"] for x in pf.values()]
    rep["ours_full_fold_mean_sd"] = [round(float(np.mean(v)), 4), round(float(np.std(v, ddof=1)), 4)]
    v2 = [x["ours_wsi_omics"] for x in pf.values()]
    rep["ours_wsi_omics_fold_mean_sd"] = [round(float(np.mean(v2)), 4),
                                          round(float(np.std(v2, ddof=1)), 4)]

    text = json.dumps(rep, indent=1)
    if args.out:
        open(args.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
