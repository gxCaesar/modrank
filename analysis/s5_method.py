#!/usr/bin/env python3
"""The method, with every selection moved INSIDE the training fold.

Pre-freeze, exploratory. NOT reportable -- but this is the arm intended for the freeze, so it is
written to be honest first and good second.

WHAT THE METHOD IS, and every line of it is a choice made from a measurement taken today rather
than from a preference:

  1. ONE SURVIVAL HEAD PER MODALITY, never a jointly-trained cross-modal model. On these folds the
     jointly-trained published models run 0.612 (SurvPath, rerun) to 0.679 (DIMAF) while a ridge
     Cox on six routine clinical fields alone reaches 0.6856. At 359 cases and 113 events the joint
     models are estimating far more parameters than the data supports.

  2. THE ESTIMATOR IS CHOSEN PER MODALITY, INSIDE THE TRAINING FOLD, from a small library. This is
     the part that has to be done in-fold and was not, in the sweep that produced it: bagged ridge
     gained the 275-dimensional omics block +0.024 and COST the 768-dimensional image block 0.056,
     so no single estimator is right for every block, and picking per block on the pooled outcome
     would carry a selection term of about sigma*sqrt(2 ln 4) = 0.012 at this task's sigma.

  3. PERCENTILE NORMALISATION WITHIN FOLD before pooling, for every arm equally. Without it a Cox
     log-hazard on one fold's scale is pooled against another's and the comparison measures scale.

  4. EQUAL-WEIGHT RANK AVERAGE, and this is a measured choice, not laziness. A 21-point weight
     sweep on the image/clinical pair peaks at exactly w = 0.50; a two-parameter fusion fitted
     IN-SAMPLE on the held-out fold reaches 0.7268 against the equal weight's 0.7191, so even with
     oracle access a fitted fusion is worth +0.008; and gating the image weight on clinical risk
     tertile scored 0.7493 against a permuted-stratum control at 0.7502 -- the control won, so that
     component is void.

WHAT IS COMPARED, and the parity problem is stated rather than avoided. Our arm uses WSI, bulk
transcriptome and routine clinical staging. The published competitors use WSI and transcriptome and
no clinical variables -- their own clinical baselines use age, sex and GRADE, never stage. So three
comparisons are reported, and the weakest of them is the honest headline:

  A. against published numbers on the same released folds        -- unpaired, no clinical parity
  B. against SurvPath given THE SAME clinical variables, paired  -- the only competitor whose
     per-case predictions we hold, so the only input-parity comparison that can be made at all
  C. our arm WITHOUT clinical variables, against the same published numbers -- what the method is
     worth on WSI+omics alone, which is the input set the competitors actually used

Comparison C is the one that decides whether this is a method result or an information result, and
it is reported whichever way it comes out.
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cindex, cpairs, cox_fit  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402
from s5_c2_estimator_slate import COARSE, DENSE, bagged_arm, ridge_arm  # noqa: E402

LIBRARY = {
    "ridge_coarse": lambda Xtr, t, e, Xte, sd: ridge_arm(Xtr, t, e, Xte, COARSE, 3, sd),
    "ridge_dense": lambda Xtr, t, e, Xte, sd: ridge_arm(Xtr, t, e, Xte, DENSE, 5, sd),
    "bagged_ridge": lambda Xtr, t, e, Xte, sd: bagged_arm(Xtr, t, e, Xte, COARSE, 40, sd),
}


def infold_estimator_arm(X, co, fi, seed=0, record=None):
    """Out-of-fold scores where the ESTIMATOR is selected by inner CV inside each training fold."""
    s = np.full(len(co.keep), np.nan)
    for tri, vai in fi:
        if len(tri) < 30 or not len(vai):
            continue
        inner = list(np.array_split(np.random.default_rng(seed + 5).permutation(len(tri)), 3))
        best, bname = -1.0, next(iter(LIBRARY))
        for nm, fn in LIBRARY.items():
            ip = np.zeros(len(tri))
            for ite in inner:
                itr = tri[np.setdiff1d(np.arange(len(tri)), ite)]
                ip[ite] = fn(X[itr], co.t[itr], co.e[itr], X[tri[ite]], seed)
            c = cindex(ip, co.t[tri], co.e[tri])
            if c == c and c > best:
                best, bname = c, nm
        if record is not None:
            record.append(bname)
        s[vai] = LIBRARY[bname](X[tri], co.t[tri], co.e[tri], X[vai], seed)
    return co.fold_pct(s)


def paired_bootstrap(a, b, t, e, reps=4000, seed=0):
    """Difference in C-index with cases -- not comparable pairs -- resampled."""
    rng = np.random.default_rng(seed)
    n = len(t)
    out = []
    for _ in range(reps):
        bs = rng.choice(n, size=n, replace=True)
        ii, jj = cpairs(t[bs], e[bs])
        if ii.size == 0:
            continue
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
    ap.add_argument("--rna", required=True)
    ap.add_argument("--sig", required=True)
    ap.add_argument("--chief", default="", help="npz of per-slide CHIEF moment features")
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    G, genes, _ = load_rna(a.rna, co.keep)
    P, pnames = pathway_matrix(G, genes, a.sig)

    blocks = {"wsi_titan": co.T, "omics_pathways": P, "clinical": co.CLIN}
    if a.chief and os.path.exists(a.chief):
        z = np.load(a.chief, allow_pickle=True)
        stem2row = {str(n): i for i, n in enumerate(z["names"])}
        import csv as _csv
        stem2case = {r["slide_id"].replace(".svs", ""): r["case_id"]
                     for r in _csv.DictReader(open(a.root + "/tcga_blca_meta.csv"))}
        percase = {}
        for stem, i in stem2row.items():
            c = stem2case.get(stem)
            if c:
                percase.setdefault(c, []).append(i)
        parts = {}
        for key in ("mean", "std", "q10", "q90"):
            M = z[key]
            parts[key] = np.vstack([
                M[percase[c]].mean(0) if c in percase else np.zeros(M.shape[1])
                for c in co.keep])
        blocks["wsi_chief_mean"] = parts["mean"]
        blocks["wsi_chief_dispersion"] = np.hstack([parts["std"], parts["q90"] - parts["q10"]])
        blocks["wsi_chief_full"] = np.hstack([parts[k] for k in ("mean", "std", "q10", "q90")])

    ii, jj = cpairs(co.t, co.e)
    gapref = None
    rep = {"artifact_type": "s5_method_infold", "reportable": False,
           "phase_of_origin": "exploratory", "n": len(co.keep), "events": int(co.e.sum()),
           "blocks": {k: int(v.shape[1]) for k, v in blocks.items()},
           "estimator_library": sorted(LIBRARY)}

    arms, picks = {}, {}
    for nm, X in blocks.items():
        rec = []
        arms[nm] = infold_estimator_arm(X, co, fi, 0, rec)
        picks[nm] = rec
        print("%-24s C=%.4f  estimators picked in-fold: %s"
              % (nm, cidx(arms[nm], ii, jj), rec), file=sys.stderr, flush=True)
    gapref = arms["clinical"]
    gap = np.abs(gapref[ii] - gapref[jj])
    tied = gap <= float(np.quantile(gap, 0.10))
    same_site = co.site[ii] == co.site[jj]

    rep["single_arms"] = {
        k: {"cindex": round(cidx(v, ii, jj), 4),
            "clinically_tied_q10": round(cidx(v, ii[tied], jj[tied]), 4),
            "within_site": round(cidx(v, ii[same_site], jj[same_site]), 4),
            "estimators_picked_in_fold": picks[k]}
        for k, v in arms.items()}
    arms["survpath_rerun"] = co.fold_pct(co.spr)
    rep["single_arms"]["survpath_rerun"] = {
        "cindex": round(cidx(arms["survpath_rerun"], ii, jj), 4),
        "clinically_tied_q10": round(cidx(arms["survpath_rerun"], ii[tied], jj[tied]), 4),
        "within_site": round(cidx(arms["survpath_rerun"], ii[same_site], jj[same_site]), 4),
        "estimators_picked_in_fold": ["official implementation, rerun"]}

    core = [k for k in blocks]
    combos = {}
    for r in range(2, len(core) + 1):
        for cmb in itertools.combinations(sorted(core), r):
            z = co.fold_pct(sum(arms[k] for k in cmb))
            combos[" + ".join(cmb)] = round(cidx(z, ii, jj), 4)
    rep["all_combinations"] = dict(sorted(combos.items(), key=lambda kv: -kv[1]))
    rep["n_combinations"] = len(combos)

    # ---- the three declared comparisons
    ours_full = co.fold_pct(arms["wsi_titan"] + arms["omics_pathways"] + arms["clinical"])
    ours_noclin = co.fold_pct(arms["wsi_titan"] + arms["omics_pathways"])
    sp_plus_clin = co.fold_pct(arms["survpath_rerun"] + arms["clinical"])

    rep["comparison_A_vs_published"] = {
        "ours_wsi_omics_clinical": round(cidx(ours_full, ii, jj), 4),
        "published_on_same_released_folds": {
            "DIMAF (MICCAI 2025, verified folds)": 0.679, "MOAD-Net (arXiv, folds unverified)": 0.691,
            "APL (arXiv/workshop)": 0.677, "PIBD (ICLR 2024)": 0.667,
            "ProtoPathway (arXiv 2026)": 0.646, "DSCASurv (BiB 2025)": 0.646,
            "OTSurv (MICCAI 2025)": 0.637, "MMP (ICML 2024)": 0.635,
            "SurvPath (CVPR 2024)": 0.625},
        "caveat": "unpaired, and the competitors had no clinical variables"}

    rep["comparison_B_input_parity_paired"] = {
        "ours": round(cidx(ours_full, ii, jj), 4),
        "survpath_plus_same_clinical": round(cidx(sp_plus_clin, ii, jj), 4),
        "paired_bootstrap": paired_bootstrap(ours_full, sp_plus_clin, co.t, co.e),
        "note": "the only competitor whose per-case predictions exist here, hence the only "
                "input-parity comparison that can be made"}

    rep["comparison_C_no_clinical"] = {
        "ours_wsi_plus_omics_only": round(cidx(ours_noclin, ii, jj), 4),
        "incumbent_DIMAF": 0.679,
        "gap": round(cidx(ours_noclin, ii, jj) - 0.679, 4),
        "bar": 0.0145,
        "decides": "whether this is a method result or an information result"}

    per_fold = {}
    for k in range(5):
        m = co.fold == k
        i2, j2 = cpairs(co.t[m], co.e[m])
        per_fold["fold_%d" % k] = {
            "n": int(m.sum()), "events": int(co.e[m].sum()),
            "ours": round(cidx(ours_full[m], i2, j2), 4),
            "survpath_plus_clinical": round(cidx(sp_plus_clin[m], i2, j2), 4)}
    rep["per_fold"] = per_fold
    rep["mean_over_folds"] = round(float(np.mean([v["ours"] for v in per_fold.values()])), 4)
    rep["sd_over_folds"] = round(float(np.std([v["ours"] for v in per_fold.values()], ddof=1)), 4)

    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
