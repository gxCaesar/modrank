#!/usr/bin/env python3
"""Leakage falsifiers, run BEFORE the protocol freeze.

Each one must be able to FAIL, and each is checked on the real pipeline rather than on a toy.
A falsifier that cannot reject is a comment.

  L1  fold integrity          every case in exactly one validation fold; train and val disjoint in
                              every fold; the fold assignment identical to the one SurvPath's own
                              rerun used, so the head-to-head is paired
  L2  no fitting across the   the standardisation moments and the ridge alpha must be computed from
      split boundary          training rows only. Tested by MUTATION: refit the same arm with the
                              moments taken from the WHOLE cohort and confirm the score MOVES. If a
                              leak makes no difference the test is blind and cannot guard anything.
  L3  planted label           a synthetic feature that is pure noise in training and a perfect copy
                              of the outcome in the validation fold. An honest fold-respecting
                              pipeline cannot use it; if the C-index rises, information is crossing
                              the boundary.
  L4  case alignment          permute which feature row belongs to which case, keeping the outcome
                              fixed. The arm must collapse to chance. If it does not, the score is
                              not coming from the features.
  L5  outcome provenance      the labels used by our arm and by SurvPath's rerun must be identical
                              case by case -- otherwise the two arms are scored on different
                              endpoints and the paired comparison is meaningless.

Exit 0 only if every falsifier reaches its expected verdict. Anything else blocks the freeze.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cindex, cpairs, cox_fit, fitapply  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402

ALPHAS = (1.0, 8.0, 64.0, 512.0, 4096.0)


def leaky_fitapply(Xtr, ttr, etr, Xte, Xall, seed=0):
    """The SAME estimator with the standardisation moments taken from the whole cohort.

    This is the mutation L2 tests against. It is deliberately wrong.
    """
    mu, sd = Xall.mean(0), Xall.std(0) + 1e-9
    best, ba = -1.0, ALPHAS[0]
    inner = np.array_split(np.random.default_rng(seed).permutation(len(ttr)), 3)
    for al in ALPHAS:
        ip = np.zeros(len(ttr))
        for ite in inner:
            itr = np.setdiff1d(np.arange(len(ttr)), ite)
            ip[ite] = ((Xtr[ite] - mu) / sd) @ cox_fit((Xtr[itr] - mu) / sd,
                                                       ttr[itr], etr[itr], al)
        c = cindex(ip, ttr, etr)
        if c == c and c > best:
            best, ba = c, al
    return ((Xte - mu) / sd) @ cox_fit((Xtr - mu) / sd, ttr, etr, ba)


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
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    ii, jj = cpairs(co.t, co.e)
    rep = {"artifact_type": "s6_leakage_falsifiers", "n": len(co.keep),
           "events": int(co.e.sum()), "falsifiers": {}}
    ok = True

    # ---------------------------------------------------------------- L1 fold integrity
    seen, dis = {}, True
    for k, (tri, vai) in enumerate(fi):
        if set(tri.tolist()) & set(vai.tolist()):
            dis = False
        for i in vai.tolist():
            seen[i] = seen.get(i, 0) + 1
    every_once = all(v == 1 for v in seen.values()) and len(seen) == len(co.keep)
    same_as_survpath = True   # Cohort.__init__ asserts this and would have raised
    v = dis and every_once and same_as_survpath
    rep["falsifiers"]["L1_fold_integrity"] = {
        "status": "ran_and_passed" if v else "ran_and_FAILED",
        "train_val_disjoint_every_fold": bool(dis),
        "each_case_in_exactly_one_val_fold": bool(every_once),
        "cases_covered": len(seen),
        "fold_assignment_identical_to_survpath_rerun": same_as_survpath,
        "command": "s6_leakage_falsifiers.py (L1)"}
    ok &= v

    # ---------------------------------------------------------------- L2 mutation test
    honest = np.full(len(co.keep), np.nan)
    leaky = np.full(len(co.keep), np.nan)
    for tri, vai in fi:
        honest[vai] = fitapply(co.T[tri], co.t[tri], co.e[tri], co.T[vai], 0)
        leaky[vai] = leaky_fitapply(co.T[tri], co.t[tri], co.e[tri], co.T[vai], co.T, 0)
    ch, cl = cidx(co.fold_pct(honest), ii, jj), cidx(co.fold_pct(leaky), ii, jj)
    moved = abs(ch - cl) > 1e-6
    rep["falsifiers"]["L2_no_fitting_across_boundary"] = {
        "status": "ran_and_passed" if moved else "ran_and_FAILED_BLIND",
        "honest_cindex": round(ch, 4), "whole_cohort_moments_cindex": round(cl, 4),
        "difference": round(cl - ch, 4),
        "reading": "the mutation must MOVE the score. If it does not, the test cannot detect a "
                   "leak of this kind and passing it would mean nothing.",
        "command": "s6_leakage_falsifiers.py (L2)"}
    ok &= moved

    # ---------------------------------------------------------------- L3 planted label
    # A FIRST VERSION OF THIS TEST WAS WRONG AND SAID SO LOUDLY: it built ONE planted column for
    # the whole cohort, writing the outcome rank into each case's row while that case sat in its
    # validation fold. But L1 has just established that every case is in exactly one validation
    # fold, so the loop overwrote all 359 rows and the column encoded the outcome EVERYWHERE --
    # including in the four folds where each case is a training row. The pipeline then learned a
    # real coefficient from real information and scored 0.8715. That is the estimator working, not
    # a leak, and reporting it as a leak would have been the day's most expensive mistake.
    #
    # A planted column that is noise in training and outcome in validation cannot be one column;
    # it has to be built INSIDE each fold, which is what this does.
    sp = np.full(len(co.keep), np.nan)
    for k, (tri, vai) in enumerate(fi):
        g = np.random.default_rng(1100 + k)
        col = g.standard_normal(len(co.keep))
        r = np.argsort(np.argsort(-co.t[vai])).astype(float)
        col[vai] = (r - r.mean()) / (r.std() + 1e-9)          # outcome ONLY on validation rows
        Xk = np.hstack([co.T, col[:, None]])
        sp[vai] = fitapply(Xk[tri], co.t[tri], co.e[tri], Xk[vai], 0)
    cp = cidx(co.fold_pct(sp), ii, jj)
    v = abs(cp - ch) < 0.05
    rep["falsifiers"]["L3_planted_label"] = {
        "status": "ran_and_rejected" if v else "ran_and_FAILED",
        "baseline_cindex": round(ch, 4), "with_planted_column": round(cp, 4),
        "difference": round(cp - ch, 4),
        "construction": "the planted column is rebuilt per fold: standard normal on every training "
                        "row, the outcome ordering on every validation row. Fitting sees only "
                        "noise in that column.",
        "reading": "a fold-respecting pipeline fits the planted coefficient on noise and cannot "
                   "exploit the validation values; a leaking one would jump.",
        "tolerance": 0.05, "command": "s6_leakage_falsifiers.py (L3)"}
    ok &= v

    # ---------------------------------------------------------------- L4 case alignment
    # Many draws, not one. A single permutation of 359 cases has a C-index standard error of
    # roughly 0.03, so one draw landing at 0.55 says nothing about whether the arm collapses.
    cms = []
    for d_ in range(30):
        perm = np.random.default_rng(23 + d_).permutation(len(co.keep))
        Tp = co.T[perm]
        sm = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            sm[vai] = fitapply(Tp[tri], co.t[tri], co.e[tri], Tp[vai], 0)
        cms.append(cidx(co.fold_pct(sm), ii, jj))
    cm, cs = float(np.mean(cms)), float(np.std(cms, ddof=1))
    v = abs(cm - 0.5) < 3 * cs / np.sqrt(len(cms))
    rep["falsifiers"]["L4_case_alignment"] = {
        "status": "ran_and_rejected" if v else "ran_and_FAILED",
        "shuffled_alignment_cindex_mean": round(cm, 4),
        "sd_over_draws": round(cs, 4), "draws": len(cms),
        "expected": 0.5, "test": "|mean - 0.5| < 3 * SE of the mean",
        "reading": "permuting which feature row belongs to which case, outcome held fixed, must "
                   "collapse the arm to chance.",
        "command": "s6_leakage_falsifiers.py (L4)"}
    ok &= v

    # ---------------------------------------------------------------- L5 outcome provenance
    import glob
    import pickle
    mism = 0
    for f in sorted(glob.glob(a.sp_dir + "/split_*_results.pkl")):
        for c, d in pickle.load(open(f, "rb")).items():
            if c in co.idx:
                i = co.idx[c]
                if abs(float(d["time"]) - co.t[i]) > 1e-9:
                    mism += 1
                if (1.0 if float(d["censorship"]) == 0.0 else 0.0) != co.e[i]:
                    mism += 1
    v = mism == 0
    rep["falsifiers"]["L5_outcome_provenance"] = {
        "status": "ran_and_passed" if v else "ran_and_FAILED",
        "label_mismatches_against_survpath_rerun": mism,
        "reading": "our arm and SurvPath's rerun must be scored on identical (time, event) per "
                   "case, or the paired comparison compares two endpoints.",
        "command": "s6_leakage_falsifiers.py (L5)"}
    ok &= v

    rep["all_passed"] = bool(ok)
    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
