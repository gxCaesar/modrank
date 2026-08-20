#!/usr/bin/env python3
"""Two measurements the cold panel earned, both of which the paper should have had already.

=========================================================================== 1. IS THE COMBINATION
=========================================================================== TRANSDUCTIVE?

The statistics-leakage reviewer read the phrase "percentile normalisation within fold" and asked
the exact right question:

    "It is not stated whether validation predictions are ranked against only a training-derived
     empirical distribution or against all predictions in the held-out fold. If held-out cases
     define one another's percentiles, the prediction for a patient depends on the other validation
     patients. That is transductive preprocessing... If validation-fold ranks were used, the authors
     should rerun the primary analysis using a transformation fitted exclusively on training
     predictions."

They did. `fold_pct` ranks each validation fold's scores against each other. No outcome is involved,
so this is not label leakage -- and it cannot change any SINGLE arm's concordance, because
concordance is rank-based and a within-fold monotone map leaves an arm's own ranking untouched. But
it can change the COMBINED arm, because the rank average depends on how the three arms are scaled
relative to one another, and that scaling is read off the held-out fold.

So it is testable, and this is the test. INDUCTIVE variant: each validation case's percentile is the
fraction of the TRAINING fold's own predicted scores that fall below it. That is deployable for a
single new patient against a fixed reference distribution, which the transductive version is not.
Everything else is identical.

If the two agree, the objection is answered with a number. If they do not, the paper has to report
the inductive one, because that is the arm a reader could actually deploy.

=========================================================================== 2. THE STAGE-VS-GRADE
=========================================================================== RESULT HAS NO INTERVAL

    "The central stage-versus-grade result has no inferential analysis."

Also correct, and worse than the first: +0.0972 is the paper's headline finding about the field and
it is reported as a bare difference while six secondary comparisons carry intervals. The two arms
are scored on the SAME cases through the SAME pipeline, so the paired case-level bootstrap already
used everywhere else applies directly.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, cox_fit, pct  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402
from s6_amend_clinical import dimaf_clinical, boot  # noqa: E402

ALPHAS = (1.0, 8.0, 64.0, 512.0, 4096.0)


def fit_train_and_val(Xtr, ttr, etr, Xte, seed=0, alphas=ALPHAS):
    """Same estimator as the primary, but returns TRAIN scores as well as validation scores.

    The training scores are the reference distribution the inductive variant needs, and they must
    come from the same fitted model -- refitting or reusing inner-CV predictions would compare two
    different things.
    """
    from blca_common import cindex
    best, ba = -1.0, alphas[0]
    inner = np.array_split(np.random.default_rng(seed).permutation(len(ttr)), 3)
    for al in alphas:
        ip = np.zeros(len(ttr))
        for ite in inner:
            itr = np.setdiff1d(np.arange(len(ttr)), ite)
            mu, sd = Xtr[itr].mean(0), Xtr[itr].std(0) + 1e-9
            ip[ite] = ((Xtr[ite] - mu) / sd) @ cox_fit((Xtr[itr] - mu) / sd, ttr[itr], etr[itr], al)
        c = cindex(ip, ttr, etr)
        if c == c and c > best:
            best, ba = c, al
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
    beta = cox_fit((Xtr - mu) / sd, ttr, etr, ba)
    return ((Xtr - mu) / sd) @ beta, ((Xte - mu) / sd) @ beta


def ecdf_against(reference, values):
    """Fraction of `reference` strictly below each value, with ties at half -- a proper ECDF.

    This is the deployable transform: a single new patient's score is placed against a distribution
    fixed at training time, with no reference to any other patient being scored at the same time.
    """
    r = np.sort(np.asarray(reference, dtype=float))
    lo = np.searchsorted(r, values, side="left")
    hi = np.searchsorted(r, values, side="right")
    return (lo + hi) / (2.0 * max(1, len(r)))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("root", "titan", "sp-dir", "omics-dir", "dimaf-dir", "out"):
        ap.add_argument("--" + f, required=True)
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    ii, jj = cpairs(co.t, co.e)
    n = len(co.keep)
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    CLIN = np.hstack([dc["age"], dc["fem"], dc["stage"]])
    GRADE = np.hstack([dc["age"], dc["fem"], dc["grade"]]) if "grade" in dc else None

    rep = {"artifact_type": "s11_cold_panel_round1_measurements",
           "phase_of_origin": "frozen_post_hoc",
           "why": "two findings from cold panel round 1, answered with measurements rather than "
                  "argument", "n": n, "events": int(co.e.sum())}

    # ================================================================= 1. transductive vs inductive
    SEEDS = (0, 1, 2, 3, 4)
    BLOCKS = {"slide": co.T, "omics": P, "clinical": CLIN}

    def arms_for_seed(seed):
        raw = {k: np.full(n, np.nan) for k in BLOCKS}
        ind = {k: np.full(n, np.nan) for k in BLOCKS}
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            for k, X in BLOCKS.items():
                str_, sva = fit_train_and_val(X[tri], co.t[tri], co.e[tri], X[vai], seed)
                raw[k][vai] = sva
                ind[k][vai] = ecdf_against(str_, sva)      # deployable: training reference only
        # transductive: rank the validation fold against itself, which is what the paper does
        tra = {k: co.fold_pct(raw[k]) for k in BLOCKS}
        return tra, ind

    def combine(d, transductive):
        s = sum(d.values())
        return co.fold_pct(s) if transductive else s

    tr_scores, in_scores, per_seed = [], [], []
    for sd in SEEDS:
        tra, ind = arms_for_seed(sd)
        ct = cidx(combine(tra, True), ii, jj)
        ci = cidx(combine(ind, False), ii, jj)
        per_seed.append({"seed": sd, "transductive": round(ct, 4), "inductive": round(ci, 4),
                         "difference": round(ci - ct, 4)})
        tr_scores.append(co.fold_pct(sum(tra.values())))
        in_scores.append(combine(ind, False))
        print("seed %d  transductive %.4f   inductive %.4f   %+.4f"
              % (sd, ct, ci, ci - ct), file=sys.stderr, flush=True)

    TR = co.fold_pct(np.vstack(tr_scores).mean(0))
    IN = np.vstack(in_scores).mean(0)
    rep["transductive_vs_inductive"] = {
        "question": "does ranking a validation fold against ITSELF, rather than against the "
                    "training fold's own score distribution, change the combined arm?",
        "per_seed": per_seed,
        "primary_transductive": round(float(np.mean([r["transductive"] for r in per_seed])), 4),
        "primary_inductive": round(float(np.mean([r["inductive"] for r in per_seed])), 4),
        "mean_difference": round(float(np.mean([r["difference"] for r in per_seed])), 4),
        "paired_case_bootstrap_inductive_minus_transductive": boot(IN, TR, co.t, co.e, 6000),
        "note": "no outcome enters either transform. A within-fold monotone map cannot change a "
                "SINGLE arm's concordance at all; only the combination can move, because the rank "
                "average depends on the relative scaling of the three arms.",
    }

    # ================================================================= 2. stage vs grade, with an interval
    if GRADE is not None:
        def arm(X, seed=0):
            s = np.full(n, np.nan)
            for tri, vai in fi:
                if len(tri) < 30 or not len(vai):
                    continue
                _, sva = fit_train_and_val(X[tri], co.t[tri], co.e[tri], X[vai], seed)
                s[vai] = sva
            return co.fold_pct(s)
        zs, zg = arm(CLIN), arm(GRADE)
        rep["stage_vs_grade_inference"] = {
            "clinical_with_stage": round(cidx(zs, ii, jj), 4),
            "clinical_with_grade": round(cidx(zg, ii, jj), 4),
            "paired_case_bootstrap": boot(zs, zg, co.t, co.e, 6000),
            "note": "the two arms are scored on the SAME cases through the SAME pipeline and differ "
                    "in one column, so the paired case-level bootstrap used elsewhere applies "
                    "directly. The paper reported this difference without an interval while six "
                    "secondary comparisons carried one.",
        }
        b = rep["stage_vs_grade_inference"]["paired_case_bootstrap"]
        print("stage - grade  %+.4f  CI [%+.4f, %+.4f]  p=%.4g"
              % (b["mean"], b["ci95"][0], b["ci95"][1], b["p_two_sided"]), file=sys.stderr)
    else:
        rep["stage_vs_grade_inference"] = {"status": "grade column absent from the clinical block"}

    open(a.out, "w").write(json.dumps(rep, indent=1) + "\n")
    print(json.dumps({k: v for k, v in rep.items() if k.startswith(("transductive", "stage"))},
                     indent=1))
    print("wrote", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
