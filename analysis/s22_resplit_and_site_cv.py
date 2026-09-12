#!/usr/bin/env python3
"""Does the ordering survive a different partition of the same patients, and a partition by hospital?

POST-FREEZE AND DESCRIPTIVE. The released five folds are one partition of 359 patients, and the
paper shows that fold-to-fold variation (0.103 across folds) is larger than any margin on this
benchmark. Two things a reviewer can ask for and that cost only CPU:

  RESPLITS. The same 24 random five-fold partitions the noise floor was measured on
  (analysis/s5_step0_headroom.py reseeded_folds, rng 1000 + r), estimator seed 0. Every construction
  is fitted on every partition: the three single arms, ModRank, ModRank with the grade arm, the
  concatenated ridge Cox and the stacked Cox, exactly as s19 defines them. Reported: each paired
  difference's distribution over partitions (mean, SD, how often it is positive) and the added-value
  inflation D's.

  SITE-GROUPED CV. TCGA cases carry a tissue source site (the second field of the barcode). Sites
  are assigned whole to five folds, largest first to the fold with fewest cases, so no hospital
  contributes to both the training and the validation side of any fold. This is the nearest thing to
  an external test this cohort allows. It is not one, and it is reported as a sensitivity analysis.

Known answer first: on the RELEASED folds this code must reproduce ModRank at seed 0 (0.7225) and
the concatenated Cox at seed 0 from s19, or it reports nothing.
"""

from __future__ import annotations

import argparse
import collections
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import ALPHAS, Cohort, cidx, cox_fit, cpairs, fitapply, pct  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix                             # noqa: E402
from s5_step0_headroom import reseeded_folds                                  # noqa: E402
from s6_amend_clinical import dimaf_clinical                                  # noqa: E402

EXT_ALPHAS = ALPHAS + (32768.0,)


def site_folds(site, k=5):
    """Whole sites to folds, largest site first into the fold holding fewest cases. Deterministic."""
    cnt = collections.Counter(site)
    order = sorted(cnt, key=lambda s: (-cnt[s], s))
    load, members = [0] * k, [[] for _ in range(k)]
    for s in order:
        j = int(np.argmin(load))
        members[j].append(s)
        load[j] += cnt[s]
    assign = np.array([next(j for j in range(k) if s in members[j]) for s in site])
    return [np.flatnonzero(assign == j) for j in range(k)], members


def run_partition(X, t, e, folds, seed=0):
    """Every construction on one partition. Returns within-fold percentile scores per construction."""
    n = len(t)
    assign = np.full(n, -1)
    for j, va in enumerate(folds):
        assign[va] = j

    def fpct(v):
        out = np.full(n, np.nan)
        for j in range(len(folds)):
            m = assign == j
            out[m] = pct(v[m])
        return out

    raw = {k: np.full(n, np.nan) for k in ("slide", "omics", "clinical", "grade", "concat")}
    stack = np.full(n, np.nan)
    for va in folds:
        tr = np.setdiff1d(np.arange(n), va)
        for k in ("slide", "omics", "clinical", "grade"):
            raw[k][va] = fitapply(X[k][tr], t[tr], e[tr], X[k][va], seed)
        raw["concat"][va] = fitapply(X["concat"][tr], t[tr], e[tr], X["concat"][va], seed,
                                     EXT_ALPHAS)
    z = {k: fpct(v) for k, v in raw.items()}
    for va in folds:                     # stacking, as s19: weights from inner out-of-fold ranks
        tr = np.setdiff1d(np.arange(n), va)
        m = len(tr)
        inner = np.array_split(np.random.default_rng(seed).permutation(m), 3)
        cols = []
        for k in ("slide", "omics", "clinical"):
            io = np.full(m, np.nan)
            for ite in inner:
                itr = np.setdiff1d(np.arange(m), ite)
                io[ite] = pct(fitapply(X[k][tr][itr], t[tr][itr], e[tr][itr], X[k][tr][ite], seed))
            cols.append(io)
        Zi = np.column_stack(cols)
        w = cox_fit(Zi - Zi.mean(0), t[tr], e[tr], 1.0)
        stack[va] = np.column_stack([z[k][va] for k in ("slide", "omics", "clinical")]) @ w
    z["ours"] = fpct(z["slide"] + z["omics"] + z["clinical"])
    z["ours_grade"] = fpct(z["slide"] + z["omics"] + z["grade"])
    z["stacked"] = fpct(stack)
    return z


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("root", "titan", "sp-dir", "omics-dir", "dimaf-dir", "out"):
        ap.add_argument("--" + f, required=True)
    ap.add_argument("--reseeds", type=int, default=24)
    a = ap.parse_args()
    t0 = time.time()
    co = Cohort(a.root, a.titan, a.sp_dir)
    n = len(co.keep)
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    CL = np.hstack([dc["age"], dc["fem"], dc["stage"]])
    X = {"slide": co.T, "omics": P, "clinical": CL,
         "grade": np.hstack([dc["age"], dc["fem"], dc["grade"]]),
         "concat": np.hstack([co.T, P, CL])}
    ii, jj = cpairs(co.t, co.e)

    # ---- known answer on the released folds
    rel = [np.array(sorted(va)) for _, va in co.fold_indices()]
    zr = run_partition(X, co.t, co.e, rel, 0)
    s19 = json.load(open(os.path.join(os.path.dirname(os.path.abspath(a.out)),
                                      "unified-fusion-and-added-value.json")))
    want = {"ours": 0.7225, "concat": s19["fusion_on_identical_inputs"]["concatenated_ridge_cox"]["per_seed"][0],
            "stacked": s19["fusion_on_identical_inputs"]["stacked_learned_weights"]["per_seed"][0]}
    got = {k: round(cidx(zr[k], ii, jj), 4) for k in want}
    print("known answer on the released folds:", got, "want", want, file=sys.stderr)
    if any(abs(got[k] - want[k]) > 5e-5 for k in want):
        print(json.dumps({"status": "error", "error_code": "released_fold_known_answer_failed"}),
              file=sys.stderr)
        return 2

    KEYS = ("slide", "omics", "clinical", "grade", "ours", "ours_grade", "concat", "stacked")
    # ---- resplits
    per = collections.defaultdict(list)
    for r in range(a.reseeds):
        f2 = reseeded_folds(n, np.random.default_rng(1000 + r))
        # reseeded_folds returns (train, validation) pairs; the partition is the validation sets
        z = run_partition(X, co.t, co.e, [np.array(sorted(va)) for _, va in f2], 0)
        for k in KEYS:
            per[k].append(cidx(z[k], ii, jj))
        print("  resplit %2d/%d  ours %.4f  concat %.4f  stacked %.4f  (%.0f s)"
              % (r + 1, a.reseeds, per["ours"][-1], per["concat"][-1], per["stacked"][-1],
                 time.time() - t0), file=sys.stderr, flush=True)
    P_ = {k: np.array(v) for k, v in per.items()}

    def dist(x):
        return {"mean": round(float(x.mean()), 4), "sd": round(float(x.std(ddof=1)), 4),
                "min": round(float(x.min()), 4), "max": round(float(x.max()), 4),
                "fraction_positive": round(float((x > 0).mean()), 4)}
    D = (P_["ours_grade"] - P_["grade"]) - (P_["ours"] - P_["clinical"])
    resplit = {"partitions": a.reseeds, "construction": "analysis/s5_step0_headroom.py "
               "reseeded_folds, rng 1000 + r; estimator seed 0; same estimator and grids as s19",
               "per_construction": {k: dist(P_[k]) for k in KEYS},
               "ours_minus_concat": dist(P_["ours"] - P_["concat"]),
               "ours_minus_stacked": dist(P_["ours"] - P_["stacked"]),
               "ours_minus_clinical": dist(P_["ours"] - P_["clinical"]),
               "clinical_minus_grade": dist(P_["clinical"] - P_["grade"]),
               "D_added_value_inflation": dist(D),
               "per_partition": {k: [round(float(x), 4) for x in P_[k]] for k in KEYS}}

    # ---- site-grouped five-fold CV
    folds_s, members = site_folds(co.site)
    zs = run_partition(X, co.t, co.e, folds_s, 0)
    per_fold = []
    for j, va in enumerate(folds_s):
        i2, j2 = cpairs(co.t[va], co.e[va])
        per_fold.append({"sites": len(members[j]), "cases": int(len(va)), "events": int(co.e[va].sum()),
                         **{k: round(cidx(zs[k][va], i2, j2), 4) for k in KEYS}})
    site = {"sites_total": len(set(co.site)), "folds": per_fold,
            "pooled": {k: round(cidx(zs[k], ii, jj), 4) for k in KEYS},
            "note": "whole tissue source sites per fold; no site on both sides of any fold"}
    rep = {"artifact_type": "s22_resplit_and_site_cv", "phase_of_origin": "post_freeze_2026-09-11",
           "reportable": True, "known_answer": {"got": got, "want": want},
           "resplits": resplit, "site_grouped_cv": site,
           "runtime_seconds": round(time.time() - t0, 1)}
    json.dump(rep, open(a.out, "w"), indent=1)
    print(json.dumps({k: resplit[k] for k in ("ours_minus_concat", "ours_minus_stacked",
                                              "ours_minus_clinical", "D_added_value_inflation")},
                     indent=1), file=sys.stderr)
    print("site-grouped pooled:", site["pooled"], file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
