#!/usr/bin/env python3
"""SurvPath at five seeds, so the comparison is seed-matched.

WHY THIS WAS OWED. Our primary metric is a mean over five inner-cross-validation seeds. SurvPath's
rerun was a single run at its default seed 1. Comparing a five-seed mean against a one-seed value
is an asymmetry a reviewer can name, and the direction it favours is not knowable in advance -- an
average is less noisy, but the single draw could have landed high or low. Seeds 2 to 5 were
therefore run with SurvPath's configuration otherwise UNCHANGED, on the same released folds.

WHAT IS REPORTED
  * SurvPath alone, per seed and averaged, against its published 0.625
  * SurvPath given the same clinical block, per seed and averaged -- the input-parity comparator
  * the paired case-level bootstrap of our arm against the seed-AVERAGED competitor, which is the
    like-for-like version of the headline comparison
  * and the competitor's own seed spread, which is a quantity this benchmark has never reported
    and which bears directly on whether its published margins are resolvable
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, fitapply, pct  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402
from s6_amend_clinical import dimaf_clinical, boot  # noqa: E402


def load_seed(dirpath, keep, idx):
    """Per-case SurvPath risk and its fold, for one seed, aligned to `keep`."""
    risk = np.full(len(keep), np.nan)
    fold = np.full(len(keep), -1)
    t = np.full(len(keep), np.nan)
    e = np.full(len(keep), np.nan)
    for k, f in enumerate(sorted(glob.glob(os.path.join(dirpath, "split_*_results.pkl")))):
        for c, v in pickle.load(open(f, "rb")).items():
            if c in idx:
                i = idx[c]
                risk[i] = float(v["risk"]); fold[i] = k
                t[i] = float(v["time"]); e[i] = 1.0 if float(v["censorship"]) == 0.0 else 0.0
    return risk, fold, t, e


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--seed-dir", required=True)
    ap.add_argument("--omics-dir", required=True)
    ap.add_argument("--dimaf-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    ii, jj = cpairs(co.t, co.e)
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    CLIN = np.hstack([dc["age"], dc["fem"], dc["stage"]])

    def arm(X, seed=0):
        s = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], seed)
        return co.fold_pct(s)

    zt, zo, zc = arm(co.T), arm(P), arm(CLIN)
    OURS = co.fold_pct(zt + zo + zc)

    rep = {"artifact_type": "s7_survpath_multiseed", "phase_of_origin": "frozen_post_hoc",
           "n": len(co.keep), "events": int(co.e.sum()),
           "what_varied": "SurvPath's --seed only; every other argument identical to the rerun that "
                          "reproduced 0.6147 against a published 0.625"}

    per_seed = {}
    zsp = {}
    # seed 1 is the original rerun already on disk
    r1, f1, t1, e1 = load_seed(a.sp_dir, co.keep, co.idx)
    assert np.allclose(t1, co.t) and np.allclose(e1, co.e), "seed-1 labels differ from the cohort's"
    zsp[1] = co.fold_pct(r1, f1)
    for s in (2, 3, 4, 5):
        d = os.path.join(a.seed_dir, "seed%d" % s)
        if not os.path.isdir(d):
            continue
        r, f, t, e = load_seed(d, co.keep, co.idx)
        bad = int(np.sum(~np.isclose(np.nan_to_num(t), np.nan_to_num(co.t))))
        assert bad == 0, "seed %d labels differ from the cohort's on %d cases" % (s, bad)
        assert (f == co.fold).all(), "seed %d fold assignment differs" % s
        zsp[s] = co.fold_pct(r, f)

    for s, z in sorted(zsp.items()):
        zc_plus = co.fold_pct(z + zc)
        per_seed[s] = {"survpath_alone": round(cidx(z, ii, jj), 4),
                       "survpath_plus_clinical": round(cidx(zc_plus, ii, jj), 4)}
        print("seed %d  survpath=%.4f  +clinical=%.4f"
              % (s, per_seed[s]["survpath_alone"], per_seed[s]["survpath_plus_clinical"]),
              file=sys.stderr, flush=True)
    rep["per_seed"] = per_seed

    alone = [v["survpath_alone"] for v in per_seed.values()]
    plus = [v["survpath_plus_clinical"] for v in per_seed.values()]
    rep["survpath_alone"] = {"mean": round(float(np.mean(alone)), 4),
                             "sd_over_seeds": round(float(np.std(alone, ddof=1)), 4),
                             "min": min(alone), "max": max(alone),
                             "published": 0.625, "published_sd_over_folds": 0.056,
                             "seeds": sorted(zsp)}
    rep["survpath_plus_clinical"] = {"mean": round(float(np.mean(plus)), 4),
                                     "sd_over_seeds": round(float(np.std(plus, ddof=1)), 4),
                                     "min": min(plus), "max": max(plus)}

    # the seed-averaged competitor: average its per-case percentile ranks across seeds, then
    # re-normalise within fold, which is exactly what our own arm does across its seeds
    Z = np.vstack([zsp[s] for s in sorted(zsp)])
    sp_avg = co.fold_pct(Z.mean(0))
    sp_avg_plus = co.fold_pct(sp_avg + zc)
    rep["seed_averaged_competitor"] = {
        "survpath_alone": round(cidx(sp_avg, ii, jj), 4),
        "survpath_plus_clinical": round(cidx(sp_avg_plus, ii, jj), 4)}

    rep["headline_like_for_like"] = {
        "ours_seed0": round(cidx(OURS, ii, jj), 4),
        "survpath_plus_clinical_seed_averaged": round(cidx(sp_avg_plus, ii, jj), 4),
        "paired": boot(OURS, sp_avg_plus, co.t, co.e, 6000),
        "note": "our arm at seed 0 against the competitor averaged over its five seeds -- the "
                "conservative direction, since averaging can only reduce the competitor's noise"}
    rep["vs_single_seed_for_comparison"] = {
        "survpath_plus_clinical_seed1_only": per_seed[1]["survpath_plus_clinical"],
        "paired": boot(OURS, co.fold_pct(zsp[1] + zc), co.t, co.e, 6000)}

    rep["reading"] = (
        "SurvPath's spread across seeds is a number this benchmark has never reported. If it is "
        "comparable to the margins the literature publishes, then those margins are not resolvable "
        "by a single run, and every published value here -- including the one this paper is "
        "measured against -- is a draw from a distribution nobody has characterised.")

    open(a.out, "w").write(json.dumps(rep, indent=1) + "\n")
    print(json.dumps({k: rep[k] for k in ("survpath_alone", "survpath_plus_clinical",
                                          "seed_averaged_competitor", "headline_like_for_like")},
                     indent=1))
    print("wrote", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
