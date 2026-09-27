#!/usr/bin/env python3
"""Does fitting the fusion help on the benchmark's other four studies, or only fail to help on one?

WHY THIS EXISTS, and it is an obligation the title created rather than a question left over. The
paper now claims that fitting the fusion adds nothing. That is shown on the development cohort,
where a concatenated Cox and a stacked Cox are both beaten in 24 of 24 re-partitions, and it is
shown nowhere else: the five-study runs scored the slide, transcriptome and clinical arms and their
equal-weight average, and never built either fusion alternative. A claim in a title that holds on
one of five studies is a claim a reviewer will ask about, correctly.

WHAT IT BUILDS, per study, on the released folds with everything else held to the five-study
protocol:

  ModRank        the equal-weight rank average of the three arms, as already reported
  concatenated   one ridge Cox on all three blocks side by side, with the penalty grid extended by
                 32,768 in its favour, because a wide model deserves the chance to regularise
  stacked        three weights learned per fold from the arms' INNER out-of-fold percentiles, ridge
                 1.0, applied to the validation fold. The inner split is what keeps the weights from
                 seeing the scores they are weighting

Both differences are the paired case-level bootstrap used everywhere else, 6,000 resamples of
patients with the pairs rebuilt inside each replicate.

WHAT WOULD FALSIFY THE TITLE. If either alternative beats the rank average by an interval excluding
zero in any study, "fitting the fusion adds nothing" is false as a general statement and has to
become a statement about the studies where it holds. This script is written to find that out rather
than to confirm the sentence, and it reports every study either way.

KNOWN ANSWERS. The rank average must reproduce the committed value for all five studies before any
alternative is scored.

Post-freeze, exploratory.
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
import blca_common                                                             # noqa: E402
from blca_common import ALPHAS, cidx, cox_fit, cpairs, fitapply, pct          # noqa: E402
from s5_stage_five_cohorts import COHORTS, fold_pct, load_cohort              # noqa: E402

REPS = 6000
EXT_ALPHAS = ALPHAS + (32768.0,)     # the concatenation's grid, extended in its favour
SEED = 0


def boot(x, y, t, e, reps=REPS, seed=0):
    """Copied from s6_confirmatory.boot, as every other interval in this work is."""
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
            "reps": int(v.size)}


def arm(X, co, alphas=ALPHAS, seed=SEED):
    """Identical to s5_stage_five_cohorts.arm, with the penalty grid exposed."""
    if X.shape[1] == 0:
        return None
    s = np.full(len(co["keep"]), np.nan)
    for tri, vai in co["fi"]:
        if len(tri) < 30 or not len(vai):
            continue
        s[vai] = fitapply(X[tri], co["t"][tri], co["e"][tri], X[vai], seed, alphas)
    return fold_pct(s, co["assign"])


def stacked(blocks, per_arm, co, seed=SEED, tuned=None):
    """Three weights per fold from the arms' inner out-of-fold percentiles, ridge 1.0.

    Copied from s19 rather than re-derived. The inner three-fold split inside each training fold is
    the part that matters: weights learned on the same scores they weight would be fitted on the
    validation fold through the back door, and the result would look like learned fusion winning.
    """
    n = len(co["keep"])
    st = np.full(n, np.nan)
    if tuned is not None:           # amendment A2's added comparator, the ridge chosen per fold
        tuned.update(st=np.full(n, np.nan), ridges=[])
    for tri, vai in co["fi"]:
        if len(tri) < 30 or not len(vai):
            continue
        m = len(tri)
        inner = np.array_split(np.random.default_rng(seed).permutation(m), 3)
        ioof = {b: np.full(m, np.nan) for b in blocks}
        for b, X in blocks.items():
            Xt = X[tri]
            for ite in inner:
                itr = np.setdiff1d(np.arange(m), ite)
                ioof[b][ite] = pct(fitapply(Xt[itr], co["t"][tri][itr], co["e"][tri][itr],
                                            Xt[ite], seed))
        Zi = np.column_stack([ioof[b] for b in blocks])
        w = cox_fit(Zi - Zi.mean(0), co["t"][tri], co["e"][tri], 1.0)
        Zv = np.column_stack([per_arm[b][vai] for b in blocks])
        st[vai] = Zv @ w
        if tuned is not None:
            al, wt = blca_common.stack_weights_tuned(Zi, co["t"][tri], co["e"][tri], inner)
            tuned["st"][vai] = Zv @ wt
            tuned["ridges"].append(al)
    if tuned is not None:
        tuned["score"] = fold_pct(tuned["st"], co["assign"])
    return fold_pct(st, co["assign"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("titan", "meta-dir", "clin-dir", "pan-dir", "sig", "frozen", "out"):
        ap.add_argument("--" + f, required=True)
    a = ap.parse_args()
    t0 = time.time()

    frozen = json.load(open(a.frozen))["cohorts"]
    d0 = pickle.load(open(a.titan, "rb"))
    E = np.asarray(d0["embeddings"])
    stem2emb = {}
    for i, f in enumerate(d0["filenames"]):
        s = str(f)
        for suf in (".svs", ".h5", ".pt"):
            if s.endswith(suf):
                s = s[:-len(suf)]
        stem2emb[s] = E[i]

    out = {"artifact_type": "s32_five_cohort_fusion_alternatives",
           "phase_of_origin": "post_freeze_2026-09-13",
           "answers": "whether 'fitting the fusion adds nothing' holds beyond the development cohort",
           "protocol": "the five-study protocol: released folds, benchmark clinical file, seed 0, "
                       "concatenation on the extended grid, stacking from inner out-of-fold "
                       "percentiles",
           "bootstrap": {"replicates": REPS, "unit": "cases", "seed": 0},
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
        CLIN = np.hstack([co["age"], co["fem"], co["stage"]])
        blocks = {"slide": co["T"], "omics": co["P"], "clinical": CLIN}
        per = {b: arm(X, co) for b, X in blocks.items()}
        if any(v is None for v in per.values()):
            out["cohorts"][c] = {"status": "an arm could not be built"}
            continue
        ours = fold_pct(per["slide"] + per["omics"] + per["clinical"], co["assign"])
        cat = arm(np.hstack([co["T"], co["P"], CLIN]), co, EXT_ALPHAS)
        tun = {} if blca_common.A2 else None
        stk = stacked(blocks, per, co, tuned=tun)

        def sc(v):
            return round(float(cidx(v, ii, jj)), 4)

        want = frozen.get(c, {}).get("OURS_wsi_omics_age_sex_stage")
        if want is not None and abs(sc(ours) - want) > 5e-5:
            drift[c] = {"got": sc(ours), "want": want}

        row = {"n": len(co["keep"]), "events": int(co["e"].sum()),
               "points": {"ModRank": sc(ours), "concatenated": sc(cat), "stacked": sc(stk)},
               "ModRank_minus_concatenated": boot(ours, cat, co["t"], co["e"]),
               "ModRank_minus_stacked": boot(ours, stk, co["t"], co["e"])}
        row["fitting_the_fusion_beats_the_rank_average"] = bool(
            row["ModRank_minus_concatenated"]["ci95"][1] < 0
            or row["ModRank_minus_stacked"]["ci95"][1] < 0)
        if tun is not None:
            row["stacked_tuned_ridge"] = {
                "point": sc(tun["score"]), "ridge_grid": list(blca_common.STACK_RIDGES),
                "ridge_per_fold": tun["ridges"],
                "ModRank_minus_stacked_tuned_ridge": boot(ours, tun["score"], co["t"], co["e"])}
        out["cohorts"][c] = row
        print("%-9s ModRank %.4f | concat %.4f d=%+.4f %s | stacked %.4f d=%+.4f %s  [%.0f s]"
              % (c, row["points"]["ModRank"], row["points"]["concatenated"],
                 row["ModRank_minus_concatenated"]["mean"],
                 row["ModRank_minus_concatenated"]["ci95"], row["points"]["stacked"],
                 row["ModRank_minus_stacked"]["mean"], row["ModRank_minus_stacked"]["ci95"],
                 time.time() - t0), file=sys.stderr, flush=True)

    scored = [v for v in out["cohorts"].values() if v.get("points")]
    out["summary"] = {
        "cohorts_scored": len(scored),
        "studies_where_a_fitted_fusion_beats_the_rank_average": [
            k for k, v in out["cohorts"].items()
            if v.get("fitting_the_fusion_beats_the_rank_average")],
        "studies_where_the_rank_average_beats_concatenation_with_an_interval_excluding_zero": [
            k for k, v in out["cohorts"].items()
            if v.get("ModRank_minus_concatenated", {}).get("ci95", [0, 0])[0] > 0],
        "studies_where_the_rank_average_beats_stacking_with_an_interval_excluding_zero": [
            k for k, v in out["cohorts"].items()
            if v.get("ModRank_minus_stacked", {}).get("ci95", [0, 0])[0] > 0],
        **({"studies_where_tuned_ridge_stacking_and_the_rank_average_separate": [
            k for k, v in out["cohorts"].items() if "stacked_tuned_ridge" in v and (
                v["stacked_tuned_ridge"]["ModRank_minus_stacked_tuned_ridge"]["ci95"][0] > 0
                or v["stacked_tuned_ridge"]["ModRank_minus_stacked_tuned_ridge"]["ci95"][1] < 0)]}
           if blca_common.A2 else {}),
        "reading": "the title's second clause needs the first list to be empty. The other two say "
                   "how much stronger than that the evidence is, and are reported whatever they are",
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
