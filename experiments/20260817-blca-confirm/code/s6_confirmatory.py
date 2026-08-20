#!/usr/bin/env python3
"""The confirmatory run, executed exactly once against a frozen protocol.

Reads `development/benchmark-protocol.json` and does what it says. Nothing here is tunable: the
arm, the estimator, the alpha grid, the combination rule, the seeds, the comparator set and the
decision rule were all fixed before this file ran, and `frozen_at` is the boundary G-M4 enforces
by mtime.

WHAT IS AND IS NOT CONFIRMATORY, stated here because it is the honest limit of this run. There is
no unopened held-out split on this benchmark -- it is five-fold cross-validation over 359 patients
and all of them were visible during S5, as they were for every published entrant. So this is not a
held-out test. What the freeze buys is that every choice is now fixed, that **seeds 1-4 have never
been executed**, and that the margin is judged against the selection-inflation term over 269 scored
candidates (0.0244) rather than the reseed bar (0.0145) alone.

If the primary metric misses, that is the result. The protocol's `if_the_frozen_run_misses` says to
return to S5 as a NEW exploratory cycle with a NEW freeze -- not to edit this one.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, fitapply  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402


def boot(x, y, t, e, reps=6000, seed=0):
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
            "reps": int(v.size), "resampled": "cases, not comparable pairs"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--protocol", required=True)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--omics-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    proto = json.load(open(a.protocol))
    assert proto["protocol_status"] == "frozen", "the protocol is not frozen"
    assert not proto["freeze_blockers"], "the protocol carries freeze blockers"
    seeds = proto["method"]["seeds"]["values"]
    alphas = tuple(proto["method"]["estimator"]["alpha_grid"])
    target = 0.679 + 0.0239                      # protocol prespecified_decision_rule
    seed_sd_bar = 0.0145

    co = Cohort(a.root, a.titan, a.sp_dir)
    assert len(co.keep) == proto["task"]["n_cases"], "cohort size differs from the frozen protocol"
    assert int(co.e.sum()) == proto["task"]["n_events"], "event count differs from the protocol"
    fi = co.fold_indices()
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    blocks = {"wsi_titan": co.T, "omics_combine": P, "clinical": co.CLIN}
    ii, jj = cpairs(co.t, co.e)
    assert ii.size == proto["task"]["n_comparable_pairs"], "pair count differs from the protocol"

    per_seed = []
    for sd in seeds:
        arms = {}
        for nm, X in blocks.items():
            s = np.full(len(co.keep), np.nan)
            for tri, vai in fi:
                if len(tri) < 30 or not len(vai):
                    continue
                s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], sd, alphas)
            arms[nm] = co.fold_pct(s)
        ours = co.fold_pct(sum(arms.values()))
        row = {"seed": sd, "OURS": round(cidx(ours, ii, jj), 4),
               **{k: round(cidx(v, ii, jj), 4) for k, v in arms.items()}}
        per_seed.append(row)
        print("seed %d  OURS=%.4f  titan=%.4f omics=%.4f clinical=%.4f"
              % (sd, row["OURS"], row["wsi_titan"], row["omics_combine"], row["clinical"]),
              file=sys.stderr, flush=True)

    vals = [r["OURS"] for r in per_seed]
    primary, primary_sd = float(np.mean(vals)), float(np.std(vals, ddof=1))

    # Comparators are computed at seed 0. The primary metric is seed-averaged; a PAIRED bootstrap
    # needs one pair of per-case vectors, not five, and averaging percentile ranks across seeds
    # before pairing would compare an ensemble against a single fit.
    sd0 = {}
    for nm, X in blocks.items():
        s = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], 0, alphas)
        sd0[nm] = co.fold_pct(s)
    ours0 = co.fold_pct(sum(sd0.values()))
    noclin0 = co.fold_pct(sd0["wsi_titan"] + sd0["omics_combine"])
    spc0 = co.fold_pct(co.fold_pct(co.spr) + sd0["clinical"])

    gap = np.abs(sd0["clinical"][ii] - sd0["clinical"][jj])
    tied = gap <= float(np.quantile(gap, 0.10))
    same_site = co.site[ii] == co.site[jj]

    rep = {
        "artifact_type": "s6_confirmatory_run",
        "protocol": {"path": os.path.relpath(a.protocol), "frozen_at": proto["frozen_at"]},
        "executed_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "phase_of_origin": "frozen",
        "n": len(co.keep), "events": int(co.e.sum()), "comparable_pairs": int(ii.size),
        "per_seed": per_seed,
        "primary": {"metric": "harrell_c_index", "value": round(primary, 4),
                    "sd_over_seeds": round(primary_sd, 4), "seeds": seeds,
                    "aggregation": "mean over seeds of the pooled out-of-fold C-index",
                    "note": "seeds 1-4 had never been executed before this run"},
        "component_arms_seed0": {k: round(cidx(v, ii, jj), 4) for k, v in sd0.items()},
        "secondary": {
            "OURS_clinically_tied_q10": round(cidx(ours0, ii[tied], jj[tied]), 4),
            "survpath_plus_clinical_clinically_tied_q10": round(cidx(spc0, ii[tied], jj[tied]), 4),
            "OURS_within_site": round(cidx(ours0, ii[same_site], jj[same_site]), 4),
            "within_site_pairs": int(same_site.sum()),
            "within_site_caveat": "1,765 of 24,219 pairs share a site; differences of ~0.006 in "
                                  "that column are not readable"},
        "comparisons": {
            "vs_clinical_cheap_baseline": {
                "ours_seed0": round(cidx(ours0, ii, jj), 4),
                "baseline": round(cidx(sd0["clinical"], ii, jj), 4),
                "paired": boot(ours0, sd0["clinical"], co.t, co.e)},
            "vs_titan_alone": {
                "ours_seed0": round(cidx(ours0, ii, jj), 4),
                "baseline": round(cidx(sd0["wsi_titan"], ii, jj), 4),
                "paired": boot(ours0, sd0["wsi_titan"], co.t, co.e)},
            "vs_survpath_plus_same_clinical_INPUT_PARITY": {
                "ours_seed0": round(cidx(ours0, ii, jj), 4),
                "baseline": round(cidx(spc0, ii, jj), 4),
                "paired": boot(ours0, spc0, co.t, co.e)},
            "vs_DIMAF_published": {
                "ours_primary": round(primary, 4), "incumbent": 0.679, "incumbent_sd": 0.043,
                "gap": round(primary - 0.679, 4),
                "paired": "NOT POSSIBLE -- no per-case predictions released",
                "fold_identity": "CONTESTED between two independent literature passes"},
            "ours_without_clinical_vs_DIMAF": {
                "ours_seed0": round(cidx(noclin0, ii, jj), 4), "incumbent": 0.679,
                "gap": round(cidx(noclin0, ii, jj) - 0.679, 4),
                "declared_in_advance_as_not_claimed": True}},
    }

    per_fold = {}
    for k in range(5):
        m = co.fold == k
        i2, j2 = cpairs(co.t[m], co.e[m])
        per_fold["fold_%d" % k] = {
            "n": int(m.sum()), "events": int(co.e[m].sum()),
            "OURS": round(cidx(ours0[m], i2, j2), 4),
            "survpath_plus_clinical": round(cidx(spc0[m], i2, j2), 4)}
    rep["per_fold_seed0"] = per_fold
    rep["folds_won_vs_survpath_plus_clinical"] = sum(
        1 for v in per_fold.values() if v["OURS"] > v["survpath_plus_clinical"])

    rep["decision"] = {
        "rule": "primary > 0.679 + 0.0239 = %.4f AND sd over seeds <= %.4f" % (target, seed_sd_bar),
        "target": round(target, 4),
        "primary": round(primary, 4),
        "primary_clears_target": bool(primary > target),
        "seed_sd": round(primary_sd, 4),
        "seed_sd_within_bar": bool(primary_sd <= seed_sd_bar),
        "VERDICT": "PASS" if (primary > target and primary_sd <= seed_sd_bar) else "MISS",
        "if_miss": "return to S5 Step 1 as a NEW exploratory cycle with a NEW freeze; do not edit "
                   "this protocol",
    }

    os.makedirs(a.out_dir, exist_ok=True)
    out = os.path.join(a.out_dir, "confirmatory.json")
    open(out, "w").write(json.dumps(rep, indent=1) + "\n")
    print(json.dumps(rep["decision"], indent=1))
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
