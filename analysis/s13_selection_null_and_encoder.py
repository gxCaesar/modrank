#!/usr/bin/env python3
"""Two experiments the round-2 panel named, both against the frozen protocol's own machinery.

=========================================================================== A. THE SELECTION NULL

Round 2, reviewer 1: the pre-registered bar was 0.679 + sigma*sqrt(2 ln N) with sigma = 0.0073.
The FROZEN protocol evaluates that term at N = 211, the ledger's size when it was written, giving
0.0239 and a target of 0.7029; the campaign finished at 269 and the same term over that gives
0.0244 and 0.7034. (This file's "as_preregistered" block below carries the SECOND pair, which is a
mislabel found by the 2026-08-20 carpet sweep after the manuscript had inherited it. The key is
left alone because the artifact it names cannot be regenerated without the cohort, and a source
that silently disagreed with its own shipped output would be worse than a name that is wrong in a
documented way. "as_frozen_in_the_protocol" beside it is read from the protocol and is the one to
use.) Neither N is the defect. 0.0073 is the SD of a PAIRED BETWEEN-ARM DELTA under re-splitting, while
E[max of N] - mu ~ sigma*sqrt(2 ln N) needs the SD of a candidate's OWN ABSOLUTE score. A paired
delta between correlated arms is smaller, because common fold noise cancels. The substitution is a
category error and it errs anti-conservatively. The reviewer named the correct experiment:

    "Report the empirical distribution of the candidate scores and the empirical max under a
     permuted-outcome null over the same configurations. That replaces a Gaussian max heuristic
     (which additionally assumes 269 INDEPENDENT looks -- these candidates share folds, data and
     components) with the actual selection distribution."

That is what this does. Per permutation the (time, event) PAIR is permuted across cases, which
destroys the association while leaving the censoring pattern and the folds exactly as they were.
Every candidate in the family is then scored out of fold and the maximum is recorded. The output is
three things the paper did not have: the null SD of a single candidate's absolute score, which is
the sigma the formula actually wanted; the empirical max distribution, which needs no independence
assumption; and the bar each implies.

=========================================================================== B. THE ENCODER CONFOUND

Round 2, reviewer 1 again: ModRank reads TITAN while the comparators read their own tile features,
and the seven-encoder sweep spans 0.137 in the slide arm alone against a claimed +0.0422 lead, so
"method" and "2025 encoder" are not separated. The decisive experiment is to hold the method fixed
and swap the encoder. If the full method's lead survives every encoder, the lead is the method. If
it collapses on six of seven, the lead is TITAN and the paper must say so.

Neither experiment changes a reported number. Both are post-freeze and descriptive, and both can
only make the paper's claims weaker or better supported.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, fitapply                    # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix                         # noqa: E402
from s6_amend_clinical import dimaf_clinical                              # noqa: E402

ALPHAS = (1.0, 8.0, 64.0, 512.0, 4096.0)


def _frozen_bar():
    """The pre-registered selection term, read from the frozen protocol and never hardcoded.

    The protocol states sigma, N and the resulting inflation in one field, so the arithmetic is
    regenerated here rather than copied: a copy is what put the ledger-final N into the manuscript.
    """
    import re
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    pr = json.load(open(os.path.join(root, "development", "benchmark-protocol.json")))
    unc = pr["metrics"]["uncertainty"]
    n = int(re.search(r"N\s*=\s*(\d+)", unc["selection_inflation_basis"]).group(1))
    sig = unc["paired_split_reseed_sd"]
    infl = round(sig * math.sqrt(2 * math.log(n)), 4)
    assert abs(infl - unc["selection_inflation"]) < 5e-4, "protocol inflation does not regenerate"
    return {"sigma": sig, "N": n, "inflation": unc["selection_inflation"],
            "bar": round(0.679 + unc["selection_inflation"], 4),
            "note": "read from development/benchmark-protocol.json; this is what was frozen"}


def oof(X, t, e, fi, n, seed=0):
    """Out-of-fold score for one block, exactly as the frozen arm builds it."""
    s = np.full(n, np.nan)
    for tri, vai in fi:
        if len(tri) < 30 or not len(vai):
            continue
        s[vai] = fitapply(X[tri], t[tri], e[tri], X[vai], seed, ALPHAS)
    return s


def family(co, arms, ii, jj):
    """Every candidate the selection could have taken from a set of scored blocks.

    Seven non-empty subsets of the three modalities, each under two combination rules the campaign
    actually used: the equal-weight rank average that was frozen, and the plain percentile sum
    before re-ranking. Fourteen candidates from three fits, which is why this is affordable.
    """
    keys = list(arms)
    out = {}
    for m in range(1, 1 << len(keys)):
        sel = [keys[k] for k in range(len(keys)) if m >> k & 1]
        tot = sum(co.fold_pct(arms[k]) for k in sel)
        out["rank_avg:" + "+".join(sel)] = cidx(co.fold_pct(tot), ii, jj)
        out["pct_sum:" + "+".join(sel)] = cidx(tot, ii, jj)
    return out


def percase(npz_path, keep, stem2case):
    z = np.load(npz_path, allow_pickle=True)
    names = [str(x) for x in z["names"]]
    X = np.asarray(z["X"], dtype=float)
    bag = {}
    for i, s in enumerate(names):
        c = stem2case.get(s)
        if c:
            bag.setdefault(c, []).append(i)
    have = [c for c in keep if c in bag]
    if len(have) < 0.9 * len(keep):
        return None, len(have)
    mu = X.mean(0)
    return np.vstack([X[bag[c]].mean(0) if c in bag else mu for c in keep]), len(have)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("protocol", "root", "titan", "sp-dir", "omics-dir", "dimaf-dir", "fm-dir", "out"):
        ap.add_argument("--" + f, required=True)
    ap.add_argument("--perms", type=int, default=200)
    a = ap.parse_args()

    proto = json.load(open(a.protocol))
    assert proto["protocol_status"] == "frozen"
    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    n = len(co.keep)
    ii, jj = cpairs(co.t, co.e)
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    CLIN = np.hstack([dc["age"], dc["fem"], dc["stage"]])
    BLOCKS = {"slide": co.T, "omics": P, "clinical": CLIN}

    rep = {"artifact_type": "s13_selection_null_and_encoder",
           "phase_of_origin": "frozen_post_hoc",
           "why": "two experiments named by cold panel round 2; neither changes a reported number",
           "n": n, "events": int(co.e.sum()), "comparable_pairs": int(ii.size)}

    # ================================================================= A. the selection null
    t0 = time.time()
    real = {k: oof(X, co.t, co.e, fi, n, 0) for k, X in BLOCKS.items()}
    real_fam = family(co, real, ii, jj)
    per_fit = (time.time() - t0) / 3.0
    print("  one block fit takes %.1f s; %d permutations will take about %.0f min"
          % (per_fit, a.perms, a.perms * 3 * per_fit / 60.0), file=sys.stderr, flush=True)

    rng = np.random.default_rng(20260818)
    null_max, null_single, null_head = [], [], []
    for b in range(a.perms):
        # permute the (time, event) PAIR, which destroys association and keeps censoring intact
        p = rng.permutation(n)
        tb, eb = co.t[p], co.e[p]
        iib, jjb = cpairs(tb, eb)
        arms = {k: oof(X, tb, eb, fi, n, 0) for k, X in BLOCKS.items()}
        fam = family(co, arms, iib, jjb)
        null_max.append(max(fam.values()))
        null_single.append(fam["rank_avg:slide"])
        null_head.append(fam["rank_avg:slide+omics+clinical"])
        if (b + 1) % 10 == 0:
            print("    perm %3d/%d  max %.4f" % (b + 1, a.perms, null_max[-1]),
                  file=sys.stderr, flush=True)

    nm, ns = np.array(null_max), np.array(null_single)
    N_HONEST = 206                      # the ledger stratum a headline arm actually competes in
    sigma_null = float(ns.std(ddof=1))
    gauss = sigma_null * math.sqrt(2 * math.log(N_HONEST))
    rep["A_selection_null"] = {
        "permutations": a.perms,
        "candidates_per_permutation": len(real_fam),
        "what_was_permuted": "the (time, event) pair across cases; folds, features and the "
                             "estimator are untouched",
        "null_single_candidate": {"mean": round(float(ns.mean()), 4),
                                  "sd": round(sigma_null, 4),
                                  "note": "THIS is the sigma the max-of-N formula wanted: the SD "
                                          "of one candidate's own absolute out-of-fold score with "
                                          "no signal present"},
        "null_max_over_the_family": {
            "mean": round(float(nm.mean()), 4), "sd": round(float(nm.std(ddof=1)), 4),
            "q95": round(float(np.quantile(nm, 0.95)), 4),
            "max": round(float(nm.max()), 4),
            "note": "needs no independence assumption, because the candidates were correlated "
                    "in exactly the way the real search was"},
        "bars": {
            # MISLABELLED, see the module docstring: these are the ledger-final figures, not
            # the frozen ones. Kept under its original name so this source cannot drift away from
            # the artifact it already emitted; the correct pair is the next entry.
            "as_preregistered": {"sigma": 0.0073, "N": 269, "inflation": 0.0244, "bar": 0.7034},
            "as_frozen_in_the_protocol": _frozen_bar(),
            "gaussian_with_the_null_sigma": {"sigma": round(sigma_null, 4), "N": N_HONEST,
                                             "inflation": round(gauss, 4),
                                             "bar": round(0.679 + gauss, 4)},
            "empirical_max_q95": {"inflation": round(float(np.quantile(nm, 0.95)) - 0.5, 4),
                                  "bar": round(0.679 + float(np.quantile(nm, 0.95)) - 0.5, 4),
                                  "note": "the null max is centred on 0.5, so its excess over "
                                          "chance is the inflation a search of this family buys"},
        },
        "real_family_for_reference": {k: round(v, 4) for k, v in sorted(
            real_fam.items(), key=lambda r: -r[1])[:6]},
    }
    for name, d in rep["A_selection_null"]["bars"].items():
        print("  bar %-32s %.4f   0.7212 %s"
              % (name, d["bar"], "CLEARS" if 0.7212 > d["bar"] else "FAILS"), file=sys.stderr)

    # ================================================================= B. the encoder confound
    import csv as _csv
    stem2case = {r["slide_id"].replace(".svs", ""): r["case_id"]
                 for r in _csv.DictReader(open(a.root + "/tcga_blca_meta.csv"))}
    enc = {"titan": co.T}
    for p in sorted(glob.glob(os.path.join(a.fm_dir, "*.npz"))):
        M, cov = percase(p, co.keep, stem2case)
        if M is not None:
            enc[os.path.basename(p)[:-4]] = M
    om = oof(P, co.t, co.e, fi, n, 0)
    cl = oof(CLIN, co.t, co.e, fi, n, 0)
    swap = {}
    for tag, X in enc.items():
        sl = oof(X, co.t, co.e, fi, n, 0)
        full = co.fold_pct(co.fold_pct(sl) + co.fold_pct(om) + co.fold_pct(cl))
        swap[tag] = {"slide_alone": round(cidx(co.fold_pct(sl), ii, jj), 4),
                     "full_method": round(cidx(full, ii, jj), 4)}
        print("  encoder %-22s slide %.4f  ModRank %.4f"
              % (tag, swap[tag]["slide_alone"], swap[tag]["full_method"]),
              file=sys.stderr, flush=True)
    fulls = [v["full_method"] for v in swap.values()]
    alones = [v["slide_alone"] for v in swap.values()]
    rep["B_encoder_swap"] = {
        "question": "is the lead the METHOD or the 2025 encoder? Hold the method fixed, swap the "
                    "slide block, keep omics and clinical identical.",
        "encoders": swap,
        "slide_arm_span": round(max(alones) - min(alones), 4),
        "full_method_span": round(max(fulls) - min(fulls), 4),
        "full_method_min": round(min(fulls), 4), "full_method_max": round(max(fulls), 4),
        "incumbent_published": 0.679,
        "encoders_whose_full_method_still_beats_the_incumbent":
            int(sum(1 for v in fulls if v > 0.679)),
        "reading": "if the full method clears the incumbent on every encoder, the lead is not "
                   "TITAN. If it clears on one, the lead is TITAN and the paper must say so.",
    }
    print("  slide-arm span %.4f  ->  full-method span %.4f  (%d of %d encoders still beat 0.679)"
          % (rep["B_encoder_swap"]["slide_arm_span"], rep["B_encoder_swap"]["full_method_span"],
             rep["B_encoder_swap"]["encoders_whose_full_method_still_beats_the_incumbent"],
             len(fulls)), file=sys.stderr)

    json.dump(rep, open(a.out, "w"), indent=1)
    print("  wrote %s" % a.out, file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
