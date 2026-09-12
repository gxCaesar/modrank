#!/usr/bin/env python3
"""GEO amendment 01: the evidence-gated variant on GSE32894 and GSE48075 (not GSE31684).

The definition is the pre-declared secondary of the external protocol frozen at 9e60ada, copied
unchanged: inside each outer
training fold, an arm enters the equal-weight rank average only if its inner out-of-fold concordance
exceeds the 95th percentile of 1,000 permutations of the training outcomes against that fixed score
(rng 10000 + 100*repeat + fold); if none passes, all enter. Same splits, estimator, comparators and
floor as s24, recomputed here so this file never reads s24's results for these cohorts.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import re
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s24_geo_external as g                                     # noqa: E402
from blca_common import cidx, cpairs, fitapply, pct              # noqa: E402

PERMS = 1000


def gated_two_arm(blocks, t, e, folds, seed, log, tag):
    n = len(t)
    z = [g.fpct(g.oof(X, t, e, folds, seed), folds) for X in blocks]
    out = np.full(n, np.nan)
    for k, va in enumerate(folds):
        tr = np.setdiff1d(np.arange(n), va)
        m = len(tr)
        i2, j2 = cpairs(t[tr], e[tr])
        rng = np.random.default_rng(10000 + 100 * seed + k)
        keep = []
        for b, X in enumerate(blocks):
            io = np.full(m, np.nan)
            for ite in np.array_split(np.random.default_rng(seed).permutation(m), 3):
                itr = np.setdiff1d(np.arange(m), ite)
                io[ite] = pct(fitapply(X[tr][itr], t[tr][itr], e[tr][itr], X[tr][ite], seed))
            c_obs = cidx(io, i2, j2)
            null = []
            for _ in range(PERMS):
                pr = rng.permutation(m)
                null.append(cidx(io, *cpairs(t[tr][pr], e[tr][pr])))
            q95 = float(np.quantile(null, 0.95))
            log.append({"cohort": tag, "repeat": seed, "fold": k, "arm": ["transcriptome", "clinical_stage"][b],
                        "inner_oof_c": round(c_obs, 4), "null_q95": round(q95, 4), "passes": bool(c_obs > q95)})
            if c_obs > q95:
                keep.append(b)
        keep = keep or [0, 1]
        acc = np.zeros(len(va))
        for b in (0, 1):
            if b in keep:
                acc = acc + z[b][va]
        out[va] = acc
    return g.fpct(out, folds)


def run(name, P, CS, t, e, log):
    res = g.evaluate(name, P, CS, None, t, e)             # ModRank, comparators, floor, as s24
    ii, jj = cpairs(t, e)
    rng = np.random.default_rng(0)
    splits = [g.strat_folds(e, g.K, rng) for _ in range(g.REPEATS)]
    per = [cidx(gated_two_arm([P, CS], t, e, folds, r, log, name), ii, jj) for r, folds in enumerate(splits)]
    res["cv"]["gated"] = {"mean": round(float(np.mean(per)), 4), "sd": round(float(np.std(per)), 4)}
    comps = ("clinical_stage", "transcriptome", "concatenated", "stacked", "modrank")
    margins = {c: round(res["cv"]["gated"]["mean"] - res["cv"][c]["mean"], 4) for c in comps}
    res["gated_decision"] = {"margins": margins, "bar": res["floor"]["bar"],
                             "succeeds": all(v > res["floor"]["bar"] for v in margins.values())}
    arms = collections.Counter((x["arm"], x["passes"]) for x in log if x["cohort"] == name)
    res["gate_pass_counts"] = {"%s_%s" % (a_, "pass" if p_ else "fail"): c for (a_, p_), c in arms.items()}
    print("  %s gated %.4f | margins %s | bar %.4f" % (name, res["cv"]["gated"]["mean"], margins,
                                                      res["floor"]["bar"]), file=sys.stderr, flush=True)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("geo-dir", "signatures", "out"):
        ap.add_argument("--" + f, required=True)
    a = ap.parse_args()
    t0 = time.time()
    ann = {"GPL6947": g.read_annot(os.path.join(a.geo_dir, "GPL6947.annot.gz"))}
    log, rep = [], {"artifact_type": "s24b_geo_gated_amendment_01", "phase_of_origin": "post_freeze_2026-09-11",
                    "amendment": "experiments/20260911-geo-external/amendment-01.md (0d89861)", "cohorts": {}}

    def num(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return np.nan

    s, ch, pr, X = g.read_series(os.path.join(a.geo_dir, "GSE32894_series_matrix.txt.gz"))
    G, genes = g.gene_matrix(pr, X, ann["GPL6947"])
    P, _ = g.pathway_matrix(G, genes, a.signatures)
    t = np.array([num(ch["time_to_dod_(months)"].get(x)) for x in s])
    ev = [ch["dod_event_(yes/no)"].get(x, "") for x in s]
    exc = np.array([ch["dod_excluded_(yes/no)"].get(x, "") == "yes" for x in s])
    e = np.array([1.0 if v == "yes" else 0.0 for v in ev])
    age = np.array([num(ch["age"].get(x)) for x in s])[:, None]
    fem = np.array([1.0 if ch["gender"].get(x, "") == "F" else 0.0 for x in s])[:, None]
    stage, _ = g.onehot([ch["tumor_stage"].get(x, "") for x in s])
    ok = np.isfinite(t) & (t > 0) & ~exc & np.array([v in ("yes", "no") for v in ev]) & np.isfinite(age[:, 0])
    rep["cohorts"]["GSE32894"] = run("GSE32894", P[ok], np.hstack([age, fem, stage])[ok], t[ok], e[ok], log)

    s, ch, pr, X = g.read_series(os.path.join(a.geo_dir, "GSE48075_series_matrix.txt.gz"))
    G, genes = g.gene_matrix(pr, X, ann["GPL6947"])
    P, _ = g.pathway_matrix(G, genes, a.signatures)
    t = np.array([num(ch["survival (mo)"].get(x)) for x in s])
    cen = [ch["dss censor"].get(x, "") for x in s]
    e = np.array([1.0 if c == "uncensored" else 0.0 for c in cen])
    age = np.array([num(ch["age (at specimen collection)"].get(x)) for x in s])[:, None]
    fem = np.array([1.0 if ch["gender"].get(x, "").lower() == "female" else 0.0 for x in s])[:, None]

    def tcat(v):
        m = re.match(r"^c?T(is|a|\d)", v.strip())
        return ("T" + m.group(1)) if m else ""
    tst, _ = g.onehot([tcat(ch["cstage"].get(x, "")) for x in s])
    ok = np.isfinite(t) & (t > 0) & np.array([c in ("censored", "uncensored") for c in cen]) & np.isfinite(age[:, 0])
    rep["cohorts"]["GSE48075"] = run("GSE48075", P[ok], np.hstack([age, fem, tst])[ok], t[ok], e[ok], log)
    rep["gate_log"] = log
    rep["runtime_seconds"] = round(time.time() - t0, 1)
    json.dump(rep, open(a.out, "w"), indent=1)
    print("wrote %s (%.0f s)" % (a.out, time.time() - t0), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
