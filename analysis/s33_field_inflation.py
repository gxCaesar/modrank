#!/usr/bin/env python3
"""The added-value inflation D in further multimodal architectures on TCGA-BLCA (protocol A).

Implements experiments/20260927-field-inflation/protocol.md, section A, which was committed before
the runs. Post-freeze and specified after the confirmatory run; nothing here touches ModRank.

For each architecture, its five retrainings are percentile-ranked within fold and averaged, and the
average is ranked again (M). M is combined with each clinical reference by the paper's rank rule,
M + R = fold percentile of (M + R), with the canonical grade and stage arms read from the committed
per-case vectors. Added value over R is C(M + R) - C(R), and D is the added value over grade minus
the added value over stage. All four concordances of each D are recomputed on each of 6,000 case
resamples drawn with seed 20260911, the draws analysis/s19 used for the paper's existing D.

KNOWN ANSWER FIRST. From the committed SurvPath vectors, through this code path, the script must
reproduce SurvPath's committed block in analysis/s19's output (added values +0.0670 and +0.0253, D
+0.0417, every concordance and all three bootstrap intervals) before it reports anything; otherwise
it exits 2 and writes nothing. The per-case table stores percentiles to five decimals, so each one is
first restored to its exact value, rank / (fold size - 1), which also checks that every committed
vector is a within-fold permutation of ranks.
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
from blca_common import cidx, cpairs, pct                                      # noqa: E402

RNG_SEED = 20260911          # s19's draws
KNOWN_FIELDS = ("C_weak", "C_stage", "C_M_plus_weak", "C_M_plus_stage", "added_over_weak",
                "added_over_stage", "D", "D_boot", "added_over_weak_boot", "added_over_stage_boot")


def fold_pct(v, fold):
    """blca_common.Cohort.fold_pct, for vectors indexed like the per-case table."""
    out = np.full(len(v), np.nan)
    for k in np.unique(fold):
        m = fold == k
        out[m] = pct(v[m])
    return out


def exact(v, fold):
    """Undo the table's five-decimal rounding: a fold percentile is rank / (fold size - 1)."""
    out = np.full(len(v), np.nan)
    for k in np.unique(fold):
        m = fold == k
        r = v[m] * (m.sum() - 1)
        assert np.abs(r - np.round(r)).max() < 1e-3, "not a fold percentile"
        assert sorted(np.round(r).astype(int)) == list(range(m.sum())), "not a permutation of ranks"
        out[m] = np.round(r) / (m.sum() - 1.0)
    return out


def load_model(run_dir, idx, n, t, e):
    """One risk per case from each seed's split_k_results.pkl, as analysis/s19.load_runs reads them."""
    seeds = sorted(d for d in os.listdir(run_dir) if d.startswith("seed"))
    per_seed = []
    for sd in seeds:
        risk, fold = np.full(n, np.nan), np.full(n, -1)
        for dirpath, _, files in os.walk(os.path.join(run_dir, sd)):
            for f in files:
                if f.startswith("split_") and f.endswith("_results.pkl"):
                    k = int(f.split("_")[1])
                    for c, v in pickle.load(open(os.path.join(dirpath, f), "rb")).items():
                        if c in idx:
                            i = idx[c]
                            risk[i], fold[i] = float(v["risk"]), k
                            assert abs(float(v["time"]) - t[i]) < 1e-4, "%s: survival time differs" % c
                            assert (float(v["censorship"]) == 0.0) == (e[i] == 1), "%s: event differs" % c
        per_seed.append((sd, risk, fold))
    return per_seed


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--percase", required=True, help="percase-canonical-vectors.jsonl")
    ap.add_argument("--committed", required=True, help="unified-fusion-and-added-value.json (s19)")
    ap.add_argument("--runs", required=True, help="directory holding one subdirectory per model")
    ap.add_argument("--models", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=6000)
    a = ap.parse_args()
    t0 = time.time()

    rows = [json.loads(x) for x in open(a.percase) if x.strip()]
    n = len(rows)
    idx = {r["case_id"]: i for i, r in enumerate(rows)}
    t = np.array([r["months"] for r in rows], float)
    e = np.array([r["event"] for r in rows], float)
    fold = np.array([r["fold"] for r in rows], int)
    zg = exact(np.array([r["clinical_grade"] for r in rows], float), fold)
    zc = exact(np.array([r["clinical_stage"] for r in rows], float), fold)
    ii, jj = cpairs(t, e)
    rng = np.random.default_rng(RNG_SEED)
    draws = [rng.choice(n, size=n, replace=True) for _ in range(a.reps)]

    def measure(alone):
        weak, stage = fold_pct(alone + zg, fold), fold_pct(alone + zc, fold)
        C = lambda v: cidx(v, ii, jj)
        pt = {"C_alone": C(alone), "C_weak": C(zg), "C_stage": C(zc),
              "C_M_plus_weak": C(weak), "C_M_plus_stage": C(stage)}
        pt["added_over_weak"] = pt["C_M_plus_weak"] - pt["C_weak"]
        pt["added_over_stage"] = pt["C_M_plus_stage"] - pt["C_stage"]
        pt["D"] = pt["added_over_weak"] - pt["added_over_stage"]
        bw, bs, bd = [], [], []
        for b in draws:
            i2, j2 = cpairs(t[b], e[b])
            cw, cs = cidx(zg[b], i2, j2), cidx(zc[b], i2, j2)
            mw, ms = cidx(weak[b], i2, j2), cidx(stage[b], i2, j2)
            bw.append(mw - cw); bs.append(ms - cs); bd.append((mw - cw) - (ms - cs))

        def ci(x):
            x = np.asarray(x)
            return {"ci95": [round(float(np.quantile(x, .025)), 4), round(float(np.quantile(x, .975)), 4)],
                    "p_two_sided": round(float(2 * min((x <= 0).mean(), (x >= 0).mean())), 4)}
        return {**{k: round(float(v), 4) for k, v in pt.items()},
                "added_over_weak_boot": ci(bw), "added_over_stage_boot": ci(bs), "D_boot": ci(bd)}

    # ---- known answer: SurvPath, from the committed vectors, through this code path
    sp = exact(np.array([r["survpath"] for r in rows], float), fold)
    k = measure(sp)
    want = json.load(open(a.committed))["added_value_inflation"]["constructions"]["survpath"]
    bad = {q: (k[q], want[q]) for q in KNOWN_FIELDS if k[q] != want[q]}
    if k["C_alone"] != want["C_M_alone"]:
        bad["C_alone"] = (k["C_alone"], want["C_M_alone"])
    if bad:
        print("s33: SurvPath known answer not reproduced: %s" % bad, file=sys.stderr)
        return 2
    print("known answer reproduced: SurvPath over grade %+.4f, over stage %+.4f, D %+.4f"
          % (k["added_over_weak"], k["added_over_stage"], k["D"]), file=sys.stderr)

    out = {"artifact_type": "s33_field_inflation", "phase_of_origin": "post_freeze_2026-09-27",
           "protocol": "experiments/20260927-field-inflation/protocol.md, section A",
           "bootstrap": {"replicates": a.reps, "unit": "cases", "seed": RNG_SEED},
           "known_answer_survpath": k, "models": {}}
    for m in a.models:
        d = os.path.join(a.runs, m)
        if not os.path.isdir(d):
            out["models"][m] = {"status": "not run"}
            continue
        per_seed = load_model(d, idx, n, t, e)
        complete = [(sd, r, f) for sd, r, f in per_seed if np.isfinite(r).all()]
        for sd, r, f in complete:
            assert (f == fold).all(), "%s %s: fold assignment differs from the released folds" % (m, sd)
        if len(complete) < len(per_seed) or not complete:
            out["models"][m] = {"status": "incomplete",
                                "seeds_complete": [sd for sd, _, _ in complete],
                                "seeds_found": [sd for sd, _, _ in per_seed]}
            continue
        M = fold_pct(np.vstack([fold_pct(r, fold) for _, r, _ in complete]).mean(0), fold)
        res = measure(M)
        res["seeds"] = [sd for sd, _, _ in complete]
        res["per_seed_C_alone"] = [round(float(cidx(fold_pct(r, fold), ii, jj)), 4)
                                   for _, r, _ in complete]
        out["models"][m] = res
        print("%-24s alone %.4f | over grade %+.4f %s | over stage %+.4f %s | D %+.4f %s"
              % (m, res["C_alone"], res["added_over_weak"], res["added_over_weak_boot"]["ci95"],
                 res["added_over_stage"], res["added_over_stage_boot"]["ci95"], res["D"],
                 res["D_boot"]["ci95"]), file=sys.stderr, flush=True)

    # Holm across the new D tests (protocol A)
    done = {m: v for m, v in out["models"].items() if "D_boot" in v}
    order = sorted(done, key=lambda m: done[m]["D_boot"]["p_two_sided"])
    run = 0.0
    for i, m in enumerate(order):
        run = max(run, min(1.0, (len(order) - i) * done[m]["D_boot"]["p_two_sided"]))
        done[m]["D_p_holm"] = round(run, 4)
    out["summary"] = {
        "models_scored": len(done),
        "D_positive": sorted(m for m in done if done[m]["D"] > 0),
        "D_interval_excludes_zero": sorted(m for m in done if done[m]["D_boot"]["ci95"][0] > 0),
        "added_over_stage_interval_covers_zero": sorted(
            m for m in done if done[m]["added_over_stage_boot"]["ci95"][0] <= 0 <= done[m]["added_over_stage_boot"]["ci95"][1]),
        "not_run_or_incomplete": sorted(m for m in out["models"] if m not in done),
    }
    out["runtime_seconds"] = round(time.time() - t0, 1)
    json.dump(out, open(a.out, "w"), indent=1)
    print("wrote %s (%.0f s)" % (a.out, time.time() - t0), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
