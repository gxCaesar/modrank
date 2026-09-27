#!/usr/bin/env python3
"""The inflation D with intervals in three studies (protocol B) and three further fitted fusions in all
five (protocol D).

Implements experiments/20260927-field-inflation/protocol.md, sections B and D, committed before the
run. Post-freeze and exploratory; nothing here changes ModRank or any committed comparison.

B. For bladder, head and neck and stomach, the studies whose clinical files carry a usable grade,
   D = [C(ModRank with grade) - C(age, sex, grade)] - [C(ModRank) - C(age, sex, stage)], with all
   four concordances recomputed on each of 6,000 case resamples (seed 0) drawn once per study. The
   four point concordances are already committed in the five-study reproduction and must be
   reproduced first.

D. Three fitted fusions of the same three arms (slide, transcriptome, age + sex + stage), each fitted
   only inside the training fold on the arms' inner out-of-fold percentiles, the inputs the
   committed stacking uses:
     simplex      three non-negative weights summing to one, the 66 points of a 0.1 grid, chosen by
                  concordance on the inner out-of-fold percentiles (first maximum in grid order)
     interaction  the stacked Cox with the three pairwise products of arm percentiles added, ridge 1.0
     gated        separate simplex weights inside each tertile of the clinical arm's inner
                  out-of-fold percentile, applied by the validation patient's clinical tertile
   Each is compared with the rank average by the paired case bootstrap (6,000, seed 0), 15
   comparisons. Before any is scored, the rank average must reproduce its committed value in every
   study, and plain stacking rebuilt from the same inner out-of-fold percentiles must reproduce the
   committed stacked value, which checks the shared path.
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
from blca_common import cidx, cox_fit, cpairs, fitapply, pct                  # noqa: E402
from s5_stage_five_cohorts import COHORTS, fold_pct, load_cohort              # noqa: E402
from s32_five_cohort_fusion import SEED, arm, boot                            # noqa: E402

REPS = 6000
B_STUDIES = ("blca", "hnsc", "stad")
GRID = [(i / 10.0, j / 10.0, (10 - i - j) / 10.0) for i in range(11) for j in range(11 - i)]
assert len(GRID) == 66


def best_weights(Z, t, e):
    """The grid point with the highest concordance on (Z, t, e); the first in grid order on a tie."""
    ii, jj = cpairs(t, e)
    cs = np.array([cidx(Z @ np.array(w), ii, jj) for w in GRID])
    cs = np.where(np.isfinite(cs), cs, -np.inf)
    return GRID[int(np.argmax(cs))], bool(np.isfinite(cs).any())


def tertile(x, cut):
    return (x > cut[0]).astype(int) + (x > cut[1]).astype(int)


def fusions(blocks, per_arm, co, seed=SEED):
    """Plain stacking (the committed known answer) and the three new fusions, from one inner split."""
    n = len(co["keep"])
    names = list(blocks)
    assert names == ["slide", "omics", "clinical"]
    out = {k: np.full(n, np.nan) for k in ("stacked", "simplex", "interaction", "gated")}
    record = []
    for k, (tri, vai) in enumerate(co["fi"]):
        if len(tri) < 30 or not len(vai):
            continue
        m = len(tri)
        ttr, etr = co["t"][tri], co["e"][tri]
        # s32.stacked's inner out-of-fold percentiles, line for line
        inner = np.array_split(np.random.default_rng(seed).permutation(m), 3)
        ioof = {b: np.full(m, np.nan) for b in names}
        for b in names:
            Xt = blocks[b][tri]
            for ite in inner:
                itr = np.setdiff1d(np.arange(m), ite)
                ioof[b][ite] = pct(fitapply(Xt[itr], ttr[itr], etr[itr], Xt[ite], seed))
        Zi = np.column_stack([ioof[b] for b in names])
        Zv = np.column_stack([per_arm[b][vai] for b in names])

        w = cox_fit(Zi - Zi.mean(0), ttr, etr, 1.0)
        out["stacked"][vai] = Zv @ w

        ws, _ = best_weights(Zi, ttr, etr)
        out["simplex"][vai] = Zv @ np.array(ws)

        def inter(Z):
            return np.column_stack([Z, Z[:, 0] * Z[:, 1], Z[:, 0] * Z[:, 2], Z[:, 1] * Z[:, 2]])
        Fi = inter(Zi)
        wi = cox_fit(Fi - Fi.mean(0), ttr, etr, 1.0)
        out["interaction"][vai] = inter(Zv) @ wi

        cut = np.quantile(Zi[:, 2], [1.0 / 3.0, 2.0 / 3.0])
        ti, tv = tertile(Zi[:, 2], cut), tertile(Zv[:, 2], cut)
        wg, informative = [], []
        for g in range(3):
            mg = ti == g
            wgt, ok = best_weights(Zi[mg], ttr[mg], etr[mg])
            wg.append(wgt)
            informative.append(ok)
        sv = np.full(len(vai), np.nan)
        for g in range(3):
            sv[tv == g] = Zv[tv == g] @ np.array(wg[g])
        out["gated"][vai] = sv
        record.append({"fold": k, "simplex_weights": ws,
                       "interaction_coefficients": [round(float(x), 4) for x in wi],
                       "gated_weights_by_clinical_tertile": wg,
                       "gated_tertile_had_comparable_pairs": informative,
                       "validation_cases_per_tertile": [int((tv == g).sum()) for g in range(3)]})
    return {k: fold_pct(v, co["assign"]) for k, v in out.items()}, record


def d_interval(arms4, t, e, reps=REPS, seed=0):
    """D and both added values, all four concordances recomputed on each resample of cases."""
    fg, asg, fs, ass = arms4
    rng = np.random.default_rng(seed)
    n = len(t)
    bw, bs, bd = [], [], []
    for _ in range(reps):
        b = rng.choice(n, size=n, replace=True)
        ii, jj = cpairs(t[b], e[b])
        if not ii.size:
            continue
        aw = cidx(fg[b], ii, jj) - cidx(asg[b], ii, jj)
        as_ = cidx(fs[b], ii, jj) - cidx(ass[b], ii, jj)
        bw.append(aw); bs.append(as_); bd.append(aw - as_)

    def ci(x):
        x = np.asarray(x)
        return {"ci95": [round(float(np.quantile(x, .025)), 4), round(float(np.quantile(x, .975)), 4)],
                "p_two_sided": round(float(2 * min((x <= 0).mean(), (x >= 0).mean())), 4),
                "reps": int(x.size)}
    return {"added_over_grade_boot": ci(bw), "added_over_stage_boot": ci(bs), "D_boot": ci(bd)}


def holm(ps):
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    adj, run = [0.0] * len(ps), 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (len(ps) - r) * ps[i]))
        adj[i] = round(run, 4)
    return adj


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("titan", "meta-dir", "clin-dir", "pan-dir", "sig", "frozen", "frozen-fusion", "out"):
        ap.add_argument("--" + f, required=True)
    a = ap.parse_args()
    t0 = time.time()

    frozen = json.load(open(a.frozen))["cohorts"]
    frozen_fusion = json.load(open(a.frozen_fusion))["cohorts"]
    d0 = pickle.load(open(a.titan, "rb"))
    E = np.asarray(d0["embeddings"])
    stem2emb = {}
    for i, f in enumerate(d0["filenames"]):
        s = str(f)
        for suf in (".svs", ".h5", ".pt"):
            if s.endswith(suf):
                s = s[:-len(suf)]
        stem2emb[s] = E[i]

    out = {"artifact_type": "s34_five_study_extensions", "phase_of_origin": "post_freeze_2026-09-27",
           "protocol": "experiments/20260927-field-inflation/protocol.md, sections B and D",
           "bootstrap": {"replicates": REPS, "unit": "cases", "seed": 0},
           "grid": "66 simplex points, step 0.1, first maximum in grid order",
           "inflation_B": {}, "fusions_D": {}}
    drift = {}
    for c in COHORTS:
        paths = (os.path.join(a.meta_dir, "tcga_%s.csv" % c),
                 os.path.join(a.clin_dir, "tcga_%s_clinical.csv" % c),
                 os.path.join(a.pan_dir, "splits", c),
                 os.path.join(a.pan_dir, "rna_%s.csv" % c))
        if not all(os.path.exists(p) for p in paths):
            out["fusions_D"][c] = {"status": "inputs missing"}
            continue
        co = load_cohort(paths[0], paths[1], paths[2], paths[3], a.sig, stem2emb)
        ii, jj = cpairs(co["t"], co["e"])

        def sc(v):
            return round(float(cidx(v, ii, jj)), 4)

        AS = np.hstack([co["age"], co["fem"]])
        CLIN = np.hstack([AS, co["stage"]])
        blocks = {"slide": co["T"], "omics": co["P"], "clinical": CLIN}
        per = {b: arm(X, co) for b, X in blocks.items()}
        ours = fold_pct(per["slide"] + per["omics"] + per["clinical"], co["assign"])
        fz = frozen[c]
        if sc(ours) != fz["OURS_wsi_omics_age_sex_stage"]:
            drift[c + "/ModRank"] = (sc(ours), fz["OURS_wsi_omics_age_sex_stage"])

        # ---------------------------------------------------------------- B
        if c in B_STUDIES:
            asg = arm(np.hstack([AS, co["grade"]]), co)
            full_grade = fold_pct(per["slide"] + per["omics"] + asg, co["assign"])
            pts = {"C_ModRank_with_grade": sc(full_grade), "C_age_sex_grade": sc(asg),
                   "C_ModRank": sc(ours), "C_age_sex_stage": sc(per["clinical"])}
            want = {"C_ModRank_with_grade": fz["control_same_arm_with_GRADE_instead"],
                    "C_age_sex_grade": fz["arms"]["age_sex_grade"],
                    "C_ModRank": fz["OURS_wsi_omics_age_sex_stage"],
                    "C_age_sex_stage": fz["arms"]["age_sex_stage"]}
            for q in pts:
                if pts[q] != want[q]:
                    drift["%s/%s" % (c, q)] = (pts[q], want[q])
            aw = cidx(full_grade, ii, jj) - cidx(asg, ii, jj)
            as_ = cidx(ours, ii, jj) - cidx(per["clinical"], ii, jj)
            row = {"n": len(co["keep"]), "events": int(co["e"].sum()),
                   "grade_levels": co["grade_levels"], **pts,
                   "added_over_grade": round(float(aw), 4), "added_over_stage": round(float(as_), 4),
                   "D": round(float(aw - as_), 4)}
            row.update(d_interval((full_grade, asg, ours, per["clinical"]), co["t"], co["e"]))
            out["inflation_B"][c] = row
            print("B %-9s D %+.4f %s | over grade %+.4f | over stage %+.4f  [%.0f s]"
                  % (c, row["D"], row["D_boot"]["ci95"], row["added_over_grade"],
                     row["added_over_stage"], time.time() - t0), file=sys.stderr, flush=True)

        # ---------------------------------------------------------------- D
        fus, record = fusions(blocks, per, co)
        if sc(fus["stacked"]) != frozen_fusion[c]["points"]["stacked"]:
            drift[c + "/stacked"] = (sc(fus["stacked"]), frozen_fusion[c]["points"]["stacked"])
        row = {"n": len(co["keep"]), "events": int(co["e"].sum()),
               "points": {"ModRank": sc(ours), "stacked_known_answer": sc(fus["stacked"]),
                          **{k: sc(fus[k]) for k in ("simplex", "interaction", "gated")}},
               "per_fold": record}
        for k in ("simplex", "interaction", "gated"):
            row["ModRank_minus_" + k] = boot(ours, fus[k], co["t"], co["e"])
        out["fusions_D"][c] = row
        print("D %-9s ModRank %.4f | simplex %.4f %s | interaction %.4f %s | gated %.4f %s  [%.0f s]"
              % (c, row["points"]["ModRank"], row["points"]["simplex"],
                 row["ModRank_minus_simplex"]["ci95"], row["points"]["interaction"],
                 row["ModRank_minus_interaction"]["ci95"], row["points"]["gated"],
                 row["ModRank_minus_gated"]["ci95"], time.time() - t0), file=sys.stderr, flush=True)

    if drift:
        print("KNOWN ANSWER DRIFT, refusing to write: %s" % drift, file=sys.stderr)
        return 2

    keys = [(c, k) for c in out["fusions_D"] for k in ("simplex", "interaction", "gated")
            if "points" in out["fusions_D"][c]]
    adj = holm([out["fusions_D"][c]["ModRank_minus_" + k]["p_two_sided"] for c, k in keys])
    for (c, k), p in zip(keys, adj):
        out["fusions_D"][c]["ModRank_minus_" + k]["p_holm_over_15"] = p
    out["summary"] = {
        "B_D_by_study": {c: out["inflation_B"][c]["D"] for c in out["inflation_B"]},
        "B_D_interval_excludes_zero": sorted(c for c, v in out["inflation_B"].items()
                                             if v["D_boot"]["ci95"][0] > 0),
        "D_comparisons": len(keys),
        "D_fitted_fusion_ahead_interval_excluding_zero": [
            "%s/%s" % (c, k) for c, k in keys
            if out["fusions_D"][c]["ModRank_minus_" + k]["ci95"][1] < 0],
        "D_rank_average_ahead_interval_excluding_zero": [
            "%s/%s" % (c, k) for c, k in keys
            if out["fusions_D"][c]["ModRank_minus_" + k]["ci95"][0] > 0],
    }
    out["runtime_seconds"] = round(time.time() - t0, 1)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1)
    print("wrote %s (%.0f s)" % (a.out, time.time() - t0), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
