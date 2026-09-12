#!/usr/bin/env python3
"""Post-freeze analyses on TCGA-BLCA that a Nature-portfolio reviewer asks for first.

POST-FREEZE AND DESCRIPTIVE. The frozen protocol's primary (0.7212) is not touched. Everything
here is computed on the same 359 cases, the same released folds and the same estimator, after the
2026-08-17 freeze, and is labelled that way in its output.

FOUR QUESTIONS, EACH OF WHICH HAD NO ANSWER ON DISK.

  1. FUSION ON IDENTICAL INPUTS. Does fitting one model per modality and averaging ranks beat the
     two obvious alternatives given exactly the same features, folds, seeds and estimator? The two
     are one penalised Cox on all three blocks concatenated (1,049 columns; its penalty grid is
     extended upward, which is generous to it), and learned fusion weights, fitted on inner
     out-of-fold scores inside each training fold (stacking). No three-modality concatenation
     existed anywhere in this project before this file.

  2. THE EVIDENCE-GATED VARIANT, on the cohort the method was frozen on. An external-validation
     protocol frozen on 2026-09-11 (commit 9e60ada) pre-declares a variant that admits an
     arm only when its inner out-of-fold concordance clears the 95th percentile of a permutation null
     inside the training fold. On TCGA-BLCA every arm should pass, in which case the variant IS the
     frozen method and must reproduce 0.7212. If an arm fails anywhere, the variant is a different
     method here too, and that is reported rather than smoothed over.

  3. HOW MUCH A WEAK CLINICAL REFERENCE INFLATES THE ADDED VALUE of the other modalities:
        D = [C(M + weak) - C(weak)] - [C(M + stage) - C(stage)]
     for three multimodal constructions (ours, SurvPath and PIBD, each given the clinical block the
     way the paper gives it), with the four concordances of each D recomputed on one case resample
     per replicate. Weak = age + sex + grade and stage = age + sex + stage, both from the
     incumbent's own released split files.

  4. THE MODALITY PARAGRAPH'S NUMBERS, recomputed with the AMENDED clinical block. The rank
     correlations and the tightest-tie restriction quoted in the Results came from a development
     atlas built on the pre-amendment 23-column block. They are recomputed here, seed 0, with only
     the clinical block changed.

KNOWN ANSWERS FIRST. The script exits non-zero and writes nothing unless it first reproduces, from
the same inputs and code path, eight values the manuscript already reports: the two clinical
constructions (0.5666, 0.6638), ours with grade and with stage (0.6856, 0.7225), the five-seed
primary (0.7212), the canonical seed-averaged ours (0.7237), and SurvPath and PIBD given the
clinical block (0.6897, 0.6989).

It also writes every canonical per-case vector it uses, because the paired comparisons in the paper
were computed from vectors that existed only in memory until now.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import pickle
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import ALPHAS, Cohort, cidx, cox_fit, cpairs, fitapply, pct  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix                             # noqa: E402
from s6_amend_clinical import boot, dimaf_clinical                            # noqa: E402

SEEDS = (0, 1, 2, 3, 4)
EXT_ALPHAS = ALPHAS + (32768.0,)       # the concatenation's grid, extended in its favour
KNOWN = {"clinical_grade_seed0": 0.5666, "clinical_stage_seed0": 0.6638,
         "ours_grade_seed0": 0.6856, "ours_seed0": 0.7225, "primary_five_seed_mean": 0.7212,
         "ours_canonical": 0.7237, "survpath_plus_clinical_canonical": 0.6897,
         "pibd_plus_clinical_best_val_canonical": 0.6989}
TOL = 5e-5


def load_runs(dirs, keep, idx, suffix=""):
    """One risk per case from split_k_results{suffix}.pkl files spread over several directories."""
    risk, fold = np.full(len(keep), np.nan), np.full(len(keep), -1)
    t, e = np.full(len(keep), np.nan), np.full(len(keep), np.nan)
    for d in dirs:
        for k in range(5):
            f = os.path.join(d, "split_%d_results%s.pkl" % (k, suffix))
            if not os.path.exists(f):
                continue
            for c, v in pickle.load(open(f, "rb")).items():
                if c in idx:
                    i = idx[c]
                    risk[i], fold[i] = float(v["risk"]), k
                    t[i] = float(v["time"])
                    e[i] = 1.0 if float(v["censorship"]) == 0.0 else 0.0
    return risk, fold, t, e


def spearman(a, b):
    ra, rb = pct(a), pct(b)
    return float(np.corrcoef(ra, rb)[0, 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("protocol", "root", "titan", "sp-dir", "omics-dir", "dimaf-dir", "out-dir"):
        ap.add_argument("--" + f, required=True)
    ap.add_argument("--sp-seed-dirs", nargs=4, required=True, help="SurvPath seeds 2-5")
    ap.add_argument("--pibd-dirs", nargs="+", required=True)
    ap.add_argument("--reps", type=int, default=6000)
    ap.add_argument("--perms", type=int, default=1000)
    a = ap.parse_args()

    proto = json.load(open(a.protocol))
    assert tuple(proto["method"]["estimator"]["alpha_grid"]) == ALPHAS, "alpha grid is not the frozen one"
    t0 = time.time()
    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    n, ii, jj = len(co.keep), *cpairs(co.t, co.e)
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    CLIN = np.hstack([dc["age"], dc["fem"], dc["stage"]])
    CLING = np.hstack([dc["age"], dc["fem"], dc["grade"]])
    BLOCKS = {"slide": co.T, "omics": P, "clinical": CLIN}
    print("cohort %d cases, %d events; blocks slide %s omics %s clinical %s"
          % (n, int(co.e.sum()), co.T.shape[1], P.shape[1], CLIN.shape[1]), file=sys.stderr)

    def arm(X, seed, alphas=ALPHAS):
        s = np.full(n, np.nan)
        for tri, vai in fi:
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], seed, alphas)
        return co.fold_pct(s)

    def C(v):
        return cidx(v, ii, jj)

    # ================================================================ per-seed arms, shared by all
    Z = {s: {k: arm(X, s) for k, X in BLOCKS.items()} for s in SEEDS}
    ZG = {s: arm(CLING, s) for s in SEEDS}
    for s in SEEDS:
        print("  seed %d arms fitted (%.0f s)" % (s, time.time() - t0), file=sys.stderr, flush=True)
    ours_s = {s: co.fold_pct(Z[s]["slide"] + Z[s]["omics"] + Z[s]["clinical"]) for s in SEEDS}
    oursg_s = {s: co.fold_pct(Z[s]["slide"] + Z[s]["omics"] + ZG[s]) for s in SEEDS}

    def canon(vs):
        return co.fold_pct(np.vstack(vs).mean(0))

    zt = canon([Z[s]["slide"] for s in SEEDS])
    zo = canon([Z[s]["omics"] for s in SEEDS])
    zc = canon([Z[s]["clinical"] for s in SEEDS])
    zg = canon([ZG[s] for s in SEEDS])
    OURS, OURSG = canon([ours_s[s] for s in SEEDS]), canon([oursg_s[s] for s in SEEDS])
    OURS_NOCLIN = co.fold_pct(zt + zo)

    # the competitors, exactly as s9 builds them
    sp1, f1, _, _ = load_runs([a.sp_dir], co.keep, co.idx)
    sp = [co.fold_pct(sp1, f1)]
    for d in a.sp_seed_dirs:
        r, f, _, _ = load_runs([d], co.keep, co.idx)
        assert np.isfinite(r).all() and (f == co.fold).all(), "SurvPath seed run %s is incomplete" % d
        sp.append(co.fold_pct(r, f))
    SP = co.fold_pct(np.vstack(sp).mean(0))
    SPC, SPG = co.fold_pct(SP + zc), co.fold_pct(SP + zg)
    PIBD = {}
    for suffix, lab in (("", "best_val"), ("_final", "final")):
        r, f, t, e = load_runs(a.pibd_dirs, co.keep, co.idx, suffix)
        assert np.isfinite(r).all(), "PIBD %s covers %d of %d" % (lab, np.isfinite(r).sum(), n)
        assert (f == co.fold).all() and np.allclose(t, co.t) and np.array_equal(e, co.e), (
            "PIBD %s labels or folds differ from the cohort" % lab)
        z = co.fold_pct(r, f)
        PIBD[lab] = {"alone": z, "stage": co.fold_pct(z + zc), "grade": co.fold_pct(z + zg)}

    # ================================================================ known answers, or nothing
    got = {"clinical_grade_seed0": C(ZG[0]), "clinical_stage_seed0": C(Z[0]["clinical"]),
           "ours_grade_seed0": C(oursg_s[0]), "ours_seed0": C(ours_s[0]),
           "primary_five_seed_mean": float(np.mean([C(ours_s[s]) for s in SEEDS])),
           "ours_canonical": C(OURS), "survpath_plus_clinical_canonical": C(SPC),
           "pibd_plus_clinical_best_val_canonical": C(PIBD["best_val"]["stage"])}
    bad = {k: (round(v, 4), KNOWN[k]) for k, v in got.items() if abs(round(v, 4) - KNOWN[k]) > TOL}
    for k in KNOWN:
        print("  known answer %-40s %.4f (want %.4f)" % (k, got[k], KNOWN[k]), file=sys.stderr)
    if bad:
        print(json.dumps({"status": "error", "error_code": "known_answers_not_reproduced",
                          "which": bad}), file=sys.stderr)
        return 2

    # ================================================================ 1-2. fusion on identical inputs
    XCAT = np.hstack([co.T, P, CLIN])
    concat_s = {s: arm(XCAT, s, EXT_ALPHAS) for s in SEEDS}
    print("  concatenation fitted (%.0f s)" % (time.time() - t0), file=sys.stderr, flush=True)

    stacked_s, gated_s, weights, gate_log = {}, {}, {}, []
    for s in SEEDS:
        st, gt = np.full(n, np.nan), np.full(n, np.nan)
        for k, (tri, vai) in enumerate(fi):
            m = len(tri)
            inner = np.array_split(np.random.default_rng(s).permutation(m), 3)
            ioof = {b: np.full(m, np.nan) for b in BLOCKS}
            for b, X in BLOCKS.items():
                Xt = X[tri]
                for ite in inner:
                    itr = np.setdiff1d(np.arange(m), ite)
                    ioof[b][ite] = pct(fitapply(Xt[itr], co.t[tri][itr], co.e[tri][itr],
                                                Xt[ite], s))
            # stacking: three weights from the inner out-of-fold percentiles, ridge 1.0
            Zi = np.column_stack([ioof[b] for b in BLOCKS])
            w = cox_fit(Zi - Zi.mean(0), co.t[tri], co.e[tri], 1.0)
            weights.setdefault(s, []).append([round(float(x), 4) for x in w])
            Zv = np.column_stack([Z[s][b][vai] for b in BLOCKS])
            st[vai] = Zv @ w
            # the gate: each arm against a permutation null of its own inner out-of-fold score
            ti, ei = co.t[tri], co.e[tri]
            i2, j2 = cpairs(ti, ei)
            keep_b = []
            rng = np.random.default_rng(10000 + 100 * s + k)
            for b in BLOCKS:
                c_obs = cidx(ioof[b], i2, j2)
                null = []
                for _ in range(a.perms):
                    pr = rng.permutation(m)
                    i3, j3 = cpairs(ti[pr], ei[pr])
                    null.append(cidx(ioof[b], i3, j3))
                q95 = float(np.quantile(null, 0.95))
                passed = bool(c_obs > q95)
                gate_log.append({"seed": s, "fold": k, "arm": b, "inner_oof_c": round(c_obs, 4),
                                 "null_q95": round(q95, 4), "passes": passed})
                if passed:
                    keep_b.append(b)
            if not keep_b:
                keep_b = list(BLOCKS)
            # summed in the frozen method's own order (slide + omics + clinical), so that when every
            # arm passes the variant is bit-identical to it rather than equal up to rounding
            acc = np.zeros(len(vai))
            for b in BLOCKS:
                if b in keep_b:
                    acc = acc + Z[s][b][vai]
            gt[vai] = acc
        stacked_s[s] = co.fold_pct(st)
        gated_s[s] = co.fold_pct(gt)
        print("  seed %d stacking and gate done (%.0f s)" % (s, time.time() - t0),
              file=sys.stderr, flush=True)

    def summary(per_seed):
        v = [C(per_seed[s]) for s in SEEDS]
        return {"per_seed": [round(x, 4) for x in v], "mean": round(float(np.mean(v)), 4),
                "sd": round(float(np.std(v, ddof=1)), 4)}

    CONCAT, STACK, GATED = (canon([d[s] for s in SEEDS]) for d in (concat_s, stacked_s, gated_s))
    all_pass = all(g["passes"] for g in gate_log)
    fusion = {
        "construction": "identical cases, released folds, seeds 0-4, estimator and standardisation; "
                        "per-seed concordance pooled over out-of-fold within-fold percentiles, and "
                        "a canonical seed-averaged score vector for paired comparisons",
        "ours_equal_weight_rank_average": {**summary(ours_s), "canonical": round(C(OURS), 4)},
        "concatenated_ridge_cox": {**summary(concat_s), "canonical": round(C(CONCAT), 4),
                                   "columns": int(XCAT.shape[1]), "alpha_grid": list(EXT_ALPHAS)},
        "stacked_learned_weights": {**summary(stacked_s), "canonical": round(C(STACK), 4),
                                    "weights_order": list(BLOCKS),
                                    "weights_per_seed_fold": {str(s): weights[s] for s in SEEDS}},
        "clinical_alone": {"canonical": round(C(zc), 4)},
        "paired": {"ours_vs_concatenated": boot(OURS, CONCAT, co.t, co.e, a.reps),
                   "ours_vs_stacked": boot(OURS, STACK, co.t, co.e, a.reps)},
    }
    gated = {"every_arm_passes_every_fold_and_seed": all_pass,
             "fraction_of_arm_folds_passing": round(float(np.mean([g["passes"] for g in gate_log])), 4),
             **summary(gated_s), "canonical": round(C(GATED), 4),
             "identical_to_the_frozen_method": bool(all(np.array_equal(gated_s[s], ours_s[s])
                                                        for s in SEEDS)),
             "log": gate_log}

    # ================================================================ 3. added value and its inflation
    rng_seed = 20260911
    constructions = {
        "ours": {"alone": OURS_NOCLIN, "weak": OURSG, "stage": OURS},
        "survpath": {"alone": SP, "weak": SPG, "stage": SPC},
        "pibd_best_val": {"alone": PIBD["best_val"]["alone"], "weak": PIBD["best_val"]["grade"],
                          "stage": PIBD["best_val"]["stage"]},
        "pibd_final": {"alone": PIBD["final"]["alone"], "weak": PIBD["final"]["grade"],
                       "stage": PIBD["final"]["stage"]},
    }
    rng = np.random.default_rng(rng_seed)
    draws = [rng.choice(n, size=n, replace=True) for _ in range(a.reps)]
    added = {}
    for name, v in constructions.items():
        pt = {"C_weak": C(zg), "C_stage": C(zc), "C_M_plus_weak": C(v["weak"]),
              "C_M_plus_stage": C(v["stage"]), "C_M_alone": C(v["alone"])}
        pt["added_over_weak"] = pt["C_M_plus_weak"] - pt["C_weak"]
        pt["added_over_stage"] = pt["C_M_plus_stage"] - pt["C_stage"]
        pt["D"] = pt["added_over_weak"] - pt["added_over_stage"]
        pt["alone_minus_weak"] = pt["C_M_alone"] - pt["C_weak"]
        pt["alone_minus_stage"] = pt["C_M_alone"] - pt["C_stage"]
        bd, bw, bs, ba = [], [], [], []
        for bs_ in draws:
            i2, j2 = cpairs(co.t[bs_], co.e[bs_])
            cw, cs = cidx(zg[bs_], i2, j2), cidx(zc[bs_], i2, j2)
            mw, ms = cidx(v["weak"][bs_], i2, j2), cidx(v["stage"][bs_], i2, j2)
            bw.append(mw - cw); bs.append(ms - cs); bd.append((mw - cw) - (ms - cs))
            ba.append(cidx(v["alone"][bs_], i2, j2) - cs)

        def ci(x):
            x = np.asarray(x)
            return {"ci95": [round(float(np.quantile(x, .025)), 4), round(float(np.quantile(x, .975)), 4)],
                    "p_two_sided": round(float(2 * min((x <= 0).mean(), (x >= 0).mean())), 4)}
        added[name] = {**{k: round(float(val), 4) for k, val in pt.items()},
                       "D_boot": ci(bd), "added_over_weak_boot": ci(bw),
                       "added_over_stage_boot": ci(bs), "alone_minus_stage_boot": ci(ba)}
        print("  added value %-14s over weak %+.4f, over stage %+.4f, D %+.4f"
              % (name, pt["added_over_weak"], pt["added_over_stage"], pt["D"]), file=sys.stderr)

    # ================================================================ 4. the modality paragraph, amended
    z0 = Z[0]
    sp_seed1 = co.fold_pct(sp1, f1)
    gap = np.abs(z0["clinical"][ii] - z0["clinical"][jj])
    thr = float(np.quantile(gap, 0.05))
    tie = gap <= thr
    it, jt = ii[tie], jj[tie]
    ours0 = ours_s[0]
    modality = {
        "clinical_block": "amended: age, sex, AJCC pathologic stage from DIMAF's split files, seed 0",
        "spearman_slide_vs_clinical": round(spearman(z0["slide"], z0["clinical"]), 4),
        "spearman_survpath_vs_clinical": round(spearman(sp_seed1, z0["clinical"]), 4),
        "spearman_omics_vs_clinical": round(spearman(z0["omics"], z0["clinical"]), 4),
        "spearman_slide_vs_omics": round(spearman(z0["slide"], z0["omics"]), 4),
        "tightest_5pct_by_clinical_gap": {
            "pairs": int(tie.sum()), "of": int(len(ii)),
            "clinical": round(cidx(z0["clinical"], it, jt), 4),
            "slide": round(cidx(z0["slide"], it, jt), 4),
            "omics": round(cidx(z0["omics"], it, jt), 4),
            "ours": round(cidx(ours0, it, jt), 4)},
        "single_arms_seed0": {k: round(C(z0[k]), 4) for k in BLOCKS},
        "the_values_this_replaces": {"spearman_slide_vs_clinical": 0.213,
                                     "spearman_survpath_vs_clinical": 0.166,
                                     "tightest_5pct_clinical": 0.5012, "tightest_5pct_slide": 0.6120,
                                     "pairs": 1281, "source": "development atlas, pre-amendment "
                                     "23-column clinical block"}}

    # ================================================================ write
    os.makedirs(a.out_dir, exist_ok=True)
    rep = {"artifact_type": "s19_unified_fusion_and_added_value",
           "phase_of_origin": "post_freeze_2026-09-11",
           "reportable": True,
           "primary_unchanged": "the frozen primary 0.7212 is reproduced below as a known answer and is "
                                "not replaced by anything in this file",
           "executed_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "n": n, "events": int(co.e.sum()), "comparable_pairs": int(len(ii)),
           "known_answers": {k: {"recomputed": round(got[k], 4), "reported": KNOWN[k]} for k in KNOWN},
           "fusion_on_identical_inputs": fusion,
           "evidence_gated_variant": gated,
           "added_value_inflation": {"weak_reference": "age + sex + grade (DIMAF files)",
                                     "stage_reference": "age + sex + stage (DIMAF files)",
                                     "bootstrap": {"replicates": a.reps, "unit": "case",
                                                   "rng_seed": rng_seed},
                                     "constructions": added},
           "modality_paragraph_amended": modality,
           "runtime_seconds": round(time.time() - t0, 1)}
    json.dump(rep, open(os.path.join(a.out_dir, "unified-fusion-and-added-value.json"), "w"),
              indent=1)
    percase = {"note": "canonical seed-averaged, within-fold percentile scores; seed-0 vectors "
                       "under seed0_*. Every paired comparison in the paper is computable from these.",
               "cases": [{"case_id": c, "fold": int(co.fold[i]), "months": float(co.t[i]),
                          "event": int(co.e[i]),
                          **{k: round(float(v[i]), 5) for k, v in (
                              ("slide", zt), ("omics", zo), ("clinical_stage", zc),
                              ("clinical_grade", zg), ("ours", OURS), ("ours_grade", OURSG),
                              ("ours_without_clinical", OURS_NOCLIN), ("concatenated", CONCAT),
                              ("stacked", STACK), ("gated", GATED), ("survpath", SP),
                              ("survpath_plus_stage", SPC), ("survpath_plus_grade", SPG),
                              ("pibd_best_val", PIBD["best_val"]["alone"]),
                              ("pibd_best_val_plus_stage", PIBD["best_val"]["stage"]),
                              ("pibd_best_val_plus_grade", PIBD["best_val"]["grade"]),
                              ("pibd_final", PIBD["final"]["alone"]),
                              ("pibd_final_plus_stage", PIBD["final"]["stage"]),
                              ("pibd_final_plus_grade", PIBD["final"]["grade"]),
                              ("seed0_slide", z0["slide"]), ("seed0_omics", z0["omics"]),
                              ("seed0_clinical_stage", z0["clinical"]),
                              ("seed0_clinical_grade", ZG[0]), ("seed0_ours", ours0))}}
                         for i, c in enumerate(co.keep)]}
    json.dump(percase, open(os.path.join(a.out_dir, "percase-canonical-vectors.json"), "w"))
    print("wrote %s (%.0f s)" % (a.out_dir, time.time() - t0), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
