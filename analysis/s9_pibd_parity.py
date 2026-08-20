#!/usr/bin/env python3
"""PIBD scored by our code, and PIBD given the same clinical block: the second parity comparator.

WHY THIS MATTERS MORE THAN ANOTHER COHORT. Nine of the ten published entries on this benchmark
cannot be placed at input parity, because none of them releases per-case predictions. Until now
SurvPath was the only one that could, so the paper rested its input-parity claim on a single
competitor. PIBD's implementation runs end to end, its split files were verified byte-identical to
the released ones on all five studies, and it is the strongest entrant we can execute. Adding it
takes the paired comparisons from one competitor to two.

TWO CHECKPOINT CONVENTIONS, AND THE DIFFERENCE IS ITSELF A FINDING. PIBD writes both
`split_k_results.pkl` at its best VALIDATION epoch and `split_k_results_final.pkl` at the last one.
Its published number is a "Best Val c-index", so the headline it reports is selected on the very
fold it is scored on -- an optimism this benchmark's tables do not flag. Both are computed here and
both are reported; the gap between them is how much that selection is worth.

WHAT IS CHECKED BEFORE ANYTHING IS SCORED. The labels and the fold assignment must be identical to
ours case by case, exactly as falsifier L5 required for SurvPath. Two arms scored on different
endpoints are not a paired comparison, whatever the arithmetic says.
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
from blca_common import Cohort, cidx, cpairs, fitapply  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402
from s6_amend_clinical import dimaf_clinical, boot  # noqa: E402
from s7_gaps import holm  # noqa: E402
from s7_survival_metrics import (chi2_sf, g_at, ibs, km_censoring, logrank,  # noqa: E402
                                 td_auc)


def load_pibd(dirpath, keep, idx, suffix):
    risk = np.full(len(keep), np.nan)
    fold = np.full(len(keep), -1)
    t = np.full(len(keep), np.nan)
    e = np.full(len(keep), np.nan)
    n = 0
    for k in range(5):
        f = os.path.join(dirpath, "split_%d_results%s.pkl" % (k, suffix))
        if not os.path.exists(f):
            continue
        for c, v in pickle.load(open(f, "rb")).items():
            if c in idx:
                i = idx[c]
                risk[i] = float(v["risk"]); fold[i] = k
                t[i] = float(v["time"])
                e[i] = 1.0 if float(v["censorship"]) == 0.0 else 0.0
                n += 1
    return risk, fold, t, e, n


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("root", "titan", "sp-dir", "seed-dir", "pibd-dir", "omics-dir", "dimaf-dir", "out"):
        ap.add_argument("--" + f, required=True)
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

    # EVERY arm in the comparison family is averaged over the SAME five seeds that define the
    # paper's primary metric. The first version of this file used seed 0 for our arms and the
    # five-seed mean for the competitor, which put +0.0334 in the Holm table against +0.0318 in the
    # text -- the same asymmetry, in miniature, that the seed-matching section of the paper is
    # about. A paired bootstrap needs one score vector per arm, so the averaging is over SCORES;
    # that is also how the seed-matched headline was constructed, and it must stay the same one.
    SEEDS = (0, 1, 2, 3, 4)
    zt = co.fold_pct(np.vstack([arm(co.T, s) for s in SEEDS]).mean(0))
    zo = co.fold_pct(np.vstack([arm(P, s) for s in SEEDS]).mean(0))
    zc = co.fold_pct(np.vstack([arm(CLIN, s) for s in SEEDS]).mean(0))
    OURS = co.fold_pct(np.vstack([co.fold_pct(arm(co.T, s) + arm(P, s) + arm(CLIN, s))
                                  for s in SEEDS]).mean(0))

    rep = {"artifact_type": "s9_pibd_input_parity", "arms_averaged_over_seeds": [0, 1, 2, 3, 4], "phase_of_origin": "frozen_post_hoc",
           "n": len(co.keep), "events": int(co.e.sum())}

    # ---------------------------------------------------------------- integrity, before scoring
    variants = {}
    for suffix, label in (("", "best_validation_epoch"), ("_final", "final_epoch")):
        r, f, t, e, n = load_pibd(a.pibd_dir, co.keep, co.idx, suffix)
        cov = int(np.isfinite(r).sum())
        lab_ok = bool(np.allclose(np.nan_to_num(t), np.nan_to_num(co.t), atol=1e-6))
        ev_ok = bool(np.array_equal(np.nan_to_num(e), np.nan_to_num(co.e)))
        fold_ok = bool((f == co.fold).all())
        variants[label] = {"cases_covered": cov, "labels_identical": lab_ok,
                           "events_identical": ev_ok, "fold_assignment_identical": fold_ok,
                           "risk": r}
        assert cov == len(co.keep), "%s covers %d of %d cases" % (label, cov, len(co.keep))
        assert lab_ok and ev_ok and fold_ok, (
            "%s: labels/events/folds differ from the cohort -- this is not a paired comparison"
            % label)
    rep["integrity"] = {k: {kk: vv for kk, vv in v.items() if kk != "risk"}
                        for k, v in variants.items()}

    # ---------------------------------------------------------------- scoring
    uniq, g = km_censoring(co.t, co.e)
    ev_t = np.sort(co.t[co.e == 1])
    hs = [float(np.quantile(ev_t, q)) for q in (0.25, 0.5, 0.75)]
    grid = [float(x) for x in np.linspace(ev_t[0], np.quantile(ev_t, 0.9), 12)]
    kt, ks = km_censoring(co.t, 1.0 - co.e)
    base = g_at(kt, ks, np.array(grid))

    def metrics(z):
        curves = np.clip(base[None, :] ** np.exp(2.0 * (z[:, None] - 0.5)), 1e-6, 1 - 1e-6)
        q = np.quantile(z, [1 / 3, 2 / 3])
        c2, df, _ = logrank(co.t, co.e, np.digitize(z, q))
        return {"cindex": round(cidx(z, ii, jj), 4),
                **{"td_auc_%dm" % round(h): round(td_auc(z, co.t, co.e, h, uniq, g), 4) for h in hs},
                "ibs": round(ibs(curves, grid, co.t, co.e, uniq, g), 4),
                "logrank_p": float("%.3g" % chi2_sf(c2, df))}

    out = {}
    for label, v in variants.items():
        z = co.fold_pct(v["risk"], co.fold)
        zpc = co.fold_pct(z + zc)
        per_fold = []
        for k in range(5):
            m = co.fold == k
            i2, j2 = cpairs(co.t[m], co.e[m])
            per_fold.append(round(cidx(z[m], i2, j2), 4))
        out[label] = {
            "pibd_alone": metrics(z),
            "pibd_plus_same_clinical": metrics(zpc),
            "per_fold_cindex": per_fold,
            "fold_mean": round(float(np.mean(per_fold)), 4),
            "fold_sd": round(float(np.std(per_fold, ddof=1)), 4),
            "_z": z, "_zpc": zpc}
        print("%-24s pooled=%.4f  fold-mean=%.4f+/-%.4f  +clinical=%.4f"
              % (label, out[label]["pibd_alone"]["cindex"], out[label]["fold_mean"],
                 out[label]["fold_sd"], out[label]["pibd_plus_same_clinical"]["cindex"]),
              file=sys.stderr, flush=True)

    rep["pibd"] = {k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")}
                   for k, v in out.items()}
    rep["published"] = {"value": 0.667, "sd_over_folds": 0.061,
                        "reported_as": "Best Val c-index, i.e. the checkpoint is chosen on the "
                                       "fold it is scored on",
                        "our_best_val_epoch_fold_mean": out["best_validation_epoch"]["fold_mean"],
                        "our_final_epoch_fold_mean": out["final_epoch"]["fold_mean"],
                        "checkpoint_selection_is_worth": round(
                            out["best_validation_epoch"]["fold_mean"]
                            - out["final_epoch"]["fold_mean"], 4)}

    # ---------------------------------------------------------------- the parity comparison
    z_best = out["best_validation_epoch"]["_zpc"]
    z_final = out["final_epoch"]["_zpc"]
    rep["comparison_ours_vs_pibd_plus_clinical"] = {
        "ours": round(cidx(OURS, ii, jj), 4),
        "pibd_plus_clinical_best_val_checkpoint": round(cidx(z_best, ii, jj), 4),
        "pibd_plus_clinical_final_epoch": round(cidx(z_final, ii, jj), 4),
        "paired_vs_best_val": boot(OURS, z_best, co.t, co.e, 6000),
        "paired_vs_final": boot(OURS, z_final, co.t, co.e, 6000),
        "note": "the best-validation-checkpoint variant is the harder comparison for us, because "
                "that checkpoint was chosen on the fold it is scored on. It is the one to quote."}

    # ---------------------------------------------------------------- Holm, family now six
    import pickle as _p
    zs = {}
    for k, fpath in enumerate(sorted(glob.glob(os.path.join(a.sp_dir, "split_*_results.pkl")))):
        for c, v in _p.load(open(fpath, "rb")).items():
            if c in co.idx:
                zs.setdefault("r", np.full(len(co.keep), np.nan))
                zs.setdefault("f", np.full(len(co.keep), -1))
                zs["r"][co.idx[c]] = float(v["risk"]); zs["f"][co.idx[c]] = k
    sp_seeds = [co.fold_pct(zs["r"], zs["f"])]
    for s in (2, 3, 4, 5):
        d = os.path.join(a.seed_dir, "seed%d" % s)
        if os.path.isdir(d):
            rr = np.full(len(co.keep), np.nan); ff = np.full(len(co.keep), -1)
            for k, fpath in enumerate(sorted(glob.glob(os.path.join(d, "split_*_results.pkl")))):
                for c, v in _p.load(open(fpath, "rb")).items():
                    if c in co.idx:
                        rr[co.idx[c]] = float(v["risk"]); ff[co.idx[c]] = k
            sp_seeds.append(co.fold_pct(rr, ff))
    SPC = co.fold_pct(co.fold_pct(np.vstack(sp_seeds).mean(0)) + zc)

    fam = {
        "ours vs omics alone": boot(OURS, zo, co.t, co.e, 6000),
        "ours vs slide alone": boot(OURS, zt, co.t, co.e, 6000),
        "ours vs clinical alone": boot(OURS, zc, co.t, co.e, 6000),
        "ours vs ours-without-clinical": boot(OURS, co.fold_pct(zt + zo), co.t, co.e, 6000),
        "ours vs SurvPath+clinical (seed-matched)": boot(OURS, SPC, co.t, co.e, 6000),
        "ours vs PIBD+clinical (best-val checkpoint)": rep[
            "comparison_ours_vs_pibd_plus_clinical"]["paired_vs_best_val"],
    }
    # THE CANONICAL ARM TABLE. Every comparison above is drawn from these and only these, because
    # having two constructions of "our five-seed arm" in one manuscript is how +0.0318 and +0.0334
    # came to sit in the same document. The paper's PRIMARY metric is a different object -- the mean
    # of five per-seed concordances, 0.7212 -- and it stays what the decision rule tested. A paired
    # bootstrap cannot be computed from a mean of concordances; it needs one score per case. So the
    # two coexist by design, and the manuscript says which is which.
    rep["canonical_arms_seed_averaged_scores"] = {
        "construction": "per seed 0-4: fit, score out-of-fold, percentile-normalise within fold. "
                        "Then average the five score vectors and re-normalise within fold. Both "
                        "sides of every comparison get this identical treatment.",
        "ours": round(cidx(OURS, ii, jj), 4),
        "slide_alone": round(cidx(zt, ii, jj), 4),
        "omics_alone": round(cidx(zo, ii, jj), 4),
        "clinical_alone": round(cidx(zc, ii, jj), 4),
        "ours_without_clinical": round(cidx(co.fold_pct(zt + zo), ii, jj), 4),
        "survpath_plus_clinical": round(cidx(SPC, ii, jj), 4),
        "pibd_plus_clinical_best_val": round(cidx(z_best, ii, jj), 4),
        "pibd_plus_clinical_final_epoch": round(cidx(z_final, ii, jj), 4),
        "not_the_primary_metric": "the pre-registered primary is 0.7212, the MEAN OF FIVE PER-SEED "
                                  "CONCORDANCES. This table's 'ours' is the concordance of the "
                                  "AVERAGED SCORE, which is a different and slightly higher "
                                  "quantity. Averaging favours the competitor if anything: its "
                                  "seeds are full retrainings, ours are inner-CV partitions."}
    raw = {k: v["p_two_sided"] for k, v in fam.items()}
    adj = holm(raw)
    rep["holm_family_of_six"] = {
        "note": "the family grew from five to six when PIBD became scoreable; Holm is recomputed "
                "over the whole family rather than the five it was first applied to, because "
                "adding a comparison and keeping the old adjustment would be the cheapest possible "
                "way to keep a result",
        "comparisons": {k: {"delta": fam[k]["mean"], "ci95": fam[k]["ci95"],
                            "p_raw": raw[k], "p_holm": adj[k],
                            "survives_at_0.05": bool(adj[k] < 0.05)} for k in fam}}

    open(a.out, "w").write(json.dumps(rep, indent=1) + "\n")
    print(json.dumps({"published_check": rep["published"],
                      "parity": {k: v for k, v in
                                 rep["comparison_ours_vs_pibd_plus_clinical"].items()
                                 if k != "note"},
                      "holm": rep["holm_family_of_six"]["comparisons"]}, indent=1))
    print("wrote", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
