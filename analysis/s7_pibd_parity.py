#!/usr/bin/env python3
"""PIBD at input parity -- the second paired comparator this benchmark permits.

WHY IT MATTERS. Until now exactly one competitor could be compared to our arm PAIRED, on the same
patients, with the same clinical variables added: SurvPath, because it is the only entrant whose
implementation we could rerun end to end. Every other value in the field's table is a published
number that cannot be given the clinical block and cannot be paired. PIBD's rerun makes it two --
and PIBD is the stronger of the pair by its published value (0.667 against 0.625) and the one whose
split files were verified byte-identical to the released ones on all five studies, 25 of 25 files.

TWO CHECKPOINT CONVENTIONS, and the choice is declared rather than made silently. PIBD writes
`split_k_results.pkl` at the best validation epoch and `split_k_results_final.pkl` at the last one.
Its own logged metric is "Best Val c-index", so the best-epoch file is what corresponds to the
number it publishes. Both are scored here; the best-epoch convention is used for the comparison,
and the final-epoch value is reported beside it so the choice is visible.

The reproduction itself is the first thing to check: five folds at 0.6038, 0.6983, 0.6250, 0.6147
and 0.7626 give a mean of 0.6609 against a published 0.667 +/- 0.061. If the pooled score computed
from the recovered per-case risks does not agree with those per-fold numbers, the recovered
predictions are not usable and nothing below is either.
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
from s7_survival_metrics import km_censoring, g_at, td_auc, ibs, logrank, chi2_sf  # noqa: E402
from s7_gaps import holm  # noqa: E402


def load_pibd(d, keep, idx, final):
    pat = "split_*_results_final.pkl" if final else "split_*_results.pkl"
    files = sorted(f for f in glob.glob(os.path.join(d, pat))
                   if final or not f.endswith("_final.pkl"))
    risk = np.full(len(keep), np.nan)
    fold = np.full(len(keep), -1)
    t = np.full(len(keep), np.nan)
    e = np.full(len(keep), np.nan)
    for k, f in enumerate(files):
        for c, v in pickle.load(open(f, "rb")).items():
            if c in idx:
                i = idx[c]
                risk[i] = float(v["risk"]); fold[i] = k
                t[i] = float(v["time"]); e[i] = 1.0 if float(v["censorship"]) == 0.0 else 0.0
    return risk, fold, t, e, len(files)


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

    zt, zo, zc = arm(co.T), arm(P), arm(CLIN)
    OURS = co.fold_pct(zt + zo + zc)

    rep = {"artifact_type": "s7_pibd_input_parity", "phase_of_origin": "frozen_post_hoc",
           "n": len(co.keep), "events": int(co.e.sum())}

    # ---- recover PIBD under both checkpoint conventions and CHECK the reproduction first
    conv = {}
    for final in (False, True):
        r, f, t, e, nf = load_pibd(a.pibd_dir, co.keep, co.idx, final)
        cov = int(np.isfinite(r).sum())
        name = "final_epoch" if final else "best_epoch"
        if cov < len(co.keep):
            conv[name] = {"status": "incomplete", "cases_recovered": cov, "files": nf}
            continue
        assert np.allclose(t, co.t), "%s: labels differ from the cohort's" % name
        assert (f == co.fold).all(), "%s: fold assignment differs from the released folds" % name
        z = co.fold_pct(r, f)
        per_fold = []
        for k in range(5):
            m = co.fold == k
            i2, j2 = cpairs(co.t[m], co.e[m])
            per_fold.append(round(cidx(z[m], i2, j2), 4))
        conv[name] = {"status": "ok", "files": nf, "cases_recovered": cov,
                      "pooled_cindex": round(cidx(z, ii, jj), 4),
                      "per_fold_from_recovered_risks": per_fold,
                      "per_fold_mean": round(float(np.mean(per_fold)), 4),
                      "per_fold_sd": round(float(np.std(per_fold, ddof=1)), 4),
                      "_z": z}
        print("%-12s pooled=%.4f  per-fold %s  mean=%.4f"
              % (name, conv[name]["pooled_cindex"], per_fold, conv[name]["per_fold_mean"]),
              file=sys.stderr, flush=True)

    rep["reproduction"] = {
        "published": {"value": 0.667, "sd_over_folds": 0.061},
        "logged_best_val_per_fold": [0.6038, 0.6983, 0.6250, 0.6147, 0.7626],
        "logged_mean": 0.6609,
        "conventions": {k: {kk: vv for kk, vv in v.items() if kk != "_z"} for k, v in conv.items()},
        "convention_used": "best_epoch, because PIBD's own logged metric is 'Best Val c-index' and "
                           "that is what corresponds to its published value"}

    use = conv.get("best_epoch")
    if not use or use.get("status") != "ok":
        rep["blocked"] = "the best-epoch predictions could not be recovered for all cases"
        open(a.out, "w").write(json.dumps(rep, indent=1) + "\n")
        print("BLOCKED", file=sys.stderr)
        return 2
    z_pibd = use["_z"]
    PIBD_C = co.fold_pct(z_pibd + zc)

    # ---- SurvPath, seed-averaged, for the two-competitor table
    from s7_survpath_multiseed import load_seed
    zs = {}
    r, f, _, _ = load_seed(a.sp_dir, co.keep, co.idx); zs[1] = co.fold_pct(r, f)
    for s in (2, 3, 4, 5):
        d = os.path.join(a.seed_dir, "seed%d" % s)
        if os.path.isdir(d):
            r, f, _, _ = load_seed(d, co.keep, co.idx); zs[s] = co.fold_pct(r, f)
    SP_avg = co.fold_pct(np.vstack([zs[s] for s in sorted(zs)]).mean(0))
    SPC_avg = co.fold_pct(SP_avg + zc)

    uniq, g = km_censoring(co.t, co.e)
    ev_t = np.sort(co.t[co.e == 1])
    hs = [float(np.quantile(ev_t, q)) for q in (0.25, 0.5, 0.75)]
    grid = [float(x) for x in np.linspace(ev_t[0], np.quantile(ev_t, 0.9), 12)]
    kt, ks = km_censoring(co.t, 1.0 - co.e)
    base = g_at(kt, ks, np.array(grid))
    gap = np.abs(zc[ii] - zc[jj])
    tied = gap <= float(np.quantile(gap, 0.10))

    def row(z):
        curves = np.clip(base[None, :] ** np.exp(2.0 * (z[:, None] - 0.5)), 1e-6, 1 - 1e-6)
        q = np.quantile(z, [1 / 3, 2 / 3]); grp = np.digitize(z, q)
        c2, df, _ = logrank(co.t, co.e, grp)
        return {"cindex": round(cidx(z, ii, jj), 4),
                "tied_q10": round(cidx(z, ii[tied], jj[tied]), 4),
                **{"td_auc_%dm" % round(h): round(td_auc(z, co.t, co.e, h, uniq, g), 4) for h in hs},
                "ibs": round(ibs(curves, grid, co.t, co.e, uniq, g), 4),
                "logrank_p": float("%.3g" % chi2_sf(c2, df))}

    rep["arms"] = {
        "PIBD rerun": row(z_pibd),
        "PIBD + the same clinical": row(PIBD_C),
        "SurvPath rerun, seed-averaged": row(SP_avg),
        "SurvPath + the same clinical, seed-averaged": row(SPC_avg),
        "OURS": row(OURS)}

    fam = {
        "ours vs PIBD + same clinical": boot(OURS, PIBD_C, co.t, co.e, 6000),
        "ours vs SurvPath + same clinical (seed-averaged)": boot(OURS, SPC_avg, co.t, co.e, 6000),
        "ours vs PIBD alone": boot(OURS, z_pibd, co.t, co.e, 6000),
    }
    raw = {k: v["p_two_sided"] for k, v in fam.items()}
    adj = holm(raw)
    rep["input_parity_comparisons"] = {
        "family": "the paired comparisons the benchmark permits -- both competitors whose "
                  "implementations could be rerun end to end, each given the same clinical block",
        "multiplicity": "Holm, as the frozen protocol declares",
        "comparisons": {k: {"delta": fam[k]["mean"], "se": fam[k]["se"], "ci95": fam[k]["ci95"],
                            "p_raw": raw[k], "p_holm": adj[k],
                            "survives_holm_at_0.05": bool(adj[k] < 0.05)} for k in fam}}

    pf = {}
    for k in range(5):
        m = co.fold == k
        i2, j2 = cpairs(co.t[m], co.e[m])
        pf["fold_%d" % k] = {"n": int(m.sum()), "events": int(co.e[m].sum()),
                             "OURS": round(cidx(OURS[m], i2, j2), 4),
                             "PIBD+clinical": round(cidx(PIBD_C[m], i2, j2), 4),
                             "SurvPath+clinical": round(cidx(SPC_avg[m], i2, j2), 4)}
    rep["per_fold"] = pf
    rep["folds_won"] = {
        "vs PIBD+clinical": sum(1 for v in pf.values() if v["OURS"] > v["PIBD+clinical"]),
        "vs SurvPath+clinical": sum(1 for v in pf.values() if v["OURS"] > v["SurvPath+clinical"])}

    open(a.out, "w").write(json.dumps(rep, indent=1) + "\n")
    print(json.dumps({"arms": rep["arms"], "comparisons": rep["input_parity_comparisons"]["comparisons"],
                      "folds_won": rep["folds_won"]}, indent=1))
    print("wrote", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
