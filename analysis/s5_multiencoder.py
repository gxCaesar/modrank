#!/usr/bin/env python3
"""Six independently pretrained slide encoders, one omics view, one clinical view.

Pre-freeze, exploratory. NOT reportable.

WHY THIS IS THE RUN THE EVIDENCE ASKED FOR. Everything measured today says the same thing about
where the remaining room is. The combination rule is at its ceiling: the best fixed weight on the
image/clinical pair is exactly 0.50, an oracle two-parameter fusion buys +0.008 over it, and gating
came back void against its own permuted control. The estimator is not the deficit either -- the
in-sample "headroom" that suggested it was turned out to be overfitting once the null was fixed.
Seven candidate families died on their falsifiers. What is left is the arms themselves, and the
image arm has been one encoder all day: TITAN at 0.6596, with CHIEF's mean pooling at 0.5836 the
only alternative and a weak one.

Six public, ungated, already-computed pan-TCGA slide representations are now on disk, every one
covering all 359 bladder cases:

    TITAN                768-d   CONCH-family slide encoder (already in use)
    UNI                 1024-d   mean-pooled patch features -- THE ENCODER MMP AND DIMAF USE
    Prov-GigaPath      10752-d   native slide encoder, all 14 layers concatenated
    GigaSSL/Phikon       512-d   Phikon tiles, GigaSSL aggregation
    GigaSSL/H-Optimus-0  512-d   H-Optimus-0 tiles, same aggregator
    GigaSSL/GigaPath     512-d   Prov-GigaPath tiles, same aggregator
    GigaSSL/CTransPath   512-d   CTransPath tiles, same aggregator

The four GigaSSL variants share an aggregator, so they are NOT four independent encoders -- they
are four tile encoders under one slide-level construction, and their errors will be more
correlated with each other than with TITAN, UNI or native GigaPath. That is stated here rather
than discovered later, and the measured correlation matrix is reported.

TWO COMBINATION RULES, BOTH FIXED BEFORE THIS RUN, BOTH REPORTED.

  primary   skill-weighted: w_v = max(0, C_inner(v) - 0.5), with C_inner estimated by inner CV
            inside the training fold only. Equal weighting is its special case. It is primary
            because equal weighting has now been damaged three separate times today by admitting a
            weak view at full weight -- CHIEF on bladder, hallmarks on bladder, and omics on BRCA.
  secondary equal weight over everything.

Neither rule looks at a held-out outcome, and no view is dropped by inspecting its score.

THE COMPARISON THAT DECIDES THE PAPER, unchanged from this morning: our arm on **WSI + omics only**
against DIMAF's 0.679 on the same released folds. It has failed under three combination rules with
one image encoder. Six encoders is the last untried input to it.
"""

from __future__ import annotations

import argparse
import csv as _csv
import glob
import itertools
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cindex, cpairs, fitapply, pct  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402


def spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = float(np.sqrt((ra @ ra) * (rb @ rb)))
    return float(ra @ rb / d) if d else float("nan")


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
    M = np.vstack([X[bag[c]].mean(0) if c in bag else mu for c in keep])
    return M, len(have)


def inner_skill(X, t, e, tri, seed=0):
    inner = np.array_split(np.random.default_rng(seed + 21).permutation(len(tri)), 3)
    ip = np.zeros(len(tri))
    for ite in inner:
        itr = np.setdiff1d(np.arange(len(tri)), ite)
        ip[ite] = fitapply(X[tri][itr], t[tri][itr], e[tri][itr], X[tri][ite], seed)
    c = cindex(ip, t[tri], e[tri])
    return float(c) if c == c else 0.5


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--omics-dir", required=True)
    ap.add_argument("--fm-dir", required=True)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    stem2case = {r["slide_id"].replace(".svs", ""): r["case_id"]
                 for r in _csv.DictReader(open(a.root + "/tcga_blca_meta.csv"))}

    blocks = {"wsi_titan": co.T}
    coverage = {"wsi_titan": len(co.keep)}
    for p in sorted(glob.glob(os.path.join(a.fm_dir, "*.npz"))):
        tag = "wsi_" + os.path.basename(p)[:-4]
        M, n = percase(p, co.keep, stem2case)
        coverage[tag] = n
        if M is not None:
            blocks[tag] = M
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    blocks["omics_combine"] = P
    blocks["clinical"] = co.CLIN

    ii, jj = cpairs(co.t, co.e)
    arms, skills = {}, {}
    for nm, X in blocks.items():
        s = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], 0)
        arms[nm] = co.fold_pct(s)
        print("%-28s %6d-d  C=%.4f" % (nm, X.shape[1], cidx(arms[nm], ii, jj)),
              file=sys.stderr, flush=True)
    arms["survpath_rerun"] = co.fold_pct(co.spr)

    rep = {"artifact_type": "s5_multiencoder", "reportable": False,
           "phase_of_origin": "exploratory", "n": len(co.keep), "events": int(co.e.sum()),
           "blocks": {k: int(v.shape[1]) for k, v in blocks.items()},
           "cases_covered_per_encoder": coverage,
           "single_arms": {k: round(cidx(v, ii, jj), 4) for k, v in arms.items()}}

    wsi = [k for k in blocks if k.startswith("wsi_")]
    rep["image_arm_rank_correlations"] = {
        "%s~%s" % (a1, b1): round(spearman(arms[a1], arms[b1]), 3)
        for a1, b1 in itertools.combinations(sorted(wsi), 2)}

    def combine(pool, mode):
        out = np.full(len(co.keep), np.nan)
        wlog = {}
        for k, (tri, vai) in enumerate(fi):
            if len(tri) < 30 or not len(vai):
                continue
            if mode == "equal":
                w = {nm: 1.0 for nm in pool}
            else:
                w = {nm: max(0.0, inner_skill(blocks[nm], co.t, co.e, tri) - 0.5) for nm in pool}
            tot = sum(w.values()) or float(len(pool))
            if sum(w.values()) <= 0:
                w = {nm: 1.0 for nm in pool}
            wlog["fold_%d" % k] = {nm: round(w[nm] / tot, 4) for nm in pool}
            acc = np.zeros(len(vai))
            for nm in pool:
                acc += w[nm] * pct(arms[nm][vai])
            out[vai] = acc / tot
        return co.fold_pct(out), wlog

    POOLS = {"wsi_only": wsi,
             "wsi_plus_omics": wsi + ["omics_combine"],
             "wsi_plus_omics_plus_clinical": wsi + ["omics_combine", "clinical"],
             "titan_only_reference": ["wsi_titan", "omics_combine", "clinical"]}
    res = {}
    for pn, pool in POOLS.items():
        res[pn] = {}
        for mode in ("skill", "equal"):
            z, wl = combine(pool, mode)
            res[pn][mode] = {"cindex": round(cidx(z, ii, jj), 4)}
            if mode == "skill":
                res[pn][mode]["weights_fold_0"] = wl.get("fold_0")
                res[pn]["_z"] = z
            print("%-32s %-6s C=%.4f" % (pn, mode, res[pn][mode]["cindex"]),
                  file=sys.stderr, flush=True)
    rep["pools"] = {k: {m: v[m] for m in ("skill", "equal")} for k, v in res.items()}

    spc = co.fold_pct(arms["survpath_rerun"] + arms["clinical"])
    ours = res["wsi_plus_omics_plus_clinical"]["_z"]
    noclin = res["wsi_plus_omics"]["_z"]

    def boot(x, y, reps=6000, seed=0):
        rng = np.random.default_rng(seed)
        n = len(co.keep)
        v = []
        for _ in range(reps):
            bs = rng.choice(n, size=n, replace=True)
            i2, j2 = cpairs(co.t[bs], co.e[bs])
            if i2.size:
                v.append(cidx(x[bs], i2, j2) - cidx(y[bs], i2, j2))
        v = np.asarray(v)
        return {"mean": round(float(v.mean()), 4), "se": round(float(v.std(ddof=1)), 4),
                "ci95": [round(float(np.quantile(v, .025)), 4),
                         round(float(np.quantile(v, .975)), 4)],
                "p_two_sided": round(float(2 * min((v <= 0).mean(), (v >= 0).mean())), 4)}

    rep["comparison_C_THE_DECIDER"] = {
        "ours_wsi_plus_omics_skill_weighted": round(cidx(noclin, ii, jj), 4),
        "ours_wsi_plus_omics_equal_weight": res["wsi_plus_omics"]["equal"]["cindex"],
        "incumbent_DIMAF_same_folds": 0.679,
        "gap_primary_rule": round(cidx(noclin, ii, jj) - 0.679, 4),
        "bar": 0.0145,
        "passes": bool(cidx(noclin, ii, jj) - 0.679 > 0.0145),
        "previously_failed_with_one_encoder": {"flat_7_views": 0.6770, "skill_5_views": 0.6767,
                                               "threshold_5_views": 0.6776, "titan_plus_omics": 0.6841}}
    rep["comparison_B_input_parity"] = {
        "ours_full": round(cidx(ours, ii, jj), 4),
        "survpath_plus_same_clinical": round(cidx(spc, ii, jj), 4),
        "paired_bootstrap": boot(ours, spc)}
    rep["comparison_A_vs_published"] = {
        "ours_full": round(cidx(ours, ii, jj), 4),
        "best_published_blca_any_protocol": {"MOAD-FNet": 0.691, "DIMAF": 0.679, "APL": 0.677},
        "gap_vs_best": round(cidx(ours, ii, jj) - 0.691, 4)}
    rep["vs_single_encoder_baseline"] = {
        "titan_only_arm": res["titan_only_reference"]["skill"]["cindex"],
        "six_encoder_arm": round(cidx(ours, ii, jj), 4),
        "paired_bootstrap": boot(ours, res["titan_only_reference"]["_z"])}

    pf = {}
    for k in range(5):
        m = co.fold == k
        i2, j2 = cpairs(co.t[m], co.e[m])
        pf["fold_%d" % k] = {"n": int(m.sum()), "events": int(co.e[m].sum()),
                             "ours": round(cidx(ours[m], i2, j2), 4),
                             "ours_no_clinical": round(cidx(noclin[m], i2, j2), 4)}
    rep["per_fold"] = pf
    v = [x["ours"] for x in pf.values()]
    rep["ours_fold_mean_sd"] = [round(float(np.mean(v)), 4), round(float(np.std(v, ddof=1)), 4)]

    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
