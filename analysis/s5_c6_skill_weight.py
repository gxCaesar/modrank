#!/usr/bin/env python3
"""C6 -- weight each view by the skill it demonstrates INSIDE the training fold.

Pre-freeze, exploratory. NOT reportable.

THE DEFICIT, and it is one this campaign created for itself an hour ago. The combination rule was
pre-specified as "use every available view, equal weight", on the strength of a weight sweep whose
maximum sat at exactly w = 0.50 and a fusion oracle worth only +0.008 over that. Then three more
views were added and the rule collapsed:

    wsi_titan 0.6596 · omics_xena 0.6749 · omics_combine 0.6510 · clinical 0.6856
    wsi_chief_mean 0.5836 · wsi_chief_dispersion 0.5554 · omics_hallmarks 0.6061

    equal weight over the three image views      0.6214   (worse than wsi_titan alone)
    equal weight over all seven                  0.6941   (worse than the three-view 0.7281)

Both declared comparisons failed on that rule: WSI+omics alone reached 0.6770 against DIMAF's
0.679, and the full arm reached 0.6941 against SurvPath-plus-the-same-clinical at 0.6963.

MECHANISM. Equal weighting is optimal among forecasters of COMPARABLE accuracy -- that is the whole
content of the small-sample combination result this campaign has been leaning on, and it is exactly
where the sweep's flat optimum came from, since the image and clinical arms differ by 0.026. It
says nothing about a pool whose members range from 0.5554 to 0.6856. A view at 0.5554 given the
same weight as one at 0.6856 is not a modest cost; it is one seventh of the ensemble spent on
something barely above chance.

THE COMPONENT, and it deliberately contains no arbitrary constant. Each view's weight is

    w_v = max(0, C_inner(v) - 0.5)

where C_inner(v) is that view's concordance estimated by inner cross-validation using ONLY the
training fold's patients. Equal weighting is the special case where every view has the same inner
skill, so the component degenerates into its own ablation rather than being a separate code path.
Nothing about the held-out fold enters the weights.

FALSIFIER, FIXED BEFORE THE RUN. C6 must beat the flat equal-weight rule over the SAME view pool by
more than 0.0145, the 2-SD paired reseed bar. It is not enough to beat it on the seven-view pool
where the flat rule is obviously handicapped: it must also not LOSE on the three-view pool, where
the flat rule is at its best.

REMOVAL ABLATION. `uniform weights` -- the identical code path with w_v = 1.

MATCHED CONTROL. Weights computed from inner CV on a PERMUTED training outcome. The permuted skills
have the same shape and the same spread of magnitudes; only which view earns which weight is
destroyed. If the control moves as much as the real weights do, the gain is the weighting
structure and not the skill estimate, and the run is void.

A THRESHOLD VARIANT is reported alongside for sensitivity, not as the arm: admit only views whose
inner skill clears 0.55 and equal-weight those. It carries an arbitrary constant, which is why it
is not the component.
"""

from __future__ import annotations

import argparse
import csv as _csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cindex, cpairs, fitapply, pct  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402


def inner_skill(X, t, e, tri, seed=0, folds=3, permute=False):
    """C-index of a view estimated by inner CV inside the training fold only."""
    rng = np.random.default_rng(seed + 21)
    tt, ee = t[tri], e[tri]
    if permute:
        q = rng.permutation(len(tri))
        tt, ee = tt[q], ee[q]
    inner = np.array_split(np.random.default_rng(seed + 21).permutation(len(tri)), folds)
    ip = np.zeros(len(tri))
    for ite in inner:
        itr = np.setdiff1d(np.arange(len(tri)), ite)
        ip[ite] = fitapply(X[tri][itr], tt[itr], ee[itr], X[tri][ite], seed)
    c = cindex(ip, tt, ee)
    return float(c) if c == c else 0.5


def build_arms(co, blocks, fi, seed=0):
    arms = {}
    for nm, X in blocks.items():
        s = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], seed)
        arms[nm] = co.fold_pct(s)
    return arms


def weighted_arm(co, blocks, arms, fi, pool, mode, seed=0):
    """Combine `pool` per fold. mode: uniform | skill | skill_permuted | threshold55."""
    out = np.full(len(co.keep), np.nan)
    used = {}
    for k, (tri, vai) in enumerate(fi):
        if len(tri) < 30 or not len(vai):
            continue
        w = {}
        for nm in pool:
            if mode == "uniform":
                w[nm] = 1.0
            else:
                c = inner_skill(blocks[nm], co.t, co.e, tri, seed,
                                permute=(mode == "skill_permuted"))
                if mode == "threshold55":
                    w[nm] = 1.0 if c >= 0.55 else 0.0
                else:
                    w[nm] = max(0.0, c - 0.5)
        tot = sum(w.values())
        if tot <= 0:
            w = {nm: 1.0 for nm in pool}
            tot = float(len(pool))
        used["fold_%d" % k] = {nm: round(w[nm] / tot, 4) for nm in pool}
        acc = np.zeros(len(vai))
        for nm in pool:
            acc += w[nm] * pct(arms[nm][vai])
        out[vai] = acc / tot
    return co.fold_pct(out), used


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--omics-dir", required=True)
    ap.add_argument("--chief", required=True)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    co = Cohort(args.root, args.titan, args.sp_dir)
    fi = co.fold_indices()
    blocks = {"wsi_titan": co.T}
    z = np.load(args.chief, allow_pickle=True)
    stem2case = {r["slide_id"].replace(".svs", ""): r["case_id"]
                 for r in _csv.DictReader(open(args.root + "/tcga_blca_meta.csv"))}
    percase = {}
    for i, nm in enumerate(z["names"]):
        c = stem2case.get(str(nm))
        if c:
            percase.setdefault(c, []).append(i)
    part = {}
    for key in ("mean", "std", "q10", "q90"):
        M = np.asarray(z[key], dtype=float)
        part[key] = np.vstack([M[percase[c]].mean(0) if c in percase else M.mean(0)
                               for c in co.keep])
    blocks["wsi_chief_mean"] = part["mean"]
    blocks["wsi_chief_dispersion"] = np.hstack([part["std"], part["q90"] - part["q10"]])
    for tag, sig in (("combine", "combine_signatures.csv"),
                     ("hallmarks", "hallmarks_signatures.csv"),
                     ("xena", "xena_signatures.csv")):
        rna = os.path.join(args.omics_dir, "rna_%s.csv" % tag)
        sp = os.path.join(args.omics_dir, sig)
        if os.path.exists(rna) and os.path.exists(sp):
            G, genes, _ = load_rna(rna, co.keep)
            P, _pn = pathway_matrix(G, genes, sp)
            if P.shape[1] >= 10:
                blocks["omics_%s" % tag] = P
    blocks["clinical"] = co.CLIN

    arms = build_arms(co, blocks, fi)
    ii, jj = cpairs(co.t, co.e)
    rep = {"artifact_type": "s5_c6_skill_weight", "reportable": False,
           "phase_of_origin": "exploratory", "n": len(co.keep), "events": int(co.e.sum()),
           "single_arms": {k: round(cidx(v, ii, jj), 4) for k, v in arms.items()},
           "bar": 0.0145,
           "falsifier": "beat the uniform rule on the SAME pool by >0.0145, on the seven-view pool "
                        "AND without losing on the three-view pool where uniform is at its best"}

    POOLS = {
        "all_seven_views": sorted(blocks),
        "three_view_titan_combine_clinical": ["wsi_titan", "omics_combine", "clinical"],
        "wsi_plus_omics_only_five": [k for k in sorted(blocks) if not k.startswith("clinical")],
    }
    res = {}
    for pname, pool in POOLS.items():
        res[pname] = {}
        for mode in ("uniform", "skill", "skill_permuted", "threshold55"):
            zz, used = weighted_arm(co, blocks, arms, fi, pool, mode)
            res[pname][mode] = {"cindex": round(cidx(zz, ii, jj), 4)}
            if mode in ("skill", "skill_permuted"):
                res[pname][mode]["weights_per_fold"] = used
            print("%-38s %-16s C=%.4f" % (pname, mode, res[pname][mode]["cindex"]),
                  file=sys.stderr, flush=True)
        u = res[pname]["uniform"]["cindex"]
        res[pname]["_verdict"] = {
            "skill_minus_uniform": round(res[pname]["skill"]["cindex"] - u, 4),
            "CONTROL_permuted_minus_uniform": round(
                res[pname]["skill_permuted"]["cindex"] - u, 4),
            "threshold55_minus_uniform": round(res[pname]["threshold55"]["cindex"] - u, 4)}
    rep["pools"] = res
    rep["read_the_control_first"] = ("if CONTROL_permuted_minus_uniform is comparable to "
                                     "skill_minus_uniform on a pool, that pool's result is VOID")

    # the two declared comparisons, recomputed under the surviving rule
    zz_all, _ = weighted_arm(co, blocks, arms, fi, POOLS["all_seven_views"], "skill")
    zz_no, _ = weighted_arm(co, blocks, arms, fi, POOLS["wsi_plus_omics_only_five"], "skill")
    spc = co.fold_pct(co.fold_pct(co.spr) + arms["clinical"])
    rep["comparison_C_THE_DECIDER"] = {
        "ours_wsi_plus_omics_skill_weighted": round(cidx(zz_no, ii, jj), 4),
        "incumbent_DIMAF": 0.679,
        "gap": round(cidx(zz_no, ii, jj) - 0.679, 4),
        "passes": bool(cidx(zz_no, ii, jj) - 0.679 > 0.0145)}
    rep["comparison_B_input_parity"] = {
        "ours_full_skill_weighted": round(cidx(zz_all, ii, jj), 4),
        "survpath_plus_same_clinical": round(cidx(spc, ii, jj), 4),
        "gap": round(cidx(zz_all, ii, jj) - cidx(spc, ii, jj), 4)}

    text = json.dumps(rep, indent=1)
    if args.out:
        open(args.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
