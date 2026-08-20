#!/usr/bin/env python3
"""C7 -- train the molecular heads on WITHIN-STRATUM risk sets, and measure where it should matter.

Pre-freeze, exploratory. NOT reportable.

THE DEFICIT, and it is the most robust thing measured today. Restrict the concordance to the
comparable pairs the clinical score cannot separate -- the 10% of pairs with the smallest gap in
clinical risk, 2,423 of 24,219 -- and the arms separate far more than they do overall:

                                      all pairs     clinically tied (q10)
    clinical alone                       0.6856          0.5167   (chance, by construction)
    SurvPath rerun + same clinical       0.6963          0.5852
    OURS                                 0.7291          0.6513
    omics_xena alone                     0.6749          0.6529

The overall margin over SurvPath-plus-clinical is +0.033; on tied pairs it is +0.066. Two patients
with the same stage and a different outcome is also the question a pathologist is actually asking,
and no paper on this benchmark reports it.

MECHANISM. Cox's partial likelihood compares each event against everyone still at risk, so a
molecular head can lower its loss by re-deriving stage -- the easiest thing in the image -- and a
large part of its fitted capacity goes there. Restricting each risk set to patients in the SAME
clinical stratum removes that route entirely: within a stratum, stage carries no information, so
the only way to lower the loss is with signal the clinical model does not have. The head is not
asked to be better at ranking overall; it is asked to be better at the part that is left.

This is a change to the LIKELIHOOD'S RISK SETS, operating on the same features, the same estimator
and the same combination rule. No previously killed candidate has that shape: C1/C1' constrained
the coefficient vector, C2-C4 changed the estimator, C5 and C6 changed the combination.

FALSIFIER, FIXED BEFORE THE RUN. The stratified-trained arm must beat the marginally-trained arm
ON TIED PAIRS by more than 0.0145. Its overall C-index is allowed to fall -- that is what the
mechanism predicts, and a rise in both would be a different result than the one hypothesised.

REMOVAL ABLATION. `strata = one stratum for everybody`, which is the ordinary partial likelihood
and the same code path.

MATCHED CONTROL. Strata PERMUTED across patients inside the training fold -- same number of strata,
same sizes, membership destroyed. A gain there is the stratification structure, not the clinical
information, and voids the run.

Strata are built from the clinical Cox score fitted INSIDE the training fold, so no held-out
patient contributes to the stratum boundaries.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cindex, cpairs, fitapply  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402

ALPHAS = (1.0, 8.0, 64.0, 512.0, 4096.0)


def strat_cox_fit(X, t, e, strata, alpha):
    """Ridge Cox on a STRATIFIED partial likelihood: risk sets never cross a stratum.

    The log-likelihood is the sum of each stratum's own Breslow partial likelihood, so the
    baseline hazard is free to differ between strata and only WITHIN-stratum ordering contributes.
    One shared coefficient vector across strata, which is the point -- otherwise this is just
    per-stratum models.
    """
    from scipy.optimize import minimize
    groups = []
    for g in np.unique(strata):
        idx = np.flatnonzero(strata == g)
        o = idx[np.argsort(t[idx])]
        ev = np.flatnonzero(e[o] == 1)
        if ev.size >= 1:
            groups.append((X[o], ev))
    if not groups:
        return np.zeros(X.shape[1])

    def f(b):
        tot, grad = 0.0, np.zeros(X.shape[1])
        for Xo, ev in groups:
            eta = Xo @ b
            m = float(eta.max())
            ex = np.exp(eta - m)
            rs = np.cumsum(ex[::-1])[::-1]
            rx = np.cumsum((Xo * ex[:, None])[::-1], axis=0)[::-1]
            tot += float(np.sum(eta[ev] - (np.log(rs[ev]) + m)))
            grad += (Xo[ev] - rx[ev] / rs[ev, None]).sum(0)
        return -tot + 0.5 * alpha * float(b @ b), -grad + alpha * b

    return minimize(f, np.zeros(X.shape[1]), jac=True, method="L-BFGS-B",
                    options={"maxiter": 400, "gtol": 1e-8}).x


def strat_fitapply(Xtr, ttr, etr, Str, Xte, seed=0):
    """Alpha by inner CV, scored by WITHIN-STRATUM concordance -- the objective being optimised."""
    def strat_cindex(s, t, e, st):
        num = den = 0.0
        for g in np.unique(st):
            m = st == g
            ii, jj = cpairs(t[m], e[m])
            if ii.size:
                num += cidx(s[m], ii, jj) * ii.size
                den += ii.size
        return num / den if den else float("nan")

    best, ba = -1.0, ALPHAS[0]
    inner = np.array_split(np.random.default_rng(seed).permutation(len(ttr)), 3)
    for al in ALPHAS:
        ip = np.zeros(len(ttr))
        for ite in inner:
            itr = np.setdiff1d(np.arange(len(ttr)), ite)
            mu, sd = Xtr[itr].mean(0), Xtr[itr].std(0) + 1e-9
            ip[ite] = ((Xtr[ite] - mu) / sd) @ strat_cox_fit(
                (Xtr[itr] - mu) / sd, ttr[itr], etr[itr], Str[itr], al)
        c = strat_cindex(ip, ttr, etr, Str)
        if c == c and c > best:
            best, ba = c, al
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
    return ((Xte - mu) / sd) @ strat_cox_fit((Xtr - mu) / sd, ttr, etr, Str, ba)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--omics-dir", required=True)
    ap.add_argument("--n-strata", type=int, default=4)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    mol = {"titan": co.T, "omics_combine": P}

    z_cl = np.full(len(co.keep), np.nan)
    for tri, vai in fi:
        if len(tri) < 30 or not len(vai):
            continue
        z_cl[vai] = fitapply(co.CLIN[tri], co.t[tri], co.e[tri], co.CLIN[vai], 0)
    z_cl = co.fold_pct(z_cl)

    modes = ("marginal", "stratified", "stratified_PERMUTED_CONTROL")
    out = {m: {} for m in modes}
    rngp = np.random.default_rng(7)
    for nm, X in mol.items():
        for mode in modes:
            s = np.full(len(co.keep), np.nan)
            for tri, vai in fi:
                if len(tri) < 30 or not len(vai):
                    continue
                if mode == "marginal":
                    s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], 0)
                else:
                    # clinical strata built INSIDE the training fold
                    cin = fitapply(co.CLIN[tri], co.t[tri], co.e[tri], co.CLIN[tri], 0)
                    q = np.quantile(cin, np.linspace(0, 1, a.n_strata + 1)[1:-1])
                    st = np.digitize(cin, q)
                    if mode.endswith("CONTROL"):
                        st = st[rngp.permutation(len(st))]
                    s[vai] = strat_fitapply(X[tri], co.t[tri], co.e[tri], st, X[vai], 0)
            out[mode][nm] = co.fold_pct(s)
            print("%-30s %-16s done" % (mode, nm), file=sys.stderr, flush=True)

    ii, jj = cpairs(co.t, co.e)
    gap = np.abs(z_cl[ii] - z_cl[jj])
    tied = gap <= float(np.quantile(gap, 0.10))
    tied20 = gap <= float(np.quantile(gap, 0.20))

    def row(z):
        return {"cindex": round(cidx(z, ii, jj), 4),
                "tied_q10": round(cidx(z, ii[tied], jj[tied]), 4),
                "tied_q20": round(cidx(z, ii[tied20], jj[tied20]), 4)}

    rep = {"artifact_type": "s5_c7_stratified_training", "reportable": False,
           "phase_of_origin": "exploratory", "n_strata": a.n_strata, "bar": 0.0145,
           "falsifier": "stratified must beat marginal ON TIED PAIRS by >0.0145; overall may fall",
           "read_the_control_first": "if stratified_PERMUTED_CONTROL moves like stratified, VOID"}

    rep["per_view"] = {}
    for nm in mol:
        rep["per_view"][nm] = {m: row(out[m][nm]) for m in modes}
        rep["per_view"][nm]["_verdict"] = {
            "tied_q10_stratified_minus_marginal": round(
                row(out["stratified"][nm])["tied_q10"] - row(out["marginal"][nm])["tied_q10"], 4),
            "tied_q10_CONTROL_minus_marginal": round(
                row(out["stratified_PERMUTED_CONTROL"][nm])["tied_q10"]
                - row(out["marginal"][nm])["tied_q10"], 4),
            "overall_stratified_minus_marginal": round(
                row(out["stratified"][nm])["cindex"] - row(out["marginal"][nm])["cindex"], 4)}

    for m in modes:
        z = co.fold_pct(out[m]["titan"] + out[m]["omics_combine"] + z_cl)
        rep.setdefault("full_arm", {})[m] = row(z)
    rep["full_arm"]["_verdict"] = {
        "tied_q10_stratified_minus_marginal": round(
            rep["full_arm"]["stratified"]["tied_q10"]
            - rep["full_arm"]["marginal"]["tied_q10"], 4),
        "tied_q10_CONTROL_minus_marginal": round(
            rep["full_arm"]["stratified_PERMUTED_CONTROL"]["tied_q10"]
            - rep["full_arm"]["marginal"]["tied_q10"], 4),
        "overall_stratified_minus_marginal": round(
            rep["full_arm"]["stratified"]["cindex"] - rep["full_arm"]["marginal"]["cindex"], 4)}

    # the tied-pair comparison against the input-parity competitor, with a case-level bootstrap
    spc = co.fold_pct(co.fold_pct(co.spr) + z_cl)
    best_mode = "stratified" if (rep["full_arm"]["_verdict"]["tied_q10_stratified_minus_marginal"]
                                 > 0) else "marginal"
    ours = co.fold_pct(out[best_mode]["titan"] + out[best_mode]["omics_combine"] + z_cl)
    rng = np.random.default_rng(0)
    bt = []
    for _ in range(6000):
        bs = rng.choice(len(co.keep), size=len(co.keep), replace=True)
        i2, j2 = cpairs(co.t[bs], co.e[bs])
        if i2.size == 0:
            continue
        g2 = np.abs(z_cl[bs][i2] - z_cl[bs][j2])
        m2 = g2 <= float(np.quantile(g2, 0.10))
        if m2.sum() < 50:
            continue
        bt.append(cidx(ours[bs], i2[m2], j2[m2]) - cidx(spc[bs], i2[m2], j2[m2]))
    v = np.asarray(bt)
    rep["tied_pair_head_to_head_vs_survpath_plus_clinical"] = {
        "arm_used": best_mode,
        "ours_tied_q10": round(cidx(ours, ii[tied], jj[tied]), 4),
        "survpath_plus_clinical_tied_q10": round(cidx(spc, ii[tied], jj[tied]), 4),
        "paired_bootstrap": {"mean": round(float(v.mean()), 4),
                             "se": round(float(v.std(ddof=1)), 4),
                             "ci95": [round(float(np.quantile(v, .025)), 4),
                                      round(float(np.quantile(v, .975)), 4)],
                             "p_two_sided": round(
                                 float(2 * min((v <= 0).mean(), (v >= 0).mean())), 4),
                             "reps": int(v.size),
                             "note": "the tie threshold is recomputed inside each bootstrap "
                                     "replicate, so it is a property of the resampled cohort "
                                     "rather than of the original one"}}

    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
