#!/usr/bin/env python3
"""Is transcriptomics a THIRD orthogonal signal on these folds, and how strong is it?

Pre-freeze, exploratory. NOT reportable.

WHY IT MATTERS NOW. Three measurements taken today agree that the two-arm ceiling is close:
the best FIXED linear fusion of the image and clinical scores, fitted in-sample on the held-out
fold, reaches 0.7268 against the 0.7191 the deployed equal-weight rule already gets; the weight
sweep is flat from w=0.30 to w=0.60 with its maximum at exactly 0.50; and a simulation calibrated
to the two arms' own C-indices puts their rank average at 0.7421 +/- 0.0238. Nothing about the
COMBINATION RULE is leaving much on the table. More information is.

WHAT HAS NEVER BEEN MEASURED HERE. Every competitor on this benchmark -- SurvPath, MCAT, MOTCat,
MMP, PIBD, DIMAF -- consumes the TCGA-BLCA transcriptome, and this project's arm never has. That
is two separate problems in one: a reviewer's parity objection, and an unexamined information
source. Neither can be settled by argument.

THREE REPRESENTATIONS, because the choice is not obvious and picking one after seeing the answers
would be the selection this campaign keeps measuring in other people's work:

  genes        ridge Cox over all 4,999 genes SurvPath ships
  pathways     the same genes averaged into the pathway groups SurvPath's own signature file
               defines, then ridge Cox -- this is close to what SurvPath's omics branch consumes
  screened     top-k genes by the univariate Cox score test, SELECTED STRICTLY INSIDE THE
               TRAINING FOLD

That last one carries a warning this project has already paid for once. Screening features against
the outcome over the whole cohort before cross-validation inflated an AUROC by +0.352 in an
earlier experiment here, and a related leak of 44 sign bits was worth +0.127 on a C-index. The
screen below refits inside every training fold and never sees a validation patient; the point of
including it at all is that a leaky version of the same arm is the single easiest way to fake this
result, so the honest version needs to be on the record.

AND THE QUESTION BEHIND ALL OF IT: not "is omics predictive" -- it is -- but whether it is
predictive of something the image and the stage do not already say. So every arm is reported with
its rank correlation to the other two and with its C-index restricted to the pairs the clinical
score cannot separate.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, cox_fit, fitapply  # noqa: E402


def spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    d = float(np.sqrt((ra @ ra) * (rb @ rb)))
    return float(ra @ rb / d) if d else float("nan")


def load_rna(path, keep):
    """Genes x cases aligned to `keep`; returns (matrix, gene names, mask of cases present)."""
    with open(path) as fh:
        rd = csv.reader(fh)
        header = next(rd)
        genes = header[1:]
        rows, ids = [], []
        for r in rd:
            ids.append(r[0])
            rows.append([float(x) if x else 0.0 for x in r[1:]])
    M = np.asarray(rows, dtype=float)
    pos = {c: i for i, c in enumerate(ids)}
    have = np.array([c in pos for c in keep])
    X = np.zeros((len(keep), M.shape[1]))
    for i, c in enumerate(keep):
        if c in pos:
            X[i] = M[pos[c]]
    return X, genes, have


def pathway_matrix(X, genes, sig_path):
    """Average the gene columns inside each signature into one column per pathway."""
    gi = {g: i for i, g in enumerate(genes)}
    with open(sig_path) as fh:
        rd = csv.reader(fh)
        names = next(rd)
        cols = [[] for _ in names]
        for r in rd:
            for j, g in enumerate(r):
                if j < len(cols) and g and g in gi:
                    cols[j].append(gi[g])
    keep_p, out = [], []
    for j, idxs in enumerate(cols):
        u = sorted(set(idxs))
        if len(u) >= 3:
            keep_p.append(names[j])
            out.append(X[:, u].mean(1))
    return (np.column_stack(out) if out else np.zeros((X.shape[0], 0))), keep_p


def score_test(X, t, e):
    """Univariate Cox score statistic at beta=0, per column. Censoring-aware by construction.

    A previous version of a screen in this project ranked genes by the difference in mean between
    patients who had an event and patients who did not, which throws the censoring away and is not
    a survival statistic at all. This is the exact score test: U = sum over events of (x_i - xbar
    over the risk set), standardised by its variance under the null.
    """
    o = np.argsort(t)
    Xo, eo = X[o], e[o]
    n = len(t)
    ev = np.flatnonzero(eo == 1)
    csum = np.cumsum(Xo[::-1], axis=0)[::-1]
    csq = np.cumsum((Xo ** 2)[::-1], axis=0)[::-1]
    sz = np.arange(n, 0, -1).astype(float)[:, None]
    mean_rs = csum / sz
    var_rs = np.maximum(csq / sz - mean_rs ** 2, 0.0)
    U = (Xo[ev] - mean_rs[ev]).sum(0)
    V = var_rs[ev].sum(0)
    return np.abs(U) / np.sqrt(np.maximum(V, 1e-12))


def screened_arm(X, t, e, fi, topk, seed=0):
    """Ridge Cox on the top-k genes, screened inside each TRAINING fold only."""
    s = np.full(X.shape[0], np.nan)
    for tri, vai in fi:
        if len(tri) < 30 or not len(vai):
            continue
        mu, sd = X[tri].mean(0), X[tri].std(0) + 1e-9
        Z = (X - mu) / sd
        st = score_test(Z[tri], t[tri], e[tri])
        sel = np.argsort(-st)[:topk]
        s[vai] = fitapply(Z[tri][:, sel], t[tri], e[tri], Z[vai][:, sel], seed)
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--rna", required=True)
    ap.add_argument("--sig", required=True)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    G, genes, have = load_rna(a.rna, co.keep)
    P, pnames = pathway_matrix(G, genes, a.sig)

    rep = {"artifact_type": "s5_omics_arm", "reportable": False,
           "phase_of_origin": "exploratory",
           "n_cases": len(co.keep), "events": int(co.e.sum()),
           "cases_with_rna": int(have.sum()),
           "cases_without_rna_get_the_training_fold_mean": int((~have).sum()),
           "n_genes": len(genes), "n_pathways": len(pnames)}

    # baseline arms on the same folds
    p_ti = np.full(len(co.keep), np.nan)
    p_cl = np.full(len(co.keep), np.nan)
    for tri, vai in fi:
        if len(tri) < 30 or not len(vai):
            continue
        p_ti[vai] = fitapply(co.T[tri], co.t[tri], co.e[tri], co.T[vai], 0)
        p_cl[vai] = fitapply(co.CLIN[tri], co.t[tri], co.e[tri], co.CLIN[vai], 0)
    z_ti, z_cl = co.fold_pct(p_ti), co.fold_pct(p_cl)

    omic = {}
    for nm, X in (("genes_all", G), ("pathways", P)):
        s = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], 0)
        omic[nm] = co.fold_pct(s)
        print("done %s" % nm, file=sys.stderr, flush=True)
    for k in (32, 128, 512):
        omic["screened_top%d" % k] = co.fold_pct(screened_arm(G, co.t, co.e, fi, k))
        print("done screened_top%d" % k, file=sys.stderr, flush=True)

    ii, jj = cpairs(co.t, co.e)
    gap = np.abs(z_cl[ii] - z_cl[jj])
    tied = gap <= float(np.quantile(gap, 0.10))

    rows = {}
    for nm, z in list(omic.items()) + [("titan", z_ti), ("clinical", z_cl),
                                       ("survpath_rerun", co.fold_pct(co.spr))]:
        rows[nm] = {
            "cindex": round(cidx(z, ii, jj), 4),
            "cindex_on_clinically_tied_pairs_q10": round(cidx(z, ii[tied], jj[tied]), 4),
            "spearman_vs_clinical": round(spearman(z, z_cl), 4),
            "spearman_vs_titan": round(spearman(z, z_ti), 4),
        }
    rep["arms"] = rows

    best = max(omic, key=lambda k: cidx(omic[k], ii, jj))
    rep["best_omics_representation"] = {
        "name": best, "cindex": round(cidx(omic[best], ii, jj), 4),
        "caveat": "chosen by looking at the pooled outcome across five representations, so it "
                  "carries a selection term of roughly sd*sqrt(2*ln 5); any arm built on it must "
                  "reselect inside the training fold before it is comparable"}

    z_om = omic[best]
    combos = {
        "titan + clinical": co.fold_pct(z_ti + z_cl),
        "omics + clinical": co.fold_pct(z_om + z_cl),
        "titan + omics": co.fold_pct(z_ti + z_om),
        "titan + omics + clinical": co.fold_pct(z_ti + z_om + z_cl),
    }
    rep["rank_average_combinations"] = {
        k: {"cindex": round(cidx(v, ii, jj), 4),
            "on_clinically_tied_pairs_q10": round(cidx(v, ii[tied], jj[tied]), 4)}
        for k, v in combos.items()}
    rep["three_way_minus_two_way"] = round(
        cidx(combos["titan + omics + clinical"], ii, jj) - cidx(combos["titan + clinical"], ii, jj), 4)

    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
