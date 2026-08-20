#!/usr/bin/env python3
"""What is the model actually reading? Bladder-cancer biology behind the two molecular arms.

`phase_of_origin: frozen_post_hoc`. Descriptive. None of it feeds back into the method, and none
of it is a discovery claim -- these are associations in 359 patients, at a sample size that the
power analysis in this same campaign shows cannot resolve the margins the method is judged on.

=========================== WHAT IS ASKED, AND WHY EACH QUESTION IS ANSWERABLE

B1  SIGNATURE SCORES. Eight axes of bladder biology, each a z-scored mean over the marker genes
    present in SurvPath's own 4,999-gene matrix: basal and luminal (the molecular subtypes that
    organise this disease), EMT, immune infiltration, proliferation, p53/cell-cycle, the FGFR3
    axis, and stroma. Gene membership is listed in the output so a reader can check it; panels
    are reported with how many of their members were available, because a panel at 6 of 13 is a
    weaker instrument than one at 12 of 13 and the difference should not be hidden in a mean.

B2  DO THEY PREDICT SURVIVAL? Univariate Cox score test per axis and per pathway, Benjamini-
    Hochberg across the 275 pathways. Descriptive: computed on the whole cohort, which is correct
    for an association and would be leakage if it selected anything. It selects nothing.

B3  WHAT DOES THE IMAGE ARM TRACK? The question a whole-slide paper has to answer. TITAN's
    dimensions are not interpretable, but its out-of-fold SCORE is, and correlating it with each
    biological axis and each pathway says what morphology is a proxy for. A slide score that were
    merely a stage detector, or merely a proliferation readout, would show it here.

B4  AND THE OMICS ARM? The same, which turns the earlier "these two arms are only rho=0.27
    correlated" into a statement about which biology each one carries.

B5  DOES THE METHOD ADD WITHIN A MOLECULAR SUBTYPE? Basal and luminal tumours differ in prognosis
    and in treatment. Within a subtype, subtype carries no information -- the same construction as
    the per-stage analysis -- so concordance there is what the model adds beyond the molecular
    class.

B6  CASE STUDIES, biologically framed. The selection rule is FROZEN BELOW, before any output was
    inspected, and every case it admits is reported.

=========================== THE CASE-STUDY RULE, FIXED BEFORE THE RUN

  POPULATION  patients in the middle tertile of the clinical arm's out-of-fold risk -- the
              stage-ambiguous band, where an image-and-molecule model is being asked to add
              something.
  SELECTION   the 4 with the highest and the 4 with the lowest full-method score, and in addition
              the 2 patients with the largest DISAGREEMENT between the image arm and the omics arm
              in each direction. Twelve cases, no discretion.
  REPORTED    the three component scores, the eight biological axis scores, stage, observed time
              and event. Every case the rule admits, including the ones that embarrass the model.
"""

from __future__ import annotations

import argparse
import collections
import math
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, fitapply  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix, score_test  # noqa: E402
from s6_amend_clinical import dimaf_clinical  # noqa: E402
from s7_survival_metrics import logrank, chi2_sf  # noqa: E402

PANELS = {
    "basal": ["KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT14", "KRT1", "CD44", "COL17A1", "DSC3",
              "GSDMC", "SDC1", "MSN"],
    "luminal": ["GATA3", "FOXA1", "PPARG", "KRT20", "UPK1A", "UPK1B", "UPK2", "UPK3A", "XBP1",
                "ERBB2", "ERBB3", "CYP2J2", "SNX31"],
    "EMT": ["ZEB1", "ZEB2", "VIM", "SNAI1", "SNAI2", "TWIST1", "CDH2", "CLDN3", "CLDN4", "CLDN7"],
    "immune": ["CD8A", "CD3D", "CD3E", "GZMB", "PRF1", "IFNG", "CXCL9", "CXCL10", "PDCD1", "CD274",
               "CTLA4", "LAG3", "ITGAE"],
    "proliferation": ["MKI67", "TOP2A", "CCNB1", "BIRC5", "AURKA", "RRM2", "TYMS", "MCM2"],
    "p53_cell_cycle": ["TP53", "RB1", "CDKN2A", "E2F3", "CCND1", "MDM2"],
    "FGFR3_axis": ["FGFR3", "CCND1", "SHH", "TP63"],
    "stroma": ["FAP", "PDGFRB", "ACTA2", "COL1A1", "COL3A1", "THBS2", "POSTN"],
}
# CLDN3/4/7 are epithelial-tight-junction genes that go DOWN in claudin-low/EMT tumours, so their
# sign is flipped before averaging. Stated here rather than silently applied.
INVERTED = {"EMT": {"CLDN3", "CLDN4", "CLDN7"}}


def spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    d = float(np.sqrt((ra @ ra) * (rb @ rb)))
    return float(ra @ rb / d) if d else float("nan")


def spearman_p(r, n):
    """Two-sided p for a Spearman correlation via the t approximation."""
    if not np.isfinite(r) or n < 5 or abs(r) >= 1:
        return float("nan")
    t = r * np.sqrt((n - 2) / max(1e-12, 1 - r * r))
    try:
        from scipy.stats import t as _t
        return float(2 * _t.sf(abs(t), n - 2))
    except Exception:
        return float("nan")


def bh(pvals):
    """Benjamini-Hochberg q-values, order preserved."""
    ok = [(k, p) for k, p in pvals.items() if np.isfinite(p)]
    m = len(ok)
    out = {k: float("nan") for k in pvals}
    if not m:
        return out
    ok.sort(key=lambda kv: kv[1])
    prev = 1.0
    for i in range(m - 1, -1, -1):
        k, p = ok[i]
        q = min(prev, p * m / (i + 1))
        out[k] = round(float(q), 4)
        prev = q
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("root", "titan", "sp-dir", "omics-dir", "dimaf-dir", "out"):
        ap.add_argument("--" + f, required=True)
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    ii, jj = cpairs(co.t, co.e)
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, pnames = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    CLIN = np.hstack([dc["age"], dc["fem"], dc["stage"]])
    gi = {g: k for k, g in enumerate(genes)}
    n = len(co.keep)

    def arm(X, seed=0):
        s = np.full(n, np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], seed)
        return co.fold_pct(s)

    zt, zo, zc = arm(co.T), arm(P), arm(CLIN)
    OURS = co.fold_pct(zt + zo + zc)
    rep = {"artifact_type": "s8_biology", "phase_of_origin": "frozen_post_hoc",
           "n": n, "events": int(co.e.sum()),
           "status": "DESCRIPTIVE. Associations in 359 patients. Nothing here feeds back into the "
                     "method and nothing here is a discovery claim."}

    # ============================================================ B1 signature scores
    Z = (G - G.mean(0)) / (G.std(0) + 1e-9)
    sig, meta = {}, {}
    for name, members in PANELS.items():
        idx = [gi[g] for g in members if g in gi]
        if len(idx) < 3:
            meta[name] = {"available": len(idx), "of": len(members), "skipped": "fewer than 3"}
            continue
        sgn = np.array([-1.0 if g in INVERTED.get(name, set()) else 1.0
                        for g in members if g in gi])
        sig[name] = (Z[:, idx] * sgn[None, :]).mean(1)
        meta[name] = {"available": len(idx), "of": len(members),
                      "genes_used": [g for g in members if g in gi],
                      "inverted": sorted(INVERTED.get(name, set()) & set(genes))}
    rep["B1_signatures"] = {"panels": meta,
                            "construction": "z-scored mean over available member genes; the "
                                            "tight-junction genes in the EMT panel are sign-flipped "
                                            "because they fall rather than rise with EMT"}

    # ============================================================ B2 do they predict survival
    b2 = {}
    for name, v in sig.items():
        s = float(score_test(v[:, None], co.t, co.e)[0])
        c = cidx(co.fold_pct(v), ii, jj)
        b2[name] = {"cox_score_statistic": round(s, 3),
                    "p": float("%.3g" % (2 * (1 - 0.5 * (1 + math.erf(s / math.sqrt(2))))
                                         if s == s else float("nan"))),
                    "cindex_as_a_raw_score": round(c, 4)}
    # pathways, BH across all 275
    st = score_test((P - P.mean(0)) / (P.std(0) + 1e-9), co.t, co.e)
    from math import erf, sqrt
    praw = {pnames[k]: float(2 * (1 - 0.5 * (1 + erf(st[k] / sqrt(2))))) for k in range(len(pnames))}
    q = bh(praw)
    top = sorted(praw, key=lambda k: praw[k])[:15]
    rep["B2_survival_association"] = {
        "axes": b2,
        "pathways_tested": len(pnames),
        "pathways_with_BH_q_below_0.10": int(sum(1 for k in q if q[k] < 0.10)),
        "top_15_pathways": [{"pathway": k, "p": float("%.3g" % praw[k]), "q_BH": q[k]} for k in top],
        "note": "computed on the whole cohort, which is correct for an ASSOCIATION and would be "
                "leakage if it selected anything. It selects nothing -- the frozen arm uses all 275."}
    print("B2 done: %d of %d pathways at q<0.10" % (rep["B2_survival_association"]["pathways_with_BH_q_below_0.10"], len(pnames)),
          file=sys.stderr, flush=True)

    # ============================================================ B3/B4 what each arm tracks
    def tracks(score, label):
        ax = {}
        for name, v in sig.items():
            r = spearman(score, v)
            ax[name] = {"spearman": round(r, 4), "p": float("%.3g" % spearman_p(r, n))}
        pr = {pnames[k]: spearman(score, P[:, k]) for k in range(len(pnames))}
        pp = {k: spearman_p(v, n) for k, v in pr.items()}
        pq = bh(pp)
        order = sorted(pr, key=lambda k: -abs(pr[k]))[:12]
        return {"axes": ax,
                "pathways_with_BH_q_below_0.05": int(sum(1 for k in pq if pq[k] < 0.05)),
                "top_12_pathways_by_absolute_correlation":
                    [{"pathway": k, "spearman": round(pr[k], 3), "q_BH": pq[k]} for k in order],
                "label": label}

    rep["B3_what_the_image_arm_tracks"] = tracks(zt, "TITAN slide-embedding arm, out-of-fold score")
    rep["B4_what_the_omics_arm_tracks"] = tracks(zo, "pathway-expression arm, out-of-fold score")
    rep["B4b_contrast"] = {
        "spearman_image_vs_omics_arm": round(spearman(zt, zo), 4),
        "axes_side_by_side": {
            k: {"image": rep["B3_what_the_image_arm_tracks"]["axes"][k]["spearman"],
                "omics": rep["B4_what_the_omics_arm_tracks"]["axes"][k]["spearman"]}
            for k in sig},
        "reading": "the two arms are only weakly correlated with each other; this says whether that "
                   "is because they read DIFFERENT biology or because one of them reads none"}
    print("B3/B4 done", file=sys.stderr, flush=True)

    # ============================================================ B5 within a molecular subtype
    bl, lu = sig.get("basal"), sig.get("luminal")
    b5 = {}
    if bl is not None and lu is not None:
        d = bl - lu                                   # basal-minus-luminal axis
        grp = (d > np.median(d)).astype(int)          # 0 = luminal-leaning, 1 = basal-leaning
        b5["axis"] = ("basal-minus-luminal signature difference, split at its median. A continuous "
                      "split, not a published classifier: the released files carry no subtype call "
                      "and we do not have the training data to reproduce one.")
        c2, df, _ = logrank(co.t, co.e, grp)
        b5["subtype_axis_itself"] = {
            "n_luminal_leaning": int((grp == 0).sum()), "n_basal_leaning": int((grp == 1).sum()),
            "events": [int(co.e[grp == 0].sum()), int(co.e[grp == 1].sum())],
            "logrank_p": float("%.3g" % chi2_sf(c2, df)),
            "cindex_of_the_axis_as_a_risk_score": round(cidx(co.fold_pct(d), ii, jj), 4)}
        for k, lab in ((0, "luminal_leaning"), (1, "basal_leaning")):
            m = grp == k
            i2, j2 = cpairs(co.t[m], co.e[m])
            if i2.size < 100:
                b5[lab] = {"n": int(m.sum()), "skipped": "fewer than 100 comparable pairs"}
                continue
            b5[lab] = {"n": int(m.sum()), "events": int(co.e[m].sum()), "pairs": int(i2.size),
                       "OURS": round(cidx(OURS[m], i2, j2), 4),
                       "image": round(cidx(zt[m], i2, j2), 4),
                       "omics": round(cidx(zo[m], i2, j2), 4),
                       "clinical": round(cidx(zc[m], i2, j2), 4)}
    rep["B5_within_molecular_subtype"] = b5
    print("B5 done", file=sys.stderr, flush=True)

    # ============================================================ B6 case studies, frozen rule
    qc = np.quantile(zc, [1 / 3, 2 / 3])
    band = (zc > qc[0]) & (zc <= qc[1])
    idx = np.flatnonzero(band)
    by_score = idx[np.argsort(OURS[idx])]
    disagree = idx[np.argsort(zt[idx] - zo[idx])]
    picked, why = [], {}
    for i in list(by_score[-4:][::-1]):
        picked.append(i); why[i] = "highest full-method score in the band"
    for i in list(by_score[:4]):
        picked.append(i); why[i] = "lowest full-method score in the band"
    for i in list(disagree[-2:][::-1]):
        if i not in why:
            picked.append(i); why[i] = "image arm far ABOVE omics arm"
    for i in list(disagree[:2]):
        if i not in why:
            picked.append(i); why[i] = "omics arm far ABOVE image arm"
    stage_lab = {k: s for k, s in enumerate(dc["stage_levels"])}
    sn = np.array([np.argmax(r) if r.sum() else -1 for r in dc["stage"]])
    cases = []
    for i in picked:
        cases.append({
            "case_id": co.keep[i], "selected_as": why[i],
            "ours": round(float(OURS[i]), 3), "image": round(float(zt[i]), 3),
            "omics": round(float(zo[i]), 3), "clinical": round(float(zc[i]), 3),
            "stage": stage_lab.get(int(sn[i])) if sn[i] >= 0 else None,
            **{"sig_" + k: round(float(v[i]), 2) for k, v in sig.items()},
            "observed_months": round(float(co.t[i]), 1), "event": bool(co.e[i] == 1)})
    hi = [c for c in cases if "highest" in c["selected_as"]]
    lo = [c for c in cases if "lowest" in c["selected_as"]]
    rep["B6_case_studies"] = {
        "selection_rule": "FIXED IN THE FILE BEFORE THE RUN: within the middle tertile of the "
                          "clinical arm's out-of-fold risk, the 4 highest and 4 lowest full-method "
                          "scores, plus the 2 largest image-over-omics and 2 largest "
                          "omics-over-image disagreements. All are reported.",
        "band_size": int(band.sum()), "cases": cases,
        "summary": {"high": {"n": len(hi), "events": sum(c["event"] for c in hi),
                             "median_months": round(float(np.median([c["observed_months"] for c in hi])), 1)},
                    "low": {"n": len(lo), "events": sum(c["event"] for c in lo),
                            "median_months": round(float(np.median([c["observed_months"] for c in lo])), 1)}},
        "caution": "eight-to-twelve patients cannot support an inference. These are illustrations "
                   "of what the model reads, chosen by a rule fixed in advance, and the ones that "
                   "contradict it are in the table."}
    print("B6 done", file=sys.stderr, flush=True)

    open(a.out, "w").write(json.dumps(rep, indent=1) + "\n")
    print("wrote", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
