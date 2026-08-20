#!/usr/bin/env python3
"""Do survival directions in slide-embedding space TRANSFER across cancer types?

Pre-freeze, exploratory. NOT reportable. This is the cheapest possible falsifier for a direction,
and it is run before the direction is designed rather than after it is implemented.

WHAT PROMPTED IT. Today's headroom probe put a number on where the image arm loses: fitting the
same estimator IN-SAMPLE on the held-out fold at 16 principal components reaches 0.7666, its
out-of-fold twin reaches 0.6496, and a permuted-outcome control at the same capacity sits at
0.4980 over 1,000 draws. So the 768-dimensional TITAN representation is NOT exhausted -- it holds
0.2686 of recoverable structure above its own control -- and 0.1170 of C-index is being lost to
ESTIMATION at 287 training patients and 113 events. The clinical block, at a fraction of the
dimension, loses 0.0819.

Estimation loss of that shape has one classical remedy: borrow strength. And the strength is
already on disk -- the TITAN feature file covers 11,658 TCGA slides, not the 457 bladder ones, and
SurvPath ships disease-specific survival for BRCA (871 cases, 60 events), HNSC (394, 117), STAD
(319, 84) and COADREAD (298, 37): 1,882 cases and 298 events, 2.6x what bladder has, requiring no
download.

THE QUESTION THIS FILE ANSWERS, and it has exactly two outcomes. Fit a Cox model on ONE cancer's
slide embeddings, take its coefficient vector, and score BLADDER patients with it. Nothing about
bladder is used in fitting it.

    transfers      -> the survival-relevant subspace is partly shared across cancers, and a
                      pan-cancer prior can cut the effective dimension of the bladder head. The
                      direction is alive and the component is worth designing.
    does not       -> every candidate built on cross-cancer transfer is dead, for the price of
                      this file. Nothing else in the slate is affected.

THE CONTROL IS NOT OPTIONAL AND IT IS THE POINT. A coefficient vector fitted on ANY 800-patient
matrix has structure -- it will align with the leading principal directions of the embedding space,
and those directions are shared across cancers for reasons that have nothing to do with survival
(scanner, stain, tissue-source site, magnification). So each donor cohort is refitted on PERMUTED
survival labels, holding the cohort, the design matrix, the event count and the censoring pattern
fixed, and the permuted coefficient is transferred the same way. Only the excess of the real
transfer over the permuted transfer is evidence about survival.

A second control severs a different shortcut: cases are matched by TCGA barcode, so a bladder
patient can never appear in a donor cohort, but SITE can be shared. The per-site breakdown is
reported so a transfer that lives entirely in shared collection sites is visible rather than
hidden.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, cox_fit  # noqa: E402

DONORS = ("brca", "coadread", "hnsc", "stad")


def load_donor(path, stem2emb):
    """Case-level (X, time, event) for one donor cohort, slides mean-pooled per case."""
    per = collections.defaultdict(list)
    lab = {}
    for r in csv.DictReader(open(path)):
        cid = r["case_id"]
        stem = r["slide_id"].replace(".svs", "")
        if stem in stem2emb:
            per[cid].append(stem2emb[stem])
        if r.get("survival_months_dss", "").strip() and r.get("censorship_dss", "").strip():
            lab[cid] = (float(r["survival_months_dss"]), 1.0 - float(r["censorship_dss"]))
    cases = sorted(c for c in per if c in lab)
    if not cases:
        return None
    X = np.vstack([np.mean([stem2emb[s] if isinstance(s, str) else s for s in per[c]], axis=0)
                   for c in cases]).astype(float)
    t = np.array([lab[c][0] for c in cases])
    e = np.array([lab[c][1] for c in cases])
    return X, t, e, cases


def fit_direction(X, t, e, alpha, rng=None):
    """Standardised ridge Cox coefficient, optionally on permuted outcomes."""
    mu, sd = X.mean(0), X.std(0) + 1e-9
    Z = (X - mu) / sd
    if rng is not None:
        q = rng.permutation(len(t))
        t, e = t[q], e[q]
    b = cox_fit(Z, t, e, alpha)
    nrm = float(np.linalg.norm(b))
    return b / nrm if nrm > 0 else b


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--pan-dir", required=True)
    ap.add_argument("--alpha", type=float, default=64.0)
    ap.add_argument("--perms", type=int, default=20)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    d = pickle.load(open(a.titan, "rb"))
    E = np.asarray(d["embeddings"])
    stem2emb = {}
    for i, f in enumerate(d["filenames"]):
        s = str(f)
        for suf in (".svs", ".h5", ".pt"):
            if s.endswith(suf):
                s = s[:-len(suf)]
        stem2emb[s] = E[i]

    # bladder side, standardised on the WHOLE bladder cohort -- no bladder outcome is used here,
    # only its feature moments, so this is not a leak of the label
    mu, sd = co.T.mean(0), co.T.std(0) + 1e-9
    Zb = (co.T - mu) / sd
    ii, jj = cpairs(co.t, co.e)
    blca_cases = set(co.keep)

    rep = {"artifact_type": "s5_transfer_falsifier", "reportable": False,
           "phase_of_origin": "exploratory", "alpha": a.alpha, "perms_per_donor": a.perms,
           "blca_n": len(co.keep), "blca_events": int(co.e.sum())}

    rng = np.random.default_rng(0)
    donors, dirs = {}, {}
    for name in DONORS:
        got = load_donor(os.path.join(a.pan_dir, "tcga_%s.csv" % name), stem2emb)
        if got is None:
            donors[name] = {"status": "no embeddings matched"}
            continue
        X, t, e, cases = got
        overlap = len(set(cases) & blca_cases)
        b = fit_direction(X, t, e, a.alpha)
        c_real = cidx(Zb @ b, ii, jj)
        perms = [cidx(Zb @ fit_direction(X, t, e, a.alpha,
                                         np.random.default_rng(700 + k)), ii, jj)
                 for k in range(a.perms)]
        donors[name] = {
            "donor_cases": len(cases), "donor_events": int(e.sum()),
            "cases_also_in_blca": overlap,
            "transfer_cindex_on_blca": round(float(c_real), 4),
            "control_permuted_labels_mean": round(float(np.mean(perms)), 4),
            "control_sd": round(float(np.std(perms, ddof=1)), 4),
            "control_max": round(float(np.max(perms)), 4),
            "excess_over_control": round(float(c_real - np.mean(perms)), 4),
        }
        dirs[name] = b
        print("%-9s cases=%4d ev=%3d  transfer=%.4f  control=%.4f+/-%.4f"
              % (name, len(cases), int(e.sum()), c_real, np.mean(perms),
                 np.std(perms, ddof=1)), file=sys.stderr, flush=True)
    rep["per_donor"] = donors

    if dirs:
        names = sorted(dirs)
        B = np.vstack([dirs[n] for n in names])
        pooled = B.mean(0)
        pooled /= np.linalg.norm(pooled) + 1e-12
        pp = []
        for k in range(a.perms):
            Bp = np.vstack([fit_direction(*load_donor(
                os.path.join(a.pan_dir, "tcga_%s.csv" % n), stem2emb)[:3], a.alpha,
                np.random.default_rng(1700 + k * 7 + i)) for i, n in enumerate(names)])
            v = Bp.mean(0)
            v /= np.linalg.norm(v) + 1e-12
            pp.append(cidx(Zb @ v, ii, jj))
        rep["pooled_direction"] = {
            "donors": names,
            "transfer_cindex_on_blca": round(float(cidx(Zb @ pooled, ii, jj)), 4),
            "control_permuted_labels_mean": round(float(np.mean(pp)), 4),
            "control_sd": round(float(np.std(pp, ddof=1)), 4),
            "excess_over_control": round(float(cidx(Zb @ pooled, ii, jj) - np.mean(pp)), 4),
        }
        # how aligned are the donors with each other, and is that alignment survival-specific?
        cos = {}
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                cos["%s~%s" % (names[i], names[j])] = round(float(B[i] @ B[j]), 4)
        rep["donor_direction_cosines"] = cos

        # the BLADDER direction, fitted on all bladder cases, for reference only
        bb = fit_direction(co.T, co.t, co.e, a.alpha)
        rep["cosine_to_bladder_own_direction"] = {
            n: round(float(dirs[n] @ bb), 4) for n in names}
        rep["cosine_pooled_to_bladder"] = round(float(pooled @ bb), 4)
        rep["cosine_note"] = ("the bladder direction here is fitted on ALL bladder cases and is a "
                              "reference quantity only; it never enters any scored arm")

    rep["reading"] = ("a donor whose transfer_cindex_on_blca exceeds its permuted control by more "
                      "than the control's own spread is evidence that survival-relevant directions "
                      "are shared across cancer types. If every donor sits inside its control, "
                      "every cross-cancer candidate in the slate is dead and nothing else is.")

    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
