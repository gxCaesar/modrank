#!/usr/bin/env python3
"""C1' -- the pan-cancer direction as an ADDITIONAL arm, not as a constraint on the bladder head.

Pre-freeze, exploratory. NOT reportable.

WHY THIS IS NOT C1 RELABELLED. C1 was killed an hour ago: constraining the bladder head to a
donor-estimated subspace gave 0.6554 against a plain-ridge baseline of 0.6596, missing its
pre-registered bar of +0.0145 and losing outright, while a matched control on permuted donor
labels reached 0.6455 and an unsupervised-subspace ablation reached 0.6638. Under the slate's
shape rule a new candidate inherits that kill unless it can name the mechanical difference that
makes it not apply. Here it is, in one sentence: C1 REPLACED the bladder head's own direction with
the donor subspace, and C1' KEEPS BOTH and combines them at the score level.

The measurement that says this is the right correction rather than a rationalisation was already
in C1's own output. The donor-pooled direction transfers to bladder at 0.6492, the bladder head
reaches 0.6596, and the cosine between them is 0.1748 -- two nearly orthogonal directions of
nearly equal strength. Projecting one onto the other throws away whichever is not kept. That is
what C1 did.

WHAT IS DIFFERENT ABOUT THIS ARM, and it is the part worth a paper if it survives: it requires NO
BLADDER TRAINING DATA AT ALL. The coefficients come from breast, colorectal, head-and-neck and
stomach cohorts; the bladder patients supply only their embeddings. Every other arm on this
benchmark is fitted on bladder.

DISCIPLINE. The donor fits never touch a bladder patient, so the transfer score itself is
fold-independent -- but the standardisation of the bladder features is not, and it is computed
from TRAINING-FOLD moments only. Percentile normalisation is applied within fold to every arm
equally, as everywhere else in this campaign.

FALSIFIER, FIXED BEFORE THE RUN. Adding the zero-shot arm to the bladder image arm must improve on
the bladder image arm alone by more than 0.0145, the 2-SD paired reseed bar measured on this code.
And its matched control -- the same construction with every donor's survival labels permuted --
must not move.
"""

from __future__ import annotations

import argparse
import collections
import csv
import itertools
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, cox_fit, fitapply  # noqa: E402
from s5_c1_pancancer_prior import bootstrap_directions, load_donor  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402

DONORS = ("brca", "coadread", "hnsc", "stad")


def donor_direction(pan_dir, stem2emb, permute, boots, seed):
    """One unit direction per donor cohort, averaged over bootstraps, then averaged over donors."""
    vs = []
    for i, nm in enumerate(DONORS):
        got = load_donor(os.path.join(pan_dir, "tcga_%s.csv" % nm), stem2emb)
        if got is None:
            continue
        X, t, e, _ = got
        B = bootstrap_directions(X, t, e, boots, np.random.default_rng(seed + i), 64.0,
                                 permute=permute)
        if B.shape[0]:
            v = B.mean(0)
            n = float(np.linalg.norm(v))
            if n > 0:
                vs.append(v / n)
    v = np.mean(vs, axis=0)
    return v / (np.linalg.norm(v) + 1e-12)


def foldwise_transfer(co, fi, direction):
    """Score bladder patients with a fixed direction, standardising by TRAINING-fold moments."""
    s = np.full(len(co.keep), np.nan)
    for tri, vai in fi:
        if len(tri) < 30 or not len(vai):
            continue
        mu, sd = co.T[tri].mean(0), co.T[tri].std(0) + 1e-9
        s[vai] = ((co.T[vai] - mu) / sd) @ direction
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--pan-dir", required=True)
    ap.add_argument("--rna", default="")
    ap.add_argument("--sig", default="")
    ap.add_argument("--boots", type=int, default=40)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    d = pickle.load(open(a.titan, "rb"))
    E = np.asarray(d["embeddings"])
    stem2emb = {}
    for i, f in enumerate(d["filenames"]):
        s = str(f)
        for suf in (".svs", ".h5", ".pt"):
            if s.endswith(suf):
                s = s[:-len(suf)]
        stem2emb[s] = E[i]

    ii, jj = cpairs(co.t, co.e)
    same_site = co.site[ii] == co.site[jj]
    gap_ref = None

    arms = {}
    # bladder-trained image arm
    p = np.full(len(co.keep), np.nan)
    pc = np.full(len(co.keep), np.nan)
    for tri, vai in fi:
        if len(tri) < 30 or not len(vai):
            continue
        p[vai] = fitapply(co.T[tri], co.t[tri], co.e[tri], co.T[vai], 0)
        pc[vai] = fitapply(co.CLIN[tri], co.t[tri], co.e[tri], co.CLIN[vai], 0)
    arms["titan_bladder_trained"] = co.fold_pct(p)
    arms["clinical"] = co.fold_pct(pc)
    gap_ref = arms["clinical"]

    # zero-shot pan-cancer arm and its matched control
    v_real = donor_direction(a.pan_dir, stem2emb, False, a.boots, 100)
    v_perm = donor_direction(a.pan_dir, stem2emb, True, a.boots, 200)
    arms["pancancer_zeroshot"] = co.fold_pct(foldwise_transfer(co, fi, v_real))
    arms["pancancer_zeroshot_PERMUTED_CONTROL"] = co.fold_pct(
        foldwise_transfer(co, fi, v_perm))

    if a.rna and a.sig:
        G, genes, _ = load_rna(a.rna, co.keep)
        P, pnames = pathway_matrix(G, genes, a.sig)
        po = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            po[vai] = fitapply(P[tri], co.t[tri], co.e[tri], P[vai], 0)
        arms["omics_pathways"] = co.fold_pct(po)

    arms["survpath_rerun"] = co.fold_pct(co.spr)

    rep = {"artifact_type": "s5_c1prime_zeroshot_arm", "reportable": False,
           "phase_of_origin": "exploratory", "n": len(co.keep), "events": int(co.e.sum()),
           "donor_bootstraps": a.boots,
           "cosine_zeroshot_to_permuted_control": round(float(v_real @ v_perm), 4),
           "falsifier": "titan + zeroshot must beat titan alone by more than 0.0145, and the "
                        "permuted-donor control must not move"}

    gap = np.abs(gap_ref[ii] - gap_ref[jj])
    tied = gap <= float(np.quantile(gap, 0.10))
    rep["single_arms"] = {
        k: {"cindex": round(cidx(v, ii, jj), 4),
            "within_site": round(cidx(v, ii[same_site], jj[same_site]), 4),
            "clinically_tied_q10": round(cidx(v, ii[tied], jj[tied]), 4)}
        for k, v in arms.items()}

    # every rank-average combination of the real arms, so the best is not cherry-picked silently
    real = [k for k in arms if "CONTROL" not in k and k != "survpath_rerun"]
    combos = {}
    for r in range(2, len(real) + 1):
        for cmb in itertools.combinations(sorted(real), r):
            z = co.fold_pct(sum(arms[k] for k in cmb))
            combos[" + ".join(cmb)] = {
                "cindex": round(cidx(z, ii, jj), 4),
                "within_site": round(cidx(z, ii[same_site], jj[same_site]), 4),
                "clinically_tied_q10": round(cidx(z, ii[tied], jj[tied]), 4)}
    rep["all_rank_average_combinations"] = dict(
        sorted(combos.items(), key=lambda kv: -kv[1]["cindex"]))
    rep["n_combinations_scored"] = len(combos)
    rep["selection_warning"] = (
        "%d combinations were scored on the pooled outcome; picking the best of them carries a "
        "selection term of about sd*sqrt(2*ln N). Any arm taken forward must reselect INSIDE the "
        "training fold before it means anything." % len(combos))

    base = cidx(arms["titan_bladder_trained"], ii, jj)
    z2 = co.fold_pct(arms["titan_bladder_trained"] + arms["pancancer_zeroshot"])
    z2c = co.fold_pct(arms["titan_bladder_trained"]
                      + arms["pancancer_zeroshot_PERMUTED_CONTROL"])
    rep["verdict"] = {
        "titan_alone": round(base, 4),
        "titan_plus_zeroshot": round(cidx(z2, ii, jj), 4),
        "gain": round(cidx(z2, ii, jj) - base, 4),
        "titan_plus_PERMUTED_control": round(cidx(z2c, ii, jj), 4),
        "control_gain": round(cidx(z2c, ii, jj) - base, 4),
        "bar": 0.0145,
        "read_the_control_first": "if control_gain is comparable to gain, this run is VOID",
    }

    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
