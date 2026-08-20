#!/usr/bin/env python3
"""S5 Step 1. The error per slice, and the one probe that decides what to build.

Pre-freeze, exploratory. NOT reportable.

THE AGGREGATE THAT SENT US HERE. Pooled over all 359 cases: clinical alone 0.6754, TITAN alone
0.6511, the two rank-averaged 0.7068. So the entire measured contribution of a whole-slide image,
through a 768-dimensional foundation-model embedding, is +0.031 over what T, N, M, stage, age and
sex already say -- and +0.031 is the size of this task's split-reseed noise. One number cannot say
whether that is because the image carries nothing extra, or because it carries something extra in
a subpopulation and is redundant everywhere else. Those have opposite next moves.

WHAT IS MEASURED.

  1. PER-SLICE C-INDEX for four arms, over slices that could differ mechanistically: collection
     site, T stage, node status, metastasis, sex, age tertile, slides per case, and the tertile of
     the clinical model's own out-of-fold risk. The reported quantity per slice is not just each
     arm's score but IMAGE GAIN = C(clinical + image) - C(clinical), because that is the thing the
     method has to make bigger.

  2. THE CONDITIONAL PROBE, which is the point of this file. A C-index is an average over
     comparable pairs, and a pair whose two patients differ in stage is a pair the clinical model
     already orders correctly. Restrict the comparable pairs to those the clinical score CANNOT
     separate -- |z_clin_i - z_clin_j| below a quantile of the observed gaps -- and score each arm
     on that restricted set.

        image at ~0.5 on clinically-tied pairs -> the image is a stage detector, and every fusion
                                                  rule in the world will not create information
                                                  that is not there. The next component must
                                                  change the inputs or the target.
        image clearly above 0.5 there          -> the information IS present and orthogonal, and
                                                  the deficit is that mean-pooled slide embeddings
                                                  plus a linear head cannot get at it. The next
                                                  component is representation or estimator.

     This is the whole fork in the campaign, and it costs one pass over the pairs.

  3. VARIATION BEFORE BELIEF. A slice a rule points at is not a bottleneck yet. For every slice
     the deviation of its image gain from the overall image gain is reported next to the
     FOLD-TO-FOLD standard deviation of that same gain, so a slice selected at 0.4 of a standard
     deviation cannot be written up as the bottleneck. A deficit that turns out to be UNIFORM is
     the more useful finding: it rules out, in one measurement, every component that special-cases
     a subpopulation.
"""

from __future__ import annotations

import argparse
import collections
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, fitapply  # noqa: E402

MIN_PAIRS = 100     # a slice thinner than this is reported but never ranked
MIN_EVENTS = 8


def spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    d = float(np.sqrt((ra @ ra) * (rb @ rb)))
    return float(ra @ rb / d) if d else float("nan")


def build_slices(co, z_clin):
    """Every slice as (group_name, per-case label). Labels are strings; '' means excluded."""
    def tert(v, names):
        q = np.quantile(v, [1 / 3, 2 / 3])
        return np.array([names[0] if x <= q[0] else (names[1] if x <= q[1] else names[2])
                         for x in v])

    def g(c, k):
        x = (co.tnm.get(c) or {}).get(k)
        return "" if x in (None, "NX", "MX", "TX") else str(x)

    tstage = np.array([g(c, "ajcc_pathologic_t")[:2] for c in co.keep])
    tgrp = np.array([{"T0": "T0-T2", "T1": "T0-T2", "T2": "T0-T2",
                      "T3": "T3", "T4": "T4"}.get(x, "") for x in tstage])
    nstat = np.array(["N0" if g(c, "ajcc_pathologic_n") == "N0"
                      else ("N+" if g(c, "ajcc_pathologic_n") else "") for c in co.keep])
    mstat = np.array(["M0" if g(c, "ajcc_pathologic_m") == "M0"
                      else ("M1" if g(c, "ajcc_pathologic_m") else "") for c in co.keep])

    cnt = collections.Counter(co.site)
    site = np.array([s if cnt[s] >= 20 else "other" for s in co.site])

    return {
        "collection_site": site,
        "T_stage": tgrp,
        "node_status": nstat,
        "metastasis": mstat,
        "sex": np.array(["female" if x else "male" for x in co.fem]),
        "age_tertile": tert(co.age, ("age_low", "age_mid", "age_high")),
        "slides_per_case": np.array(["one_slide" if x == 1 else "multi_slide"
                                     for x in co.n_slides]),
        "clinical_risk_tertile": tert(z_clin, ("clin_low", "clin_mid", "clin_high")),
        "fold": np.array(["fold_%d" % k for k in co.fold]),
    }


def conditional_probe(co, arms, z_clin, qs=(0.05, 0.10, 0.20, 0.40)):
    """Each arm's C-index restricted to comparable pairs the clinical score cannot separate."""
    ii, jj = cpairs(co.t, co.e)
    gap = np.abs(z_clin[ii] - z_clin[jj])
    out = {}
    for q in qs:
        thr = float(np.quantile(gap, q))
        m = gap <= thr
        if m.sum() < MIN_PAIRS:
            continue
        out["clinically_tied_q%02d" % int(q * 100)] = {
            "pairs": int(m.sum()), "pairs_total": int(ii.size),
            "clinical_gap_threshold": round(thr, 4),
            "arms": {k: round(cidx(v, ii[m], jj[m]), 4) for k, v in arms.items()},
        }
    out["all_pairs"] = {"pairs": int(ii.size), "pairs_total": int(ii.size),
                        "clinical_gap_threshold": None,
                        "arms": {k: round(cidx(v, ii, jj), 4) for k, v in arms.items()}}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()

    # ---- out-of-fold component scores on the released folds
    p_ti = np.full(len(co.keep), np.nan)
    p_cl = np.full(len(co.keep), np.nan)
    for tri, vai in fi:
        if len(tri) < 30 or not len(vai):
            continue
        p_ti[vai] = fitapply(co.T[tri], co.t[tri], co.e[tri], co.T[vai], 0)
        p_cl[vai] = fitapply(co.CLIN[tri], co.t[tri], co.e[tri], co.CLIN[vai], 0)

    z_ti, z_cl, z_sp = co.fold_pct(p_ti), co.fold_pct(p_cl), co.fold_pct(co.spr)
    arms = {"clinical": z_cl, "titan": z_ti, "survpath": z_sp,
            "titan+clinical": co.fold_pct(z_ti + z_cl)}

    rep = {"artifact_type": "s5_step1_error_atlas", "reportable": False,
           "phase_of_origin": "exploratory",
           "independent_unit": "one TCGA case", "metric": "Harrell C-index",
           "direction": "higher is better",
           "n": len(co.keep), "events": int(co.e.sum()),
           "min_pairs_to_rank": MIN_PAIRS, "min_events_to_rank": MIN_EVENTS}

    # ---- how much of the image score is just the clinical score
    rep["redundancy"] = {
        "spearman_titan_vs_clinical": round(spearman(z_ti, z_cl), 4),
        "spearman_survpath_vs_clinical": round(spearman(z_sp, z_cl), 4),
        "spearman_titan_vs_survpath": round(spearman(z_ti, z_sp), 4),
        "reading": "a high correlation between the image score and the clinical score is the "
                   "mechanism behind a small image gain; a low one means the two are orthogonal "
                   "and the fusion rule, not the representation, is what is failing",
    }

    # ---- the conditional probe
    rep["conditional_probe"] = conditional_probe(co, arms, z_cl)

    # ---- per-slice atlas
    ii_all, jj_all = cpairs(co.t, co.e)
    overall = {k: cidx(v, ii_all, jj_all) for k, v in arms.items()}
    overall_gain = overall["titan+clinical"] - overall["clinical"]
    rep["overall"] = {k: round(v, 4) for k, v in overall.items()}
    rep["overall_image_gain"] = round(overall_gain, 4)

    slices = build_slices(co, z_cl)
    atlas = {}
    for gname, lab in slices.items():
        rows = {}
        for lv in sorted(set(lab) - {""}):
            m = lab == lv
            i2, j2 = cpairs(co.t[m], co.e[m])
            row = {"n": int(m.sum()), "events": int(co.e[m].sum()), "pairs": int(i2.size),
                   "rankable": bool(i2.size >= MIN_PAIRS and co.e[m].sum() >= MIN_EVENTS)}
            if i2.size:
                for k, v in arms.items():
                    row[k] = round(cidx(v[m], i2, j2), 4)
                row["image_gain"] = round(row["titan+clinical"] - row["clinical"], 4)
                row["gain_minus_overall"] = round(row["image_gain"] - overall_gain, 4)
            rows[lv] = row
        atlas[gname] = rows
    rep["slices"] = atlas

    # ---- variation before belief: fold-to-fold spread of the SAME quantity
    per_fold_gain = []
    for k in range(5):
        m = co.fold == k
        i2, j2 = cpairs(co.t[m], co.e[m])
        if i2.size:
            per_fold_gain.append(cidx(arms["titan+clinical"][m], i2, j2)
                                 - cidx(arms["clinical"][m], i2, j2))
    sd_gain = float(np.std(per_fold_gain, ddof=1))
    rep["gain_variation"] = {
        "per_fold_image_gain": [round(x, 4) for x in per_fold_gain],
        "fold_to_fold_sd": round(sd_gain, 4),
        "rule": "a slice whose gain_minus_overall is under one fold-to-fold SD is not a "
                "bottleneck; it is where the noise happened to land this time",
    }

    ranked = []
    for gname, rows in atlas.items():
        for lv, row in rows.items():
            if row.get("rankable") and "image_gain" in row:
                ranked.append({"slice": "%s=%s" % (gname, lv), "n": row["n"],
                               "image_gain": row["image_gain"],
                               "gain_minus_overall": row["gain_minus_overall"],
                               "sds_from_overall": round(row["gain_minus_overall"] / sd_gain, 2)
                               if sd_gain else None})
    ranked.sort(key=lambda r: r["image_gain"])
    rep["worst_slices_by_image_gain"] = ranked[:8]
    rep["best_slices_by_image_gain"] = ranked[-5:][::-1]
    rep["any_slice_beyond_one_sd"] = [r["slice"] for r in ranked
                                      if r["sds_from_overall"] is not None
                                      and abs(r["sds_from_overall"]) >= 1.0]

    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
