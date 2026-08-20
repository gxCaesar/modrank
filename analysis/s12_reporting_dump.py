#!/usr/bin/env python3
"""Dump the per-case and per-pathway tables the reporting figures need, from the frozen run.

WHY THIS IS NOT AN UNBLINDING EVENT, stated first because it is the thing to get wrong. Nothing
here fits, selects, tunes or compares anything new. It re-executes the frozen arm construction
verbatim and writes out quantities the confirmatory script computed and then summarised away: the
per-case out-of-fold score of every arm, and the survival association of every one of the 275
pathways rather than the fifteen strongest. A reader who wants calibration, per-fold Kaplan-Meier
or the full pathway landscape cannot get them from a summary, and a paper that reports only the
top fifteen of 275 has shown its reader the tail it liked.

THE GUARD. The recomputed seed-averaged concordance must equal the frozen primary to four decimal
places or this script fails. If the environment, the inputs or the code path have drifted since the
freeze, that number moves, and every figure built on this dump would be built on a different run
than the one the paper reports.

Inputs are the same three the confirmatory run read. The interpreter must be the frozen one:
Python 3.8.20 with numpy 1.24.4, scipy 1.10.1, scikit-learn 1.3.2.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, fitapply                    # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix                         # noqa: E402
from s6_amend_clinical import dimaf_clinical                              # noqa: E402
from s8_biology import bh, score_test, spearman                           # noqa: E402

FROZEN_PRIMARY = 0.7212


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("protocol", "root", "titan", "sp-dir", "omics-dir", "dimaf-dir", "out"):
        ap.add_argument("--" + f, required=True)
    a = ap.parse_args()

    proto = json.load(open(a.protocol))
    assert proto["protocol_status"] == "frozen", "the protocol is not frozen"
    seeds = proto["method"]["seeds"]["values"]
    alphas = tuple(proto["method"]["estimator"]["alpha_grid"])

    co = Cohort(a.root, a.titan, a.sp_dir)
    assert len(co.keep) == proto["task"]["n_cases"]
    assert int(co.e.sum()) == proto["task"]["n_events"]
    fi = co.fold_indices()
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, pnames = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    # AMENDMENT A1. The frozen run's clinical block took stage from a cached query that
    # selected the first diagnosis record per case, and 26 cases disagreed with the
    # incumbent's own file. The reported primary is the AMENDED one, built from DIMAF's
    # released split files, so this dump must use the same block or it reproduces a
    # number the paper does not report. The guard below is what caught the difference.
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    CLIN = np.hstack([dc["age"], dc["fem"], dc["stage"]])
    blocks = {"wsi_titan": co.T, "omics_combine": P, "clinical": CLIN}
    ii, jj = cpairs(co.t, co.e)
    n = len(co.keep)

    # --- fold membership per case, which the summary never recorded
    fold_of = np.full(n, -1, dtype=int)
    for k, (tri, vai) in enumerate(fi):
        fold_of[vai] = k
    assert (fold_of >= 0).all(), "every case must sit in exactly one validation fold"

    # --- the frozen arm construction, verbatim
    per_seed, arms_by_seed = [], []
    for sd in seeds:
        arms = {}
        for nm, X in blocks.items():
            s = np.full(n, np.nan)
            for tri, vai in fi:
                if len(tri) < 30 or not len(vai):
                    continue
                s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], sd, alphas)
            arms[nm] = co.fold_pct(s)
        ours = co.fold_pct(sum(arms.values()))
        arms["OURS"] = ours
        arms_by_seed.append(arms)
        per_seed.append(round(cidx(ours, ii, jj), 4))
        print("  seed %d  OURS=%.4f" % (sd, per_seed[-1]), file=sys.stderr, flush=True)

    primary = float(np.mean(per_seed))
    assert abs(primary - FROZEN_PRIMARY) < 5e-5, (
        "recomputed primary %.4f does not match the frozen %.4f; the run has drifted and no "
        "figure may be built on this dump" % (primary, FROZEN_PRIMARY))
    print("  primary %.4f matches the frozen value" % primary, file=sys.stderr)

    mean_arm = {k: np.vstack([a[k] for a in arms_by_seed]).mean(0)
                for k in ("wsi_titan", "omics_combine", "clinical", "OURS")}

    # --- per-case table. Seed 0 is what the paired comparisons use; the seed mean is what the
    #     primary reports. Both are written so a figure never has to choose silently.
    cases = []
    for i in range(n):
        cases.append({
            "case_id": str(co.keep[i]), "fold": int(fold_of[i]),
            "months": round(float(co.t[i]), 3), "event": int(co.e[i]),
            "seed0": {k: round(float(arms_by_seed[0][k][i]), 5)
                      for k in ("wsi_titan", "omics_combine", "clinical", "OURS")},
            "seed_mean": {k: round(float(mean_arm[k][i]), 5) for k in mean_arm}})

    # --- every pathway, not the fifteen strongest
    st = score_test((P - P.mean(0)) / (P.std(0) + 1e-9), co.t, co.e)
    praw = {pnames[k]: float(2 * (1 - 0.5 * (1 + math.erf(abs(st[k]) / math.sqrt(2)))))
            for k in range(len(pnames))}
    q = bh(praw)
    slide0, om0 = arms_by_seed[0]["wsi_titan"], arms_by_seed[0]["omics_combine"]
    # score_test returns |U|/sqrt(V) and is UNSIGNED by construction. The project's convention,
    # set in the biology figure, recovers direction from whether the variable scores above or
    # below chance as a raw risk score, so the same convention is used here rather than a second
    # one invented for this dump.
    paths = [{"pathway": pnames[k], "cox_score_abs": round(float(st[k]), 4),
              "cindex_raw": round(float(cidx(co.fold_pct(P[:, k]), ii, jj)), 4),
              "p": float("%.4g" % praw[pnames[k]]), "q_BH": round(float(q[pnames[k]]), 4),
              "rho_slide": round(float(spearman(slide0, P[:, k])), 4),
              "rho_omics": round(float(spearman(om0, P[:, k])), 4)}
             for k in range(len(pnames))]
    assert len(paths) == 275, "the benchmark's own grouping is 275 pathways; got %d" % len(paths)

    # --- the transcriptome and the clinical block as DATA, so a figure can show what the model
    #     reads rather than only what it concluded. The pathway matrix is z-scored per pathway
    #     across the cohort, which is a display transform and changes no reported number.
    Z = (P - P.mean(0)) / (P.std(0) + 1e-9)
    order_by_assoc = np.argsort(-np.abs(st))
    TOPK = 40
    keep = order_by_assoc[:TOPK]
    stage_lab = dc.get("stage_label")
    if stage_lab is None:
        # the one-hot columns carry the level; recover the label from the column that is set
        lv = dc.get("stage_levels") or ["Stage I", "Stage II", "Stage III", "Stage IV"]
        oh = dc["stage"]
        stage_lab = [lv[int(np.argmax(oh[i]))] if oh[i].max() > 0 else "unknown"
                     for i in range(len(co.keep))]
    for i, c in enumerate(cases):
        c["stage"] = str(stage_lab[i])
        c["age"] = round(float(dc["age"][i][0]), 1)
        c["female"] = int(dc["fem"][i][0])
        c["pathway_z"] = [round(float(Z[i, k]), 3) for k in keep]

    rep = {"artifact_type": "s12_reporting_dump",
           "pathway_z_columns": [pnames[k] for k in keep],
           "pathway_z_note": "the %d pathways with the strongest univariate association, z-scored "
                             "across the cohort for display; the arm itself uses all 275" % TOPK,
           "phase_of_origin": "frozen_post_hoc_reporting",
           "reportable": "descriptive only; nothing here fits, selects or compares anything new",
           "n": n, "events": int(co.e.sum()), "comparable_pairs": int(ii.size),
           "seeds": list(seeds), "per_seed_OURS": per_seed, "primary": round(primary, 4),
           "frozen_primary_matched": True,
           "fold_sizes": [int((fold_of == k).sum()) for k in range(len(fi))],
           "cases": cases, "pathways": paths}
    with open(a.out, "w") as fh:
        json.dump(rep, fh, indent=1)
    print("  wrote %s  (%d cases, %d pathways)" % (a.out, len(cases), len(paths)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
