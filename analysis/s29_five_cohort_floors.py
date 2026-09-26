#!/usr/bin/env python3
"""Each of the benchmark's five studies, given its own noise floor.

WHY. The five-study validation reads each margin against its bootstrap interval and nothing else.
The bladder result is held to a second and stricter standard: a margin must also exceed twice the
standard deviation of the same paired difference under re-partitioning, because an interval says how
far the estimate would move if these patients were resampled while a re-partition says how far it
would move if somebody re-ran the study. The other four studies had no such number, so a reviewer
asking "how much of that margin is partition luck" had no answer for them.

WHAT IT MEASURES, per study. Twenty-four re-partitions built the way the bladder floor was built
(`s5_step0_headroom.reseeded_folds`, rng 1000 + r, estimator seed 0), each one refitting the slide,
transcriptome and clinical arms from scratch and recombining them by the same equal-weight rank
rule. The floor is the standard deviation across those 24 of the PAIRED difference
C(ModRank) - C(clinical arm), and the bar is twice it. The same 24 also give the spread of ModRank's
own concordance, which is reported beside it because the two answer different questions.

WHAT THIS NUMBER IS NOT. The bladder floor in the main text is 0.0073, measured with the incumbent's
clinical block and five seeds. This run uses the benchmark's released clinical file and one seed,
which is the five-study protocol. The bladder row here is therefore NOT 0.0073 and is not meant to
be: it is the same quantity measured under the other protocol, and having it lets the four
non-bladder floors be read as comparable to something.

KNOWN ANSWERS. Partition 0 is the benchmark's own released folds, so before any re-partition is
drawn the run must reproduce the committed point estimates for all five studies. It refuses to write
if it does not.

Post-freeze, exploratory. About forty minutes on one CPU host, dominated by the 768-dimensional
slide arm on the 871-patient breast study.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from blca_common import cidx, cpairs                                          # noqa: E402
from s5_step0_headroom import reseeded_folds                                  # noqa: E402
from s5_stage_five_cohorts import COHORTS, arm, fold_pct, load_cohort         # noqa: E402

RESEEDS = 24


def repartition(co, folds):
    """A copy of the cohort dict carrying a different partition, so `arm` can be reused verbatim."""
    assign = np.full(len(co["keep"]), -1)
    fi = []
    for k, (tri, vai) in enumerate(folds):
        fi.append((np.asarray(tri), np.asarray(vai)))
        assign[np.asarray(vai)] = k
    d = dict(co)
    d["fi"], d["assign"] = fi, assign
    return d


def score(co, ii, jj):
    """ModRank and the clinical arm on whatever partition `co` carries."""
    AS = np.hstack([co["age"], co["fem"]])
    ASS = np.hstack([AS, co["stage"]])
    z_ass, z_tit, z_omi = arm(ASS, co), arm(co["T"], co), arm(co["P"], co)
    if z_ass is None:
        return None, None
    ours = fold_pct(z_tit + z_omi + z_ass, co["assign"])
    return float(cidx(ours, ii, jj)), float(cidx(z_ass, ii, jj))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("titan", "meta-dir", "clin-dir", "pan-dir", "sig", "frozen", "out"):
        ap.add_argument("--" + f, required=True)
    ap.add_argument("--reseeds", type=int, default=RESEEDS)
    a = ap.parse_args()
    t0 = time.time()

    frozen = json.load(open(a.frozen))["cohorts"]
    d0 = pickle.load(open(a.titan, "rb"))
    E = np.asarray(d0["embeddings"])
    stem2emb = {}
    for i, f in enumerate(d0["filenames"]):
        s = str(f)
        for suf in (".svs", ".h5", ".pt"):
            if s.endswith(suf):
                s = s[:-len(suf)]
        stem2emb[s] = E[i]

    out = {"artifact_type": "s29_five_cohort_noise_floors",
           "phase_of_origin": "post_freeze_2026-09-12",
           "construction": "s5_step0_headroom.reseeded_folds, rng 1000 + r, estimator seed 0, the "
                           "five-study protocol (benchmark clinical file, one seed)",
           "floor_is": "the standard deviation across re-partitions of C(ModRank) - C(clinical arm)",
           "not_the_same_as": "the main text's bladder floor of 0.0073, which uses the incumbent's "
                              "clinical block and five seeds. The bladder row here measures the same "
                              "quantity under the other protocol and is not expected to equal it",
           "reseeds": a.reseeds, "cohorts": {}}

    drift = {}
    for c in COHORTS:
        paths = (os.path.join(a.meta_dir, "tcga_%s.csv" % c),
                 os.path.join(a.clin_dir, "tcga_%s_clinical.csv" % c),
                 os.path.join(a.pan_dir, "splits", c),
                 os.path.join(a.pan_dir, "rna_%s.csv" % c))
        if not all(os.path.exists(p) for p in paths):
            out["cohorts"][c] = {"status": "inputs missing"}
            continue
        co = load_cohort(paths[0], paths[1], paths[2], paths[3], a.sig, stem2emb)
        ii, jj = cpairs(co["t"], co["e"])

        # partition 0: the released folds, which must reproduce the committed values
        ours0, clin0 = score(co, ii, jj)
        want = frozen.get(c, {}).get("OURS_wsi_omics_age_sex_stage")
        if want is not None and abs(round(ours0, 4) - want) > 5e-5:
            drift[c] = {"got": round(ours0, 4), "want": want}

        deltas, ours_v = [], []
        n = len(co["keep"])
        for r in range(a.reseeds):
            co_r = repartition(co, reseeded_folds(n, np.random.default_rng(1000 + r)))
            o, cl = score(co_r, ii, jj)
            deltas.append(o - cl)
            ours_v.append(o)
            if (r + 1) % 8 == 0:
                print("  %-9s re-partition %2d/%d  [%.0f s]" % (c, r + 1, a.reseeds, time.time() - t0),
                      file=sys.stderr, flush=True)
        d = np.asarray(deltas)
        margin = ours0 - clin0
        out["cohorts"][c] = {
            "n": n, "events": int(co["e"].sum()),
            "released_folds": {"ModRank": round(ours0, 4), "clinical": round(clin0, 4),
                               "margin": round(margin, 4)},
            "floor": {"sd_of_paired_delta": round(float(d.std(ddof=1)), 4),
                      "bar_2sd": round(2 * float(d.std(ddof=1)), 4),
                      "mean_delta": round(float(d.mean()), 4),
                      "min_delta": round(float(d.min()), 4), "max_delta": round(float(d.max()), 4),
                      "partitions_with_a_positive_delta": int((d > 0).sum())},
            "ModRank_spread": {"sd": round(float(np.std(ours_v, ddof=1)), 4),
                               "min": round(float(np.min(ours_v)), 4),
                               "max": round(float(np.max(ours_v)), 4)},
            "margin_clears_its_own_bar": bool(margin > 2 * float(d.std(ddof=1)))}
        row = out["cohorts"][c]
        print("%-9s n=%4d ev=%3d | margin %+.4f | floor %.4f bar %.4f | clears %s | %d/%d partitions positive  [%.0f s]"
              % (c, n, row["events"], margin, row["floor"]["sd_of_paired_delta"],
                 row["floor"]["bar_2sd"], row["margin_clears_its_own_bar"],
                 row["floor"]["partitions_with_a_positive_delta"], a.reseeds, time.time() - t0),
              file=sys.stderr, flush=True)

    scored = [v for v in out["cohorts"].values() if v.get("floor")]
    out["summary"] = {
        "cohorts_scored": len(scored),
        "cohorts_whose_margin_clears_twice_their_own_re-partition_sd": [
            k for k, v in out["cohorts"].items() if v.get("margin_clears_its_own_bar")],
        "reading": "a margin below its own bar is partition luck at this size and is reported as "
                   "such, whatever its bootstrap interval says",
    }
    out["runtime_seconds"] = round(time.time() - t0, 1)
    if drift:
        out["known_answer_drift"] = drift
        print("KNOWN ANSWER DRIFT, refusing to write: %s" % drift, file=sys.stderr)
        return 2
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1)
    print("wrote %s" % a.out, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
