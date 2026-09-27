#!/usr/bin/env python3
"""Does the model predict equally well for the groups it can be checked against?

WHY THIS EXISTS. The TRIPOD+AI checklist has two rows this work answered with "Not addressed":
item 3c, health inequalities, and item 14, fairness. For a clinical prediction model at a general
journal those are rows a referee reads. They were answerable all along: the released metadata
carries sex and age for every patient, and the canonical per-case scores are already committed, so
the analysis needs no refitting of anything.

WHAT IT CANNOT ANSWER, and this is stated in the output rather than left for a reader to discover.
The frozen inputs carry age, sex and contributing site. They do NOT carry race, ethnicity,
insurance status or any measure of socioeconomic position, so this is a subgroup-performance
analysis on the axes the data supports and not a full fairness audit. Saying "fairness: addressed"
on the strength of two axes would be worse than the honest "not addressed" it replaces.

WHAT IT MEASURES, per stratum, from the committed seed-averaged out-of-fold scores:

  the concordance of each arm inside the stratum, with the patients and events that entered
  ModRank minus the corrected clinical arm inside the stratum, paired at the patient level
  ModRank's concordance in one stratum minus the other, which is the quantity a fairness question
    actually asks: does the model serve one group better than another

The paired bootstrap is the one used everywhere else in this work: 6,000 resamples of patients,
comparable pairs rebuilt inside each replicate, seed 0. A stratum is reported only if it clears the
floors the stage-subgroup analysis used, 25 patients and 100 comparable pairs.

KNOWN ANSWERS. Pooled over all 359 patients the vectors must give the committed canonical values,
0.7237 for ModRank and 0.6644 for the clinical arm, or the file is not the one the paper reports
and nothing below means anything.

Post-freeze, exploratory. No model is refitted; this reads scores that already exist.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from blca_common import cidx, cpairs                                         # noqa: E402

REPS = 6000
MIN_CASES, MIN_PAIRS = 25, 100
KNOWN = {"ours": 0.7237, "clinical_stage": 0.6644}
ARMS = ("clinical_stage", "slide", "omics", "ours")


def boot(x, y, t, e, reps=REPS, seed=0):
    """Identical to s6_confirmatory.boot. Cases are resampled and the pairs are rebuilt inside each
    replicate; reusing the original pairs would make every interval too narrow."""
    rng = np.random.default_rng(seed)
    n = len(t)
    v = []
    for _ in range(reps):
        bs = rng.choice(n, size=n, replace=True)
        ii, jj = cpairs(t[bs], e[bs])
        if ii.size:
            v.append(cidx(x[bs], ii, jj) - cidx(y[bs], ii, jj))
    v = np.asarray(v)
    return {"mean": round(float(v.mean()), 4), "se": round(float(v.std(ddof=1)), 4),
            "ci95": [round(float(np.quantile(v, .025)), 4),
                     round(float(np.quantile(v, .975)), 4)],
            "p_two_sided": round(float(2 * min((v <= 0).mean(), (v >= 0).mean())), 4),
            "reps": int(v.size), "resampled": "cases, not comparable pairs"}


def boot_between(xa, ta, ea, xb, tb, eb, reps=REPS, seed=0):
    """One arm's concordance in stratum A minus the same arm's concordance in stratum B. The two
    strata are disjoint patients, so each is resampled independently: there is nothing to pair."""
    rng = np.random.default_rng(seed)
    na, nb = len(ta), len(tb)
    v = []
    for _ in range(reps):
        a = rng.choice(na, size=na, replace=True)
        b = rng.choice(nb, size=nb, replace=True)
        ia, ja = cpairs(ta[a], ea[a])
        ib, jb = cpairs(tb[b], eb[b])
        if ia.size and ib.size:
            v.append(cidx(xa[a], ia, ja) - cidx(xb[b], ib, jb))
    v = np.asarray(v)
    return {"mean": round(float(v.mean()), 4), "se": round(float(v.std(ddof=1)), 4),
            "ci95": [round(float(np.quantile(v, .025)), 4),
                     round(float(np.quantile(v, .975)), 4)],
            "p_two_sided": round(float(2 * min((v <= 0).mean(), (v >= 0).mean())), 4),
            "reps": int(v.size), "resampled": "the two strata independently, being disjoint"}


def boot_gap_difference(va, ta, ea, vb, tb, eb, arm1, arm2, reps=REPS, seed=0):
    """(arm1 in A minus arm1 in B) minus (arm2 in A minus arm2 in B), on the same independent
    resamples of the two disjoint strata that boot_between draws: whether one arm depends on the
    stratum more than another."""
    rng = np.random.default_rng(seed)
    na, nb = len(ta), len(tb)
    v = []
    for _ in range(reps):
        a = rng.choice(na, size=na, replace=True)
        b = rng.choice(nb, size=nb, replace=True)
        ia, ja = cpairs(ta[a], ea[a])
        ib, jb = cpairs(tb[b], eb[b])
        if ia.size and ib.size:
            v.append((cidx(va[arm1][a], ia, ja) - cidx(vb[arm1][b], ib, jb))
                     - (cidx(va[arm2][a], ia, ja) - cidx(vb[arm2][b], ib, jb)))
    v = np.asarray(v)
    return {"mean": round(float(v.mean()), 4), "se": round(float(v.std(ddof=1)), 4),
            "ci95": [round(float(np.quantile(v, .025)), 4),
                     round(float(np.quantile(v, .975)), 4)],
            "p_two_sided": round(float(2 * min((v <= 0).mean(), (v >= 0).mean())), 4),
            "reps": int(v.size), "resampled": "the two strata independently, being disjoint"}


def stratum(rows, name):
    t = np.array([r["months"] for r in rows], dtype=float)
    e = np.array([r["event"] for r in rows], dtype=float)
    ii, jj = cpairs(t, e)
    out = {"n": len(rows), "events": int(e.sum()), "comparable_pairs": int(ii.size)}
    if len(rows) < MIN_CASES or ii.size < MIN_PAIRS:
        out["skipped"] = "below the reporting floor of %d patients and %d pairs" % (MIN_CASES, MIN_PAIRS)
        return out, None
    vec = {a: np.array([r[a] for r in rows], dtype=float) for a in ARMS}
    for a in ARMS:
        out[a] = round(cidx(vec[a], ii, jj), 4)
    out["ours_minus_clinical"] = boot(vec["ours"], vec["clinical_stage"], t, e)
    return out, (vec, t, e)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("vectors", "meta", "out"):
        ap.add_argument("--" + f, required=True)
    ap.add_argument("--race", default="", help="optional patient-level race map from s30")
    a = ap.parse_args()

    cases = json.load(open(a.vectors))["cases"]
    meta = {r["case_id"]: r for r in csv.DictReader(open(a.meta))}
    race = json.load(open(a.race)) if a.race and os.path.exists(a.race) else {}
    for c in cases:
        m = meta.get(c["case_id"])
        if m is None:
            print("no metadata row for %s" % c["case_id"], file=sys.stderr)
            return 2
        c["is_female"] = int(float(m["is_female"]))
        c["age"] = float(m["age"])
        c["race"] = (race.get(c["case_id"], {}) or {}).get("race") or "not recorded"

    # known answers first: these vectors must be the ones the paper reports
    t = np.array([c["months"] for c in cases], dtype=float)
    e = np.array([c["event"] for c in cases], dtype=float)
    ii, jj = cpairs(t, e)
    drift = {}
    import blca_common
    known = KNOWN
    if blca_common.A2:            # amendment A2: the canonical values come from the regenerated s19 run
        _u = json.load(open(os.path.join(os.path.dirname(os.path.abspath(a.vectors)),
                                         "unified-fusion-and-added-value.json")))
        known = {"ours": _u["known_answers"]["ours_canonical"]["recomputed"],
                 "clinical_stage": _u["added_value_inflation"]["constructions"]["ours"]["C_stage"]}
    for arm, want in known.items():
        got = round(cidx(np.array([c[arm] for c in cases], dtype=float), ii, jj), 4)
        if abs(got - want) > 5e-5:
            drift[arm] = {"got": got, "want": want}

    med = float(np.median([c["age"] for c in cases]))
    SPLITS = {
        "sex": {"men": [c for c in cases if c["is_female"] == 0],
                "women": [c for c in cases if c["is_female"] == 1]},
        "age": {"at or below the median (%.0f years)" % med: [c for c in cases if c["age"] <= med],
                "above the median": [c for c in cases if c["age"] > med]},
    }
    if race:
        # Every recorded category, including the ones that will not clear the floors. A race axis
        # that silently drops its small groups reports the largest group's concordance and calls it
        # a fairness result; the counts are the finding here, so they all appear.
        groups = sorted({c["race"] for c in cases})
        SPLITS["race"] = {g: [c for c in cases if c["race"] == g] for g in groups}

    out = {"artifact_type": "s28_subgroup_fairness",
           "phase_of_origin": "post_freeze_2026-09-12",
           "answers": "TRIPOD+AI items 3c (health inequalities) and 14 (fairness)",
           "refits_nothing": "reads the committed seed-averaged out-of-fold scores",
           "axes_available": ["sex", "age", "contributing site (reported separately as grouped CV)"],
           "axes_not_available": "race, ethnicity, insurance status and socioeconomic position are "
                                 "absent from the released files, so this is subgroup performance "
                                 "on the axes the data supports, not a full fairness audit",
           "bootstrap": {"replicates": REPS, "unit": "cases", "seed": 0},
           "floors": {"min_cases": MIN_CASES, "min_pairs": MIN_PAIRS},
           "pooled": {arm: round(cidx(np.array([c[arm] for c in cases], dtype=float), ii, jj), 4)
                      for arm in ARMS},
           "splits": {}}

    for axis, groups in SPLITS.items():
        rec, keep = {}, {}
        for label, rows in groups.items():
            rec[label], vec = stratum(rows, label)
            if vec is not None:
                keep[label] = vec
        if len(keep) == 2:
            (la, (va, ta, ea)), (lb, (vb, tb, eb)) = list(keep.items())
            rec["ours_%s_minus_%s" % (la.split()[0], lb.split()[0])] = boot_between(
                va["ours"], ta, ea, vb["ours"], tb, eb)
            # the clinical arm's stratum gap minus ModRank's (SI Note 7), emitted here since
            # amendment A2; it had been computed once and quoted from the README
            rec["clinical_gap_minus_ours_gap"] = boot_gap_difference(
                va, ta, ea, vb, tb, eb, "clinical_stage", "ours")
        out["splits"][axis] = rec
        for label in groups:
            r = rec[label]
            if "skipped" in r:
                print("%-8s %-34s n=%3d ev=%3d  SKIPPED" % (axis, label, r["n"], r["events"]),
                      file=sys.stderr)
            else:
                print("%-8s %-34s n=%3d ev=%3d | clinical %.4f  ours %.4f | ours-clinical %+.4f %s"
                      % (axis, label, r["n"], r["events"], r["clinical_stage"], r["ours"],
                         r["ours_minus_clinical"]["mean"], r["ours_minus_clinical"]["ci95"]),
                      file=sys.stderr)

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
