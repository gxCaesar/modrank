#!/usr/bin/env python3
"""The slide arm on pairs that neither the transcriptome arm nor the clinical model separates.

WHY. The sentence "restricting to pairs that neither the transcriptome nor stage can separate
(3,893 of 24,219), the slide arm still reaches 0.6212", with 0.6396 and 0.5604 on the pairs the
transcriptome arm separates least, has been in the manuscript since August. Found 2026-09-11 while
binding the Nature Communications build's numbers: no committed script or results file produces
any of the four values. They exist only in a markdown audit note. This is the defect the external
review found in the subtype table, and it gets the same repair: a committed emitter, a known-answer
check, and the definition written down.

WHAT THE NOTE ACTUALLY COMPUTED. A first version of this script read "stage" literally (both
patients in the same stage group) and reproduced the first two values but not the last two, under
either stage field. The note's "stage" was the CLINICAL MODEL, age, sex and pathologic stage from
the incumbent's split files (the amended block): pairs inside the 40% smallest gap on the
transcriptome arm AND inside the 40% smallest gap on the clinical arm. The development-era 23-column
block does not reproduce it (3,969 pairs, 0.6165), and neither does either stage-group reading,
which stays below as a sensitivity row.

ROUNDING MATTERS HERE, AND WHY THE CLINICAL ARM IS REFITTED. The committed per-case vectors are
rounded to five decimals. Within-fold percentiles from folds of different sizes can then collide
across folds, and a restriction defined by a quantile of score GAPS is sensitive to exactly those
collisions: the rounded clinical vector gives 3,892 pairs and 0.6209 where the full-precision arm
gives 3,893 and 0.6212. The clinical arm is therefore refitted here at full precision with the
frozen estimator. The slide and transcriptome vectors reproduce every known answer as committed.

INPUTS, all committed or frozen:
  - seed-0 within-fold percentile scores of the slide and transcriptome arms, from
    experiments/20260911-blca-posthoc/results/percase-canonical-vectors.json (written by s19);
  - the amended clinical arm and the development-era 23-column arm, refitted here at seed 0 with
    the frozen estimator (blca_common.fitapply) on the frozen inputs;
  - stage groups from the incumbent's released split files (amended) and from tnm.json
    (development era), for the sensitivity row only.

KNOWN ANSWERS, refused unless reproduced: 24,219 comparable pairs; seed-0 slide 0.6596 and
transcriptome 0.6510 on all pairs; the refitted amended clinical arm at 0.6638 and the
development-era arm at 0.6856 on all pairs (amendment-A1-clinical-provenance.json, c6.json); and the
note's own values.

Post-freeze and descriptive. Nothing here changes a reported primary.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, fitapply       # noqa: E402
from s6_amend_clinical import dimaf_clinical                 # noqa: E402

KNOWN = {"pairs": 24219, "slide_all": 0.6596, "omics_all": 0.6510, "clinical_all": 0.6638,
         "clinical_dev_all": 0.6856}
NOTE = {"q40_pairs": 9688, "q40_slide": 0.6396, "q40_omics": 0.5604,
        "double_pairs": 3893, "double_slide": 0.6212}
Q = 0.40


def stage_dimaf(dimaf_dir, cases):
    rec = {}
    for k in range(5):
        for f in ("train", "test"):
            for r in csv.DictReader(open(os.path.join(dimaf_dir, "%d_%s.csv" % (k, f)))):
                rec.setdefault(r["case_id"], r)
    bad = ("[Not Available]", "[Unknown]", "", "[Not Evaluated]")
    out = []
    for c in cases:
        v = (rec.get(c) or {}).get("ajcc_pathologic_tumor_stage", "")
        out.append("" if v in bad else v)
    return out


def stage_tnm(tnm_json, cases):
    tnm = json.load(open(tnm_json))
    return [str((tnm.get(c) or {}).get("ajcc_pathologic_stage") or "") for c in cases]


def tied(z, ii, jj, q=Q):
    gap = np.abs(z[ii] - z[jj])
    return gap <= float(np.quantile(gap, q))


def arms_on(mask, arms, ii, jj):
    return {"pairs": int(mask.sum()),
            **{k: round(cidx(v, ii[mask], jj[mask]), 4) for k, v in arms.items()}}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("vectors", "root", "titan", "sp-dir", "dimaf-dir", "tnm", "out"):
        ap.add_argument("--" + f, required=True)
    a = ap.parse_args()

    cases = json.load(open(a.vectors))["cases"]
    ids = [c["case_id"] for c in cases]
    t = np.array([c["months"] for c in cases], float)
    e = np.array([c["event"] for c in cases], float)
    z_sl = np.array([c["seed0_slide"] for c in cases], float)
    z_om = np.array([c["seed0_omics"] for c in cases], float)
    # both clinical arms, refitted at full precision with the frozen estimator at seed 0
    co = Cohort(a.root, a.titan, a.sp_dir)
    if list(co.keep) != ids:
        print(json.dumps({"status": "error", "error_code": "case_order_differs"}), file=sys.stderr)
        return 2
    if not (np.allclose(co.t, t) and np.array_equal(co.e, e)):
        print(json.dumps({"status": "error", "error_code": "labels_differ"}), file=sys.stderr)
        return 2

    def arm(X):
        p = np.full(len(ids), np.nan)
        for tri, vai in co.fold_indices():
            p[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], 0)
        return co.fold_pct(p)

    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    z_cl = arm(np.hstack([dc["age"], dc["fem"], dc["stage"]]))
    z_dev = arm(co.CLIN)
    z_cl_rounded = np.array([c["seed0_clinical_stage"] for c in cases], float)

    ii, jj = cpairs(t, e)
    got = {"pairs": int(ii.size), "slide_all": round(cidx(z_sl, ii, jj), 4),
           "omics_all": round(cidx(z_om, ii, jj), 4),
           "clinical_all": round(cidx(z_cl, ii, jj), 4),
           "clinical_dev_all": round(cidx(z_dev, ii, jj), 4)}
    for k in KNOWN:
        print("  known answer %-16s %s (want %s)" % (k, got[k], KNOWN[k]), file=sys.stderr)
    import blca_common
    a2_ok = blca_common.A2 and got["pairs"] == KNOWN["pairs"] and all(
        abs(got[k] - KNOWN[k]) <= blca_common.A2_SANITY for k in KNOWN if k != "pairs")
    if got != KNOWN and not a2_ok:
        print(json.dumps({"status": "error", "error_code": "known_answer_not_reproduced"}),
              file=sys.stderr)
        return 2

    arms = {"slide": z_sl, "omics": z_om, "clinical": z_cl}
    m_om = tied(z_om, ii, jj)
    amended = m_om & tied(z_cl, ii, jj)
    dev = m_om & tied(z_dev, ii, jj)
    rounded = m_om & tied(z_cl_rounded, ii, jj)
    reproduced = {"q40_pairs": int(m_om.sum()) == NOTE["q40_pairs"],
                  "q40_slide": round(cidx(z_sl, ii[m_om], jj[m_om]), 4) == NOTE["q40_slide"],
                  "q40_omics": round(cidx(z_om, ii[m_om], jj[m_om]), 4) == NOTE["q40_omics"],
                  "double_pairs": int(amended.sum()) == NOTE["double_pairs"],
                  "double_slide": round(cidx(z_sl, ii[amended], jj[amended]), 4)
                  == NOTE["double_slide"]}
    for k, v in reproduced.items():
        print("  note value   %-16s %s" % (k, "reproduced" if v else "NOT reproduced"),
              file=sys.stderr)
    # amendment A2: the August note's values are pre-A2 and are checked exactly under BLCA_A2=0 only
    if not all(reproduced.values()) and not blca_common.A2:
        print(json.dumps({"status": "error", "error_code": "note_not_reproduced"}), file=sys.stderr)
        return 2

    sens = {}
    for name, lab in (("amended_stage_group", stage_dimaf(a.dimaf_dir, ids)),
                      ("development_era_stage_group", stage_tnm(a.tnm, ids))):
        s = np.array(lab, dtype=object)
        same = (s[ii] != "") & (s[jj] != "") & (s[ii] == s[jj])
        sens[name] = arms_on(m_om & same, arms, ii, jj)

    rep = {"artifact_type": "s19c_double_tied_probe", "phase_of_origin": "post_freeze_2026-09-11",
           "reportable": True,
           "known_answer": {"reproduced": not blca_common.A2, **got, "note_values_reproduced": reproduced,
                            "amendment": "A2" if blca_common.A2 else None},
           "definitions": {
               "scores": "seed-0 within-fold percentiles",
               "transcriptome_tied": "comparable pairs whose transcriptome-arm gap is at or below "
                                     "its 0.40 quantile",
               "clinically_tied": "the same rule on the clinical arm's gap",
               "stage_group_tied": "sensitivity only: the same non-missing stage group on both sides"},
           "all_pairs": arms_on(np.ones(ii.size, bool), arms, ii, jj),
           "transcriptome_tied_q40": arms_on(m_om, arms, ii, jj),
           "transcriptome_and_clinically_tied_q40": {
               "amended_clinical_block": arms_on(amended, arms, ii, jj),
               "development_era_clinical_block": {
                   "pairs": int(dev.sum()), "slide": round(cidx(z_sl, ii[dev], jj[dev]), 4)},
               "amended_block_from_the_rounded_committed_vector": {
                   "pairs": int(rounded.sum()),
                   "slide": round(cidx(z_sl, ii[rounded], jj[rounded]), 4),
                   "why_reported": "the committed vectors are rounded to five decimals, and a "
                                   "gap-quantile restriction is sensitive to the cross-fold "
                                   "collisions that rounding creates"}},
           "sensitivity_transcriptome_tied_and_same_stage_group": sens}
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(rep, open(a.out, "w"), indent=1)
    r = rep["transcriptome_and_clinically_tied_q40"]["amended_clinical_block"]
    print("amended block: %d pairs, slide %.4f, transcriptome %.4f, clinical %.4f | stage-group "
          "sensitivity (amended) %d pairs, slide %.4f"
          % (r["pairs"], r["slide"], r["omics"], r["clinical"],
             sens["amended_stage_group"]["pairs"], sens["amended_stage_group"]["slide"]),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
