#!/usr/bin/env python3
"""The committed emitter for why-grade-fails.json, which had none.

WHY. experiments/20260817-blca-confirm/results/why-grade-fails.json carries the numbers behind the
paper's mechanism (grade takes two levels with 94.4% of bladder patients in one of them, normalised
entropy 0.311, against 0.665 for stage) and the five-study table of the Supplementary Information.
It was written on 2026-08-17 by a script that was never committed. An external review on
2026-09-11 found the gap. This script recomputes every field from the benchmark's own released files
and REFUSES to write unless the result equals the committed file, so the committed file becomes a
reproduced artefact rather than a transcription.

WHAT.
  - level composition: the benchmark's clinical files (tcga_<study>_clinical.csv), one row per
    patient (the first), over the patients of the benchmark's metadata file. "N/A" and empty are
    missing. Normalised Shannon entropy is H / ln(K) over the K observed levels.
  - concordance of the age-sex-grade and age-sex-stage arms: read from the five-study results
    (experiments/20260911-five-study-repro/results/stage5.json), which re-derived them from restored
    inputs today with 80 of 80 fields identical to the original run. They are not refitted here.
  - the Pearson correlation between grade entropy and the stage-minus-grade gain, over the studies
    with a usable grade (a consistency check over three points, as the file already says).
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import math
import os
import sys

STUDIES = ("blca", "brca", "coadread", "hnsc", "stad")
MISSING = ("", "N/A", "nan", "NaN")


def composition(vals):
    v = [x for x in vals if x not in MISSING]
    cnt = collections.Counter(v)
    if len(cnt) < 2:
        return None
    n = len(v)
    h = -sum(k / n * math.log(k / n) for k in cnt.values())
    return {"n_with_a_value": n, "n_levels": len(cnt), "modal_share": round(max(cnt.values()) / n, 3),
            "normalised_entropy": round(h / math.log(len(cnt)), 3),
            "levels": dict(sorted(cnt.items(), key=lambda kv: -kv[1]))}


def pearson(x, y):
    mx, my = sum(x) / len(x), sum(y) / len(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--meta-dir", required=True)
    ap.add_argument("--clin-dir", required=True)
    ap.add_argument("--stage5", required=True)
    ap.add_argument("--committed", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    s5 = json.load(open(a.stage5))["cohorts"]
    committed = json.load(open(a.committed))
    cohorts = {}
    for st in STUDIES:
        cases = sorted({r["case_id"] for r in csv.DictReader(open(os.path.join(a.meta_dir, "tcga_%s.csv" % st)))})
        clin = {}
        for r in csv.DictReader(open(os.path.join(a.clin_dir, "tcga_%s_clinical.csv" % st))):
            clin.setdefault(r["case_id"], r)
        g = composition([(clin.get(c) or {}).get("grade", "") for c in cases])
        s = composition([(clin.get(c) or {}).get("stage", "") for c in cases])
        arms = s5[st]["arms"]
        cohorts[st] = {"n": len(cases), "grade": g, "stage": s,
                       "c_age_sex_grade": arms["age_sex_grade"], "c_age_sex_stage": arms["age_sex_stage"],
                       "stage_minus_grade": round(arms["age_sex_stage"] - arms["age_sex_grade"], 4)}
    with_grade = [st for st in STUDIES if cohorts[st]["grade"]]
    r = round(pearson([cohorts[st]["grade"]["normalised_entropy"] for st in with_grade],
                      [cohorts[st]["stage_minus_grade"] for st in with_grade]), 3)

    diffs = []
    for st in STUDIES:
        for key, val in cohorts[st].items():
            if val != committed["cohorts"][st].get(key):
                diffs.append((st, key, val, committed["cohorts"][st].get(key)))
    if r != committed["gain_vs_grade_entropy_pearson"]:
        diffs.append(("all", "gain_vs_grade_entropy_pearson", r, committed["gain_vs_grade_entropy_pearson"]))
    for d in diffs:
        print("  differs: %s %s recomputed %r committed %r" % d, file=sys.stderr)
    if diffs:
        print(json.dumps({"status": "error", "error_code": "does_not_reproduce_the_committed_file"}),
              file=sys.stderr)
        return 2
    rep = dict(committed)
    rep["reproduced_by"] = ("analysis/s7b_why_grade_fails.py on 2026-09-11: every field recomputed from "
                            "the benchmark's released clinical and metadata files and stage5.json, "
                            "equal to the committed file")
    json.dump(rep, open(a.out, "w"), indent=1)
    print("reproduced: %d studies, grade entropy %s, Pearson %.3f over %d studies with a grade"
          % (len(STUDIES), {st: cohorts[st]["grade"]["normalised_entropy"] for st in with_grade}, r,
             len(with_grade)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
