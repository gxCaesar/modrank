#!/usr/bin/env python3
"""Fetch the recorded race and ethnicity for this cohort, and find out whether they can be used.

WHY. The TRIPOD+AI fairness rows were answered on sex and age because the benchmark's released files
carry nothing else. Race is recorded by TCGA and is one request away, so "the data does not have it"
was a statement about the released subset rather than about the study. This asks the Genomic Data
Commons directly.

WHAT IT WRITES. Two files with different rules. The patient-level map goes to the frozen-inputs
directory, which is gitignored and not redistributed, because it is another group's per-patient
record and this repository points at sources rather than copying them. The counts go to a committed
result file, because a count is the thing the manuscript can report and a reader needs in order to
judge whether any subgroup claim was possible at all.

WHAT IT IS LIKELY TO FIND, written here before the run so the answer is not shaped by it: TCGA-BLCA
is heavily white, so the question is not what the concordance is in each group but whether any group
other than the largest clears the reporting floors this work already uses, 25 patients and 100
comparable pairs. If none does, the honest report is that the cohort cannot support the analysis,
with the counts shown.

KNOWN ANSWERS. The project returns 412 cases and all 359 of ours must match one of them, or the
identifiers have drifted and nothing below is about our patients.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import sys
import urllib.parse
import urllib.request

API = "https://api.gdc.cancer.gov/cases"
EXPECT_PROJECT_CASES = 412


def fetch(project):
    flt = {"op": "in", "content": {"field": "project.project_id", "value": [project]}}
    q = urllib.parse.urlencode({"filters": json.dumps(flt),
                                "fields": "submitter_id,demographic.race,demographic.ethnicity",
                                "format": "JSON", "size": "1000"})
    with urllib.request.urlopen("%s?%s" % (API, q), timeout=120) as r:
        return json.load(r)["data"]["hits"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--meta", required=True, help="the cohort's released metadata csv")
    ap.add_argument("--project", default="TCGA-BLCA")
    ap.add_argument("--out-map", required=True, help="patient-level, into gitignored frozen inputs")
    ap.add_argument("--out-counts", required=True, help="counts only, committed")
    a = ap.parse_args()

    hits = fetch(a.project)
    if len(hits) != EXPECT_PROJECT_CASES:
        print("%s returned %d cases, expected %d" % (a.project, len(hits), EXPECT_PROJECT_CASES),
              file=sys.stderr)
        return 2
    rec = {h["submitter_id"]: {"race": (h.get("demographic") or {}).get("race"),
                               "ethnicity": (h.get("demographic") or {}).get("ethnicity")}
           for h in hits}

    ours = sorted({r["case_id"] for r in csv.DictReader(open(a.meta))})
    missing = [c for c in ours if c not in rec]
    if missing:
        print("%d of our cases are not in the %s response: %s"
              % (len(missing), a.project, missing[:5]), file=sys.stderr)
        return 2

    os.makedirs(os.path.dirname(os.path.abspath(a.out_map)), exist_ok=True)
    json.dump({c: rec[c] for c in ours}, open(a.out_map, "w"), indent=1)

    race = collections.Counter((rec[c]["race"] or "missing") for c in ours)
    eth = collections.Counter((rec[c]["ethnicity"] or "missing") for c in ours)
    counts = {"artifact_type": "s30_gdc_demographics_counts",
              "phase_of_origin": "post_freeze_2026-09-12",
              "source": "Genomic Data Commons cases endpoint, %s, fields submitter_id, "
                        "demographic.race and demographic.ethnicity" % a.project,
              "patient_level_map": "written to the gitignored frozen-inputs tree and not "
                                   "redistributed; this file carries counts only",
              "project_cases": len(hits), "cohort_cases": len(ours),
              "race": dict(race.most_common()), "ethnicity": dict(eth.most_common())}
    os.makedirs(os.path.dirname(os.path.abspath(a.out_counts)), exist_ok=True)
    json.dump(counts, open(a.out_counts, "w"), indent=1)
    print("project %s: %d cases, %d matched to this cohort" % (a.project, len(hits), len(ours)),
          file=sys.stderr)
    for k, v in race.most_common():
        print("  %-32s %3d" % (k, v), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
