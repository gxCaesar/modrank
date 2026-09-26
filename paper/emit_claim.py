#!/usr/bin/env python3
"""Print one registered claim's value, recomputed from the rows under it where the rows are local.

Each headline claim of the paper is checked with this script. Its last line of stdout is the value
the manuscript prints. Where a per-case table exists the number is RECOMPUTED
from it with the project's own concordance (analysis/blca_common.cindex), and the recomputation is
checked against the committed summary before anything is printed: a disagreement exits 2 rather
than printing either value.

    /usr/bin/python3 paper/emit_claim.py <claim_id> [--rows <per-case .jsonl>]
    /usr/bin/python3 paper/emit_claim.py --export-jsonl     # writes the per-case table the claims read
"""
from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))

import numpy as np  # noqa: E402
from blca_common import cindex  # noqa: E402

POSTHOC = os.path.join(ROOT, "experiments", "20260911-blca-posthoc", "results")
PERCASE_JSON = os.path.join(POSTHOC, "percase-canonical-vectors.json")
PERCASE_JSONL = os.path.join(POSTHOC, "percase-canonical-vectors.jsonl")
A1 = os.path.join(ROOT, "experiments", "20260817-blca-confirm", "results",
                  "amendment-A1-clinical-provenance.json")
UNIFIED = os.path.join(POSTHOC, "unified-fusion-and-added-value.json")
GEO = os.path.join(ROOT, "experiments", "20260911-geo-external", "results", "geo-external.json")
FIVE = os.path.join(ROOT, "experiments", "20260912-five-cohort-intervals", "results",
                    "five-cohort-intervals.json")
FUSION = os.path.join(ROOT, "experiments", "20260913-five-cohort-fusion", "results",
                      "five-cohort-fusion.json")


def load(p):
    with open(p) as fh:
        return json.load(fh)


ROWS = [PERCASE_JSONL]


def rows():
    # The registered per-unit table. It is exported from the committed JSON by --export-jsonl and
    # must carry the same cases, so a stale export is caught here rather than read as data.
    with open(ROWS[0]) as fh:
        cases = [json.loads(line) for line in fh if line.strip()]
    if [c["case_id"] for c in cases] != [c["case_id"] for c in load(PERCASE_JSON)["cases"]]:
        print("emit_claim: the per-case table no longer matches its source; re-export it",
              file=sys.stderr)
        sys.exit(2)
    return cases


def conc(cases, key):
    t = np.array([c["months"] for c in cases], float)
    e = np.array([c["event"] for c in cases], int)
    return cindex(np.array([c[key] for c in cases], float), t, e)


def agree(label, recomputed, committed, tol=5e-5):
    if abs(recomputed - committed) > tol:
        print("emit_claim: %s recomputes to %.6f from the rows but the committed summary says %.6f"
              % (label, recomputed, committed), file=sys.stderr)
        sys.exit(2)


def modrank_tcga_benchmark():
    # The primary is the mean of five per-seed concordances. Only seed 0's per-case vector is
    # local, so seed 0 is recomputed from rows and must equal its committed per-seed value.
    a1 = load(A1)
    per_seed = {s["seed"]: s["OURS"] for s in a1["per_seed"]}
    agree("seed-0 ModRank", conc(rows(), "seed0_ours"), per_seed[0])
    value = float(np.mean([per_seed[s] for s in a1["primary"]["seeds"]]))
    agree("primary", value, a1["primary"]["value"])
    return "%.4f" % value


def stage_reference_correction():
    cases = rows()
    d = conc(cases, "seed0_clinical_stage") - conc(cases, "seed0_clinical_grade")
    agree("stage minus grade", d, load(A1)["arms_seed0"]["stage_minus_grade"])
    return "%+.4f" % d


def clinical_reference_inflation():
    cases = rows()
    d = ((conc(cases, "ours_grade") - conc(cases, "clinical_grade"))
         - (conc(cases, "ours") - conc(cases, "clinical_stage")))
    agree("D for ModRank", d, load(UNIFIED)["added_value_inflation"]["constructions"]["ours"]["D"])
    return "%+.4f" % d


GEO_ROWS = [os.path.join(ROOT, "experiments", "20260927-geo-percase", "results", "geo-percase.jsonl")]


def geo_inflation_replication():
    # D in GSE32894 from the repeat-averaged per-case scores the 2026-09-27 re-run of s24 wrote,
    # checked against the committed 2026-09-11 summary.
    with open(GEO_ROWS[0]) as fh:
        cases = [r for r in (json.loads(x) for x in fh if x.strip()) if r["analysis"] == "GSE32894"]
    d = ((conc(cases, "modrank_with_weak") - conc(cases, "clinical_weak"))
         - (conc(cases, "modrank") - conc(cases, "clinical_stage")))
    agree("GEO D, GSE32894", d, load(GEO)["cohorts"]["GSE32894"]["D"]["point"])
    return "%+.4f" % d


def five_study_validation():
    fcc = load(FIVE)["cohorts"]
    n = sum(1 for v in fcc.values() if v["differences"]["ours_minus_clinical_stage"]["ci95"][0] > 0)
    return str(n)


def fusion_null_five_studies():
    fus = load(FUSION)["cohorts"]
    n = sum(1 for v in fus.values()
            for k in ("ModRank_minus_concatenated", "ModRank_minus_stacked")
            if v[k]["ci95"][0] > 0 or v[k]["ci95"][1] < 0)
    return str(n)


def export_jsonl():
    cases = load(PERCASE_JSON)["cases"]
    with open(PERCASE_JSONL, "w") as fh:
        for c in cases:
            fh.write(json.dumps(c, sort_keys=True) + "\n")
    print("%d rows -> %s" % (len(cases), os.path.relpath(PERCASE_JSONL, ROOT)))


CLAIMS = {f.__name__: f for f in (modrank_tcga_benchmark, stage_reference_correction,
                                  clinical_reference_inflation, geo_inflation_replication,
                                  five_study_validation, fusion_null_five_studies)}


def main(argv):
    if argv == ["--export-jsonl"]:
        export_jsonl()
        return 0
    if len(argv) == 3 and argv[1] == "--rows":
        (GEO_ROWS if argv[0] == "geo_inflation_replication" else ROWS)[0] = os.path.join(ROOT, argv[2])
        argv = argv[:1]
    if len(argv) != 1 or argv[0] not in CLAIMS:
        print("usage: emit_claim.py <%s> [--rows <per-case .jsonl>] | --export-jsonl"
              % "|".join(CLAIMS), file=sys.stderr)
        return 2
    print(CLAIMS[argv[0]]())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
