#!/usr/bin/env python3
"""Amendment A2: the headline quantities before and after the correction, side by side.

The before value of each quantity is read from the result file as committed at the last commit
whose results were all produced with the original tie handling (--before, recorded in
experiments/20260928-amendment-a2/protocol.md); the after value is read from the working tree once
the A2 results have replaced them. Nothing is recomputed here. A quantity that exists only after
A2 (the tuned-ridge stacking comparator) has a before value of null.

The Supplementary Information's amendment table is bound to this file by paper/check_numbers.py.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONF = "experiments/20260817-blca-confirm/results/"
POST = "experiments/20260911-blca-posthoc/results/"
FI = "experiments/20260927-field-inflation/analysis-results/"
A1 = CONF + "amendment-A1-clinical-provenance.json"
UF = POST + "unified-fusion-and-added-value.json"
RS = POST + "resplit-and-site-cv.json"
CAL = POST + "calibration-and-decision-curve.json"
PP = CONF + "pibd-parity.json"
FC = "experiments/20260912-five-cohort-intervals/results/five-cohort-intervals.json"
FU = "experiments/20260913-five-cohort-fusion/results/five-cohort-fusion.json"
GX = "experiments/20260911-geo-external/results/geo-external.json"
SN = "experiments/20260818-selection-null/results/selection-null.json"
HOLM = ("holm_family_of_six", "comparisons")
FUS = ("fusion_on_identical_inputs",)
AV = ("added_value_inflation", "constructions")
NB = ("constructions",)

# (group, label, file, path)
ROWS = [
    ("TCGA-BLCA", "ModRank, mean over five seeds", A1, ("primary", "value")),
    ("TCGA-BLCA", "ModRank, SD over seeds", A1, ("primary", "sd_over_seeds")),
    ("TCGA-BLCA", "clinical reference with grade", A1, ("arms_seed0", "clinical_age_sex_GRADE_from_DIMAF_file")),
    ("TCGA-BLCA", "clinical reference with stage", A1, ("arms_seed0", "clinical_age_sex_stage_from_DIMAF_file")),
    ("TCGA-BLCA", "stage minus grade", A1, ("arms_seed0", "stage_minus_grade")),
    ("TCGA-BLCA", "ModRank minus best verified published", A1, ("comparisons", "vs_DIMAF_published", "gap")),
    ("TCGA-BLCA", "ModRank minus SurvPath + clinical", PP, HOLM + ("ours vs SurvPath+clinical (seed-matched)", "delta")),
    ("TCGA-BLCA", "ModRank minus PIBD + clinical", PP, HOLM + ("ours vs PIBD+clinical (best-val checkpoint)", "delta")),
    ("TCGA-BLCA", "ModRank minus clinical, Holm p", PP, HOLM + ("ours vs clinical alone", "p_holm")),
    ("TCGA-BLCA", "calibration slope", CAL, ("ours", "calibration_slope")),
    ("TCGA-BLCA", "IPA at two years", CAL, ("ours", "by_horizon", "24.0", "ipa")),
    ("TCGA-BLCA", "selection-null bar", SN, ("A_selection_null", "bars", "empirical_max_q95", "bar")),
    ("TCGA-BLCA", "inductive variant", CONF + "cold-panel-round1-measurements.json",
     ("transductive_vs_inductive", "primary_inductive")),
    ("Fusion, TCGA-BLCA", "concatenation, mean", UF, FUS + ("concatenated_ridge_cox", "mean")),
    ("Fusion, TCGA-BLCA", "stacking (ridge 1), mean", UF, FUS + ("stacked_learned_weights", "mean")),
    ("Fusion, TCGA-BLCA", "stacking (tuned ridge), mean", UF, FUS + ("stacked_tuned_ridge", "mean")),
    ("Fusion, TCGA-BLCA", "ModRank minus concatenation", UF, FUS + ("paired", "ours_vs_concatenated", "mean")),
    ("Fusion, TCGA-BLCA", "ModRank minus stacking (ridge 1)", UF, FUS + ("paired", "ours_vs_stacked", "mean")),
    ("Fusion, TCGA-BLCA", "ModRank minus stacking (tuned ridge)", UF,
     FUS + ("paired", "ours_vs_stacked_tuned_ridge", "mean")),
    ("Fusion, TCGA-BLCA", "re-partitions, minus stacking", RS, ("resplits", "ours_minus_stacked", "mean")),
    ("Fusion, TCGA-BLCA", "re-partitions, minus concatenation", RS, ("resplits", "ours_minus_concat", "mean")),
    ("Fusion, TCGA-BLCA", "held-out sites, ModRank", RS, ("site_grouped_cv", "pooled", "ours")),
    ("Fusion, TCGA-BLCA", "held-out sites, stacking", RS, ("site_grouped_cv", "pooled", "stacked")),
    ("Added value", "D, ModRank", UF, AV + ("ours", "D")),
    ("Added value", "D, SurvPath", UF, AV + ("survpath", "D")),
    ("Added value", "D, PIBD", UF, AV + ("pibd_best_val", "D")),
    ("Added value", "excess over no signal, ModRank", FI + "inflation-null-blca.json",
     NB + ("ModRank", "excess_over_no_signal", "mean")),
    ("Added value", "excess over no signal, SurvPath", FI + "inflation-null-blca.json",
     NB + ("SurvPath", "excess_over_no_signal", "mean")),
    ("Added value", "excess over no signal, PIBD", FI + "inflation-null-blca.json",
     NB + ("PIBD", "excess_over_no_signal", "mean")),
    ("Added value", "D, DeepMISL with pathways", FI + "field-inflation-a.json",
     ("models", "deepmisl_wsi_pathways", "D")),
    ("Added value", "GSE32894, D", GX, ("cohorts", "GSE32894", "D", "point")),
    ("Added value", "GSE32894, excess over no signal", FI + "inflation-null-geo.json",
     NB + ("GSE32894", "excess_over_no_signal", "mean")),
    ("Added value", "GSE31684, excess over no signal", FI + "inflation-null-geo.json",
     NB + ("GSE31684", "excess_over_no_signal", "mean")),
    ("Added value", "GSE19915, stage minus grade reference", FI + "reference-gap-external.json",
     ("cohorts", "GSE19915", "gap")),
]
for c, name in (("blca", "BLCA"), ("brca", "BRCA"), ("coadread", "COADREAD"), ("hnsc", "HNSC"), ("stad", "STAD")):
    ROWS.append(("Five studies", "%s, ModRank" % name, FC, ("cohorts", c, "points", "OURS_wsi_omics_age_sex_stage")))
    ROWS.append(("Five studies", "%s, ModRank minus clinical" % name, FC,
                 ("cohorts", c, "differences", "ours_minus_clinical_stage", "mean")))
    ROWS.append(("Five studies", "%s, ModRank minus stacking (tuned ridge)" % name, FU,
                 ("cohorts", c, "stacked_tuned_ridge", "ModRank_minus_stacked_tuned_ridge", "mean")))


def get(d, path):
    for k in path:
        if isinstance(d, dict) and k in d:
            d = d[k]
        elif isinstance(d, list) and isinstance(k, int) and k < len(d):
            d = d[k]
        else:
            return None
    return d


def at_commit(commit, rel):
    out = subprocess.run(["git", "-C", ROOT, "show", "%s:%s" % (commit, rel)], capture_output=True, text=True)
    if out.returncode != 0:
        raise SystemExit(json.dumps({"status": "error", "error_code": "before_file_unreadable", "file": rel}))
    return json.loads(out.stdout)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--before", required=True, help="commit whose result files are all pre-A2")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    before, after, rows, missing = {}, {}, [], []
    for group, label, rel, path in ROWS:
        if rel not in before:
            before[rel] = at_commit(a.before, rel)
            after[rel] = json.load(open(os.path.join(ROOT, rel)))
        b, v = get(before[rel], path), get(after[rel], path)
        if v is None:
            missing.append(label)
        signed = (("minus" in label and "Holm" not in label) or label.startswith("D,")
                  or label.endswith(", D") or "excess" in label)
        rows.append({"group": group, "quantity": label, "signed": signed, "before": b, "after": v,
                     "change": None if b is None or v is None else round(v - b, 4),
                     "source": rel, "field": "/".join(str(k) for k in path)})
    if missing:
        print(json.dumps({"status": "error", "error_code": "after_value_missing", "which": missing}),
              file=sys.stderr)
        return 2
    rep = {"artifact_type": "s37_amendment_a2_before_after", "phase_of_origin": "amendment_A2",
           "before_commit": a.before,
           "reads": "values as stored in each result file; nothing recomputed",
           "rows": rows}
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(rep, open(a.out, "w"), indent=1)
    for r in rows:
        print("%-20s %-45s %8s -> %8s" % (r["group"][:20], r["quantity"][:45], r["before"], r["after"]),
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
