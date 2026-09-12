#!/usr/bin/env python3
"""Source Data for the Nature Communications build: one file per display item.

The journal asks for the values behind each figure. This assembles them from the SAME committed
result files the figure builders read, so a Source Data file cannot drift from the figure: every
entry names its source file and the path inside it, and the script fails if a path is missing.

Figures 5 and 6 already have a builder-written source file beside the artwork
(fig6_robustness.source.json, fig7_added_value_utility.source.json); those are copied verbatim
rather than rebuilt, which is the one way to be sure they match what was drawn.

Run:  /usr/bin/python3 paper/nc-submission/make_source_data.py
"""

from __future__ import annotations

import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(HERE, "source_data")
CONF = os.path.join(ROOT, "experiments", "20260817-blca-confirm", "results")
POST = os.path.join(ROOT, "experiments", "20260911-blca-posthoc", "results")
S5 = os.path.join(ROOT, "development", "s5-results")
GEO = os.path.join(ROOT, "experiments", "20260911-geo-external", "results")
FIG = os.path.join(ROOT, "paper", "figures")


def load(path):
    with open(path) as fh:
        return json.load(fh)


def dig(obj, path, where):
    if not path:
        return obj                     # an empty path means the whole file is the panel's source
    cur = obj
    for part in path.split("/"):
        if isinstance(cur, list):
            cur = cur[int(part)]
        elif part in cur:
            cur = cur[part]
        else:
            raise KeyError("%s: no such field %s (stopped at %r)" % (where, path, part))
    return cur


def panel(source_rel, path, description):
    """One panel's values, with the file and field they came from."""
    obj = load(os.path.join(ROOT, source_rel))
    return {"description": description, "source_file": source_rel, "field": path,
            "values": dig(obj, path, source_rel)}


def main():
    os.makedirs(OUT, exist_ok=True)
    written = []

    fig1 = {
        "figure": "Figure 1",
        "title": "The clinical reference, the method and the protocol",
        "a_clinical_constructions": panel(
            "experiments/20260817-blca-confirm/results/amendment-A1-clinical-provenance.json",
            "arms_seed0", "the two clinical arms rebuilt from the incumbent's own split files"),
        "b_parameter_counts": panel(
            "experiments/20260817-blca-confirm/results/parameter-counts.json", "ours",
            "our coefficient count, and the competitors' own trainer output beside it"),
        "b_competitor_counts": {
            "survpath": panel("experiments/20260817-blca-confirm/results/parameter-counts.json",
                              "survpath/parameters", "printed by SurvPath's own trainer"),
            "pibd": panel("experiments/20260817-blca-confirm/results/parameter-counts.json",
                          "pibd/parameters", "printed by PIBD's own trainer")},
        "c_protocol_and_falsifiers": panel(
            "development/s5-results/leakage.json", "falsifiers",
            "the five leakage falsifiers, two of which rejected"),
        "d_grade_and_stage_composition": panel(
            "experiments/20260817-blca-confirm/results/why-grade-fails.json", "cohorts/blca",
            "level composition and normalised entropy of grade and stage"),
    }
    fig2 = {
        "figure": "Figure 2",
        "title": "The missing variable, quantified",
        "a_published_entries": panel(
            "experiments/20260817-blca-confirm/results/published-benchmark-table.json", "entries",
            "every published value on this benchmark, with whether its folds could be verified"),
        "b_census": panel("development/clinical-baseline-census.json", "papers",
                          "one row per method paper: what it compares against"),
        "c_five_study_swap": panel(
            "experiments/20260911-five-study-repro/results/stage5.json", "cohorts",
            "the same one-column swap in each of the five studies"),
        "d_gain_against_grade_entropy": panel(
            "experiments/20260817-blca-confirm/results/why-grade-fails.json", "cohorts",
            "each study's grade entropy and the gain from swapping grade for stage"),
        "e_four_clinical_constructions": dict(panel(
            "development/s5-results/decomp.json", "",
            "four constructions of the clinical block, alone and with the other two modalities"),
            provenance_note="This is the one panel in the figure set whose values have no committed "
                            "emitter: the file was written on 2026-08-17 by a script that was never "
                            "committed. Two of its four rows (the amended block, and grade in place "
                            "of stage) reproduce exactly in amendment-A1-clinical-provenance.json; "
                            "the other two use a GDC record that the amendment replaced."),
        "f_within_stage_strata": panel(
            "experiments/20260817-blca-confirm/results/biology.json", "B5_within_molecular_subtype",
            "concordance inside strata where the clinical variable carries little information"),
    }
    fig3 = {
        "figure": "Figure 3",
        "title": "What limits a method at 359 patients and 113 events",
        "a_capacity_against_noise": panel(
            "experiments/20260817-blca-confirm/results/regime-evidence.json", "capacity_sweep",
            "a ridge Cox fitted to pure noise and scored in sample"),
        "b_weight_sweep": panel(
            "experiments/20260817-blca-confirm/results/regime-evidence.json", "fusion",
            "the 21-point weight sweep, the held-out-fold oracle and equal weighting"),
        "c_gated_component_and_control": panel(
            "development/s5-results/c2.json", "C5_gated_fusion",
            "the gated component against its permuted-stratum control"),
        "d_power": panel(
            "experiments/20260817-blca-confirm/results/reviewer-gaps.json", "m_detectable_margin",
            "power against a true difference at this cohort's patient-resampling standard error"),
        "e_competitor_across_seeds": panel(
            "experiments/20260817-blca-confirm/results/pibd-parity.json",
            "canonical_arms_seed_averaged_scores",
            "the competitor arms as seed-averaged scores"),
        "f_checkpoint_convention": panel(
            "experiments/20260817-blca-confirm/results/pibd-parity.json", "pibd",
            "PIBD at its best-validation checkpoint and at its final epoch"),
    }
    fig4 = {
        "figure": "Figure 4",
        "title": "The confirmatory result and the comparisons that do not survive correction",
        "a_holm_family": panel(
            "experiments/20260817-blca-confirm/results/pibd-parity.json",
            "holm_family_of_six/comparisons",
            "all six prespecified comparisons, their intervals and their Holm-adjusted p values"),
        "e_modality_subsets": panel(
            "experiments/20260817-blca-confirm/results/ablation-generalisation.json",
            "B_ablation_bladder", "all seven modality subsets and the leave-one-out differences"),
        "f_cost_against_benefit": panel(
            "experiments/20260817-blca-confirm/results/parameter-counts.json", "",
            "the two measured parameter counts against concordance"),
    }
    fig7 = {
        "figure": "Figure 7",
        "title": "The three modalities: what each is, and what each carries",
        "b_seven_representations": panel(
            "development/s5-results/c6.json", "single_arms",
            "seven representations across the three families, on identical patients and folds"),
        "c_redundancy": panel(
            "experiments/20260911-blca-posthoc/results/modality-atlas-amended.json", "redundancy",
            "Spearman correlation between the arms' out-of-fold scores"),
        "d_conditional_probe": panel(
            "experiments/20260911-blca-posthoc/results/modality-atlas-amended.json",
            "conditional_probe", "pairs restricted by how far apart the clinical model places them"),
        "g_tertile_separation": panel(
            "experiments/20260817-blca-confirm/results/biology.json", "B5_within_molecular_subtype",
            "log-rank separation across each arm's own risk tertiles"),
    }

    for name, obj in (("Figure_1", fig1), ("Figure_2", fig2), ("Figure_3", fig3),
                      ("Figure_4", fig4), ("Figure_7", fig7)):
        with open(os.path.join(OUT, name + ".json"), "w") as fh:
            json.dump(obj, fh, indent=1)
        written.append(name + ".json")

    # the two figures whose builder already wrote the plotted values
    for src, dst, title in (("fig6_robustness.source.json", "Figure_5.json",
                             "Fusion on identical inputs, re-partitions and held-out hospitals"),
                            ("fig7_added_value_utility.source.json", "Figure_6.json",
                             "Added value, its inflation, calibration and net benefit")):
        s = os.path.join(FIG, src)
        if not os.path.exists(s):
            print("missing builder source file: %s" % src, file=sys.stderr)
            return 2
        obj = {"figure": dst.replace("_", " ").replace(".json", ""), "title": title,
               "source_file": os.path.relpath(s, ROOT), "values": load(s)}
        with open(os.path.join(OUT, dst), "w") as fh:
            json.dump(obj, fh, indent=1)
        written.append(dst)

    print("wrote %d Source Data files to %s" % (len(written), os.path.relpath(OUT, ROOT)))
    for w in sorted(written):
        print("  %s  %d bytes" % (w, os.path.getsize(os.path.join(OUT, w))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
