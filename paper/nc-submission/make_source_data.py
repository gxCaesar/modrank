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
        "d_grade_and_stage_composition": dict(panel(
            "experiments/20260817-blca-confirm/results/why-grade-fails.json", "cohorts/blca",
            "level composition and normalised entropy of grade and stage, in the clinical file "
            "released with the benchmark"),
            provenance_note="Two released files carry a stage field for these patients. This panel "
                            "draws the benchmark's own clinical file, five levels over 334 valued "
                            "cases, because it is the only clinical source common to all five "
                            "studies and so is what the five-study comparison can use. The arms in "
                            "panel a, and every concordance in the paper, read the incumbent's "
                            "split files instead: four AJCC stage groups over 357 valued cases, "
                            "whose composition is the next entry and is what Figure 7a draws."),
        "d_stage_groups_the_model_reads": panel(
            "experiments/20260817-blca-confirm/results/per-stage-subgroup.json", "strata",
            "the four AJCC stage groups of the incumbent's split files, the field every arm in "
            "this paper actually reads, with the patients and events in each"),
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
        "e_three_clinical_constructions": dict(panel(
            "experiments/20260912-decomp-emitter/results/clinical-block-decomposition.json",
            "constructions",
            "three constructions of the clinical block, alone and with the other two modalities"),
            provenance_note="Until 2026-09-12 this panel read a file with no committed emitter and "
                            "drew four constructions. Rebuilding them found that two rows, not one, "
                            "could not be reproduced: the amended block and the grade block match "
                            "exactly, the frozen block does not (0.6856 and 0.7291 rebuilt against "
                            "0.6863 and 0.7285 in that file, with the rebuild agreeing with the "
                            "frozen run's own clinical arm), and the corrected-GDC row cannot be "
                            "attempted because the rule that chose one diagnosis per case was never "
                            "recorded. The panel now draws the three a script regenerates."),
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
    # Figure 7 is the five-study validation, added 2026-09-12. The modalities figure, which was
    # Figure 7 until then, moves to 8: display items are numbered in order of first citation and the
    # validation subsection is cited before the slide-contribution one.
    fig7 = {
        "figure": "Figure 7",
        "title": "The same recipe on the benchmark's five studies, with intervals",
        "a_vs_corrected_clinical": panel(
            "experiments/20260912-five-cohort-intervals/results/five-cohort-intervals.json",
            "cohorts", "every arm, every paired difference and its interval, per study"),
        "c_best_published_same_folds": panel(
            "experiments/20260817-blca-confirm/results/ablation-generalisation.json", "A_generalisation",
            "per study, the values published on the released folds and best_published_verified_folds, "
            "the value marked in panel c (PIBD outside bladder, DIMAF in bladder)"),
    }
    fig8 = {
        "figure": "Figure 8",
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
                      ("Figure_4", fig4), ("Figure_7", fig7), ("Figure_8", fig8)):
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
    return write_workbook(sorted(written))


def flatten(obj, path=""):
    """Every scalar under obj as (path, value), lists indexed, in document order."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from flatten(v, "%s/%s" % (path, k) if path else str(k))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from flatten(v, "%s[%d]" % (path, i))
    else:
        yield path, obj


def panels(obj, name=""):
    """(panel name, source file, field, values) for every panel in a Source Data object."""
    if isinstance(obj, dict) and "values" in obj and "source_file" in obj:
        if not name:                        # a builder-written file: one panel per key of values
            for pk, pv in obj["values"].items():
                yield pk, obj["source_file"], "", pv
        else:
            yield name, obj["source_file"], obj.get("field", ""), obj["values"]
    elif isinstance(obj, dict):
        for k, v in obj.items():
            if k not in ("figure", "title", "description"):
                yield from panels(v, "%s/%s" % (name, k) if name else k)


def write_workbook(names):
    """The journal's form of the same files (added 2026-09-27): one Excel workbook, one sheet per
    figure, each panel flattened losslessly to rows of panel, source file, field and value. The JSON
    files stay the canonical copy and the workbook is rebuilt from them on every run."""
    from openpyxl import Workbook
    wb = Workbook()
    wb.remove(wb.active)
    for name in names:
        obj = load(os.path.join(OUT, name))
        ws = wb.create_sheet(name.replace(".json", ""))
        ws.append([obj.get("figure", ""), obj.get("title", "")])
        ws.append(["Each row is one plotted or reported value: the figure panel, the committed result "
                   "file it was read from, the field inside that file, and the value."])
        ws.append([])
        ws.append(["panel", "source file", "field", "value"])
        n = 0
        for pname, src, field, values in panels(obj):
            for sub, value in flatten(values):
                ws.append([pname, src, "/".join(x for x in (field, sub) if x), value])
                n += 1
        if n == 0:
            print("workbook sheet %s has no values" % name, file=sys.stderr)
            return 2
        print("  %-10s %5d rows" % (name.replace(".json", ""), n))
    path = os.path.join(OUT, "SourceData.xlsx")
    wb.save(path)
    print("  SourceData.xlsx  %d sheets, %d bytes" % (len(names), os.path.getsize(path)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
