#!/usr/bin/env python3
"""Figure 7 -- the method on the benchmark's other four studies, with the intervals.

WHY THIS FIGURE EXISTS. The benchmark ships five TCGA studies and every method paper reports all
five, so evaluating on five is the benchmark's own convention. The three-modality arm had been run
on all of them since the S5 campaign and had never carried an interval, which is why the manuscript
could only call them a replication of the stage finding. This figure is the measurement that lets
them be called validation, and it is drawn so that the two cohorts where the method does NOT
separate are as visible as the three where it does.

  a  the difference that a validation claim rests on: the three-modality arm minus the corrected
     clinical reference, per cohort, with the 95% interval from 6,000 patient resamples. Three of
     five exclude zero. The two that do not are drawn identically, not faded.
  b  the same arm minus itself with grade substituted for stage. Only COADREAD excludes zero, and
     the panel says so rather than leaning on the four positive point estimates.
  c  every arm's concordance per cohort, against the value SurvPath reports on these same released
     folds. The best published value under ANY protocol is drawn as a separate open marker, because
     on BRCA and COADREAD it is far above everything here and a reader must not be able to mistake
     the narrow comparison for the wide one.

WHAT IS DELIBERATELY NOT DRAWN. No pooled estimate across cohorts. The five studies are different
diseases with different event counts and a pooled number would invite exactly the reading the
per-cohort intervals refuse: that the method is validated as a whole. It is validated in three of
five, and the figure is the list.
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import style                                                            # noqa: E402
from schematic import panel_letter                                      # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
IV = json.load(open(os.path.join(ROOT, "experiments", "20260912-five-cohort-intervals",
                                 "results", "five-cohort-intervals.json")))

style.apply()
C, INK = style.ROLE, style.INK

ORDER = ["blca", "brca", "coadread", "hnsc", "stad"]
NAME = {"blca": "BLCA (bladder)", "brca": "BRCA", "coadread": "COADREAD", "hnsc": "HNSC",
        "stad": "STAD"}
assert set(ORDER) == set(IV["cohorts"]), "the figure and the result file disagree on the cohorts"

fig = plt.figure(figsize=(style.width(6.785), 3.05))
axa = fig.add_axes([0.145, 0.175, 0.245, 0.700])
axb = fig.add_axes([0.455, 0.175, 0.245, 0.700])
axc = fig.add_axes([0.775, 0.175, 0.205, 0.700])

yb = np.arange(len(ORDER))[::-1]


def forest(ax, key, title, xlim):
    """One difference per cohort, with its interval. A cohort whose interval covers zero is drawn
    with the same ink as one whose interval does not: the zero line is what separates them, not a
    styling choice, because fading the null results is how a figure argues instead of reporting."""
    n_excl = 0
    for yi, k in zip(yb, ORDER):
        d = IV["cohorts"][k]["differences"][key]
        lo, hi = d["ci95"]
        excludes = lo > 0
        n_excl += excludes
        ax.plot([lo, hi], [yi, yi], color="#8A93A3", lw=1.3, solid_capstyle="round", zorder=2)
        ax.scatter([d["mean"]], [yi], s=30, zorder=3, edgecolor="none",
                   color=C["ours"] if excludes else "#FFFFFF")
        if not excludes:                       # an open marker reads as "interval covers zero"
            ax.scatter([d["mean"]], [yi], s=30, zorder=4, facecolor="none",
                       edgecolor=C["ours"], linewidths=0.9)
    ax.axvline(0.0, color=INK, lw=0.8, ls=":", zorder=1)
    ax.set_yticks(yb)
    ax.set_xlim(*xlim)
    ax.set_xlabel("difference in concordance", fontsize=6.4)
    ax.set_title(title, fontsize=6.3, pad=4)
    return n_excl


panel_letter(axa, "a", x=-0.34, y=1.05)   # these axes carry two-line cohort labels reaching almost
                                          # to the axes top; -0.235 still put the letter on top of
                                          # "BLCA (bladder)". Further left and above clears it.
n_a = forest(axa, "ours_minus_clinical_stage",
             "ModRank $-$ corrected clinical arm", (-0.13, 0.24))
axa.set_yticklabels(["%s\n$n=%d$, %d ev" % (NAME[k], IV["cohorts"][k]["n_cases"],
                                            IV["cohorts"][k]["events"]) for k in ORDER],
                    fontsize=5.9)
axa.text(0.5, -0.215, "%d of 5 exclude zero" % n_a, transform=axa.transAxes, ha="center",
         va="top", fontsize=5.8, color="#6B7280")

panel_letter(axb, "b", x=-0.06, y=1.05)  # outside the plot and clear of the left spine/axis
n_b = forest(axb, "ours_minus_grade_control",
             "ModRank $-$ same arm with grade", (-0.13, 0.24))
axb.set_yticklabels([])
axb.text(0.5, -0.215, "%d of 5 exclude zero" % n_b, transform=axb.transAxes, ha="center",
         va="top", fontsize=5.8, color="#6B7280")

# --------------------------------------------------------------------------- c  the two references
panel_letter(axc, "c", x=-0.05, y=1.05)
for yi, k in zip(yb, ORDER):
    p = IV["cohorts"][k]["points"]
    rel = IV["cohorts"][k]["published_best_verified_released_folds"][1]
    any_ = IV["cohorts"][k]["published_best_any_protocol"][1]
    axc.plot([rel, p["OURS_wsi_omics_age_sex_stage"]], [yi, yi], color="#B8BEC9", lw=1.3,
             solid_capstyle="round", zorder=1)
    axc.scatter([rel], [yi], s=24, color=C["competitor"], marker=style.ROLE_MARKER["competitor"],
                zorder=3, edgecolor="none")
    axc.scatter([p["OURS_wsi_omics_age_sex_stage"]], [yi], s=26, color=C["ours"],
                marker=style.ROLE_MARKER["ours"], zorder=3, edgecolor="none")
    axc.scatter([any_], [yi], s=26, zorder=3, facecolor="none", edgecolor=INK, linewidths=0.8)
axc.set_yticks(yb)
axc.set_yticklabels([])
axc.set_xlim(0.55, 0.90)          # room on the right for the legend, which sat on STAD's ModRank
                                  # point at 0.63 in the first render. Found by reading the PDF.
axc.set_xlabel("concordance", fontsize=6.4)
axc.scatter([], [], s=24, color=C["competitor"], marker=style.ROLE_MARKER["competitor"],
            label="SurvPath, these folds")
axc.scatter([], [], s=26, color=C["ours"], marker=style.ROLE_MARKER["ours"], label="ModRank")
axc.scatter([], [], s=26, facecolor="none", edgecolor=INK, linewidths=0.8,
            label="best published, any protocol")
# "upper left" covered the BLCA row; the legend is nearly as wide as this narrow panel, so ANY
# in-axes corner sits close to a full-width row of data at some y. Above the axes, one row, in the
# gap the title used to have to itself -- the title's own words fold into the legend's title, as in
# panel c of Figure 2 -- is the placement that touches no row's markers.
axc.legend(fontsize=5.1, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1,
           handletextpad=0.3, borderpad=0.2, labelspacing=0.2,
           title="against two references", title_fontsize=6.3)

style.save(fig, style.out(HERE, "fig8_five_cohort_validation.pdf"))
print("wrote fig8_five_cohort_validation.pdf   %.3f x %.3f in" % tuple(fig.get_size_inches()))
print("  a: %d of 5 separate from the corrected clinical arm" % n_a)
print("  b: %d of 5 separate from the grade control" % n_b)
print("  %d cases, %d events" % (IV["summary"]["total_cases"], IV["summary"]["total_events"]))
