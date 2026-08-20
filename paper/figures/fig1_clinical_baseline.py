#!/usr/bin/env python3
"""Figure 1 -- the benchmark's clinical baseline, as published and as it should be.

Panel A: bladder. Every published entrant's concordance, with the two clinical constructions
underneath them -- the grade-based one the field uses, and the stage-based one built from the same
files. The point of the panel is the position of the stage bar relative to the method bars.

Panel B: the same swap on all five TCGA studies the benchmark covers, as a paired difference.

The source-data file is written beside the PDF, because a Source Data upload assembled at
submission time from numbers that live only inside a figure is the version of this that fails.
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style  # noqa: E402

style.apply()
HERE = os.path.dirname(os.path.abspath(__file__))

# ---- source data. Every value traces to a run under experiments/ or development/.
PUBLISHED = [
    ("MOAD-FNet", 0.691, "unverified"),
    ("DIMAF", 0.679, "verified"),
    ("APL", 0.677, "uncheckable"),
    ("PIBD", 0.667, "verified"),
    ("DSCASurv", 0.646, "unverified"),
    ("OTSurv", 0.637, "unverified"),
    ("MMP", 0.635, "unverified"),
    ("SurvPath", 0.625, "verified"),
    ("MOTCat", 0.596, "verified"),
]
OURS = 0.7212
CLIN_GRADE, CLIN_STAGE = 0.5666, 0.6638
SWAP = [("BLCA", 0.5658, 0.6218), ("BRCA", 0.4828, 0.5686), ("COADREAD", 0.4690, 0.6898),
        ("HNSC", 0.5027, 0.5450), ("STAD", 0.5358, 0.5819)]

fig = plt.figure(figsize=(7.0, 3.1))
gs = fig.add_gridspec(1, 2, width_ratios=[1.28, 1.0], wspace=0.24,
                      left=0.208, right=0.985, top=0.895, bottom=0.155)

# ---------------------------------------------------------------- panel A
# Reading order is best to worst, top to bottom: ours, the published field, then the two clinical
# constructions set below a gap so they read as a different kind of object rather than as two more
# methods.
axA = fig.add_subplot(gs[0, 0])
names = [n for n, _, _ in PUBLISHED]
vals = [v for _, v, _ in PUBLISHED]
pos_pub = np.arange(1, len(names) + 1)
POS_OURS, POS_STAGE, POS_GRADE = 0.0, len(names) + 1.7, len(names) + 2.7

axA.barh([POS_OURS], [OURS], height=0.62, color=style.NOMINAL[0], edgecolor="none")
axA.barh(pos_pub, vals, height=0.62, color=style.NEUTRAL, edgecolor="none")
axA.barh([POS_STAGE], [CLIN_STAGE], height=0.62, color=style.NOMINAL[3], edgecolor="none")
axA.barh([POS_GRADE], [CLIN_GRADE], height=0.62, color=style.NOMINAL[1], edgecolor="none")

axA.set_yticks([POS_OURS] + list(pos_pub) + [POS_STAGE, POS_GRADE])
axA.set_yticklabels(["ours"] + names + ["clinical: age+sex+stage",
                                        "clinical: age+sex+grade"])
for t, w in zip(axA.get_yticklabels(), ["bold"] + ["normal"] * len(names) + ["bold", "normal"]):
    t.set_fontweight(w)
axA.set_ylim(POS_GRADE + 0.75, POS_OURS - 0.75)
axA.set_xlim(0.50, 0.775)
axA.set_xlabel("concordance index, TCGA-BLCA disease-specific survival")
axA.axvline(CLIN_STAGE, color=style.INK, lw=0.7, ls=(0, (3, 2)), zorder=0)
axA.set_title("A   the field against a correct clinical arm", loc="left",
              fontweight="bold")
for yy, vv in zip([POS_OURS] + list(pos_pub) + [POS_STAGE, POS_GRADE],
                  [OURS] + vals + [CLIN_STAGE, CLIN_GRADE]):
    axA.text(vv + 0.004, yy, "%.3f" % vv, va="center", ha="left", fontsize=6.4)

# ---------------------------------------------------------------- panel B
axB = fig.add_subplot(gs[0, 1])
coh = [c for c, _, _ in SWAP]
g = np.array([a for _, a, _ in SWAP])
s_ = np.array([b for _, _, b in SWAP])
x = np.arange(len(coh))
axB.hlines(x, g, s_, color=style.INK, lw=0.9, zorder=1)
axB.scatter(g, x, s=26, color=style.NOMINAL[1], zorder=2, label="age+sex+grade")
axB.scatter(s_, x, s=26, color=style.NOMINAL[3], zorder=2, label="age+sex+stage")
axB.set_yticks(x)
axB.set_yticklabels(coh)
axB.set_ylim(len(coh) - 0.4, -1.15)     # headroom at the top for the legend
axB.set_xlim(0.44, 0.755)
axB.set_xlabel("concordance index")
axB.axvline(0.5, color=style.NEUTRAL, lw=0.7, zorder=0)
axB.set_title("B   the same swap, five studies", loc="left", fontweight="bold")
axB.legend(loc="upper right", handletextpad=0.35, borderaxespad=0.1, labelspacing=0.25)
for xx, a, b in zip(x, g, s_):
    axB.text(b + 0.007, xx, "+%.3f" % (b - a), va="center", ha="left", fontsize=6.4,
             color=style.NOMINAL[3])

style.save(fig, os.path.join(HERE, "fig1_clinical_baseline.pdf"))

src = {"panel_A": {"published": [{"method": n, "cindex": v, "fold_identity": f}
                                 for n, v, f in PUBLISHED],
                   "ours": OURS, "clinical_age_sex_stage": CLIN_STAGE,
                   "clinical_age_sex_grade": CLIN_GRADE},
       "panel_B": [{"cohort": c, "age_sex_grade": a, "age_sex_stage": b, "difference": round(b - a, 4)}
                   for c, a, b in SWAP],
       "provenance": "panel A ours and clinical arms: experiments/20260817-blca-confirm/results/"
                     "amendment-A1-clinical-provenance.json and amended-secondary.json; published "
                     "values from each primary paper's own table; panel B: the five-cohort swap "
                     "described in development/five-cohort-stage-evidence.md"}
json.dump(src, open(os.path.join(HERE, "fig1_source_data.json"), "w"), indent=1)
print("wrote fig1_clinical_baseline.pdf and fig1_source_data.json")
