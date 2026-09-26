#!/usr/bin/env python3
"""Figure 4 -- the ten cases the frozen rule admitted, and what each modality said about them.

WHAT THIS IS FOR. The paper argues about morphology for nineteen pages and, before this figure,
never showed any. Here the reader looks at the tissue directly. The ten cases are the five ModRank
scores highest and the five it scores lowest, selected by a rule fixed before the run and reported
whole, including the cases that embarrass the model. Under each slide the three modality
percentiles are drawn as a stack, so agreement and disagreement between morphology, transcriptome
and the report are visible per patient rather than only in aggregate.

READ IT DOWNWARD, NOT ACROSS. The top row is high risk and the bottom row is low, so the honest
reading is the disagreements: TCGA-XF-AAN1 sits in the bottom row and died at 31.4 months, and
TCGA-2F-A9KT sits in the top row and was alive at 110.5. A figure of ten cases cannot support an
inference and this one is labelled as descriptive, which is also why it is drawn from a
pre-registered selection rule rather than from whichever cases looked convincing.
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import style                                                             # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "experiments", "20260817-blca-confirm", "results")
CASES = json.load(open(os.path.join(RES, "interpret-calibrate-cases.json")))["I_case_studies"]
IDX = json.load(open(os.path.join(HERE, "assets", "case_series", "index.json")))

style.apply()
C, INK, LAD = style.ROLE, style.INK, style.LADDER
GREY = "#8C8C8C"
MOD = {"image": C["slide"], "omics": C["omics"], "clinical": C["clinical"]}

# Slides differ in shape by more than a factor of three. Stretching them to a common cell would
# distort the tissue, which is not a thing to do to histology, so each is PADDED to the cell's
# aspect instead and the pad is the slide's own background white.
CELL_ASPECT = 1.34


def pad_to(img, aspect):
    h, w = img.shape[:2]
    if w / h < aspect:
        need = int(round(h * aspect)) - w
        pad = np.full((h, need // 2, 3), 0.965, dtype=img.dtype)
        return np.concatenate([pad, img, np.full((h, need - need // 2, 3), 0.965,
                                                 dtype=img.dtype)], axis=1)
    need = int(round(w / aspect)) - h
    pad = np.full((need // 2, w, 3), 0.965, dtype=img.dtype)
    return np.concatenate([pad, img, np.full((need - need // 2, w, 3), 0.965,
                                             dtype=img.dtype)], axis=0)


hi = [c for c in CASES["cases"] if c["selected_as"] == "high"]
lo = [c for c in CASES["cases"] if c["selected_as"] == "low"]
assert len(hi) == 5 and len(lo) == 5, "the rule admits five of each; got %d and %d" % (len(hi),
                                                                                      len(lo))
hi.sort(key=lambda c: -c["our_score_percentile"])
lo.sort(key=lambda c: c["our_score_percentile"])

FIGW, FIGH = 6.785, 5.45   # 172.3 mm (PI correction 2026-09-27: left as is; do not widen to 180 mm)
fig = plt.figure(figsize=(FIGW, FIGH))

L, R, TOPY = 0.052, 0.988, 0.905
CW = (R - L) / 5.0
IMH, BARH = 0.268, 0.048          # slide cell and the modality stack, in figure fractions
ROWY = {0: TOPY - IMH, 1: TOPY - IMH - 0.487}

for row, group in enumerate((hi, lo)):
    for col, c in enumerate(group):
        x = L + col * CW
        y = ROWY[row]
        img = plt.imread(os.path.join(HERE, "assets", "case_series", c["case_id"] + ".png"))
        img = pad_to(img[..., :3], CELL_ASPECT)
        ax = fig.add_axes([x + 0.006, y, CW - 0.012, IMH])
        ax.imshow(img, interpolation="lanczos", aspect="auto")
        ax.set_xticks([]); ax.set_yticks([])
        edge = C["ours"] if row == 0 else "#B8BEC9"
        for sp in ax.spines.values():
            sp.set_edgecolor(edge); sp.set_linewidth(1.1)

        # the three modality percentiles, as a stack the eye can compare down a column
        for k, key in enumerate(("image", "omics", "clinical")):
            by = y - 0.026 - k * (BARH / 3.0 + 0.004)
            axb = fig.add_axes([x + 0.006, by, CW - 0.012, BARH / 3.0])
            axb.set_xlim(0, 1); axb.set_ylim(0, 1); axb.axis("off")
            axb.add_patch(Rectangle((0, 0.12), 1.0, 0.76, facecolor="#EDEFF2", edgecolor="none"))
            axb.add_patch(Rectangle((0, 0.12), c[key], 0.76, facecolor=MOD[key], edgecolor="none"))
            if col == 0:
                axb.text(-0.045, 0.5, key, transform=axb.transAxes, ha="right", va="center",
                         fontsize=5.0, color=MOD[key])

        died = bool(c["event"])
        fig.text(x + CW / 2, y + IMH + 0.014, c["case_id"].replace("TCGA-", ""), ha="center",
                 va="bottom", fontsize=5.6, color=INK, fontweight="bold")
        # The field is the AJCC stage GROUP (Stage I-IV). An earlier label rewrote it as "pTIII", a
        # pathologic-T prefix on a stage-group numeral, which is valid in neither system.
        fig.text(x + 0.008, y - 0.104, "%s" % c["stage"], ha="left",
                 va="top", fontsize=5.2, color=GREY)
        fig.text(x + CW - 0.014, y - 0.104, "%.1f mo %s" % (c["observed_months"],
                                                            "died" if died else "alive"),
                 ha="right", va="top", fontsize=5.2,
                 color=C["competitor"] if died else GREY,
                 fontweight="bold" if died else "normal")

for row, (lab, sub) in enumerate((("Five highest", "ModRank risk percentile 0.93 to 0.97"),
                                  ("Five lowest", "0.07 to 0.10"))):
    fig.text(L, ROWY[row] + IMH + 0.052, lab, ha="left", va="bottom", fontsize=6.4,
             fontweight="bold", color=INK)
    fig.text(L + 0.115, ROWY[row] + IMH + 0.053, sub, ha="left", va="bottom", fontsize=5.6,
             color=GREY)

# at 0.030 this line sat under the bottom row's own stage/outcome text (y - 0.104 with the bottom
# row's y around 0.150 puts that text's own bottom edge at almost the same height); lower still
# clears it without crowding the figure's own bottom edge.
fig.text(L, 0.006, "Descriptive. Ten cases cannot resolve an association and none of this fed "
         "back into the method.", ha="left", va="bottom", fontsize=5.4, color=GREY)

n_died_hi = sum(c["event"] for c in hi)
n_died_lo = sum(c["event"] for c in lo)
assert {r["case_id"] for r in IDX["cases"]} == {c["case_id"] for c in CASES["cases"]}, (
    "the fetched slides and the case series must be the same ten cases")

style.save(fig, os.path.join(HERE, "fig4_case_array.pdf"))
print("wrote fig4_case_array.pdf   %.3f x %.3f in" % (FIGW, FIGH))
print("  %d high (%d died), %d low (%d died)" % (len(hi), n_died_hi, len(lo), n_died_lo))
