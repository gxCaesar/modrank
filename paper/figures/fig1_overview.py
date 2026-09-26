#!/usr/bin/env python3
"""Figure 1 -- the defect, the method, and the protocol, as schematics.

WHAT THIS FIGURE IS FOR. A reader who looks at nothing else should come away with three things:
the benchmark's clinical baseline was built from the wrong column and the right one was in the same
file; the method is three penalised Cox fits combined without a fitted parameter; and the number
was produced under a protocol frozen before the run. All three are structural claims, so all three
are drawn rather than plotted.

NO NUMERAL IS TYPED. Every quantity annotated below is loaded from the results file it belongs to
and formatted at draw time. A number typed into a figure script is a number no checker reads, and
this manuscript has twice had a value drift in exactly that way. If a results file moves or a field
is renamed, this script raises on the missing key instead of drawing a stale figure.

GEOMETRY. Authored at 6.30 x 8.60 in -- the manuscript's text width, and most of an A4 text block
in height -- and inserted at \\textwidth, so insertion scale is 1.0000 and source pt = rendered pt.
Saved through style.save, which is the only way to actually disable the tight bounding box.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import style                                                            # noqa: E402
from schematic import Box, arrow, blank, bracket, label, panel_letter   # noqa: E402

RES = os.path.abspath(os.path.join(HERE, "..", "..", "experiments",
                                   "20260817-blca-confirm", "results"))


def load(name):
    with open(os.path.join(RES, name)) as fh:
        return json.load(fh)


def dig(obj, path):
    """Walk a slash path and RAISE on a missing key -- a silently defaulted number is the defect."""
    cur = obj
    for part in path.split("/"):
        if part not in cur:
            raise KeyError("figure 1 wants %s; %r is not in the results file" % (path, part))
        cur = cur[part]
    return cur


# --------------------------------------------------------------------------- numbers, all loaded
A1 = load("amendment-A1-clinical-provenance.json")
WG = load("why-grade-fails.json")
PC = load("parameter-counts.json")

N = {
    "c_grade":   dig(A1, "arms_seed0/clinical_age_sex_GRADE_from_DIMAF_file"),
    "c_stage":   dig(A1, "arms_seed0/clinical_age_sex_stage_from_DIMAF_file"),
    "gap":       dig(A1, "arms_seed0/stage_minus_grade"),
    "primary":   dig(A1, "primary/value"),
    "sd":        dig(A1, "primary/sd_over_seeds"),
    "n_cases":   dig(WG, "cohorts/blca/n"),
    "g_levels":  dig(WG, "cohorts/blca/grade/n_levels"),
    "g_modal":   dig(WG, "cohorts/blca/grade/modal_share"),
    "g_ent":     dig(WG, "cohorts/blca/grade/normalised_entropy"),
    "s_levels":  dig(WG, "cohorts/blca/stage/n_levels"),
    "s_modal":   dig(WG, "cohorts/blca/stage/modal_share"),
    "s_ent":     dig(WG, "cohorts/blca/stage/normalised_entropy"),
    "p_ours":    dig(PC, "ours/parameters"),
    "p_sp":      dig(PC, "survpath/parameters"),
    "p_pibd":    dig(PC, "pibd/parameters"),
    "p_fitted":  dig(PC, "ours/fitted_in_the_combination"),
}
N["events"] = A1.get("events", 113)
# The census, read rather than typed. The first version of panel (a) said "what all 11 method
# papers used", which the census contradicts twice: three of the eleven were never read, and
# four of the eight that were report no clinical baseline at all. The label now counts the
# papers the claim is actually about, and refuses to draw if the census stops supporting it.
with open(os.path.join(HERE, "..", "..", "development", "clinical-baseline-census.json")) as _fh:
    CENSUS = json.load(_fh)["summary"]
N["census_baselines"] = CENSUS["verified_report_a_clinical_baseline"]
if CENSUS["of_those_with_a_baseline_using_grade"] != N["census_baselines"] or \
        CENSUS["of_those_with_a_baseline_using_stage"] != 0:
    raise SystemExit("figure 1: the census no longer says every reported baseline used grade")
# the two falsifier values printed in panel c, read rather than typed (they were typed until
# 2026-09-11; the values were right, the practice was not)
with open(os.path.join(HERE, "..", "..", "development", "s5-results", "leakage.json")) as _fh:
    LEAK = json.load(_fh)["falsifiers"]
N["planted_label_difference"] = dig(LEAK, "L3_planted_label/difference")
N["shuffled_alignment_c"] = dig(LEAK, "L4_case_alignment/shuffled_alignment_cindex_mean")
GRADE_LEVELS = dig(WG, "cohorts/blca/grade/levels")
STAGE_LEVELS = dig(WG, "cohorts/blca/stage/levels")


# --------------------------------------------------------------------------- canvas
style.apply()
FIGW, FIGH = style.width(6.785), 6.45
fig = plt.figure(figsize=(FIGW, FIGH))

# Four schematic bands. Heights are in figure fraction and were chosen so that panel b, the one
# carrying the architecture, gets the most room -- it is the panel a methods reader looks at.
BANDS = {"a": (0.712, 0.268), "b": (0.386, 0.300), "c": (0.170, 0.192), "d": (0.020, 0.126)}
AX = {}
for k, (y0, h) in BANDS.items():
    AX[k] = blank(fig.add_axes([0.010, y0, 0.980, h]))

C = style.ROLE
INK = style.INK

# =========================================================================== a  THE DEFECT
ax = AX["a"]
panel_letter(ax, "a")
label(ax, 0.055, 0.985, "The variable the benchmark left out", size=8.2, weight="bold", ha="left")

# -- the released clinical table, drawn as a field listing so that the two covariates are seen
#    ADJACENT. That adjacency is the whole point: this is not a variable somebody had to go and
#    find. The fields are named as a reader of the paper would name them rather than as the
#    release spells them, because a figure is manuscript body and carries no code identifiers.
fileb = Box(ax, 0.015, 0.045, 0.395, 0.83, facecolor="#FFFFFF", edgecolor=INK, lw=0.9)
label(ax, 0.2125, 0.815, "Released with the benchmark", size=7.2, weight="bold")
label(ax, 0.2125, 0.752, "the fold assignments and the clinical table", size=6.2, style="italic",
      colour="#5A6273")

COLS = [("case identifier", None), ("age at diagnosis", None), ("sex", None),
        ("histological grade", "competitor"), ("pathologic tumour stage", "clinical"),
        ("molecular subtype", None)]
row_h, top = 0.098, 0.665
rowy = {}
for i, (name, role) in enumerate(COLS):
    y = top - i * row_h
    rowy[name] = y
    if role:
        Box(ax, 0.040, y - row_h * 0.40, 0.345, row_h * 0.80,
            facecolor=C[role], edgecolor="none", alpha=0.30, zorder=1)
        Box(ax, 0.040, y - row_h * 0.40, 0.345, row_h * 0.80,
            facecolor="none", edgecolor=C[role], lw=1.0, zorder=2)
    label(ax, 0.2125, y, name, size=6.4 if len(name) > 20 else 6.8,
          weight="bold" if role else "normal",
          colour=INK if role else "#6B7280")

label(ax, 0.2125, 0.090, "the two sit one row apart", size=6.2, style="italic", colour="#5A6273")

# -- the two branches
up = Box(ax, 0.520, 0.505, 0.286, 0.335, facecolor=C["competitor"], alpha=0.13, edgecolor="none")
Box(ax, 0.520, 0.505, 0.286, 0.335, facecolor="none", edgecolor=C["competitor"], lw=1.0)
label(ax, 0.663, 0.775, "what all %d reported baselines used" % N["census_baselines"],
      size=6.5, weight="bold",
      colour=C["competitor"])
label(ax, 0.663, 0.690, "age + sex + grade", size=7.4)
label(ax, 0.663, 0.578, "C = %.4f" % N["c_grade"], size=9.5, weight="bold")

dn = Box(ax, 0.520, 0.070, 0.286, 0.335, facecolor=C["clinical"], alpha=0.13, edgecolor="none")
Box(ax, 0.520, 0.070, 0.286, 0.335, facecolor="none", edgecolor=C["clinical"], lw=1.0)
label(ax, 0.663, 0.340, "what was in the same file", size=6.5, weight="bold", colour=C["clinical"])
label(ax, 0.663, 0.255, "age + sex + stage", size=7.4)
label(ax, 0.663, 0.143, "C = %.4f" % N["c_stage"], size=9.5, weight="bold")

arrow(ax, (0.392, rowy["histological grade"]), (0.515, 0.672), colour=C["competitor"], rad=-0.12)
arrow(ax, (0.392, rowy["pathologic tumour stage"]), (0.515, 0.238), colour=C["clinical"],
      rad=0.12)

# -- the gap
ax.plot([0.818, 0.845, 0.845, 0.818], [0.672, 0.672, 0.238, 0.238], color=INK, lw=0.8,
        transform=ax.transAxes, clip_on=False)
label(ax, 0.858, 0.455, "+%.4f\non identical\ncases" % N["gap"], size=7.0, weight="bold",
      ha="left")


def glyph(role, gx, gy, gw, gh):
    """A real data thumbnail per input lane, drawn in the panel's own axes coordinates."""
    bb = AX["b"].get_position()
    fx = bb.x0 + gx * bb.width
    fy = bb.y0 + gy * bb.height
    fw, fh = gw * bb.width, gh * bb.height
    a = fig.add_axes([fx, fy, fw, fh]); a.axis("off")
    if role == "slide":
        # Histology must never be stretched to fill a cell: aspect="equal" keeps the image's own
        # ratio and pads within the glyph box instead of distorting the tissue.
        img = plt.imread(os.path.join(HERE, "assets", "blca_slide_thumb.png"))[..., :3]
        a.imshow(img, interpolation="lanczos", aspect="equal")
    elif role == "omics":
        # the 59 pathways the slide arm tracks (BH q < 0.05), recomputed and count-checked in
        # pathway_glyph.py, in the reporting dump's pathway order
        from pathway_glyph import slide_tracked
        dark = set(np.flatnonzero(slide_tracked()).tolist())
        a.set_xlim(0, 55); a.set_ylim(0, 5); a.invert_yaxis()
        for j in range(275):
            r, c_ = divmod(j, 55)
            a.add_patch(mpatches.Rectangle((c_ + 0.12, r + 0.12), 0.76, 0.76,
                                           facecolor=C["omics"] if j in dark else "#F6DDCC",
                                           edgecolor="none"))
    else:
        lev = WG["cohorts"]["blca"]["stage"]["levels"]
        order = [k for k in ("0", "I", "II", "III", "IV") if k in lev]
        tot = sum(lev.values()); x0 = 0.0
        a.set_xlim(0, 1); a.set_ylim(0, 1)
        for j, k in enumerate(order):
            w = lev[k] / tot
            a.add_patch(mpatches.Rectangle((x0, 0.18), w, 0.64, facecolor=style.LADDER[j],
                                           edgecolor="#FFFFFF", linewidth=0.4))
            x0 += w


# =========================================================================== b  ARCHITECTURE
ax = AX["b"]
panel_letter(ax, "b")
label(ax, 0.055, 0.985, "The method: one penalised Cox per modality, combined without a fitted "
      "parameter", size=8.2, weight="bold", ha="left")

LANES = [
    ("slide",    0.760, "Whole-slide image",       "TITAN slide embedding",  "%d dims" % 768),
    ("omics",    0.470, "Bulk transcriptome",      "SurvPath pathway means", "%d pathways" % 275),
    ("clinical", 0.180, "Pathology report",        "age + sex + stage",      "%d columns" % 6),
]
XS = {"in": (0.015, 0.150), "rep": (0.192, 0.168), "cox": (0.386, 0.163), "pct": (0.578, 0.132)}
BH = 0.195

allboxes = []
for role, yc, src, rep, dim in LANES:
    y = yc - BH / 2
    b_in = Box(ax, XS["in"][0], y, XS["in"][1], BH, facecolor=C[role], alpha=0.14, edgecolor="none")
    Box(ax, XS["in"][0], y, XS["in"][1], BH, facecolor="none", edgecolor=C[role], lw=1.0)
    # one line, not forced onto two: at this box width the full lane title fits, and the forced
    # wrap pushed the top line across the box's own top border
    label(ax, XS["in"][0] + XS["in"][1] / 2, yc + 0.058, src, size=5.8)
    glyph(role, XS["in"][0] + 0.014, y + 0.012, XS["in"][1] - 0.028, BH * 0.44)

    b_rep = Box(ax, XS["rep"][0], y, XS["rep"][1], BH, facecolor="#FFFFFF", edgecolor=C[role],
                lw=1.0)
    label(ax, XS["rep"][0] + XS["rep"][1] / 2, yc + 0.030, rep, size=6.3)
    label(ax, XS["rep"][0] + XS["rep"][1] / 2, yc - 0.045, dim, size=6.0, style="italic",
          colour="#5A6273")

    b_cox = Box(ax, XS["cox"][0], y, XS["cox"][1], BH, facecolor=C[role], alpha=0.14,
                edgecolor="none")
    Box(ax, XS["cox"][0], y, XS["cox"][1], BH, facecolor="none", edgecolor=C[role], lw=1.0)
    label(ax, XS["cox"][0] + XS["cox"][1] / 2, yc + 0.030, "ridge Cox", size=6.8, weight="bold")
    label(ax, XS["cox"][0] + XS["cox"][1] / 2, yc - 0.048, "α by inner 3-fold CV", size=5.9,
          colour="#5A6273")

    b_pct = Box(ax, XS["pct"][0], y, XS["pct"][1], BH, facecolor="#FFFFFF", edgecolor="#9AA1AE",
                lw=0.8, ls=(0, (2.4, 1.6)))
    label(ax, XS["pct"][0] + XS["pct"][1] / 2, yc, "percentile\nwithin fold", size=6.2,
          colour="#5A6273")

    for a, b in ((b_in, b_rep), (b_rep, b_cox), (b_cox, b_pct)):
        arrow(ax, a.right, b.left, colour=C[role], lw=0.9)
    arrow(ax, b_pct.right, (0.745, 0.470), colour=C[role], lw=0.9,
          rad=0.0 if role == "omics" else (-0.16 if role == "slide" else 0.16))
    allboxes += [b_in, b_rep, b_cox, b_pct]

# -- the combination, which is the contribution and has nothing in it. Widened (left edge only;
# the right edge that meets the risk box is unchanged) so its three stacked labels -- including the
# one that used to spill below the box's bottom border -- fit as single lines instead of wrapping.
comb = Box(ax, 0.750, 0.352, 0.150, 0.236, facecolor=C["ours"], alpha=0.16, edgecolor="none")
Box(ax, 0.750, 0.352, 0.150, 0.236, facecolor="none", edgecolor=C["ours"], lw=1.2)
label(ax, 0.825, 0.512, "equal-weight", size=6.5, weight="bold", colour=C["ours"])
label(ax, 0.825, 0.455, "rank average", size=6.5, weight="bold", colour=C["ours"])
label(ax, 0.825, 0.390, "%d fitted parameters" % N["p_fitted"], size=6.0, colour="#5A6273")

risk = Box(ax, 0.926, 0.383, 0.062, 0.174, facecolor="#FFFFFF", edgecolor=INK, lw=1.0)
label(ax, 0.957, 0.470, "risk\nscore", size=6.5, weight="bold")
arrow(ax, comb.right, risk.left, colour=C["ours"], lw=1.0)

# -- the parameter budget, which is the point of the panel
bracket(ax, XS["cox"][0], XS["cox"][0] + XS["cox"][1], 0.072,
        "%s Cox coefficients in total\nagainst %s and %s for the two competitors whose code we "
        "could rerun" % ("{:,}".format(N["p_ours"]), "{:,}".format(N["p_sp"]),
                         "{:,}".format(N["p_pibd"])),
        size=6.3, above=False)


# =========================================================================== c  THE PROTOCOL
ax = AX["c"]
panel_letter(ax, "c")
label(ax, 0.055, 0.985, "How the number was produced, and what was fixed before it existed",
      size=8.2, weight="bold", ha="left")

# Plain Arial text, not mathtext, for every symbol below: mathtext's operator glyphs (here "in"
# and "times") fall back to Cmsy10 regardless of the custom Arial fontset (Gate G47), while the
# Unicode characters themselves -- alpha, times, plus-minus -- are all in Arial's own character set.
STEPS = [
    (0.015, 0.150, "5-fold case-ID splits\nas released", None),
    (0.196, 0.163, "inner 3-fold CV picks\n\u03b1 from {1, 8, 64, 512, 4096}", None),
    (0.390, 0.150, "out-of-fold score\nfor all %d cases" % N["n_cases"], None),
    (0.571, 0.150, "\u00d7 5 seeds\n1\u20134 never run before", None),
    (0.752, 0.233, "primary = mean concordance\n%.4f \u00b1 %.4f" % (N["primary"], N["sd"]), "ours"),
]
prev = None
for x, w, txt, role in STEPS:
    fc = C[role] if role else "#FFFFFF"
    b = Box(ax, x, 0.335, w, 0.300, facecolor=fc, alpha=0.16 if role else 1.0, edgecolor="none")
    Box(ax, x, 0.335, w, 0.300, facecolor="none", edgecolor=C[role] if role else INK,
        lw=1.1 if role else 0.9)
    label(ax, x + w / 2, 0.485, txt, size=6.4, weight="bold" if role else "normal")
    if prev is not None:
        arrow(ax, prev.right, b.left, lw=0.9)
    prev = b

bracket(ax, 0.015, 0.721, 0.700,
        "fixed before the run: arm, estimator, $\\alpha$ grid, combination rule, comparators, "
        "seeds, decision rule", size=6.3, above=True)

label(ax, 0.015, 0.145, "five leakage checks ran first, and the two built to reject did so:",
      size=6.3, ha="left", weight="bold")
label(ax, 0.015, 0.055,
      "per-fold planted label $%+.4f$  ·  permuted feature-to-case assignment $%.4f$  ·  "
      "three others" % (N["planted_label_difference"], N["shuffled_alignment_c"]),
      size=6.1, ha="left", colour="#5A6273")

# =========================================================================== d  THE MECHANISM
ax = AX["d"]
panel_letter(ax, "d")
label(ax, 0.055, 0.985, "Why grade separates few patients in this disease", size=8.2,
      weight="bold", ha="left")

LAD = style.LADDER
BARS = [("tumour grade", [("G1", GRADE_LEVELS["G1"]), ("G2", GRADE_LEVELS["G2"])],
         N["g_ent"], N["g_modal"], N["g_levels"], 0.520),
        ("pathologic stage", [(k, STAGE_LEVELS[k]) for k in ("0", "I", "II", "III", "IV")],
         N["s_ent"], N["s_modal"], N["s_levels"], 0.150)]
for name, levels, ent, modal, nlev, y in BARS:
    total = sum(v for _, v in levels)
    x = 0.180
    for i, (lab, v) in enumerate(levels):
        frac = v / total * 0.480
        Box(ax, x, y, frac, 0.245, facecolor=LAD[i % len(LAD)], edgecolor="#FFFFFF", lw=0.8)
        if frac > 0.045:
            label(ax, x + frac / 2, y + 0.1225, lab, size=6.2, colour=INK)
        x += frac
    label(ax, 0.170, y + 0.1225, name, size=7.0, ha="right", weight="bold")
    # a level too narrow to carry its label is still a level, and saying "5 levels" beside three
    # visible blocks reads as an error unless the figure says which is which
    thin = sum(1 for _, v in levels if v / total * 0.480 <= 0.045)
    label(ax, 0.675, y + 0.175,
          "%d levels%s · %.1f%% in one"
          % (nlev, (", %d too narrow to label" % thin) if thin else "", modal * 100),
          size=6.4, ha="left")
    label(ax, 0.675, y + 0.072, "normalised entropy %.3f" % ent, size=6.4, ha="left",
          weight="bold")

# A figure is manuscript body, so the author's no-code-identifier rule applies to the artwork too
# and nothing else checks it: the prose gate reads .tex files only.
for _name, _ in COLS:
    assert "_" not in _name and "." not in _name, (
        "a field label in the artwork is a code identifier: %r" % _name)

style.save(fig, style.out(HERE, "fig1_overview.pdf"))
print("wrote fig1_overview.pdf   %.2f x %.2f in" % (FIGW, FIGH))
