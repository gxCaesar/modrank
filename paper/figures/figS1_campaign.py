#!/usr/bin/env python3
"""Supplementary Figure S1 -- the development campaign, its generalisation and its representations.

Three thin two-panel figures used to carry this. Measured against the rest of the package they were
the loosest displays in it: the largest empty horizontal band ran to 14% of one figure's height and
another used only 83% of its width. Six panels on one page, three across, is the same content at
roughly twice the density, and it puts the three questions a reader asks about the campaign in one
place instead of three.

  a  the seven components against the bar each was pre-registered to clear. None cleared it.
  b  the four that carry a matched control, with the control drawn beside the effect.
  c  ModRank against the best published entry per study. Bladder wins, four of five lose.
  d  the one-column swap, grade to stage, in the same five studies. Five of five.
  e  seven slide encoders scored as single arms on identical cases and folds.
  f  TITAN alone against a pool of all seven, paired on cases. The pool is significantly worse.
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import style                                        # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "experiments", "20260817-blca-confirm", "results")
CS = json.load(open(os.path.join(RES, "component-slate.json")))
ABL = json.load(open(os.path.join(RES, "ablation-generalisation.json")))
WG = json.load(open(os.path.join(RES, "why-grade-fails.json")))
M = json.load(open(os.path.join(ROOT, "development", "s5-results", "multi.json")))

style.apply()
C, INK, LAD = style.ROLE, style.INK, style.LADDER

FIGW, FIGH = 6.30, 6.45
fig = plt.figure(figsize=(FIGW, FIGH))
# TWO across, not three. The panels carry long row labels (component names, study names, encoder
# names) and three columns on a 6.3 in canvas leaves 2.1 in each, of which the labels alone want
# 1.3. Two columns at 3.15 in is the tightest packing these six panels actually fit into.
RH = 0.208
COL = {"L": 0.245, "R": 0.735}
ROW = {0: 0.740, 1: 0.428, 2: 0.108}
axa = fig.add_axes([COL["L"] + 0.012, ROW[0], 0.218, RH])
axb = fig.add_axes([COL["R"], ROW[0], 0.180, RH])
axc = fig.add_axes([COL["L"], ROW[1], 0.230, RH])
axd = fig.add_axes([COL["R"], ROW[1], 0.230, RH])
axe = fig.add_axes([COL["L"], ROW[2], 0.230, RH])
axf = fig.add_axes([COL["R"], ROW[2], 0.120, RH])

comps_1 = list(reversed(CS["components"]))
y_1 = np.arange(len(comps_1))
BAR_1 = CS["bar"]

# =========================================================================== a  predicted vs got
for yi, c in zip(y_1, comps_1):
    axa.plot([c["predicted"], c["observed"]], [yi, yi], color="#C9CFD8", lw=1.2, zorder=1,
             solid_capstyle="round")
    axa.scatter([c["predicted"]], [yi], s=22, facecolor="#FFFFFF", edgecolor="#6B7280",
                linewidth=0.9, zorder=3)
    axa.scatter([c["observed"]], [yi], s=26, color=C["competitor"], zorder=4, edgecolor="none")
axa.axvline(0, color=INK, lw=0.8, zorder=0)
axa.axvspan(-BAR_1, BAR_1, color=C["published"], alpha=0.18, lw=0, zorder=0)
axa.set_yticks(y_1)
axa.set_yticklabels(["%s  %s" % (c["id"], c["label"]) for c in comps_1], fontsize=6.2)
axa.set_xlim(-0.075, 0.055)
axa.set_xlabel("effect on concordance", fontsize=6.4)
axa.scatter([], [], s=22, facecolor="#FFFFFF", edgecolor="#6B7280", linewidth=0.9,
            label="predicted, fixed before the run")
axa.scatter([], [], s=26, color=C["competitor"], label="observed")
axa.legend(fontsize=5.9, loc="upper center", bbox_to_anchor=(0.5, -0.155), ncol=2,
           handletextpad=0.4, columnspacing=1.4, borderpad=0.15)
axa.text(0.0, 0.35, "$\\pm$%.4f bar" % BAR_1, fontsize=5.9, color="#5A6273",
         ha="center", va="center", rotation=90)
n_clear_1 = sum(1 for c in comps_1 if c["observed"] > BAR_1)
axa.set_title("components against their bar (%d of %d)" % (n_clear_1, len(comps_1)), fontsize=6.4, pad=3)

# =========================================================================== b  the controls
withc_1 = [c for c in comps_1 if c["control"] is not None]
yb_1 = np.arange(len(withc_1))
for yi, c in zip(yb_1, withc_1):
    axb.barh([yi + 0.18], [c["observed"]], height=0.32, color=C["ours"], edgecolor="none")
    axb.barh([yi - 0.18], [c["control"]], height=0.32, color=C["competitor"], edgecolor="none")
axb.axvline(0, color=INK, lw=0.8)
axb.set_yticks(yb_1)
axb.set_yticklabels([c["id"] for c in withc_1], fontsize=6.2)
axb.set_xlabel("effect", fontsize=6.4)
axb.set_xlim(-0.048, 0.030)   # room on the right for the void annotation
import matplotlib.patches as mpatches
axb.legend(handles=[mpatches.Patch(color=C["ours"], label="component"),
                    mpatches.Patch(color=C["competitor"], label="matched control")],
           fontsize=5.8, loc="upper center", bbox_to_anchor=(0.5, -0.155), ncol=2,
           handlelength=1.0, handletextpad=0.4, columnspacing=1.0, borderpad=0.15)
# a control above a NEGATIVE effect is a component that simply failed, not a void. Only a
# positive effect that its own control beats is void, and `status` is what the gate assigned.
void_on_control = [c for c in withc_1
                   if c["status"] == "void" and c["observed"] > 0 and c["control"] > c["observed"]]
assert void_on_control, "panel b exists to show a control beating its component; none does"
for c in void_on_control:
    yi = list(withc_1).index(c)
    axb.annotate("the control\nwins $\\Rightarrow$ void", (c["control"], yi),
                 textcoords="offset points", xytext=(5, 13), fontsize=5.7, ha="left",
                 va="center", color=C["competitor"], fontweight="bold", linespacing=1.25)
axb.set_title("matched controls (%d)" % len(withc_1), fontsize=6.4, pad=3)

NAME_2 = {"blca": "BLCA (bladder)", "brca": "BRCA", "coadread": "COADREAD", "hnsc": "HNSC",
        "stad": "STAD"}
GEN_2 = ABL["A_generalisation"]["cohorts"]
order_2 = [k for k in ("blca", "brca", "coadread", "hnsc", "stad") if k in GEN_2]
assert len(order_2) == 5, "the benchmark covers five studies; found %d" % len(order_2)

# =========================================================================== a  does it travel
y_2 = np.arange(len(order_2))[::-1]
won_2 = 0
for yi, k in zip(y_2, order_2):
    g = GEN_2[k]
    ours, best = g["primary_mean_over_seeds"], g["best_published_verified_folds"]
    beat = ours > best
    won_2 += beat
    axc.plot([best, ours], [yi, yi], color="#C9CFD8", lw=1.3, zorder=1, solid_capstyle="round")
    axc.scatter([best], [yi], s=24, color=C["published"], zorder=3, edgecolor="none")
    axc.scatter([ours], [yi], s=28, zorder=4, edgecolor=C["ours"], linewidth=1.1,
                color=C["ours"] if beat else "#FFFFFF")
    axc.text(max(ours, best) + 0.006, yi, "%+.3f" % (ours - best), va="center", fontsize=5.9,
             color=INK if beat else C["competitor"], fontweight="bold" if beat else "normal")
axc.set_yticks(y_2)
axc.set_yticklabels([NAME_2[k] for k in order_2], fontsize=6.4)
axc.set_xlim(0.55, 0.79)
axc.set_xlabel("concordance", fontsize=6.4)
axc.scatter([], [], s=24, color=C["published"], label="best published, verified folds")
axc.scatter([], [], s=28, color=C["ours"], label="ours (filled = we win)")
axc.legend(fontsize=5.9, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2,
           handletextpad=0.4, columnspacing=1.2, borderpad=0.15)
axc.set_title("method, per study (%d of %d)" % (won_2, len(order_2)), fontsize=6.4, pad=3,
              color=C["competitor"] if won_2 < len(order_2) else INK)

# =========================================================================== b  does the defect
swon_2 = 0
for yi, k in zip(y_2, order_2):
    co = WG["cohorts"][k]
    g, s = co["c_age_sex_grade"], co["c_age_sex_stage"]
    swon_2 += s > g
    axd.plot([g, s], [yi, yi], color="#C9CFD8", lw=1.3, zorder=1, solid_capstyle="round")
    axd.scatter([g], [yi], s=24, color=C["competitor"], zorder=3, edgecolor="none")
    axd.scatter([s], [yi], s=24, color=C["clinical"], zorder=3, edgecolor="none")
    axd.text(s + 0.008, yi, "+%.3f" % (s - g), va="center", fontsize=5.9, color=INK)
axd.set_yticks(y_2)
axd.set_yticklabels([])
axd.set_xlim(0.43, 0.79)
axd.set_xlabel("concordance", fontsize=6.4)
axd.scatter([], [], s=24, color=C["competitor"], label="age + sex + grade")
axd.scatter([], [], s=24, color=C["clinical"], label="age + sex + stage")
axd.legend(fontsize=5.9, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2,
           handletextpad=0.4, columnspacing=1.2, borderpad=0.15)
axd.set_title("grade to stage, per study (%d of %d)" % (swon_2, len(order_2)), fontsize=6.4, pad=3)

PRETTY_3 = {"wsi_titan": "TITAN", "wsi_provgigapath": "Prov-GigaPath",
          "wsi_gigassl_gigapath": "GigaSSL / GigaPath",
          "wsi_gigassl_ctranspath": "GigaSSL / CTransPath",
          "wsi_gigassl_optimus": "GigaSSL / Optimus",
          "wsi_gigassl_phikon": "GigaSSL / Phikon",
          "wsi_uni_meanpool": "UNI, mean-pooled"}
arms_3 = [(PRETTY_3[k], v, M["cases_covered_per_encoder"][k], M["blocks"][k])
        for k, v in M["single_arms"].items() if k in PRETTY_3]
arms_3.sort(key=lambda r: r[1])
assert len(arms_3) == 7, "seven slide encoders were scored; found %d" % len(arms_3)

# =========================================================================== a  single arms_3
y_3 = np.arange(len(arms_3))
for yi, (name, v, n, dim) in enumerate(arms_3):
    lead = name == "TITAN"
    axe.barh([yi], [v], color=C["slide"] if lead else "#C9CFD8", height=0.66, edgecolor="none")
    axe.text(v + 0.004, yi, "%.4f" % v, va="center", fontsize=5.9,
             fontweight="bold" if lead else "normal", color=INK if lead else "#6B7280")
    # inside the bar only when the bar is long enough to hold it; otherwise after the value, or
    # the two shortest bars_3 overprint their own numbers
    inside = v > 0.545
    axe.text(0.503 if inside else v + 0.036, yi,
             "%d dims" % dim + ("" if n == 359 else ", %d cases" % n), va="center",
             fontsize=5.3, color=("#FFFFFF" if lead else "#6B7280") if inside else "#9AA1AE")
axe.axvline(0.5, color=INK, lw=0.8, zorder=3)
axe.set_yticks(y_3)
axe.set_yticklabels([a[0] for a in arms_3], fontsize=6.2)
axe.set_xlim(0.50, 0.70)
axe.set_xlabel("concordance, slide arm alone", fontsize=6.4)
axe.set_title("seven slide encoders", fontsize=6.4, pad=3)

# =========================================================================== b  the ensemble
VS_3 = M["vs_single_encoder_baseline"]
pb_3 = VS_3["paired_bootstrap"]
bars_3 = [("TITAN\nalone", VS_3["titan_only_arm"], C["slide"]),
        ("all seven\nencoders", VS_3["six_encoder_arm"], "#C9CFD8")]
xs_3 = np.arange(len(bars_3))
_floor_3 = min(b[1] for b in bars_3) - 0.012
for x, b in zip(xs_3, bars_3):
    axf.plot([x, x], [_floor_3, b[1]], color=b[2], lw=1.4, zorder=2, solid_capstyle="round")
    axf.scatter([x], [b[1]], s=38, color=b[2], zorder=3, edgecolor="none")
    axf.text(x, b[1] + 0.003, "%.4f" % b[1], ha="center", va="bottom", fontsize=6.0,
             fontweight="bold")
axf.set_xticks(xs_3)
axf.set_xticklabels([b[0] for b in bars_3], fontsize=6.0)
axf.set_xlim(-0.62, 1.62)
axf.set_ylim(0.672, 0.748)
axf.set_ylabel("concordance", fontsize=6.4)
assert pb_3["mean"] < 0 and pb_3["p_two_sided"] < 0.05, (
    "this panel exists to show the ensemble losing SIGNIFICANTLY; it does not")
axf.annotate("", xy=(1, VS_3["six_encoder_arm"] + 0.006), xytext=(0, VS_3["titan_only_arm"] + 0.006),
             arrowprops=dict(arrowstyle="-|>", color=C["competitor"], lw=1.0, mutation_scale=8))
# the delta label sits to the RIGHT of the arrow, not above it: va="bottom" at the top of the
# data range pushes text out of the axes and into the title.
axf.text(0.58, VS_3["titan_only_arm"] + 0.006, "%.4f\n$p$ = %.3f" % (pb_3["mean"], pb_3["p_two_sided"]),
         ha="left", va="center", fontsize=6.0, color=C["competitor"], fontweight="bold",
         linespacing=1.3)
axf.set_title("TITAN against all seven", fontsize=6.4, pad=3, color=C["competitor"])


for tag, ax, dx in (("a", axa, -0.66), ("b", axb, -0.22), ("c", axc, -0.52),
                    ("d", axd, -0.10), ("e", axe, -0.62), ("f", axf, -0.55)):
    ax.text(dx, 1.26, tag, transform=ax.transAxes, fontsize=7.5, fontweight="bold", color=INK,
            ha="left", va="top")

style.save(fig, os.path.join(HERE, "figS1_campaign.pdf"))
print("wrote figS1_campaign.pdf   %.2f x %.2f in" % (FIGW, FIGH))
print("  %d components, %d studies, %d encoders" % (len(comps_1), len(order_2), len(arms_3)))
