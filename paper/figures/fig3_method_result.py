#!/usr/bin/env python3
"""Figure 3 -- the method's result, and the two things about it that are not flattering.

Six panels, and the figure is deliberately built so that a reader cannot take the win without also
taking its limits. Panels a and f are the claim; b, c, d and e are what it rests on; and a carries
the correction that kills four of six comparisons rather than putting it in a footnote.

  a  the whole prespecified comparison family under Holm. TWO survive of six, and the two
     input-parity comparisons are among the four that do not. Failing comparisons are drawn OPEN
     and are not omitted -- a forest plot showing only the significant arms is the version of this
     panel that would be dishonest.
  b  Kaplan-Meier by risk tertile with a risk table, because a late tail resting on five patients
     looks identical to one resting on fifty unless the numbers are printed.
  c  the incumbent as the benchmark publishes it, beside the same predictions given the same
     clinical block. Its tertile medians are non-monotone on the left and monotone on the right.
  d  time-dependent AUC across horizons, with n-at-risk along the top for the same reason as b.
  e  all seven modality subsets, so "which part matters" is answered by the figure and not by a
     sentence.
  f  cost against benefit, on a log axis spanning four orders of magnitude. Two competitor points
     are MEASURED -- each printed by its own trainer -- and the other seven sit on a labelled
     "not reported" band rather than being invented.
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import style                                   # noqa: E402

RES = os.path.abspath(os.path.join(HERE, "..", "..", "experiments",
                                   "20260817-blca-confirm", "results"))
load = lambda n: json.load(open(os.path.join(RES, n)))

PAR = load("pibd-parity.json")
ABL = load("ablation-generalisation.json")
PC = load("parameter-counts.json")
PUB = load("published-benchmark-table.json")
FSD = json.load(open(os.path.join(HERE, "figure-source-data.json")))

style.apply()
C, INK, RISK = style.ROLE, style.INK, style.RISK

FIGW, FIGH = style.width(6.785), 6.45
fig = plt.figure(figsize=(FIGW, FIGH))
axa = fig.add_axes([0.300, 0.828, 0.395, 0.140])
axb = fig.add_axes([0.100, 0.632, 0.355, 0.122])
axc1 = fig.add_axes([0.600, 0.632, 0.170, 0.122])
axc2 = fig.add_axes([0.800, 0.632, 0.170, 0.122])
axd = fig.add_axes([0.100, 0.322, 0.355, 0.122])
axe = fig.add_axes([0.620, 0.322, 0.268, 0.122])
axf = fig.add_axes([0.115, 0.070, 0.630, 0.168])


# =========================================================================== a  forest
FAM = PAR["holm_family_of_six"]["comparisons"]
SHORT = {"ours vs slide alone": "vs slide alone",
         "ours vs omics alone": "vs omics alone",
         "ours vs clinical alone": "vs clinical alone",
         "ours vs ours-without-clinical": "vs ours without clinical",
         "ours vs SurvPath+clinical (seed-matched)": "vs SurvPath + clinical",
         "ours vs PIBD+clinical (best-val checkpoint)": "vs PIBD + clinical"}
items = sorted(FAM.items(), key=lambda kv: kv[1]["delta"])
assert len(items) == 6, "the family is six comparisons; got %d" % len(items)
y = np.arange(len(items))
for yi, (k, v) in zip(y, items):
    lo, hi = v["ci95"]
    surv = v["survives_at_0.05"]
    axa.plot([lo, hi], [yi, yi], color=INK, lw=1.0, solid_capstyle="round", zorder=2)
    axa.scatter([v["delta"]], [yi], s=30, zorder=3,
                color=C["ours"] if surv else "#FFFFFF",
                edgecolor=C["ours"], linewidth=1.1)
    axa.text(0.128, yi, "%.3f / %.3f" % (v["p_raw"], v["p_holm"]), fontsize=5.9, va="center",
             fontweight="bold" if surv else "normal", color=INK if surv else "#6B7280")
axa.axvline(0, color="#9AA1AE", lw=0.8, zorder=1)
axa.set_yticks(y)
axa.set_yticklabels([SHORT[k] for k, _ in items], fontsize=6.4)
axa.set_xlim(-0.03, 0.125)
axa.set_xlabel("difference in concordance, cases resampled", fontsize=6.4)
axa.set_ylim(-1.0, len(items) - 0.1)
axa.text(0.128, len(items) - 0.45, "$p$ / Holm $q$", fontsize=6.0, va="center", fontweight="bold")
n_surv = sum(1 for _, v in items if v["survives_at_0.05"])
axa.set_title("six prespecified comparisons under Holm (%d of %d survive)"
              % (n_surv, len(items)), fontsize=6.4, pad=6)
axa.text(-0.028, -0.85, "filled = survives Holm at $\\alpha$=0.05,  open = does not", fontsize=5.9,
         color="#5A6273", ha="left", va="center")

# =========================================================================== b, c  Kaplan-Meier
KM = FSD["F1_kaplan_meier"]
GRID = KM["ours"]["groups"]["low"]["at_risk_grid"]


def draw_km(ax, arm, title, risk_table=True, ylab=True):
    d = KM[arm]
    for lab in ("low", "middle", "high"):
        g = d["groups"][lab]
        ax.step(g["times"], g["survival"], where="post", color=RISK[lab], lw=1.3,
                label="%s (%d ev/%d)" % (lab, g["events"], g["n"]))
    ax.set_xlim(0, 96)
    ax.set_ylim(0, 1.02)
    ax.set_xticks([0, 24, 48, 72, 96])
    ax.set_xlabel("months", fontsize=6.4)
    if ylab:
        ax.set_ylabel("disease-specific survival", fontsize=6.4)
    else:
        ax.set_yticklabels([])
    ax.set_title(title, fontsize=6.3, pad=3)
    ax.text(0.97, 0.95, "log-rank $\\chi^2$ = %.1f\n$p$ = %s" % (d["logrank_chi2"], d["logrank_p"]),
            transform=ax.transAxes, ha="right", va="top", fontsize=6.0, linespacing=1.3)
    if risk_table:
        ax.legend(fontsize=5.9, loc="lower left", handlelength=1.1, handletextpad=0.4,
                  borderpad=0.15, labelspacing=0.25)
        for j, lab in enumerate(("low", "middle", "high")):
            ar = d["groups"][lab]["at_risk"]
            for x, v in zip(GRID, ar):
                ax.annotate(str(v), (x, 0), xycoords=("data", "axes fraction"),
                            textcoords="offset points", xytext=(0, -32 - j * 8),
                            ha="center", fontsize=5.4, color=RISK[lab], annotation_clip=False)
        ax.annotate("at risk", (0, 0), xycoords=("axes fraction", "axes fraction"),
                    textcoords="offset points", xytext=(-34, -32), ha="left", fontsize=5.4,
                    color=INK, annotation_clip=False)
    return d


d_ours = draw_km(axb, "ours", "ours, by risk tertile")
med = [KM["ours"]["groups"][l]["median_months"] for l in ("low", "middle", "high")]
axb.annotate("median %.1f mo" % med[2], (med[2], 0.5), textcoords="offset points",
             xytext=(9, 4), fontsize=5.9, color=RISK["high"])
axb.text(0.97, 0.60, "low and middle:\nmedian not reached", transform=axb.transAxes, ha="right",
         va="top", fontsize=5.7, color="#5A6273", linespacing=1.3)

for ax, arm, ttl, yl in ((axc1, "survpath_as_released", "incumbent,\nas released", True),
                         (axc2, "survpath_plus_clinical", "the same, given\nthe same stage", False)):
    dd = draw_km(ax, arm, ttl, risk_table=False, ylab=yl)
    ax.set_ylabel("disease-specific survival" if yl else "", fontsize=6.6)
    m = [dd["groups"][l]["median_months"] for l in ("low", "middle", "high")]
    ok = all(v is not None for v in m) and m[0] >= m[1] >= m[2]
    ax.text(0.5, 0.22, ("medians %s" % " / ".join("%.0f" % v if v else "n.r." for v in m))
            + ("\nmonotone" if ok else "\nNOT monotone"),
            transform=ax.transAxes, ha="center", va="top", fontsize=5.8,
            color=INK if ok else C["competitor"], fontweight="bold" if not ok else "normal",
            linespacing=1.3)


# =========================================================================== d  time-dependent AUC
AUC = FSD["F2_time_dependent_auc"]
SERIES = [("ours", "ours", C["ours"], "-", "o"),
          ("pibd_plus_clinical", "PIBD + clinical", C["competitor"], "--", "s"),
          ("survpath_plus_clinical", "SurvPath + clinical", C["competitor"], ":", "^"),
          ("clinical_alone", "clinical alone", C["clinical"], "-.", "D"),
          ("slide", "slide alone", C["slide"], (0, (4, 1, 1, 1)), "v"),
          ("omics", "omics alone", C["omics"], (0, (1, 1)), "P")]
h = AUC["ours"]["horizons"]
for key, lab, col, ls, mk in SERIES:
    axd.plot(h, AUC[key]["auc"], color=col, ls=ls, marker=mk, ms=2.6, lw=1.0, label=lab,
             markevery=3)
axd.set_xlabel("horizon (months)", fontsize=6.4)
axd.set_ylabel("IPCW time-dependent AUC", fontsize=6.4)
axd.set_ylim(0.50, 0.83)
axd.legend(fontsize=5.6, ncol=2, loc="lower left", handlelength=2.0, handletextpad=0.4,
           columnspacing=0.9, borderpad=0.15, labelspacing=0.22)
# n at risk along the top, for the same reason the KM panel carries a risk table
for x, v in list(zip(h, AUC["_at_risk_at_horizons"]))[::4]:
    axd.annotate(str(v), (x, 1.0), xycoords=("data", "axes fraction"),
                 textcoords="offset points", xytext=(0, 3), ha="center", fontsize=5.3,
                 color="#6B7280", annotation_clip=False)
axd.annotate("n at risk", (0.0, 1.0), xycoords="axes fraction", textcoords="offset points",
             xytext=(2, 9), ha="left", fontsize=5.3, color="#6B7280", annotation_clip=False)

# =========================================================================== e  ablation
AR = ABL["B_ablation_bladder"]["arms"]
PV = ABL["B_ablation_bladder"]["paired_vs_full"]
NICE = {"wsi": ("slide", ("slide",)), "omics": ("omics", ("omics",)),
        "clin": ("clinical", ("clinical",)),
        "wsi + omics": ("slide + omics", ("slide", "omics")),
        "wsi + clin": ("slide + clinical", ("slide", "clinical")),
        "omics + clin": ("omics + clinical", ("omics", "clinical")),
        "wsi + omics + clin": ("all three", ("slide", "omics", "clinical"))}
rows = sorted(AR.items(), key=lambda kv: kv[1])
ye = np.arange(len(rows))
for yi, (k, v) in zip(ye, rows):
    parts = NICE[k][1]
    full = len(parts) == 3
    # A DOT, not a bar: this axis starts at 0.548, and a bar encodes its value from zero, so a bar
    # here overstated every difference on the panel (found 2026-09-11 in the 180 mm rebuild).
    axe.plot([0.593, v], [yi, yi], color="#E3E6EB", lw=0.6, zorder=1)
    axe.scatter([v], [yi], s=24 if full else 16, zorder=5, edgecolor="none",
                color=C["ours"] if full else "#8C93A1")
    # a mark of one colour per subset would need seven colours for an unordered variable; instead
    # the modalities present are shown as dots, which is what the panel is actually about
    for j, role in enumerate(("slide", "omics", "clinical")):
        axe.scatter([0.5595 + j * 0.0095], [yi], s=9, zorder=4,
                    color=C[role] if role in parts else "#FFFFFF",
                    edgecolor=C[role] if role in parts else "#C9CFD8", linewidth=0.7)
    axe.text(v + 0.0035, yi, "%.4f" % v, va="center", ha="left", fontsize=5.8,
             color=INK, fontweight="bold" if full else "normal")
    if k in PV:
        axe.text(0.7315, yi, "$-$%.3f  ($p$=%.2f)" % (abs(PV[k]["mean"]), PV[k]["p_two_sided"]),
                 va="center", ha="left", fontsize=5.2, color="#5A6273")
axe.set_yticks(ye)
axe.set_yticklabels([NICE[k][0] for k, _ in rows], fontsize=6.2)
axe.set_xlim(0.548, 0.730)
axe.set_xticks([0.60, 0.65, 0.70])
axe.set_xlabel("concordance", fontsize=6.4)
axe.set_title("all seven modality subsets", fontsize=6.4, pad=6, loc="left")
axe.text(1.0, 1.06, "leave-one-out", transform=axe.transAxes, fontsize=5.6, ha="right",
         color="#5A6273")


# =========================================================================== f  cost vs benefit
# Two axes, because the x variable does not EXIST for seven of the ten entries. Putting them at an
# invented position, or omitting them, are both worse than giving them their own strip and saying
# what the strip means.
axfs = fig.add_axes([0.800, 0.070, 0.075, 0.168])

MEASURED = [("ours", PC["ours"]["parameters"],
             load("amendment-A1-clinical-provenance.json")["primary"]["value"], C["ours"]),
            ("SurvPath", PC["survpath"]["parameters"],
             [e for e in PUB["entries"] if e["method"] == "SurvPath"][0]["cindex"],
             C["competitor"]),
            ("PIBD", PC["pibd"]["parameters"],
             [e for e in PUB["entries"] if e["method"] == "PIBD"][0]["cindex"], C["competitor"])]
for name, p, c, col in MEASURED:
    axf.scatter([p], [c], s=46, color=col, edgecolor="none", zorder=3)
    axf.annotate("%s\n%s params" % (name, "{:,}".format(p)), (p, c), textcoords="offset points",
                 xytext=(0, 9), ha="center", fontsize=5.9, color=col, fontweight="bold",
                 linespacing=1.25)
axf.set_xscale("log")
axf.set_xlim(3e2, 1.2e8)
axf.set_ylim(0.575, 0.755)
axf.set_xlabel("trainable parameters (log scale), as printed by each method's own trainer",
               fontsize=6.4)
axf.set_ylabel("concordance", fontsize=6.4)
ratio = PC["pibd"]["ratio_to_ours"]
axf.annotate("", xy=(PC["pibd"]["parameters"], 0.605), xytext=(PC["ours"]["parameters"], 0.605),
             arrowprops=dict(arrowstyle="<->", color="#9AA1AE", lw=0.8))
axf.text(np.sqrt(PC["ours"]["parameters"] * PC["pibd"]["parameters"]), 0.590,
         r"$\times$%s" % "{:,}".format(ratio), ha="center", fontsize=6.4, color="#5A6273",
         fontweight="bold")

others = sorted([e for e in PUB["entries"] if e["method"] not in ("SurvPath", "PIBD")],
                key=lambda e: -e["cindex"])
axfs.scatter([0.35] * len(others), [e["cindex"] for e in others], s=20,
             color=C["published"], edgecolor="none", zorder=3)
# 0.679 and 0.677 are two thousandths apart and their labels are not: push them apart and draw a
# leader to the point, rather than letting two names overprint and calling it a plot
YLO, YHI, MINSEP = 0.575, 0.755, 0.0135
ylab = []
for e in others:
    yy = e["cindex"]
    if ylab and ylab[-1] - yy < MINSEP:
        yy = ylab[-1] - MINSEP
    ylab.append(yy)
assert min(ylab) > YLO, "the label stack ran off the bottom of the strip"
for e, yy in zip(others, ylab):
    axfs.plot([0.42, 0.62], [e["cindex"], yy], color="#C9CFD8", lw=0.5, zorder=2,
              clip_on=False)
    axfs.annotate("%s %.3f" % (e["method"], e["cindex"]), (0.66, yy), fontsize=5.3, va="center",
                  ha="left", color="#6B7280", annotation_clip=False)
axfs.set_xlim(0, 1)
axfs.set_ylim(YLO, YHI)
axfs.set_xticks([])
axfs.set_yticklabels([])
axfs.set_title("count not\nreported", fontsize=6.0, pad=2, color="#6B7280")
axfs.spines["left"].set_visible(False)
assert len(MEASURED) + len(others) == len(PUB["entries"]) + 1, (
    "every published entry must be on one axis or the other, plus ours")

for tag, xx, yy in (("a", 0.012, 0.995), ("b", 0.012, 0.760), ("c", 0.525, 0.760),
                    ("d", 0.012, 0.508), ("e", 0.560, 0.508), ("f", 0.012, 0.290)):
    fig.text(xx, yy, tag, fontsize=7.5, fontweight="bold", color=INK, ha="left", va="top")

style.save(fig, style.out(HERE, "fig3_method_result.pdf"))
print("wrote fig3_method_result.pdf   %.2f x %.2f in" % (FIGW, FIGH))
print("  a: %d of %d comparisons survive Holm" % (n_surv, len(items)))
print("  f: %d measured points, %d on the not-reported strip" % (len(MEASURED), len(others)))
