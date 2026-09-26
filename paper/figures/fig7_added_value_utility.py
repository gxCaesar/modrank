#!/usr/bin/env python3
"""Figure 7 (NC build): how much the clinical reference moves the added value, and what the model is
worth at the bedside.

  a  added value over the grade-based and the stage-based reference, three constructions, 95% CI
  b  the inflation D, TCGA-BLCA (three constructions) and two independent GEO cohorts
  c  grouped calibration at 24 months, fitted inside the training folds, ModRank and clinical
  d  net benefit at 24 months: ModRank, the clinical model, treat all, treat none

NO NUMERAL IS TYPED. Every value is read from a committed result file:
experiments/20260911-blca-posthoc/results/{unified-fusion-and-added-value,calibration-and-decision-curve}.json
and experiments/20260911-geo-external/results/geo-external.json. The plotted values are written to
fig7_added_value_utility.source.json beside the figure.

Authored at 180 mm (7.087 in), the Nature double-column width; dot and interval marks only, since
none of these axes can start at zero without hiding the differences; every categorical mark carries
a marker shape as well as a colour.
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

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PH = os.path.join(ROOT, "experiments", "20260911-blca-posthoc", "results")
AV = json.load(open(os.path.join(PH, "unified-fusion-and-added-value.json")))["added_value_inflation"]["constructions"]
CAL = json.load(open(os.path.join(PH, "calibration-and-decision-curve.json")))
GEO = json.load(open(os.path.join(ROOT, "experiments", "20260911-geo-external", "results",
                                  "geo-external.json")))["cohorts"]

style.apply()
C, INK = style.ROLE, style.INK
W, H = 7.087, 5.20
fig = plt.figure(figsize=(W, H))
src = {}


def panel(letter, x, y):
    fig.text(x, y, letter, fontsize=10, fontweight="bold", va="top", ha="left")


# ---------------------------------------------------------------- a  added value, two references
axa = fig.add_axes([0.19, 0.58, 0.29, 0.36])
CONS = [("ours", "ModRank", C["ours"], "o"), ("survpath", "SurvPath + clinical", C["competitor"], "s"),
        ("pibd_best_val", "PIBD + clinical", C["published"], "D")]
src["a"] = {}
for k, (key, lab, col, mk) in enumerate(CONS):
    v = AV[key]
    y = len(CONS) - 1 - k
    for ref, off, fill in (("weak", 0.13, "white"), ("stage", -0.13, col)):
        pt = v["added_over_weak"] if ref == "weak" else v["added_over_stage"]
        ci = v["added_over_weak_boot"]["ci95"] if ref == "weak" else v["added_over_stage_boot"]["ci95"]
        axa.plot(ci, [y + off, y + off], color=col, lw=1.1, zorder=2)
        axa.scatter([pt], [y + off], s=34, marker=mk, facecolor=fill, edgecolor=col, linewidths=1.2,
                    zorder=3)
        src["a"]["%s_%s" % (key, ref)] = {"point": pt, "ci95": ci}
axa.axvline(0, color=INK, lw=0.7, ls=":")
axa.set_yticks(range(len(CONS)))
axa.set_yticklabels([c[1] for c in CONS][::-1])
axa.set_xlabel("added concordance over the clinical reference")
axa.scatter([], [], s=34, marker="o", facecolor="white", edgecolor=INK, label="over grade (open)")
axa.scatter([], [], s=34, marker="o", facecolor=INK, edgecolor=INK, label="over stage (filled)")
axa.legend(loc="center right", bbox_to_anchor=(1.03, 0.40), fontsize=7, handletextpad=0.3)
axa.set_title("added value depends on the reference", fontsize=8, fontweight="bold")
panel("a", 0.005, 0.985)

# ---------------------------------------------------------------- b  the inflation D
axb = fig.add_axes([0.70, 0.58, 0.28, 0.36])
ROWS = [("ModRank", AV["ours"]["D"], AV["ours"]["D_boot"]["ci95"], C["ours"], "o"),
        ("SurvPath + clinical", AV["survpath"]["D"], AV["survpath"]["D_boot"]["ci95"], C["competitor"], "s"),
        ("PIBD + clinical", AV["pibd_best_val"]["D"], AV["pibd_best_val"]["D_boot"]["ci95"], C["published"], "D"),
        ("GSE32894, transcriptome", GEO["GSE32894"]["D"]["point"], GEO["GSE32894"]["D"]["ci95"], C["omics"], "^"),
        ("GSE31684, transcriptome", GEO["GSE31684"]["D"]["point"], GEO["GSE31684"]["D"]["ci95"], C["omics"], "v")]
src["b"] = {}
for k, (lab, d, ci, col, mk) in enumerate(ROWS):
    y = len(ROWS) - 1 - k
    axb.plot(ci, [y, y], color=col, lw=1.1)
    axb.scatter([d], [y], s=34, marker=mk, color=col, zorder=3)
    src["b"][lab] = {"D": d, "ci95": ci}
axb.axvline(0, color=INK, lw=0.7, ls=":")
axb.axhline(1.5, color="#C9CFD8", lw=0.6)
axb.text(axb.get_xlim()[1], 1.55, "external cohorts", fontsize=6.8, ha="right", va="bottom", color="#5A6273")
axb.set_yticks(range(len(ROWS)))
axb.set_yticklabels([r[0] for r in ROWS][::-1])
axb.set_xlabel("inflation D from a grade-based reference")
axb.set_title("the inflation, measured", fontsize=8, fontweight="bold")
panel("b", 0.515, 0.985)

# ---------------------------------------------------------------- c  grouped calibration, 24 months
axc = fig.add_axes([0.19, 0.09, 0.29, 0.36])
src["c"] = {}
for key, lab, col, mk in (("ours", "ModRank", C["ours"], "o"), ("clinical", "clinical", C["clinical"], "s")):
    g = CAL[key]["grouped_24m"]
    x = [q["predicted"] for q in g]
    y = [q["observed_km"] for q in g]
    axc.plot(x, y, color=col, lw=1.0, marker=mk, linestyle=style.ROLE_LINESTYLE[key], ms=4.5,
             label=lab)
    src["c"][key] = g
lim = [0.0, 0.75]
axc.plot(lim, lim, color=INK, lw=0.7, ls=":")
axc.set_xlim(lim); axc.set_ylim(lim)
axc.set_xlabel("predicted 2-year risk, by quintile")
axc.set_ylabel("observed (Kaplan-Meier)")
axc.legend(loc="upper left", fontsize=7)
axc.text(0.72, 0.05, "slope %.2f (ModRank)\n%.2f (clinical)" % (CAL["ours"]["calibration_slope"],
                                                               CAL["clinical"]["calibration_slope"]),
         fontsize=6.8, ha="right", va="bottom", color="#5A6273")
axc.set_title("calibration inside the training folds", fontsize=8, fontweight="bold")
panel("c", 0.005, 0.49)

# ---------------------------------------------------------------- d  decision curve, 24 months
axd = fig.add_axes([0.70, 0.09, 0.28, 0.36])
th = sorted(float(k) for k in CAL["ours"]["net_benefit_24m"])
series = [("ModRank", [CAL["ours"]["net_benefit_24m"][str(round(x, 2))] for x in th], C["ours"], "o",
           style.ROLE_LINESTYLE["ours"]),
          ("clinical", [CAL["clinical"]["net_benefit_24m"][str(round(x, 2))] for x in th],
           C["clinical"], "s", style.ROLE_LINESTYLE["clinical"]),
          ("treat all", [CAL["treat_all_net_benefit_24m"][str(round(x, 2))] for x in th],
           C["published"], "^", style.ROLE_LINESTYLE["published"])]
src["d"] = {"thresholds": th}
for lab, v, col, mk, ls in series:
    axd.plot(th, v, color=col, marker=mk, linestyle=ls, ms=3.8, lw=1.1, label=lab)
    src["d"][lab] = v
axd.axhline(0, color=INK, lw=0.7, ls=":", label="treat none")
axd.set_ylim(-0.12, 0.33)
axd.set_xlabel("threshold 2-year risk")
axd.set_ylabel("net benefit")
axd.legend(loc="upper right", fontsize=7)
axd.set_title("decision curve, 2 years", fontsize=8, fontweight="bold")
panel("d", 0.515, 0.49)

out = os.path.join(HERE, "fig7_added_value_utility.pdf")
style.save(fig, out)
json.dump(src, open(os.path.join(HERE, "fig7_added_value_utility.source.json"), "w"), indent=1)
print("wrote %s  %.3f x %.3f in" % (out, W, H))
