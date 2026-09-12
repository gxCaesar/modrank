#!/usr/bin/env python3
"""Figure 4 -- what limits a method at 359 patients and 113 events.

This is why the method is parameter-free, and it is the figure a methods referee reads first. Each
panel is a measurement that was taken BEFORE the method was designed, and each one closed off a
direction that would otherwise have looked worth taking.

  a  capacity. A ridge Cox fitted to PURE NOISE and scored in-sample reaches 0.7474 at sixteen
     parameters. The shaded band between the real in-sample fit and that null is the entire
     evidence a capacity increase can produce, and it SHRINKS from d=8 to d=16.
  b  the fusion rule is already at its ceiling: a 21-point sweep peaks exactly at equal weight, and
     an oracle allowed to fit two parameters on the held-out outcomes buys +0.008.
  c  the gated component, drawn with its control ABOVE it. The control wins, so by the
     pre-registered rule the component is void. Drawing the control taller than the result is the
     only honest way to show a void.
  d  what margin this cohort could have detected at all.
  e  two sources of variation the benchmark's tables do not report, both measured by rerunning
     competitors, both larger than the median gap between consecutive published entries.

Panels d and e make the same point by three independent routes -- a power calculation, a measured
seed spread and a measured checkpoint effect -- and being adjacent is the point.
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import style                                    # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "experiments", "20260817-blca-confirm", "results")
S5 = os.path.join(ROOT, "development", "s5-results")
load = lambda d, n: json.load(open(os.path.join(d, n)))

REG = load(RES, "regime-evidence.json")
PAR = load(RES, "pibd-parity.json")
MS = load(RES, "survpath-multiseed.json")
PUB = load(RES, "published-benchmark-table.json")
FSD = json.load(open(os.path.join(HERE, "figure-source-data.json")))
C2 = load(S5, "c2.json")
FUS = load(S5, "fusion.json")

style.apply()
C, INK, LAD = style.ROLE, style.INK, style.LADDER

FIGW, FIGH = style.width(6.785), 6.45
fig = plt.figure(figsize=(FIGW, FIGH))
axa = fig.add_axes([0.100, 0.690, 0.360, 0.280])
axb = fig.add_axes([0.610, 0.690, 0.355, 0.280])
axc = fig.add_axes([0.100, 0.390, 0.250, 0.195])
axd = fig.add_axes([0.470, 0.390, 0.230, 0.195])
axe = fig.add_axes([0.815, 0.390, 0.150, 0.195])
axf = fig.add_axes([0.100, 0.075, 0.865, 0.215])

# =========================================================================== a  capacity
titan = [r for r in REG["capacity_sweep"]["arms"] if r["arm"] == "TITAN"]
titan.sort(key=lambda r: r["d"])
dd = [r["d"] for r in titan]
real = [r["in_sample"] for r in titan]
null = [r["corrected_null"] for r in titan]
oof = [r["out_of_fold"] for r in titan]
axa.fill_between(dd, null, real, color=C["ours"], alpha=0.18, zorder=1, lw=0,
                 label="all the evidence a capacity increase can give")
axa.plot(dd, real, color=INK, lw=1.2, marker="o", ms=3.2, label="fitted in-sample")
axa.plot(dd, null, color=C["competitor"], lw=1.2, ls="--", marker="s", ms=3.2,
         label="fitted to PURE NOISE, scored in-sample")
axa.plot(dd, oof, color=C["ours"], lw=1.2, ls=":", marker="^", ms=3.2, label="out of fold")
axa.set_xscale("log", base=2)
axa.set_xticks(dd)
axa.set_xticklabels([str(x) for x in dd])
axa.set_xlabel("capacity $d$ (principal components fitted)", fontsize=6.4)
axa.set_ylabel("concordance", fontsize=6.4)
axa.set_ylim(0.48, 0.82)
axa.legend(fontsize=5.6, loc="upper left", handlelength=1.8, borderpad=0.15, labelspacing=0.28)
top = titan[-1]
axa.annotate("pure noise reaches\n%.4f at $d$=%d" % (top["corrected_null"], top["d"]),
             (top["d"], top["corrected_null"]), textcoords="offset points", xytext=(-8, -40),
             ha="right", fontsize=6.0, color=C["competitor"], fontweight="bold", linespacing=1.3)
for r in titan[-2:]:
    axa.annotate("%+.4f" % r["excess"], (r["d"], (r["in_sample"] + r["corrected_null"]) / 2),
                 textcoords="offset points", xytext=(6, 6), fontsize=5.9, va="center",
                 color=C["ours"], fontweight="bold")
axa.set_title("in-sample inflation against capacity", fontsize=6.4, pad=3)


# =========================================================================== b  weight sweep
sw = FUS["A_weight_sweep"]["grid"]
ws = sorted(float(k) for k in sw)
vs = [sw["%.2f" % w] for w in ws]
axb.plot(ws, vs, color=C["ours"], lw=1.3, marker="o", ms=2.6)
best_w = float(FUS["A_weight_sweep"]["best_w_on_image"])
axb.scatter([best_w], [sw["%.2f" % best_w]], s=46, color=C["ours"], zorder=4, edgecolor="#FFFFFF",
            linewidth=0.9)
axb.axvline(0.5, color="#9AA1AE", lw=0.7, ls=(0, (2, 2)), zorder=0)
for key, lab, col in (("B_linear_fusion_oracle", "linear fusion oracle, 2 params fitted in-sample",
                       "#5A6273"),
                      ("C_gated_fusion_oracle", "gated fusion oracle", C["competitor"])):
    v = FUS[key]["oracle_in_sample"]
    axb.axhline(v, color=col, lw=0.9, ls="--", zorder=1)
    axb.text(0.02, v + 0.0022, "%s  %.4f" % (lab, v), fontsize=5.6, color=col, va="bottom")
axb.set_xlabel("weight on the image arm", fontsize=6.4)
axb.set_ylabel("concordance", fontsize=6.4)
axb.set_xlim(-0.02, 1.02)
axb.set_title("weight sweep on the image arm", fontsize=6.4, pad=3)
axb.annotate("equal weight,\nnothing fitted", (best_w, sw["%.2f" % best_w]),
             textcoords="offset points", xytext=(0, -30), ha="center", fontsize=6.0,
             color=C["ours"], fontweight="bold", linespacing=1.3)

# =========================================================================== c  the void
c5 = C2["C5_gated_fusion"]
bars = [("flat rank\naverage", c5["flat_rank_average"], "#C9CFD8"),
        ("gated on\nthe real\ntertile", c5["gated_on_clinical_tertile"], C["ours"]),
        ("gated on a\nPERMUTED\ntertile", c5["gated_on_PERMUTED_tertile_CONTROL"],
         C["competitor"])]
xs = np.arange(len(bars))
base = min(b[1] for b in bars)
for x, b in zip(xs, bars):
    axc.plot([x, x], [base - 0.0012, b[1]], color=b[2], lw=1.2, zorder=2, solid_capstyle="round")
    axc.scatter([x], [b[1]], s=34, color=b[2], zorder=3, edgecolor="none")
    axc.text(x, b[1] + 0.0006, "%.4f" % b[1], ha="center", va="bottom", fontsize=5.9,
             fontweight="bold" if b[2] != "#C9CFD8" else "normal")
axc.set_xticks(xs)
axc.set_xticklabels([b[0] for b in bars], fontsize=5.5)
axc.set_xlim(-0.62, 2.62)     # the leftmost value label ran into the y tick labels
axc.set_ylim(0.7425, 0.7555)
axc.set_ylabel("concordance", fontsize=6.4)
assert c5["gated_on_PERMUTED_tertile_CONTROL"] > c5["gated_on_clinical_tertile"], (
    "this panel exists to show the control WINNING; it does not")
axc.set_title("control beats component", fontsize=6.4, pad=3, color=C["competitor"])

# =========================================================================== d  power
P = FSD["F3_power"]
axd.plot(P["margins"], P["power"], color=INK, lw=1.3)
lo, hi = P["published_consecutive_gap_range"]
axd.axvspan(lo, hi, color=C["published"], alpha=0.22, lw=0, zorder=0)
axd.axhline(0.80, color="#9AA1AE", lw=0.7, ls=(0, (2, 2)))
axd.axvline(P["detectable_at_80pc"], color="#9AA1AE", lw=0.7, ls=(0, (2, 2)))
for key, lab, col in (("ours_vs_pibd", "vs PIBD", C["competitor"]),
                      ("ours_vs_survpath", "vs SurvPath", C["competitor"]),
                      ("ours_vs_incumbent_table", "vs table incumbent", C["ours"])):
    m = P["marked"][key]
    axd.scatter([m["margin"]], [m["power"]], s=26, color=col, zorder=4, edgecolor="none")
axd.text(P["marked"]["ours_vs_incumbent_table"]["margin"] + 0.003,
         P["marked"]["ours_vs_incumbent_table"]["power"],
         "our three margins\ncarry %.2f, %.2f, %.2f" % tuple(
             P["marked"][k]["power"] for k in ("ours_vs_pibd", "ours_vs_survpath",
                                               "ours_vs_incumbent_table")),
         fontsize=5.8, va="center", linespacing=1.3)
axd.text((lo + hi) / 2, 0.93, "all eight published\nconsecutive gaps", fontsize=5.7, ha="center",
         va="top", color="#5A6273", linespacing=1.3)
axd.text(P["detectable_at_80pc"] + 0.004, 0.845, "80%% power\nat %.4f" % P["detectable_at_80pc"],
         fontsize=5.8, va="bottom", color="#5A6273", linespacing=1.3)
axd.set_xlim(0, 0.09)
axd.set_ylim(0, 1.02)
axd.set_xlabel("true difference in concordance", fontsize=6.4)
axd.set_ylabel("power", fontsize=6.4)
axd.set_title("power at SE %.4f" % P["se"], fontsize=6.4, pad=3)


# =========================================================================== e  seed spread
seeds = sorted(MS["per_seed"].items(), key=lambda kv: int(kv[0]))
vals = [v["survpath_alone"] for _, v in seeds]
pubv = [e["cindex"] for e in PUB["entries"] if e["method"] == "SurvPath"][0]
axe.scatter([0.5] * len(vals), vals, s=26, color=C["competitor"], edgecolor="none", zorder=3)
axe.plot([0.5, 0.5], [min(vals), max(vals)], color=C["competitor"], lw=0.9, zorder=2)
axe.axhline(pubv, color=INK, lw=0.9, ls="--")
axe.text(0.95, pubv, " published\n %.3f" % pubv, fontsize=5.8, va="center", ha="right",
         linespacing=1.3)
axe.text(0.05, max(vals) + 0.0018, "SD %.4f" % MS["survpath_alone"]["sd_over_seeds"], fontsize=6.0,
         fontweight="bold", color=C["competitor"], va="bottom")
axe.set_xlim(0, 1)
axe.set_xticks([])
axe.set_ylim(0.605, 0.638)
axe.set_ylabel("concordance", fontsize=6.4)
axe.set_title("one competitor,\nfive seeds", fontsize=6.4, pad=3)

# =========================================================================== f  checkpoint
PB = PAR["pibd"]
bv, fe = PB["best_validation_epoch"], PB["final_epoch"]
yb = np.arange(5)
for k in range(5):
    a, b = fe["per_fold_cindex"][k], bv["per_fold_cindex"][k]
    axf.plot([a, b], [k, k], color="#B8BEC9", lw=1.4, zorder=1, solid_capstyle="round")
    axf.scatter([a], [k], s=24, color="#9AA1AE", zorder=3, edgecolor="none")
    axf.scatter([b], [k], s=24, color=C["competitor"], zorder=3, edgecolor="none")
axf.set_yticks(list(yb) + [5.6])
axf.set_yticklabels(["fold %d" % (k + 1) for k in range(5)] + ["mean"], fontsize=6.2)
axf.scatter([fe["fold_mean"]], [5.6], s=44, color="#9AA1AE", zorder=3, edgecolor="none")
axf.scatter([bv["fold_mean"]], [5.6], s=44, color=C["competitor"], zorder=3, edgecolor="none")
axf.plot([fe["fold_mean"], bv["fold_mean"]], [5.6, 5.6], color=INK, lw=1.6, zorder=2)
worth = PAR["published"]["checkpoint_selection_is_worth"]
axf.annotate("+%.4f" % worth, ((fe["fold_mean"] + bv["fold_mean"]) / 2, 5.6),
             textcoords="offset points", xytext=(0, 9), ha="center", fontsize=6.4,
             fontweight="bold", color=INK)
pubp = [e["cindex"] for e in PUB["entries"] if e["method"] == "PIBD"][0]
axf.axvline(pubp, color=INK, lw=0.8, ls="--", zorder=0)
axf.text(pubp + 0.003, 4.35, "published\n%.3f" % pubp, fontsize=5.9, va="center",
         linespacing=1.3)
# the published median gap, drawn at the same scale, is what makes the size of the effect land
med = FSD["F3_power"]["published_consecutive_gap_median"]
x0 = 0.545
axf.plot([x0, x0 + med], [-0.85, -0.85], color=C["published"], lw=3.4, solid_capstyle="butt",
         zorder=3)
axf.text(x0 + med + 0.005, -0.85,
         "for scale: the median gap between consecutive published entries (%.3f)" % med,
         fontsize=5.9, va="center", color="#5A6273")
axf.set_ylim(-1.6, 6.4)
axf.set_xlim(0.535, 0.80)
axf.set_xlabel("concordance", fontsize=6.4)
axf.scatter([], [], s=24, color="#9AA1AE", label="final epoch")
axf.scatter([], [], s=24, color=C["competitor"], label="best-validation epoch, the published convention")
axf.legend(fontsize=5.9, loc="lower right", handletextpad=0.4, borderpad=0.2)
axf.set_title("checkpoint selection, per fold, worth "
              "%.1f$\\times$ the median published gap" % (worth / med), fontsize=6.4, pad=4)

for tag, xx, yy in (("a", 0.010, 0.988), ("b", 0.520, 0.988), ("f", 0.010, 0.322)):
    fig.text(xx, yy, tag, fontsize=7.5, fontweight="bold", color=INK, ha="left", va="top")
# c, d and e sit in the row whose titles run to two lines. Anchoring the letter to the AXES rather
# than the figure puts it above its own title instead of level with it.
for tag, ax, dx in (("c", axc, -0.235), ("d", axd, -0.255), ("e", axe, -0.390)):
    ax.text(dx, 1.22, tag, transform=ax.transAxes, fontsize=7.5, fontweight="bold", color=INK,
            ha="left", va="top")

style.save(fig, style.out(HERE, "fig4_regime.pdf"))
print("wrote fig4_regime.pdf   %.2f x %.2f in" % (FIGW, FIGH))
print("  checkpoint effect %.4f = %.1f x the median published gap %.3f" % (worth, worth / med, med))
