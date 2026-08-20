#!/usr/bin/env python3
"""Figure 5 -- what the two molecular-facing arms are reading.

DESCRIPTIVE, and the figure says so. These are associations in 359 patients, at a sample size the
power analysis in Figure 4 shows cannot resolve the margins the method is judged on. Nothing here
fed back into the method and nothing here is a discovery claim.

The figure is built so its central claim could have failed. If the slide arm were a molecular
proxy -- a roundabout way of reading the transcriptome -- panel b would show it, and the paper's
reason for carrying two molecular-facing arms instead of one would collapse. It does not.

  a  are the eight axes prognostic at all. This has to be established BEFORE asking what an arm
     tracks, or a correlation with an unprognostic signature reads as a finding.
  b  the contrast, arm by arm. The omics arm is a molecular-subtype reader; the slide arm is not.
  c  does the method add WITHIN a molecular class, where the class itself carries no information.
  d  the strongest pathway associations of the slide arm, of 275 tested.
  e  the twelve case studies the frozen rule admitted, drawn as a swimmer plot because it shows
     censoring natively -- and every case the rule admitted is here, including the two that
     embarrass the model.
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
BIO = json.load(open(os.path.join(ROOT, "experiments", "20260817-blca-confirm", "results",
                                  "biology.json")))
style.apply()
C, INK, LAD = style.ROLE, style.INK, style.LADDER

FIGW, FIGH = 6.785, 6.45
fig = plt.figure(figsize=(FIGW, FIGH))
axa = fig.add_axes([0.175, 0.735, 0.310, 0.215])
axb = fig.add_axes([0.640, 0.735, 0.325, 0.215])
axc = fig.add_axes([0.115, 0.455, 0.300, 0.185])
axd = fig.add_axes([0.660, 0.455, 0.230, 0.185])
axe = fig.add_axes([0.175, 0.108, 0.640, 0.262])

NICE = {"basal": "basal", "luminal": "luminal", "EMT": "EMT / claudin-low",
        "immune": "immune T-cell", "proliferation": "proliferation",
        "p53_cell_cycle": "p53 / cell cycle", "FGFR3_axis": "FGFR3 axis",
        "stroma": "stroma / fibroblast"}

# =========================================================================== a  prognostic?
B2 = BIO["B2_survival_association"]["axes"]
rows = []
for k, v in B2.items():
    # the score statistic is unsigned in the file; its DIRECTION is recoverable from whether the
    # axis scores above or below chance as a raw risk score, and the sign is the whole point here
    sign = 1.0 if v["cindex_as_a_raw_score"] > 0.5 else -1.0
    rows.append((k, sign * v["cox_score_statistic"], v["bh_q"]))
rows.sort(key=lambda r: r[1])
y = np.arange(len(rows))
for yi, (k, z, q) in zip(y, rows):
    sig = q < 0.05
    axa.barh([yi], [z], color=C["ours"] if z > 0 else C["omics"], height=0.62,
             alpha=1.0 if sig else 0.30, edgecolor="none")
    axa.text(z + (0.12 if z > 0 else -0.12), yi, "%.3f" % q, va="center",
             ha="left" if z > 0 else "right", fontsize=5.6,
             fontweight="bold" if sig else "normal", color=INK if sig else "#6B7280")
axa.axvline(0, color=INK, lw=0.8)
for v in (-1.96, 1.96):
    axa.axvline(v, color="#9AA1AE", lw=0.7, ls=(0, (2, 2)))
axa.set_yticks(y)
axa.set_yticklabels([NICE[k] for k, _, _ in rows], fontsize=6.3)
axa.set_xlim(-4.6, 4.6)
axa.set_xlabel("Cox score statistic (signed)", fontsize=6.4)
n_sig = sum(1 for _, _, q in rows if q < 0.05)
axa.set_title("molecular axes, signed Cox score (%d of %d after BH)"
              % (n_sig, len(rows)), fontsize=6.3, pad=3)
axa.text(0.99, 0.02, "solid = BH $q<0.05$", transform=axa.transAxes, ha="right", fontsize=5.7,
         color="#5A6273")


# =========================================================================== b  the contrast
SBS = BIO["B4b_contrast"]["axes_side_by_side"]
keys = [k for k, _, _ in rows]                       # same order as panel a, so the eye carries over
yb = np.arange(len(keys))
for yi, k in zip(yb, keys):
    im, om = SBS[k]["image"], SBS[k]["omics"]
    axb.plot([im, om], [yi, yi], color="#B8BEC9", lw=1.2, zorder=1, solid_capstyle="round")
    axb.scatter([im], [yi], s=24, color=C["slide"], zorder=3, edgecolor="none")
    axb.scatter([om], [yi], s=24, color=C["omics"], zorder=3, edgecolor="none")
axb.axvline(0, color=INK, lw=0.8)
axb.set_yticks(yb)
axb.set_yticklabels([])
axb.text(1.0, -0.30, "rows as in a", transform=axb.transAxes, fontsize=5.4,
         color="#6B7280", ha="right", va="top")
axb.set_xlim(-0.46, 0.46)
axb.set_xlabel("Spearman $\\rho$ with the arm's out-of-fold score", fontsize=6.4)
axb.scatter([], [], s=24, color=C["slide"], label="slide arm")
axb.scatter([], [], s=24, color=C["omics"], label="omics arm")
axb.legend(fontsize=6.0, loc="lower right", handletextpad=0.35, borderpad=0.2)
rho = BIO["B4b_contrast"]["spearman_image_vs_omics_arm"]
axb.set_title("subtype association per arm ($\\rho$ "
              "%.2f)" % rho, fontsize=6.3, pad=3)

# =========================================================================== c  within subtype
B5 = BIO["B5_within_molecular_subtype"]
groups = [("luminal-leaning", B5["luminal_leaning"]), ("basal-leaning", B5["basal_leaning"])]
ARMS = [("clinical", "clinical"), ("image", "slide"), ("omics", "omics"), ("OURS", "ours")]
xc = np.arange(len(groups))
w = 0.20
for j, (key, role) in enumerate(ARMS):
    vals = [g[1][key] for g in groups]
    # a bar encodes its value from zero, so a y-axis starting at 0.5 overstates these gaps.
    # The fourth panel in this package to carry that defect, and the last.
    xx = xc + (j - 1.5) * w
    for x, v in zip(xx, vals):
        axc.plot([x, x], [0.508, v], color=C[role], lw=1.1, solid_capstyle="round", zorder=1)
        axc.text(x, v + 0.006, "%.2f" % v, ha="center", va="bottom", fontsize=5.0, rotation=90,
                 color=C[role], fontweight="bold")
    axc.plot(xx, vals, "o", color=C[role], ms=3.2, markeredgecolor="none", zorder=3,
             label={"OURS": "ours"}.get(key, key))
axc.axhline(0.5, color="#9AA1AE", lw=0.7, ls=(0, (2, 2)))
axc.set_xticks(xc)
axc.set_xticklabels(["%s\n$n$=%d, %d events" % (g[0], g[1]["n"], g[1]["events"]) for g in groups],
                    fontsize=6.0)
axc.set_ylim(0.502, 0.79)
axc.set_ylabel("concordance within the class", fontsize=6.4)
axc.legend(fontsize=5.8, ncol=2, loc="upper left", handlelength=1.0, handletextpad=0.4,
           columnspacing=0.8, borderpad=0.2)
axc.set_title("concordance within a molecular class",
              fontsize=6.3, pad=3)

# =========================================================================== d  pathways
TOP = BIO["B3_what_the_image_arm_tracks"]["top_12_pathways_by_absolute_correlation"][:10]
TOP = sorted(TOP, key=lambda r: r["spearman"])
yd = np.arange(len(TOP))
for yi, r in zip(yd, TOP):
    axd.plot([0, r["spearman"]], [yi, yi], color="#B8BEC9", lw=1.0, zorder=1)
    axd.scatter([r["spearman"]], [yi], s=20, zorder=3, edgecolor="none",
                color=C["slide"] if r["spearman"] > 0 else C["omics"])
axd.axvline(0, color=INK, lw=0.8)
axd.set_yticks(yd)
def short(name, n=34):
    t = name.replace("_", " ").replace("HALLMARK ", "")
    return t if len(t) <= n else t[:n - 1].rstrip() + "\u2026"


axd.set_yticklabels([short(r["pathway"]) for r in TOP], fontsize=5.5)
axd.set_xlabel("Spearman $\\rho$ with the slide arm", fontsize=6.4)
n_q = BIO["B3_what_the_image_arm_tracks"]["pathways_with_BH_q_below_0.05"]
n_tested = BIO["B2_survival_association"]["pathways_tested"]
axd.set_title("the %d of %d pathways at BH $q<0.05$, ten strongest" % (n_q, n_tested),
              fontsize=6.3, pad=3)


# =========================================================================== e  case studies
CS = BIO["B6_case_studies"]
cases = CS["cases"]
GROUP = [("highest full-method score in the band", "revised UP by the model", C["competitor"]),
         ("lowest full-method score in the band", "revised DOWN by the model", C["ours"]),
         ("image arm far ABOVE omics arm", "slide arm far above omics", C["slide"]),
         ("omics arm far ABOVE image arm", "omics arm far above slide", C["omics"])]
ordered, bands = [], []
for sel, lab, col in GROUP:
    grp = [c for c in cases if c["selected_as"] == sel]
    bands.append((len(ordered), len(grp), lab, col))
    ordered += [(c, col) for c in grp]
assert len(ordered) == len(cases), "the grouping dropped %d of %d cases the rule admitted" % (
    len(cases) - len(ordered), len(cases))

ye = np.arange(len(ordered))[::-1]
for yi, (c, col) in zip(ye, ordered):
    t = c["observed_months"]
    axe.plot([0, t], [yi, yi], color=col, lw=2.6, solid_capstyle="butt", alpha=0.85, zorder=2)
    axe.scatter([t], [yi], s=30, zorder=4, color=col if c["event"] else "#FFFFFF",
                edgecolor=col, linewidth=1.1, marker="o" if c["event"] else "o")
    axe.text(-1.6, yi, "%s   %s" % (c["case_id"], c["stage"]), ha="right", va="center",
             fontsize=5.5, color=INK)
    axe.text(t + 1.6, yi, "%.1f mo" % t, ha="left", va="center", fontsize=5.4, color="#5A6273")
seen = []
for start, n, lab, col in bands:
    if lab not in [t for t, _ in seen]:
        seen.append((lab, col))
handles = [plt.Line2D([], [], color=c, lw=2.6, label=t) for t, c in seen]
axe.set_yticks([])
axe.set_xlim(0, 78)
axe.set_ylim(-0.8, len(ordered) - 0.2)
axe.set_xlabel("observed time (months)", fontsize=6.4)
axe.scatter([], [], s=30, color=INK, label="died of disease")
axe.scatter([], [], s=30, facecolor="#FFFFFF", edgecolor=INK, linewidth=1.1, label="censored")
marks = axe.legend(fontsize=6.0, loc="lower right", handletextpad=0.4, borderpad=0.25)
axe.add_artist(marks)
axe.legend(handles=handles, fontsize=5.9, loc="upper right", handletextpad=0.5, borderpad=0.25,
           labelspacing=0.35, title="selected as", title_fontsize=5.9)
n_ev = sum(1 for c, _ in ordered if c["event"])
axe.set_title("the twelve cases the frozen rule admitted (%d events)" % n_ev, fontsize=6.3, pad=4)


for tag, xx, yy in (("a", 0.010, 0.992), ("b", 0.505, 0.992), ("c", 0.010, 0.678),
                    ("d", 0.545, 0.678), ("e", 0.010, 0.400)):
    fig.text(xx, yy, tag, fontsize=7.5, fontweight="bold", color=INK, ha="left", va="top")

style.save(fig, os.path.join(HERE, "fig5_biology.pdf"))
print("wrote fig5_biology.pdf   %.2f x %.2f in" % (FIGW, FIGH))
print("  a: %d of %d axes prognostic after BH" % (n_sig, len(rows)))
print("  e: %d cases, %d events, %d groups" % (len(ordered), n_ev, len(bands)))
