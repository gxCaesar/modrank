#!/usr/bin/env python3
"""Figure 2 -- the missing variable, quantified.

Figure 1 asserted the defect and drew its mechanism. This figure is where it is measured, and the
three panels answer three different objections a referee will raise in order:

  a  "so what" -- the corrected clinical baseline is not merely better than the grade one, it sits
     above two of the four published entries on verified folds and five of all nine. The dot plot
     is ordered and the unverifiable entries are open, so the reader counts both.
  b  "bladder is a special case" -- the same one-column swap in all five studies the benchmark
     covers, each on its own cases.
  c  "then why does it help so much more in some" -- because the swap buys most where grade is
     least informative, which is the mechanism rather than a restatement of it. THREE POINTS. The
     panel says so inside itself and not only in the caption, because a correlation coefficient
     printed beside three dots is read as an estimate unless it is labelled as not being one.

Panel c is the one that could mislead, so it carries its own disclaimer at full size.
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import style                                        # noqa: E402
from schematic import label, panel_letter           # noqa: E402

RES = os.path.abspath(os.path.join(HERE, "..", "..", "experiments",
                                   "20260817-blca-confirm", "results"))
load = lambda n: json.load(open(os.path.join(RES, n)))

PUB = load("published-benchmark-table.json")
A1 = load("amendment-A1-clinical-provenance.json")
WG = load("why-grade-fails.json")

style.apply()
C, INK, LAD = style.ROLE, style.INK, style.LADDER

CENSUS = json.load(open(os.path.abspath(os.path.join(HERE, "..", "..", "development",
                                                     "clinical-baseline-census.json"))))
# Panel e used to read development/s5-results/decomp.json, which had no committed producer. Two of
# its four rows turned out to be unprovenanced rather than one: B and D reproduce exactly from
# committed code, A does not (0.6856/0.7291 rebuilt against 0.6863/0.7285 in that file, with the
# rebuild agreeing with the frozen run's own clinical arm), and C cannot be attempted because the
# rule that chose one diagnosis per case was never recorded. The panel now draws the three
# constructions a script regenerates and drops the one it cannot.
DECOMP = json.load(open(os.path.abspath(os.path.join(
    HERE, "..", "..", "experiments", "20260912-decomp-emitter", "results",
    "clinical-block-decomposition.json"))))["constructions"]
STRATA = load("per-stage-subgroup.json")["strata"]
CPM = load("cold-panel-round1-measurements.json")
# The five-study swap carried only point estimates until 2026-09-12. Panel c now draws the interval
# beside each one, because the main text was corrected the same day to say that the swap is positive
# in all five studies and separable from zero in two, and a figure showing five bare gaps would be
# making the stronger claim the text had just given up.
IV = json.load(open(os.path.abspath(os.path.join(
    HERE, "..", "..", "experiments", "20260912-five-cohort-intervals", "results",
    "five-cohort-intervals.json"))))["cohorts"]

fig = plt.figure(figsize=(style.width(6.785), 6.45))
# THREE ROWS OF TWO rather than four rows of one-and-two. Panels a and d are both eleven-row
# lists and both were full width, which is what made this figure tall; side by side they cost one
# row instead of two. The venue column that used to sit beside d moves into the caption, since at
# half width it would have squeezed the marks the panel exists to show.
axa = fig.add_axes([0.178, 0.618, 0.300, 0.330])
axd = fig.add_axes([0.672, 0.606, 0.190, 0.322])
axb = fig.add_axes([0.135, 0.352, 0.325, 0.152])
axc = fig.add_axes([0.660, 0.352, 0.295, 0.152])
axe = fig.add_axes([0.150, 0.088, 0.300, 0.138])
axf = fig.add_axes([0.660, 0.088, 0.290, 0.138])

# =========================================================================== a  the ranking
# A DOT PLOT, not bars. The axis cannot start at zero without flattening every difference the panel
# exists to show, and a bar drawn from 0.50 overstates each gap by its truncated length -- the defect
# this lab's figure rules name explicitly. Entries whose folds could not be verified are drawn OPEN,
# so the count the caption states ("two of the four on verified folds, five of all nine") can be read
# off the panel rather than taken on trust.
def _verified(e):
    return e["folds"].startswith("verified") or e["folds"] in ("definitional",
                                                                "as reported in SurvPath")


rows = [(e["method"], e["cindex"], "published", _verified(e)) for e in PUB["entries"]]
rows.append(("clinical, age+sex+stage", A1["arms_seed0"]["clinical_age_sex_stage_from_DIMAF_file"],
             "clinical", True))
rows.append(("clinical, age+sex+grade", A1["arms_seed0"]["clinical_age_sex_GRADE_from_DIMAF_file"],
             "competitor", True))
rows.append(("ours", A1["primary"]["value"], "ours", True))
rows.sort(key=lambda r: r[1])

y = np.arange(len(rows))
for yi, r in zip(y, rows):
    axa.plot([0.50, r[1]], [yi, yi], color="#E3E6EC", lw=0.8, zorder=1)
    col = C[r[2]]
    axa.scatter([r[1]], [yi], s=30, zorder=3, linewidths=1.1,
                facecolor=col if r[3] else "white", edgecolor=col)
axa.set_yticks(y)
axa.set_yticklabels([r[0] for r in rows], fontsize=6.4)
for t_, r in zip(axa.get_yticklabels(), rows):
    if r[2] != "published":
        t_.set_fontweight("bold")
    elif not r[3]:
        t_.set_color("#6B7280")
axa.set_xlim(0.50, 0.76)
axa.set_xlabel("concordance on TCGA-BLCA disease-specific survival, $n=%d$" % PUB["cohort"]["n"],
               fontsize=6.4)
for yi, r in zip(y, rows):
    axa.text(r[1] + 0.006, yi, "%.4f" % r[1], va="center", fontsize=6.1,
             fontweight="bold" if r[2] != "published" else "normal",
             color=INK if r[2] != "published" else "#6B7280")

# the rule that makes the claim countable, on both denominators
stage_c = A1["arms_seed0"]["clinical_age_sex_stage_from_DIMAF_file"]
below = sum(1 for e in PUB["entries"] if e["cindex"] < stage_c)
ver = [e for e in PUB["entries"] if _verified(e)]
below_v = sum(1 for e in ver if e["cindex"] < stage_c)
axa.axvline(stage_c, color=C["clinical"], lw=0.9, ls=(0, (3, 2)), zorder=0)
# The count itself ("two of four on verified folds, five of all nine") is stated in the caption and
# can be read off the markers; an in-panel sentence collided with the value labels.
axa.scatter([], [], s=30, facecolor="#9AA1AD", edgecolor="#9AA1AD", label="released folds")
axa.scatter([], [], s=30, facecolor="white", edgecolor="#9AA1AD", linewidths=1.1,
            label="folds unverifiable")
_lg = axa.legend(loc="upper left", fontsize=5.8, frameon=True, handletextpad=0.3,
                 borderaxespad=0.3, borderpad=0.3)
_lg.get_frame().set_facecolor("white")
_lg.get_frame().set_edgecolor("none")

# =========================================================================== c  five studies (lettered b before 2026-09-11)
order = ["blca", "brca", "coadread", "hnsc", "stad"]
NAME = {"blca": "BLCA (bladder)", "brca": "BRCA", "coadread": "COADREAD", "hnsc": "HNSC",
        "stad": "STAD"}
yb = np.arange(len(order))[::-1]
for yi, k in zip(yb, order):
    co = WG["cohorts"][k]
    g, s = co["c_age_sex_grade"], co["c_age_sex_stage"]
    axb.plot([g, s], [yi, yi], color="#B8BEC9", lw=1.4, zorder=1, solid_capstyle="round")
    axb.scatter([g], [yi], s=26, color=C["competitor"], marker=style.ROLE_MARKER["competitor"],
                zorder=3, edgecolor="none")
    axb.scatter([s], [yi], s=26, color=C["clinical"], marker=style.ROLE_MARKER["clinical"],
                zorder=3, edgecolor="none")
    d = IV[k]["differences"]["stage_minus_grade"]
    excl = d["ci95"][0] > 0
    axb.text(s + 0.010, yi, "+%.3f%s" % (s - g, "" if excl else " (n.s.)"), va="center",
             fontsize=6.0, color=INK if excl else "#6B7280")
axb.set_yticks(yb)
axb.set_yticklabels([NAME[k] for k in order], fontsize=6.5)
axb.set_xlim(0.43, 0.79)
axb.set_xlabel("concordance", fontsize=6.4)
axb.scatter([], [], s=26, color=C["competitor"], marker=style.ROLE_MARKER["competitor"],
            label="age + sex + grade")
axb.scatter([], [], s=26, color=C["clinical"], marker=style.ROLE_MARKER["clinical"],
            label="age + sex + stage")
# an in-axes legend at "upper right" sat on top of the BLCA row's markers and its "+0.056" label.
# One row above the axes, carrying the panel's own title as the legend title, collides with nothing.
_lgb = axb.legend(fontsize=6.1, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2,
                  handletextpad=0.35, borderpad=0.2, columnspacing=1.1,
                  title="grade $\\rightarrow$ stage, five studies", title_fontsize=6.4)
_lgb.get_title().set_fontweight("normal")


# =========================================================================== d  the mechanism (lettered c before 2026-09-11)
pts = [(k, WG["cohorts"][k]["grade"]["normalised_entropy"],
        WG["cohorts"][k]["stage_minus_grade"])
       for k in order if isinstance(WG["cohorts"][k].get("grade"), dict)]
missing = [k for k in order if not isinstance(WG["cohorts"][k].get("grade"), dict)]
assert len(pts) + len(missing) == len(order), "a cohort vanished between the two groups"

axc.scatter([p[1] for p in pts], [p[2] for p in pts], s=34, color=C["clinical"],
            edgecolor="none", zorder=3)
for k, e, gn in pts:
    axc.annotate(k.upper(), (e, gn), textcoords="offset points", xytext=(0, 7), ha="center",
                 fontsize=6.0, color=INK)
axc.set_xlabel("normalised entropy of grade", fontsize=6.4)
axc.set_ylabel("gain from the swap", fontsize=6.4)
axc.set_xlim(0.22, 0.76)
axc.text(0.97, 0.99, "Pearson $-%.3f$" % abs(WG["gain_vs_grade_entropy_pearson"]),
         transform=axc.transAxes, fontsize=6.3, fontweight="bold", va="top", ha="right")
# the disclaimer is INSIDE the panel. A coefficient printed beside three dots reads as an estimate.
axc.text(0.97, 0.84, "$n=3$ cohorts.\nA consistency check,\nnot an estimate.",
         transform=axc.transAxes, fontsize=6.0, va="top", ha="right", color="#5A6273",
         linespacing=1.35)
axc.text(0.03, 0.05, "%s:\nno usable grade\nin the shipped file"
         % ", ".join(m.upper() for m in missing),
         transform=axc.transAxes, fontsize=5.9, ha="left", va="bottom", color="#6B7280",
         linespacing=1.3)


# =========================================================================== b  the census itself (lettered d before 2026-09-11)
# The paper's central factual claim is about ELEVEN papers, and until now it appeared only as prose.
# One row per paper, one mark per question, and the three that could not be read say so rather than
# being filled in from what the others do.
PAPERS = CENSUS["papers"]
assert len(PAPERS) == 11, "the census covers eleven papers; found %d" % len(PAPERS)
yd = np.arange(len(PAPERS))[::-1]
COLS = [("reports a\nbaseline", "reports_clinical_baseline"),
        ("from\ngrade", "_grade"), ("uses\nstage", "uses_stage")]
GRADE_WORDS = ("grade",)
n_grade = n_stage = n_none = n_unver = 0
for yi, p in zip(yd, PAPERS):
    unver = p["status"] != "verified"
    cov = p.get("covariates") or []
    vals = [None if unver else bool(p["reports_clinical_baseline"]),
            None if unver else any(any(w in str(c).lower() for w in GRADE_WORDS) for c in cov),
            None if unver else bool(p["uses_stage"])]
    if unver:
        n_unver += 1
    else:
        n_none += (not vals[0]); n_grade += bool(vals[1]); n_stage += bool(vals[2])
    for j, v in enumerate(vals):
        x = j
        if v is None:
            axd.plot([x - 0.13, x + 0.13], [yi, yi], color="#B8BEC9", lw=1.2,
                     solid_capstyle="round", zorder=3)
        else:
            face = (C["competitor"] if j == 1 else C["clinical"]) if v else "#FFFFFF"
            axd.scatter([x], [yi], s=42, facecolor=face, zorder=3,
                        edgecolor=C["competitor"] if j == 1 else C["clinical"], linewidth=0.9)
    axd.text(-0.62, yi, p["name"], ha="right", va="center", fontsize=6.1,
             color="#8C8C8C" if unver else INK)
axd.set_xlim(-0.62, 2.62)
axd.set_ylim(-0.8, len(PAPERS) - 0.2)
axd.set_xticks(range(3))
axd.set_xticklabels([c[0] for c in COLS], fontsize=6.1)
axd.xaxis.set_ticks_position("top")
axd.tick_params(axis="x", length=0, pad=2)
axd.set_yticks([])
for sp in axd.spines.values():
    sp.set_visible(False)
axd.text(1.0, -0.95, "filled = yes,  open = no,  bar = the paper's own text could not be read",
         ha="center", va="top", fontsize=5.7, color="#6B7280")
axd.set_title("the eleven-paper census", fontsize=6.4, pad=19, fontweight="bold")
assert n_stage == 0, "the claim of this panel is that NO paper uses stage; %d do" % n_stage

# =========================================================================== e  four constructions
KEYS = ["A", "B", "D"]
assert set(KEYS) == set(DECOMP), "the emitter scored %s; the panel draws %s" % (set(DECOMP), KEYS)
SHORT = {"A": "before amendment", "B": "amended", "D": "grade for stage"}
ye = np.arange(len(KEYS))[::-1]
for yi, k in zip(ye, KEYS):
    alone, combined = DECOMP[k]["clinical_alone"], DECOMP[k]["with_both_modalities"]
    role_alone = "competitor" if k == "D" else "clinical"
    col = C[role_alone]
    axe.plot([alone, combined], [yi, yi], color="#B8BEC9", lw=1.4, zorder=1, solid_capstyle="round")
    axe.scatter([alone], [yi], s=26, color=col, marker=style.ROLE_MARKER[role_alone], zorder=3,
                edgecolor="none")
    axe.scatter([combined], [yi], s=26, color=C["ours"], marker=style.ROLE_MARKER["ours"], zorder=3,
                edgecolor="none")
axe.set_yticks(ye)
axe.set_yticklabels([SHORT[k] for k in KEYS], fontsize=6.2)
axe.set_xlim(0.54, 0.76)
axe.set_xlabel("concordance", fontsize=6.3)
axe.scatter([], [], s=26, color=C["clinical"], marker=style.ROLE_MARKER["clinical"],
            label="clinical block alone")
axe.scatter([], [], s=26, color=C["ours"], marker=style.ROLE_MARKER["ours"],
            label="with slide and transcriptome")
axe.legend(fontsize=5.7, loc="upper center", bbox_to_anchor=(0.5, -0.42), ncol=2,
           handletextpad=0.35, columnspacing=1.1, borderpad=0.15)
axe.set_title("three clinical constructions",
              fontsize=6.3, pad=4)

# =========================================================================== f  inside one stage
SHOWN = [(k, v) for k, v in STRATA.items() if "skipped" not in v]
assert len(SHOWN) == 3, "three strata clear the reporting floor; found %d" % len(SHOWN)
ARMS_F = [("clinical", "clinical"), ("wsi_titan", "slide"), ("omics", "omics"), ("OURS", "ours")]
MODF = {"clinical": C["clinical"], "slide": C["slide"], "omics": C["omics"], "ours": C["ours"]}
w = 0.20
xf = np.arange(len(SHOWN))
for i, (key, fam) in enumerate(ARMS_F):
    xx = xf + (i - 1.5) * w
    vv = [v[key] for _, v in SHOWN]
    for x_, v_ in zip(xx, vv):
        axf.plot([x_, x_], [0.405, v_], color=MODF[fam], lw=1.1, linestyle=style.ROLE_LINESTYLE[fam],
                 solid_capstyle="round", zorder=1)
    axf.plot(xx, vv, marker=style.ROLE_MARKER[fam], linestyle="none", color=MODF[fam], ms=3.4,
             markeredgecolor="none", zorder=3)
axf.axhline(0.5, color=INK, lw=0.7, ls=":", zorder=3)
# the four stems were identified only by colour, marker and linestyle, with no legend anywhere in
# the figure. Proxy handles carry all three, in the same corner every stage's tallest stem clears
# (Stage IV's cluster, nearest this corner, tops out well below 0.65).
NAME_F = {"clinical": "clinical", "slide": "slide", "omics": "transcriptome", "ours": "ModRank"}
for key, fam in ARMS_F:
    axf.plot([], [], color=MODF[fam], marker=style.ROLE_MARKER[fam],
             linestyle=style.ROLE_LINESTYLE[fam], lw=1.1, ms=3.4, markeredgecolor="none",
             label=NAME_F[fam])
# one row above the axes: inside, it sat on the Stage III stems (2026-09-27)
axf.legend(fontsize=5.3, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, handlelength=1.6,
           handletextpad=0.35, columnspacing=0.8, borderpad=0.15, frameon=False)
axf.set_xticks(xf)
axf.set_xticklabels(["%s\n$n=%d$, %d ev" % (k.replace("Stage ", "Stage\u2009"), v["n"],
                                            v["events"]) for k, v in SHOWN], fontsize=5.9)
axf.set_ylim(0.40, 0.76)
axf.set_ylabel("concordance", fontsize=6.3)
axf.set_title("within the incumbent's stage groups", fontsize=6.3, pad=13)
worst = min(v["clinical"] for _, v in SHOWN)
axf.text(0.5, -0.40, "clinical falls to %.3f, below chance, in Stage III" % worst,
         transform=axf.transAxes, ha="center", va="top", fontsize=5.6, color="#6B7280")

# letters in reading order (a b / c d / e f); the census was "d" and the five-study swap "b" until
# the 180 mm rebuild on 2026-09-11 showed the top row reading a, d
for tag, xx, yy in (("a", 0.012, 0.986), ("b", 0.545, 0.986), ("c", 0.012, 0.548),
                    ("d", 0.545, 0.548), ("e", 0.012, 0.268), ("f", 0.545, 0.268)):
    fig.text(xx, yy, tag, fontsize=7.5, fontweight="bold", color=INK, ha="left", va="top")

style.save(fig, style.out(HERE, "fig2_missing_variable.pdf"))
print("wrote fig2_missing_variable.pdf   %.3f x %.3f in" % tuple(fig.get_size_inches()))
print("  panel a: %d rows, corrected baseline clears %d of %d verified, %d of %d published"
      % (len(rows), below_v, len(ver), below, len(PUB["entries"])))
print("  panel d: %d cohorts plotted, %d excluded for no usable grade" % (len(pts), len(missing)))
