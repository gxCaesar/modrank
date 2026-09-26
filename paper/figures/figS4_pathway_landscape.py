#!/usr/bin/env python3
"""Supplementary Figure S4 -- all 275 pathways, not the fifteen the main text had room for.

Three things a reader cannot get from a top-fifteen table.

  a  the survival landscape. Every pathway's signed Cox score against its p, with the
     Benjamini-Hochberg threshold drawn. NOT ONE of the 275 clears q<0.10 as a univariate
     predictor. That is the honest headline of this panel and it is why the omics arm's value
     comes from the ensemble rather than from any single programme.
  b  what the two molecular-facing arms track, pathway by pathway, against each other. The two
     ARMS correlate at rho 0.28 as scores (computed below from the dump's seed-averaged scores),
     yet their pathway-association profiles correlate at r 0.82: they disagree about patients far
     more than they disagree about biology. Until 2026-09-11 the artwork TYPED rho=0.48, which is
     the slide arm against the incumbent's joint score, not against the omics arm.
  c  the distribution of |rho| per arm, which is the same fact as b without the reader having to
     judge a cloud by eye.

Built from the frozen reporting dump, which asserts its own primary against the frozen 0.7212
before writing anything.
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import style                                                             # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
D = json.load(open(os.path.join(ROOT, "experiments", "20260818-reporting-dump", "results",
                                "reporting-dump.json")))
assert D["frozen_primary_matched"], "the dump did not reproduce the frozen primary"

style.apply()
C, INK = style.ROLE, style.INK
GREY, FAINT = "#8C8C8C", "#C9CFD8"

P = D["pathways"]
assert len(P) == 275, "the benchmark's grouping is 275 pathways; the dump has %d" % len(P)
# the statistic is unsigned by construction; direction comes from whether the pathway scores
# above or below chance as a raw risk score, which is the convention the main biology figure
# already uses. Inventing a second convention for a supplement would be worse than either.
sgn = np.array([1.0 if p["cindex_raw"] > 0.5 else -1.0 for p in P])
score = sgn * np.array([p["cox_score_abs"] for p in P])
pval = np.array([p["p"] for p in P])
qval = np.array([p["q_BH"] for p in P])
rs = np.array([p["rho_slide"] for p in P])
ro = np.array([p["rho_omics"] for p in P])
names = [p["pathway"].replace("_", " ") for p in P]
n_bh = int((qval < 0.10).sum())

FIGW, FIGH = 6.30, 4.60   # 160 mm (PI correction 2026-09-27: this is the SI text block width,
# 453.6 pt, and every figure must insert at scale 1.0000 -- 180 mm is a valid Nature artwork width
# in the abstract but wrong for this document)
fig = plt.figure(figsize=(FIGW, FIGH))
axa = fig.add_axes([0.098, 0.575, 0.520, 0.345])
axb = fig.add_axes([0.740, 0.575, 0.235, 0.345])
axc = fig.add_axes([0.098, 0.105, 0.520, 0.305])
axd = fig.add_axes([0.740, 0.105, 0.235, 0.305])


def panel(tag, x, y):
    fig.text(x, y, tag, fontsize=7.5, fontweight="bold", color=INK, ha="left", va="top")


# =========================================================================== a  the landscape
y = -np.log10(np.clip(pval, 1e-12, None))
axa.scatter(score, y, s=7, color=C["omics"], alpha=0.55, edgecolor="none", zorder=3)
# the smallest p that would have cleared BH at 0.10, drawn where it actually falls
k = np.argsort(pval)
bh_line = 0.10 * 1.0 / len(P)
axa.axhline(-np.log10(bh_line), color=C["competitor"], lw=0.9, ls=(0, (3, 2)), zorder=2)
axa.text(0.02, -np.log10(bh_line), "BH $q<0.10$ needs $p<%.1g$" % bh_line,
         transform=axa.get_yaxis_transform(), fontsize=5.2, color=C["competitor"],
         va="bottom", ha="left")
idx = k[0]
# 34 was calibrated for the 160 mm build and, for this pathway's name, truncated mid-parenthesis
# ("Epidermal Growth Factor Receptor ("). The 180 mm build has more room; rstrip so a future run's
# longest name cannot land on a dangling "(" or trailing space either.
label_name = names[idx][:40].rstrip(" (")
axa.annotate(label_name, (score[idx], y[idx]), textcoords="offset points",
             xytext=(-7, -1), fontsize=5.2, color=INK, ha="right", va="center")
axa.text(0.98, 0.06, "strongest of 275, and it still\ndoes not clear the correction",
         transform=axa.transAxes, ha="right", va="bottom", fontsize=5.2, color=GREY,
         linespacing=1.35)
axa.axvline(0, color=INK, lw=0.7)
axa.set_xlabel("Cox score statistic, signed by direction", fontsize=6.4)
axa.set_ylabel("$-\\log_{10} p$", fontsize=6.4)
axa.set_title("275 pathways against survival (%d clear BH $q<0.10$)" % n_bh, fontsize=6.4, pad=3)
panel("a", 0.010, 0.980)

# =========================================================================== b  do they agree
axb.scatter(rs, ro, s=6, color=GREY, alpha=0.6, edgecolor="none", zorder=3)
lim = float(max(np.abs(rs).max(), np.abs(ro).max())) * 1.12
axb.plot([-lim, lim], [-lim, lim], color=C["competitor"], lw=0.8, ls=(0, (3, 2)), zorder=2)
axb.set_xlim(-lim, lim); axb.set_ylim(-lim, lim)
axb.set_xlabel("$\\rho$ with the slide arm", fontsize=6.4)
axb.set_ylabel("$\\rho$ with the omics arm", fontsize=6.4)
r_between = float(np.corrcoef(rs, ro)[0, 1])
_sl = np.array([c_["seed_mean"]["wsi_titan"] for c_ in D["cases"]])
_om = np.array([c_["seed_mean"]["omics_combine"] for c_ in D["cases"]])
rho_arms = float(np.corrcoef(np.argsort(np.argsort(_sl)), np.argsort(np.argsort(_om)))[0, 1])
axb.set_title("association profiles,\narm against arm", fontsize=6.4, pad=3)
axb.text(0.04, 0.96, "$r=%.2f$. The two arms\nrank pathways alike even\nthough their own scores\nshare only $\\rho=%.2f$" % (r_between, rho_arms), transform=axb.transAxes, fontsize=5.2, va="top",
         color=GREY, linespacing=1.35)
panel("b", 0.655, 0.980)

# =========================================================================== c  ranked |rho|
order = np.argsort(-np.abs(rs))
axc.plot(np.arange(len(P)), np.abs(rs)[order], lw=1.1, color=C["slide"],
         linestyle=style.ROLE_LINESTYLE["slide"], marker=style.ROLE_MARKER["slide"], markevery=20,
         ms=3.0, label="slide arm")
axc.plot(np.arange(len(P)), np.sort(np.abs(ro))[::-1], lw=1.1, color=C["omics"],
         linestyle=style.ROLE_LINESTYLE["omics"], marker=style.ROLE_MARKER["omics"], markevery=20,
         ms=3.0, label="omics arm")
axc.set_xlabel("pathways, ranked by $|\\rho|$ with that arm", fontsize=6.4)
axc.set_ylabel("$|\\rho|$", fontsize=6.4)
axc.set_xlim(0, len(P))
axc.legend(fontsize=5.6, loc="upper right", handletextpad=0.4, borderpad=0.15)
axc.set_title("how far the association extends", fontsize=6.4, pad=3)
panel("c", 0.010, 0.470)

# =========================================================================== d  the same, binned
bins = np.linspace(0, max(np.abs(rs).max(), np.abs(ro).max()) * 1.02, 13)
_, _, _hist_patches = axd.hist([np.abs(rs), np.abs(ro)], bins=bins, color=[C["slide"], C["omics"]])
for _patches, _role in zip(_hist_patches, ("slide", "omics")):
    for _p in _patches:
        _p.set_hatch(style.ROLE_HATCH[_role])
        _p.set_edgecolor(style.ROLE_EDGE[_role])
        _p.set_linewidth(0.4)
axd.set_xlabel("$|\\rho|$", fontsize=6.4)
axd.set_ylabel("pathways", fontsize=6.4)
axd.set_title("distribution of $|\\rho|$", fontsize=6.4, pad=3)
panel("d", 0.655, 0.470)

assert n_bh == 0, ("panel a is written around NO pathway clearing BH; %d do, so the caption and "
                   "the Results sentence both need rewriting" % n_bh)

style.save(fig, os.path.join(HERE, "figS4_pathway_landscape.pdf"))
print("wrote figS4_pathway_landscape.pdf   %.2f x %.2f in" % (FIGW, FIGH))
print("  %d pathways, %d clear BH q<0.10, r between the two rho columns %.3f"
      % (len(P), n_bh, r_between))
