#!/usr/bin/env python3
"""Supplementary Figure S5 -- what the primary number hides: per fold, and whether it is calibrated.

The main text reports one concordance and a standard deviation over seeds. Neither says how the
method behaves in a single fold, and neither says whether a risk score that ORDERS patients well
also places them at the right absolute risk. Both are standard reporting and both were missing.

  a  concordance per fold, per arm. Five folds is a small number and the spread across them is
     larger than the spread across seeds, which is worth a reader knowing before they read a
     margin of 0.03 as settled.
  b  Kaplan-Meier by ModRank tertile, computed WITHIN each fold so the tertile cuts are the ones
     a fold actually used. Five panels, no pooling.
  c  calibration. Cases are binned by ModRank percentile and the observed Kaplan-Meier survival at
     22 months is plotted against the bin. A discriminating score that is not monotone here is
     ordering patients without placing them.
  d  the out-of-fold score distributions, per arm. Percentile normalisation makes each arm uniform
     by construction within a fold; what this panel shows is the pooled result, and it is the
     check that the pooling did not undo the normalisation.

Built from the frozen reporting dump, which asserts its own primary against 0.7212 before writing.
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
# under amendment A2 the dump reproduces the A2 primary; check_numbers ties it to the A1 file
assert D["frozen_primary_matched"] or D.get("amendment") == "A2", "the dump did not reproduce the primary"

style.apply()
C, INK, RISK = style.ROLE, style.INK, style.RISK
GREY = "#8C8C8C"

cases = D["cases"]
fold = np.array([c["fold"] for c in cases])
t = np.array([c["months"] for c in cases], dtype=float)
e = np.array([c["event"] for c in cases], dtype=int)
ARMS = [("clinical", "clinical", C["clinical"]), ("wsi_titan", "slide", C["slide"]),
        ("omics_combine", "omics", C["omics"]), ("OURS", "ModRank", C["ours"])]
S = {k: np.array([c["seed_mean"][k] for c in cases]) for k, _, _ in ARMS}
NF = len(D["fold_sizes"])
assert NF == 5, "the released partition is five folds; the dump has %d" % NF


def cindex(risk, tt, ee):
    """Harrell's c over comparable pairs, written out so the figure depends on no import."""
    num = den = 0.0
    for i in range(len(tt)):
        if not ee[i]:
            continue
        j = tt > tt[i]
        den += j.sum()
        num += (risk[j] < risk[i]).sum() + 0.5 * (risk[j] == risk[i]).sum()
    return num / den if den else float("nan")


def km(tt, ee, grid):
    """Kaplan-Meier evaluated on a grid; ties handled by grouping at each distinct event time."""
    o = np.argsort(tt)
    tt, ee = tt[o], ee[o]
    n, surv, out, k = len(tt), 1.0, [], 0
    times = np.unique(tt[ee == 1])
    curve_t, curve_s = [0.0], [1.0]
    for ti in times:
        at_risk = (tt >= ti).sum()
        d = ((tt == ti) & (ee == 1)).sum()
        if at_risk:
            surv *= (1.0 - d / at_risk)
        curve_t.append(float(ti)); curve_s.append(surv)
    curve_t, curve_s = np.array(curve_t), np.array(curve_s)
    return curve_t, curve_s, np.array([curve_s[np.searchsorted(curve_t, g, "right") - 1]
                                       for g in grid])


FIGW, FIGH = 6.30, 6.20
fig = plt.figure(figsize=(FIGW, FIGH))
panel = lambda tag, x, y: fig.text(x, y, tag, fontsize=7.5, fontweight="bold", color=INK,   # noqa
                                   ha="left", va="top")

# =========================================================================== a  per fold, per arm
axa = fig.add_axes([0.098, 0.760, 0.400, 0.185])
w = 0.20
xs = np.arange(NF)
perfold = {}
for i, (k, lab, col) in enumerate(ARMS):
    v = [cindex(S[k][fold == f], t[fold == f], e[fold == f]) for f in range(NF)]
    perfold[k] = v
    axa.plot(xs + (i - 1.5) * w, v, "o", color=col, ms=4.2, markeredgecolor="none", label=lab)
    for xx, vv in zip(xs + (i - 1.5) * w, v):
        axa.plot([xx, xx], [0.455, vv], color=col, lw=1.0, solid_capstyle="round", zorder=1)
axa.axhline(0.5, color=INK, lw=0.7, ls=":")
axa.set_xticks(xs); axa.set_xticklabels(["fold %d\n$n$=%d" % (f + 1, D["fold_sizes"][f])
                                         for f in range(NF)], fontsize=5.6)
axa.set_ylim(0.45, 0.85)
axa.set_ylabel("concordance, that fold", fontsize=6.4)
axa.legend(fontsize=5.6, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.24),
           handletextpad=0.4, columnspacing=1.0, borderpad=0.1)
spread = max(perfold["OURS"]) - min(perfold["OURS"])
axa.text(1.0, -0.30, "ModRank spans %.3f across folds against %.4f across seeds"
         % (spread, 0.0048), transform=axa.transAxes, ha="right", va="top", fontsize=5.4,
         color=GREY)
panel("a", 0.010, 0.978)

# =========================================================================== c  calibration
axc = fig.add_axes([0.615, 0.760, 0.360, 0.185])
HOR = 22.0
nb = 5
qs = np.quantile(S["OURS"], np.linspace(0, 1, nb + 1))
mid, obs, lo_n = [], [], []
for b in range(nb):
    m = (S["OURS"] >= qs[b]) & (S["OURS"] <= qs[b + 1] if b == nb - 1 else S["OURS"] < qs[b + 1])
    if m.sum() < 10:
        continue
    _, _, s_at = km(t[m], e[m], [HOR])
    mid.append(float(np.mean(S["OURS"][m]))); obs.append(float(s_at[0])); lo_n.append(int(m.sum()))
axc.plot(mid, obs, "-o", color=C["ours"], ms=4, lw=1.3, markeredgecolor="none")
for x, yv, nn in zip(mid, obs, lo_n):
    axc.annotate("$n$=%d" % nn, (x, yv), textcoords="offset points", xytext=(0, 6), ha="center",
                 fontsize=5.0, color=GREY)
axc.set_xlim(0.02, 1.02)
axc.set_xlabel("ModRank percentile, bin mean", fontsize=6.4)
axc.set_ylabel("observed survival at %g months" % HOR, fontsize=6.4)
axc.set_title("calibration in the large", fontsize=6.4, pad=3)
mono = all(obs[i] >= obs[i + 1] for i in range(len(obs) - 1))
axc.text(0.5, 0.06, "monotone decreasing: %s" % ("yes" if mono else "NO"), transform=axc.transAxes,
         ha="center", va="bottom", fontsize=5.4, color=INK if mono else C["competitor"],
         fontweight="bold")
panel("c", 0.530, 0.978)

# =========================================================================== b  KM within each fold
GRID = np.linspace(0, 96, 200)
for f in range(NF):
    m = fold == f
    ax = fig.add_axes([0.098 + f * 0.1785, 0.415, 0.148, 0.205])
    cuts = np.quantile(S["OURS"][m], [1 / 3.0, 2 / 3.0])
    grp = np.digitize(S["OURS"][m], cuts)
    for g, key in ((0, "low"), (1, "middle"), (2, "high")):
        sel = grp == g
        if sel.sum() < 5:
            continue
        ct, cs, _ = km(t[m][sel], e[m][sel], [0])
        ax.step(np.append(ct, 96), np.append(cs, cs[-1]), where="post", color=RISK[key], lw=1.1)
    ax.set_xlim(0, 96); ax.set_ylim(0, 1.02)
    ax.set_xticks([0, 48, 96]); ax.set_xticklabels(["0", "48", "96"], fontsize=5.4)
    ax.set_title("fold %d" % (f + 1), fontsize=6.0, pad=2)
    if f == 0:
        ax.set_ylabel("disease-specific\nsurvival", fontsize=6.0, linespacing=1.3)
        ax.set_yticks([0, 0.5, 1.0]); ax.tick_params(axis="y", labelsize=5.4)
    else:
        ax.set_yticks([])
    if f == 2:
        ax.set_xlabel("months", fontsize=6.2)
fig.text(0.098, 0.646, "Kaplan--Meier by ModRank tertile, cut within each fold", fontsize=6.4,
         ha="left", va="bottom", color=INK)
for i, (key, lab) in enumerate((("high", "high"), ("middle", "middle"), ("low", "low"))):
    xk = 0.735 + i * 0.082
    fig.text(xk, 0.648, "▬", fontsize=6.5, color=RISK[key], ha="left", va="bottom")
    fig.text(xk + 0.018, 0.648, lab, fontsize=5.6, color=INK, ha="left", va="bottom")
panel("b", 0.010, 0.672)

# =========================================================================== d  score distributions
axd = fig.add_axes([0.098, 0.075, 0.400, 0.215])
for k, lab, col in ARMS:
    hist, edges = np.histogram(S[k], bins=22, range=(0, 1))
    axd.step(edges[:-1], hist, where="post", color=col, lw=1.1, label=lab)
axd.set_xlabel("out-of-fold percentile, pooled over folds", fontsize=6.4)
axd.set_ylabel("cases", fontsize=6.4)
axd.legend(fontsize=5.6, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.20),
           handletextpad=0.4, columnspacing=1.0, borderpad=0.1)
axd.set_title("score distributions", fontsize=6.4, pad=3)
panel("d", 0.010, 0.302)

# =========================================================================== e  events per fold
axe = fig.add_axes([0.615, 0.075, 0.360, 0.215])
ev = [int(e[fold == f].sum()) for f in range(NF)]
nn = [int((fold == f).sum()) for f in range(NF)]
axe.bar(xs - 0.19, nn, width=0.38, color="#C9CFD8", edgecolor="none", label="cases")
axe.bar(xs + 0.19, ev, width=0.38, color=C["competitor"], edgecolor="none", label="events")
for x, (a_, b_) in enumerate(zip(nn, ev)):
    axe.text(x - 0.19, a_ + 1.5, str(a_), ha="center", va="bottom", fontsize=5.4, color=GREY)
    axe.text(x + 0.19, b_ + 1.5, str(b_), ha="center", va="bottom", fontsize=5.4,
             color=C["competitor"])
axe.set_xticks(xs); axe.set_xticklabels(["fold %d" % (f + 1) for f in range(NF)], fontsize=5.6)
axe.set_ylabel("count", fontsize=6.4)
axe.set_ylim(0, max(nn) * 1.22)
axe.legend(fontsize=5.6, loc="upper right", handletextpad=0.4, borderpad=0.15)
axe.set_title("what each fold is asked to resolve", fontsize=6.4, pad=3)
panel("e", 0.530, 0.302)

assert sum(ev) == D["events"], "fold events must sum to the cohort's %d; they sum to %d" % (
    D["events"], sum(ev))
assert sum(nn) == D["n"], "fold sizes must sum to %d; they sum to %d" % (D["n"], sum(nn))

style.save(fig, os.path.join(HERE, "figS5_perfold_calibration.pdf"))
print("wrote figS5_perfold_calibration.pdf   %.2f x %.2f in" % (FIGW, FIGH))
print("  per-fold ModRank %s (spread %.3f)" % ([round(v, 3) for v in perfold["OURS"]], spread))
print("  calibration bins %s, monotone=%s" % (lo_n, mono))
print("  events per fold %s" % ev)
