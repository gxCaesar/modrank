#!/usr/bin/env python3
"""Figure 5 -- all three modalities across the whole cohort, on one patient axis.

WHAT MAKES THIS DIFFERENT FROM FIGURE 3. Figure 3 says what each modality is and what it is worth.
This one shows all 359 patients at once, ordered left to right by ModRank risk, with every modality
drawn against that same ordering. A reader can then see, without taking anyone's word for it, that
the ordering is not a restatement of stage and not a restatement of any single transcriptomic
programme.

  a  histology. Twenty-four slides sampled at even risk quantiles, at matched physical resolution.
  b  transcriptome. The forty pathways with the strongest univariate association, z-scored, one
     column per patient. The arm itself reads all 275; forty is what fits and they are the forty a
     reader would ask to see.
  c  pathologic stage, one cell per patient. If the model were re-deriving stage this row would be
     sorted. It is not.
  d  the three arm percentiles and where the deaths fall.

Built from the frozen reporting dump, which asserts its own primary against 0.7212 before writing.
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import style                                                             # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
D = json.load(open(os.path.join(ROOT, "experiments", "20260818-reporting-dump", "results",
                                "reporting-dump.json")))
assert D["frozen_primary_matched"], "the dump did not reproduce the frozen primary"

style.apply()
C, INK, LAD, RISK = style.ROLE, style.INK, style.LADDER, style.RISK
GREY = "#8C8C8C"

cases = sorted(D["cases"], key=lambda c: c["seed_mean"]["OURS"])
n = len(cases)
risk = np.array([c["seed_mean"]["OURS"] for c in cases])
Z = np.array([c["pathway_z"] for c in cases]).T          # pathways x cases
stage = [c["stage"] for c in cases]
event = np.array([c["event"] for c in cases])
STAGES = ["Stage I", "Stage II", "Stage III", "Stage IV", "unknown"]
SCOL = {"Stage I": LAD[5], "Stage II": LAD[4], "Stage III": LAD[1], "Stage IV": LAD[0],
        "unknown": "#E4E7EC"}

# the ladder is an ORDERED palette and stage is an ordered variable, which is the one case where
# lightness may carry the meaning. Checked here rather than assumed.
assert all(s in SCOL for s in set(stage)), "an unmapped stage level would render as nothing"

FIGW, FIGH = 6.785, 5.60   # 172.3 mm (PI correction 2026-09-27: left as is; do not widen to 180 mm)
fig = plt.figure(figsize=(FIGW, FIGH))
L, W = 0.088, 0.892
panel = lambda t, x, y: fig.text(x, y, t, fontsize=7.5, fontweight="bold", color=INK,   # noqa
                                 ha="left", va="top")

# =========================================================================== a  histology strip
STRIP = os.path.join(HERE, "assets", "risk_strip")
have_strip = os.path.isdir(STRIP) and os.path.exists(os.path.join(STRIP, "index.json"))
if have_strip:
    idx = json.load(open(os.path.join(STRIP, "index.json")))["cases"]
    pos = {c["case_id"]: i for i, c in enumerate(cases)}
    shown = [r for r in idx if r["case_id"] in pos]
    shown.sort(key=lambda r: pos[r["case_id"]])
    k = len(shown)
    cw = W / k
    for j, r in enumerate(shown):
        img = plt.imread(os.path.join(STRIP, r["case_id"] + ".png"))[..., :3]
        ax = fig.add_axes([L + j * cw + 0.0012, 0.800, cw - 0.0024, 0.135])
        ax.imshow(img, interpolation="lanczos", aspect="auto")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    fig.text(L, 0.945, "%d slides at even risk quantiles" % k, fontsize=6.0, color=GREY,
             ha="left", va="bottom")
else:
    k = 0
panel("a", 0.010, 0.972)

# =========================================================================== b  transcriptome
axb = fig.add_axes([L, 0.478, W, 0.290])
cmap = LinearSegmentedColormap.from_list("bwr", ["#4C72B0", "#FFFFFF", "#C44E52"])
axb.imshow(np.clip(Z, -2.5, 2.5), aspect="auto", cmap=cmap, vmin=-2.5, vmax=2.5,
           interpolation="nearest")
axb.set_xticks([]); axb.set_yticks([])
axb.set_ylabel("40 pathways", fontsize=6.4)
cax = fig.add_axes([L + W - 0.085, 0.786, 0.085, 0.009])
cax.imshow(np.linspace(-2.5, 2.5, 128)[None, :], aspect="auto", cmap=cmap)
cax.set_xticks([0, 127]); cax.set_xticklabels(["$-2.5$", "$+2.5$"], fontsize=5.0)
cax.set_yticks([]); cax.tick_params(length=0, pad=1)
cax.text(-0.06, 0.5, "$z$", transform=cax.transAxes, ha="right", va="center", fontsize=5.4,
         color=GREY)
panel("b", 0.010, 0.782)

# =========================================================================== c  pathologic stage
axc = fig.add_axes([L, 0.400, W, 0.045])
axc.set_xlim(0, n); axc.set_ylim(0, 1); axc.axis("off")
for i, s in enumerate(stage):
    axc.add_patch(Rectangle((i, 0), 1.0, 1.0, facecolor=SCOL[s], edgecolor="none"))
fig.text(L - 0.006, 0.4225, "stage", fontsize=6.4, color=INK, ha="right", va="center")
for j, s in enumerate([x for x in STAGES if x in set(stage)]):
    xk = L + 0.02 + j * 0.072
    fig.text(xk, 0.372, "\u25a0", fontsize=6.5, color=SCOL[s], ha="left", va="bottom")
    fig.text(xk + 0.016, 0.373, s.replace("Stage ", ""), fontsize=5.6, color=INK, ha="left",
             va="bottom")
panel("c", 0.010, 0.452)

# =========================================================================== d  the arms
axd = fig.add_axes([L, 0.140, W, 0.195])
for key, lab, role in (("clinical", "clinical", "clinical"), ("wsi_titan", "slide", "slide"),
                      ("omics_combine", "transcriptome", "omics")):
    v = np.array([c["seed_mean"][key] for c in cases])
    axd.scatter(np.arange(n), v, s=2.2, marker=style.ROLE_MARKER[role], color=C[role], alpha=0.55,
                edgecolor="none", label=lab)
axd.plot(np.arange(n), risk, color=C["ours"], lw=1.4, linestyle=style.ROLE_LINESTYLE["ours"],
         zorder=5, label="ModRank")
axd.set_xlim(0, n); axd.set_ylim(-0.02, 1.02)
axd.set_ylabel("out-of-fold percentile", fontsize=6.4)
axd.set_xticks([])
axd.legend(fontsize=5.6, ncol=4, loc="upper center", bbox_to_anchor=(0.72, 1.30),
           handletextpad=0.4, columnspacing=1.2, borderpad=0.1, markerscale=2.6)
axe = fig.add_axes([L, 0.098, W, 0.030])
axe.set_xlim(0, n); axe.set_ylim(0, 1); axe.axis("off")
for i in np.flatnonzero(event):
    axe.plot([i + 0.5, i + 0.5], [0.15, 0.95], color=C["competitor"], lw=0.55)
fig.text(L - 0.006, 0.113, "deaths", fontsize=6.4, color=C["competitor"], ha="right", va="center")
panel("d", 0.010, 0.348)

# the claim panel c exists to support: the ordering is NOT a re-derivation of stage
sr = np.array([STAGES.index(s) for s in stage], dtype=float)
ok = sr < 4
rho = float(np.corrcoef(np.argsort(np.argsort(risk[ok])), np.argsort(np.argsort(sr[ok])))[0, 1])
assert abs(rho) < 0.6, ("panel c is drawn to show that the risk ordering is not stage restated; "
                        "Spearman between them is %.3f" % rho)

fig.text(L, 0.040, "359 patients, ordered left to right by ModRank risk. Every row is drawn "
         "against that one ordering.", fontsize=6.0, color=GREY, ha="left", va="bottom")
fig.text(L + W, 0.040, "Spearman between the ordering and stage: %.2f" % rho, fontsize=6.0,
         color=GREY, ha="right", va="bottom")
fig.text(L, 0.018, "low risk", fontsize=6.0, color=GREY, ha="left", va="bottom")
fig.text(L + W, 0.018, "high risk", fontsize=6.0, color=GREY, ha="right", va="bottom")


style.save(fig, os.path.join(HERE, "fig5_cohort_modalities.pdf"))
print("wrote fig5_cohort_modalities.pdf   %.3f x %.3f in" % (FIGW, FIGH))
print("  %d cases, %d pathways, %d slides in the strip, Spearman risk vs stage %.3f"
      % (n, Z.shape[0], k, rho))
