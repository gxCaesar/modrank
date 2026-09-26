#!/usr/bin/env python3
"""Figure 3 -- the three modalities: what each one is, and what each one carries.

WHY THIS FIGURE EXISTS. The manuscript is a multimodal paper whose display items, before this one,
never showed the reader what the modalities ARE. The method schematic named them, the result
figures scored them, and nothing in between said what a whole-slide image contributes that a
pathology report does not, or how far the three scores are from being the same score written three
ways. Every panel here answers that from a measurement already in the frozen evidence, and several
of those measurements were sitting unused.

  a  the three modalities as they enter the model: what is measured, how it is represented, how
     many dimensions that representation carries, and what it is worth alone. The strip beside the
     transcriptome is the real thing -- 275 pathway groups, the 59 the slide arm tracks at
     BH q<0.05 drawn dark -- and the stage bar beside the clinical block is the composition of the
     field the model actually reads, the incumbent's four AJCC stage groups over 357 valued cases.
     The five-level composition in Figure 1d is a different released file and is labelled as one.
  b  seven representations across the three families, on identical cases and identical folds.
     Choosing WITHIN a modality moves the number as much as choosing between modalities, which is
     the fact that makes a single-representation multimodal comparison hard to interpret.
  c  how far the three risk scores are from being the same score. Spearman between arms, out of
     fold. Image and clinical share 0.2445 (amended clinical block): nearly orthogonal, which is why
     the equal-weight average
     works and why it is not the representation that limits the fusion.
  d  THE PANEL THIS FIGURE IS FOR. Restrict the comparable pairs to those the clinical model cannot
     separate and score the arms again. On the tightest 1,213 pairs of 24,219 the clinical arm is a
     coin flip at 0.4959 and the slide arm holds 0.6220. That is the whole-slide image earning its
     place, measured rather than asserted, and it is the answer to "what does the image add".
  e  discrimination over time, per modality, at three horizons plus the integrated Brier score.
  f  what capacity buys per modality, against the permuted-stratum control that voids it.
  g  how each modality separates its own tertiles: log-rank against the group medians it produces.

NOT SHOWN, AND DELIBERATELY. No slide thumbnail. TCGA diagnostic slides are open access, but this
repository holds no image bytes, only the frozen 768-dimensional embedding, and a figure asserting
a slide it did not read would be the one thing this paper cannot afford. The tissue token in a is
drawn as a schematic and is labelled as one.
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import style                                                            # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "experiments", "20260817-blca-confirm", "results")
S5 = os.path.join(ROOT, "development", "s5-results")
load = lambda d, n: json.load(open(os.path.join(d, n)))                  # noqa: E731

# Panels c and d read the AMENDED atlas (analysis/s19b_modality_atlas_amended.py), which re-runs the
# development atlas with the manuscript's clinical block. The original atlas.json used the
# pre-amendment 23-column block, and so did every clinical value this figure drew until 2026-09-11.
ATLAS = json.load(open(os.path.join(ROOT, "experiments", "20260911-blca-posthoc", "results",
                                    "modality-atlas-amended.json")))
A1 = load(RES, "amendment-A1-clinical-provenance.json")
CLIN_AMENDED = A1["arms_seed0"]["clinical_age_sex_stage_from_DIMAF_file"]
C6 = load(S5, "c6.json")
OF = load(S5, "oracle_fix.json")
SM = load(RES, "survival-metrics.json")
WG = load(RES, "why-grade-fails.json")
BIO = load(RES, "biology.json")

# Two image assets, retrieved once from the GDC tile service and committed beside this builder so
# the figure is reproducible without the network. Case TCGA-HQ-A2OF, a diagnostic slide that is in
# this study's own 359-case partition. TCGA open-access data carry no restriction on analysis or
# publication beyond the prohibition on reidentification; the acknowledgement the source asks for
# is in the manuscript's Acknowledgements, and the figure caption names the case.
SLIDE = plt.imread(os.path.join(HERE, "assets", "blca_slide_thumb.png"))[..., :3]
TILE = plt.imread(os.path.join(HERE, "assets", "blca_slide_tile.png"))[..., :3]
# The locator box is COMPUTED by the retrieval script and read here, not tuned by eye. That
# script asserts the box lands on tissue rather than on slide background.
ASSET = json.load(open(os.path.join(HERE, "assets", "blca_slide_asset.json")))
SLIDE_CASE = ASSET["case"]
_b = ASSET["locator_box_px"]
LOC = (_b["x"], _b["y"], _b["w"], _b["h"])
assert ASSET["locator_box_tissue_fraction"] > 0.5, (
    "the locator box in panel a must point at tissue; the asset records %.2f"
    % ASSET["locator_box_tissue_fraction"])

style.apply()
C, INK, LAD = style.ROLE, style.INK, style.LADDER
GREY, FAINT = "#8C8C8C", "#C9CFD8"

# One colour per modality family, held for the whole figure. A reader who learns the mapping in a
# must not have to relearn it in g.
MOD = {"slide": C["slide"], "omics": C["omics"], "clinical": C["clinical"], "ours": C["ours"]}
# the second channel, held for the whole figure alongside MOD so a role's marker and linestyle
# never drift from its colour (CIEDE2000 2.04 between ours and clinical under simulated CVD)
MK = {fam: style.ROLE_MARKER[fam] for fam in MOD}
LS = {fam: style.ROLE_LINESTYLE[fam] for fam in MOD}

FIGW, FIGH = style.width(6.785), 6.45
fig = plt.figure(figsize=(FIGW, FIGH))


def panel(tag, x, y):
    fig.text(x, y, tag, fontsize=7.5, fontweight="bold", color=INK, ha="left", va="top")


# ============================================================ a  the three modalities as they enter
axa = fig.add_axes([0.0, 0.735, 1.0, 0.235])
axa.set_xlim(0, 1)
axa.set_ylim(0, 1)
axa.axis("off")

N_PATH = BIO["B2_survival_association"]["pathways_tested"]
N_TRACKED = BIO["B3_what_the_image_arm_tracks"]["pathways_with_BH_q_below_0.05"]
STAGE = {k: v["n"] for k, v in load(RES, "per-stage-subgroup.json")["strata"].items()}
ALONE = {"slide": C6["single_arms"]["wsi_titan"], "omics": C6["single_arms"]["omics_combine"],
         "clinical": CLIN_AMENDED}

CARDS = [
    ("slide", "Whole-slide H&E", "one 768-dim slide embedding per case,\nmean-pooled over a "
     "case's slides", "768 dims"),
    ("omics", "Bulk transcriptome", "4,999 genes averaged into the\nbenchmark's own 275 pathway "
     "groups", "275 dims"),
    ("clinical", "Pathology report", "age, sex and pathologic stage,\none-hot for the categorical",
     "6 columns"),
]
CW, GAP = 0.288, 0.028
X0 = 0.055
for i, (key, title, body, dims) in enumerate(CARDS):
    x = X0 + i * (CW + GAP)
    axa.add_patch(FancyBboxPatch((x, 0.05), CW, 0.86, boxstyle="round,pad=0.006,rounding_size=0.01",
                                 facecolor="#FFFFFF", edgecolor=MOD[key], linewidth=1.0,
                                 transform=axa.transAxes, zorder=1))
    axa.add_patch(Rectangle((x, 0.80), CW, 0.11, facecolor=MOD[key], edgecolor="none",
                            transform=axa.transAxes, zorder=2))
    axa.text(x + CW / 2, 0.855, title, ha="center", va="center", fontsize=6.8, color="#FFFFFF",
             fontweight="bold", transform=axa.transAxes, zorder=3)

    tx, ty, tw, th = x + 0.018, 0.500, CW - 0.036, 0.255     # the token box, one per modality
    if key == "slide":
        # The real thing. A diagnostic slide from a case in this cohort, read through the GDC's
        # own tile service, with one tile at the magnification the benchmark's slide pipeline
        # works at. Drawing a schematic here would have been the one place this paper implied a
        # measurement it did not make.
        axw = fig.add_axes([0.098, 0.8525, 0.108, 0.0599]); axw.axis("off")
        axw.imshow(SLIDE, interpolation="lanczos")
        for sp in axw.spines.values():
            sp.set_visible(False)
        axw.add_patch(Rectangle((LOC[0], LOC[1]), LOC[2], LOC[3], facecolor="none",
                                edgecolor=C["competitor"], linewidth=0.9, zorder=5))
        axt = fig.add_axes([0.225, 0.8525, 0.0755, 0.0599]); axt.axis("off")
        axt.imshow(TILE, interpolation="lanczos")
        for sp in axt.spines.values():
            sp.set_visible(False)
        axt.add_patch(Rectangle((0, 0), TILE.shape[1] - 1, TILE.shape[0] - 1, facecolor="none",
                                edgecolor=C["competitor"], linewidth=0.9,
                                transform=axt.transData, zorder=5))
        fig.text(0.2135, 0.8815, "$\\rightarrow$", ha="center", va="center", fontsize=7,
                 color=MOD[key])
        axa.text(tx + tw / 2, ty - 0.022, "whole slide, then one tile, then the frozen encoder",
                 ha="center", va="top", fontsize=5.4, color=GREY, style="italic",
                 transform=axa.transAxes, zorder=5)
    elif key == "omics":
        # the pathways the slide arm tracks (BH q < 0.05), recomputed and count-checked in
        # pathway_glyph.py; the count printed under the glyph is the same N_TRACKED
        from pathway_glyph import slide_tracked
        dark = set(np.flatnonzero(slide_tracked()).tolist())
        per = 25
        rows = int(np.ceil(N_PATH / per))
        for j in range(N_PATH):
            r, c_ = divmod(j, per)
            gx = tx + 0.004 + c_ * (tw - 0.008) / per
            gy = ty + th * 0.90 - (r + 1) * (th * 0.86) / rows
            axa.add_patch(Rectangle((gx, gy), (tw - 0.008) / per * 0.82, (th * 0.86) / rows * 0.74,
                                    facecolor=MOD[key] if j in dark else "#F6DDCC",
                                    edgecolor="none", transform=axa.transAxes, zorder=4))
        axa.text(tx + tw / 2, ty - 0.022, "%d pathway groups, %d tracked by the slide arm"
                 % (N_PATH, N_TRACKED), ha="center", va="top", fontsize=5.4, color=GREY,
                 style="italic", transform=axa.transAxes, zorder=5)
    else:
        # The stage composition the clinical block ACTUALLY holds, which is not the one the
        # five-study entropy table uses. Two released files carry a stage field for these patients:
        # the benchmark's own clinical file (five levels including a stage 0, 334 valued) and the
        # incumbent's split files (four AJCC stage groups, 357 valued). The amendment rebuilt every
        # clinical variable from the SECOND, so this panel, which shows what enters the model, reads
        # the second as well. Drawing the first here labelled "what the clinical block holds" is the
        # defect an external review found on 2026-09-12.
        order = ["Stage I", "Stage II", "Stage III", "Stage IV"]
        tot = sum(STAGE.values())
        cx = tx
        for j, lev in enumerate([o for o in order if o in STAGE]):
            w = tw * STAGE[lev] / tot
            axa.add_patch(Rectangle((cx, ty + th * 0.40), w, th * 0.34, facecolor=LAD[j % len(LAD)],
                                    edgecolor="#FFFFFF", linewidth=0.5,
                                    transform=axa.transAxes, zorder=4))
            if w > 0.030:
                axa.text(cx + w / 2, ty + th * 0.57, lev.replace("Stage ", ""), ha="center",
                         va="center", fontsize=5.8, color=INK, transform=axa.transAxes, zorder=5)
            cx += w
        axa.text(tx + tw / 2, ty - 0.022,
                 "pathologic stage, four groups, %d%% in the largest" % round(100 * max(
                     STAGE.values()) / tot), ha="center", va="top", fontsize=5.4,
                 color=GREY, style="italic", transform=axa.transAxes, zorder=5)

    axa.text(x + CW / 2, 0.400, body, ha="center", va="top", fontsize=6.2, color=INK,
             linespacing=1.45, transform=axa.transAxes, zorder=5)
    axa.text(x + 0.018, 0.130, dims, ha="left", va="center", fontsize=6.4, color=MOD[key],
             fontweight="bold", transform=axa.transAxes, zorder=5)
    axa.text(x + CW - 0.018, 0.130, "alone  %.4f" % ALONE[key], ha="right", va="center",
             fontsize=6.4, color=INK, transform=axa.transAxes, zorder=5)

axa.annotate("", xy=(X0 + 3 * CW + 2 * GAP + 0.030, 0.115), xytext=(X0 + 3 * CW + 2 * GAP + 0.004,
             0.115), arrowprops=dict(arrowstyle="-", color="#FFFFFF", lw=0.1),
             transform=axa.transAxes)
panel("a", 0.010, 0.982)

# ============================================================ b  seven representations, three families
axb = fig.add_axes([0.175, 0.545, 0.375, 0.155])
REPS = [("wsi_titan", "TITAN", "slide"), ("wsi_chief_mean", "CHIEF, mean", "slide"),
        ("wsi_chief_dispersion", "CHIEF, dispersion", "slide"),
        ("omics_xena", "Xena expression", "omics"), ("omics_combine", "SurvPath pathways", "omics"),
        ("omics_hallmarks", "Hallmark sets", "omics"),
        ("clinical", "age + sex + stage", "clinical")]
vals = [(lab, CLIN_AMENDED if k == "clinical" else C6["single_arms"][k], fam)
        for k, lab, fam in REPS]
vals.sort(key=lambda r: r[1])
y = np.arange(len(vals))
# a dot plot: the axis starts at 0.50, and a bar drawn from there would overstate every gap
for yi, (_, v, fam) in zip(y, vals):
    axb.plot([0.50, v], [yi, yi], color="#E3E6EC", lw=0.8, zorder=1)
    axb.scatter([v], [yi], s=26, color=MOD[fam], marker=MK[fam], edgecolor="none", zorder=3)
for yi, (lab, v, fam) in zip(y, vals):
    axb.text(v + 0.005, yi, "%.4f" % v, va="center", fontsize=5.9, color=INK)
axb.axvline(0.5, color=INK, lw=0.7, zorder=3)
axb.set_yticks(y)
axb.set_yticklabels([lab for lab, _, _ in vals], fontsize=6.0)
axb.set_xlim(0.50, 0.72)
axb.set_xlabel("concordance, that representation alone", fontsize=6.3)
axb.set_title("seven representations, three families",
              fontsize=6.3, pad=11)
spread_img = max(v for _, v, f in vals if f == "slide") - min(v for _, v, f in vals if f == "slide")
spread_om = max(v for _, v, f in vals if f == "omics") - min(v for _, v, f in vals if f == "omics")
axb.text(0.5, 1.015, "spread within the image family %.3f, within the omics family %.3f"
         % (spread_img, spread_om), transform=axb.transAxes, ha="center", va="bottom",
         fontsize=5.6, color=GREY)
panel("b", 0.010, 0.722)

# ============================================================ c  are the three scores the same score
axc = fig.add_axes([0.700, 0.545, 0.175, 0.155])
R = ATLAS["redundancy"]
# NOTE the middle row is the INCUMBENT's joint score, not a transcriptome-only arm: the frozen
# redundancy measurement was taken against SurvPath, which reads slide and transcriptome together.
# Labelling it "omics" would be a quiet misstatement of what was correlated with what.
NAMES = ["slide", "SurvPath", "clinical"]
NCOL = [MOD["slide"], C["competitor"], MOD["clinical"]]
RHO = np.array([[1.0, R["spearman_titan_vs_survpath"], R["spearman_titan_vs_clinical"]],
                [R["spearman_titan_vs_survpath"], 1.0, R["spearman_survpath_vs_clinical"]],
                [R["spearman_titan_vs_clinical"], R["spearman_survpath_vs_clinical"], 1.0]])
axc.imshow(RHO, cmap="Blues", vmin=0.0, vmax=1.0, aspect="equal")
for i in range(3):
    for j in range(3):
        axc.text(j, i, "1" if i == j else "%.2f" % RHO[i, j], ha="center", va="center",
                 fontsize=6.6, color="#FFFFFF" if RHO[i, j] > 0.55 else INK,
                 fontweight="bold" if i != j else "normal")
axc.set_xticks(range(3)); axc.set_yticks(range(3))
axc.set_xticklabels(NAMES, fontsize=6.0)
axc.set_yticklabels(NAMES, fontsize=6.0)
for k in range(3):
    axc.get_xticklabels()[k].set_color(NCOL[k]); axc.get_yticklabels()[k].set_color(NCOL[k])
axc.set_title("Spearman between the arms", fontsize=6.3, pad=11)
axc.text(0.5, 1.015, "SurvPath = slide and transcriptome read jointly", transform=axc.transAxes,
         ha="center", va="bottom", fontsize=5.5, color=GREY)
axc.text(0.5, -0.20, "%.3f between image and clinical:\nnearly orthogonal"
         % R["spearman_titan_vs_clinical"], transform=axc.transAxes, ha="center", va="top",
         fontsize=5.7, color=GREY, linespacing=1.45)
panel("c", 0.600, 0.722)

# ============================================================ d  what the image adds where the
# ============================================================    clinical model cannot separate
axd = fig.add_axes([0.105, 0.328, 0.845, 0.124])
CP = ATLAS["conditional_probe"]
KEYS = ["clinically_tied_q05", "clinically_tied_q10", "clinically_tied_q20",
        "clinically_tied_q40", "all_pairs"]
xs = np.arange(len(KEYS))
SERIES = [("clinical", "clinical", "clinical alone"),
          ("titan", "slide", "slide alone"),
          ("survpath", "omics", "slide + transcriptome, SurvPath"),
          ("titan+clinical", "ours", "slide + clinical, rank-averaged")]
for arm, fam, lab in SERIES:
    v = [CP[k]["arms"][arm] for k in KEYS]
    axd.plot(xs, v, color=MOD[fam], marker=MK[fam], linestyle=LS[fam], ms=3.6, lw=1.3, label=lab,
             zorder=3, markeredgecolor="none")
axd.axhline(0.5, color=INK, lw=0.7, ls=":", zorder=1)
axd.text(xs[0] - 0.34, 0.502, "chance", fontsize=5.5, color=INK, va="bottom", ha="left")
tie = CP["clinically_tied_q05"]["arms"]
axd.annotate("", xy=(0, tie["titan"]), xytext=(0, tie["clinical"]),
             arrowprops=dict(arrowstyle="<->", color=INK, lw=0.8, mutation_scale=7), zorder=4)
axd.text(0.085, tie["clinical"] + 0.012, "+%.3f where the\nreport says nothing"
         % (tie["titan"] - tie["clinical"]), fontsize=6.0, color=INK, va="bottom", ha="left",
         fontweight="bold", linespacing=1.35)
axd.set_xticks(xs)
axd.set_xticklabels(["%s\n%s pairs" % (lab, "{:,}".format(CP[k]["pairs"]))
                     for k, lab in zip(KEYS, ["tightest 5%", "10%", "20%", "40%",
                                              "all pairs"])], fontsize=6.0)
axd.set_xlim(-0.45, len(KEYS) - 0.55)
axd.set_ylim(0.48, 0.74)
axd.set_ylabel("concordance on the\nretained pairs", fontsize=6.3, linespacing=1.3)
axd.set_xlabel("comparable pairs retained, ordered by how far apart the CLINICAL model puts the two "
               "patients", fontsize=6.3, labelpad=2)
axd.legend(fontsize=5.9, loc="upper left", ncol=2, handletextpad=0.5, columnspacing=1.4,
           borderpad=0.2, handlelength=1.8)
# pad=13 put this title level with panel b's x-axis label, in the narrow gap between the two panel
# rows -- pad=4 keeps it above panel d's own axes instead, clear of the row above.
axd.set_title("concordance on clinically tied pairs", fontsize=6.4, pad=4,
              fontweight="bold")
panel("d", 0.010, 0.486)

# --- one modality key for the bottom row, so e and f need no legend of their own
axk = fig.add_axes([0.100, 0.012, 0.850, 0.014]); axk.axis("off")
axk.set_xlim(0, 1); axk.set_ylim(0, 1)
for i, (lab, fam) in enumerate([("clinical", "clinical"), ("slide", "slide"),
                                ("transcriptome", "omics"), ("all three", "ours")]):
    xk = 0.300 + i * 0.115
    axk.plot([xk, xk + 0.020], [0.5, 0.5], color=MOD[fam], lw=1.6, linestyle=LS[fam],
             solid_capstyle="round")
    axk.plot([xk + 0.010], [0.5], marker=MK[fam], color=MOD[fam], ms=3.0, markeredgecolor="none")
    axk.text(xk + 0.026, 0.5, lab, va="center", ha="left", fontsize=6.0, color=INK)

# ============================================================ e  discrimination over time
axe = fig.add_axes([0.100, 0.088, 0.215, 0.130])
HOR = SM["horizons_months"]
ARMS_E = [("clinical", "clinical"), ("wsi_titan", "slide"), ("omics", "omics"), ("OURS", "ours")]
for key, fam in ARMS_E:
    a = SM["arms"][key]
    axe.plot(HOR, [a["td_auc_8m"], a["td_auc_14m"], a["td_auc_22m"]], color=MOD[fam],
             marker=MK[fam], linestyle=LS[fam], ms=3.0, lw=1.2, markeredgecolor="none")
axe.set_xticks(HOR)
axe.set_xticklabels(["%.0f" % h for h in HOR], fontsize=6.0)
axe.set_xlabel("horizon (months)", fontsize=6.6)
axe.set_ylabel("IPCW time-dependent AUC", fontsize=6.6)
axe.set_title("time-dependent AUC", fontsize=6.3, pad=3)
panel("e", 0.010, 0.262)

# ============================================================ f  what capacity buys, per modality
axf = fig.add_axes([0.435, 0.088, 0.195, 0.130])
DIMS = [2, 4, 8, 16]
for fam, pre in (("slide", "titan"), ("clinical", "clinical")):
    real = [OF["blocks"]["%s_d%d" % (pre, d)]["real_minus_corrected_null"] for d in DIMS]
    axf.plot(DIMS, real, color=MOD[fam], marker=MK[fam], linestyle=LS[fam], ms=3.2, lw=1.3,
             markeredgecolor="none", label=fam)
axf.axhline(0.0, color=INK, lw=0.7, ls=":")
axf.set_xscale("log", base=2)
axf.set_xticks(DIMS)
axf.set_xticklabels([str(d) for d in DIMS], fontsize=6.0)
axf.set_xlabel("principal components fitted", fontsize=6.6)
axf.set_ylabel("above the permuted control", fontsize=6.6)
axf.set_title("capacity above its own control", fontsize=6.3, pad=3)
panel("f", 0.345, 0.262)

# ============================================================ g  how each modality separates
axg = fig.add_axes([0.755, 0.088, 0.195, 0.130])
ARMS_G = [("clinical", "clinical"), ("wsi_titan", "slide"), ("omics", "omics"), ("OURS", "ours")]
yg = np.arange(len(ARMS_G))[::-1]
chi = [SM["arms"][k]["logrank_chi2"] for k, _ in ARMS_G]
axg.barh(yg, chi, color=[MOD[f] for _, f in ARMS_G], height=0.62, edgecolor="none")
for yi, (k, _) in zip(yg, ARMS_G):
    a = SM["arms"][k]
    med = a["group_median_survival_months"]
    hi = med[2] if med[2] is not None else float("nan")
    axg.text(a["logrank_chi2"] + 1.4, yi, "high-risk median\n%s mo" %
             ("%.0f" % hi if hi == hi else "n/r"), va="center", fontsize=5.4, color=INK,
             linespacing=1.3)
axg.set_yticks(yg)
axg.set_yticklabels(["clinical", "slide", "omics", "all three"], fontsize=6.0)
axg.set_xlim(0, max(chi) * 1.78)
axg.set_xlabel("log-rank $\\chi^2$, 2 df, own tertiles", fontsize=6.6)
axg.set_title("tertile separation", fontsize=6.3,
              pad=3)
panel("g", 0.665, 0.262)


# --- the structural claims this figure makes, asserted so a silent data change cannot pass.
assert len(CARDS) == 3, "three modalities enter the model; found %d" % len(CARDS)
assert len(REPS) == 7, "seven representations were scored; found %d" % len(REPS)
assert abs(tie["clinical"] - 0.5) < 0.02, (
    "panel d exists because the clinical arm is a COIN FLIP on the tightest ties; it reads "
    "%.4f" % tie["clinical"])
assert tie["titan"] - tie["clinical"] > 0.05, (
    "panel d exists because the slide arm holds where the clinical arm does not; the gap is "
    "%+.4f" % (tie["titan"] - tie["clinical"]))
assert R["spearman_titan_vs_clinical"] < 0.5, (
    "panel c exists because image and clinical are nearly orthogonal; rho is %.3f"
    % R["spearman_titan_vs_clinical"])
assert CP["clinically_tied_q05"]["pairs_total"] == CP["all_pairs"]["pairs"], (
    "the tie quantiles and the full set must be drawn from one pair universe")

style.save(fig, style.out(HERE, "fig3_modalities.pdf"))
print("wrote fig3_modalities.pdf   %.3f x %.3f in" % (FIGW, FIGH))
print("  a: 3 modality cards | b: %d representations | c: rho(slide,clinical)=%.3f"
      % (len(REPS), R["spearman_titan_vs_clinical"]))
print("  d: clinical %.4f vs slide %.4f on the tightest %d pairs"
      % (tie["clinical"], tie["titan"], CP["clinically_tied_q05"]["pairs"]))
