#!/usr/bin/env python3
"""Figure 6 (NC build): the fusion choice on identical inputs, and whether it survives the partition.

  a  the released folds, seeds 0-4: ModRank against a single ridge Cox on all three blocks
     concatenated, learned stacking weights, and the clinical arm
  b  the noise floor's 24 random re-partitions: each construction's concordance per partition, joined
     by partition so that the pairing a paired claim rests on is visible
  c  whole tissue source sites held out (33 sites, five folds): per-fold and pooled concordance

NO NUMERAL IS TYPED. Values come from
experiments/20260911-blca-posthoc/results/{unified-fusion-and-added-value,resplit-and-site-cv}.json,
and the plotted values are written to fig6_robustness.source.json. Authored at 180 mm. Dots and
lines only; every construction carries its own marker as well as its colour.
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
FU = json.load(open(os.path.join(PH, "unified-fusion-and-added-value.json")))["fusion_on_identical_inputs"]
RS = json.load(open(os.path.join(PH, "resplit-and-site-cv.json")))
if "per_partition" not in RS["resplits"]:
    raise SystemExit("figure 6 needs the per-partition values; rerun analysis/s22_resplit_and_site_cv.py")

style.apply()
C, INK = style.ROLE, style.INK
W, H = 7.087, 3.05
fig = plt.figure(figsize=(W, H))
CONS = [("ours", "ModRank", C["ours"], "o"), ("stacked", "stacked (learned weights)", C["omics"], "s"),
        ("concat", "concatenated Cox", C["competitor"], "D"), ("clinical", "clinical (stage)", C["clinical"], "^")]
src = {}


def panel(letter, x):
    fig.text(x, 0.975, letter, fontsize=10, fontweight="bold", va="top", ha="left")


# ---------------------------------------------------------------- a  released folds, seeds 0-4
axa = fig.add_axes([0.065, 0.27, 0.22, 0.60])
seeds = {"ours": FU["ours_equal_weight_rank_average"]["per_seed"],
         "stacked": FU["stacked_learned_weights"]["per_seed"],
         "concat": FU["concatenated_ridge_cox"]["per_seed"]}
src["a"] = seeds
for x, (key, lab, col, mk) in enumerate(CONS[:3]):
    v = seeds[key]
    axa.scatter(np.full(len(v), x) + np.linspace(-0.12, 0.12, len(v)), v, s=16, marker=mk, color=col,
                zorder=3)
    axa.plot([x - 0.22, x + 0.22], [np.mean(v)] * 2, color=INK, lw=1.0)
axa.set_xticks(range(3))
axa.set_xticklabels(["ModRank", "stacked", "concat."], fontsize=7)
axa.set_ylabel("concordance")
axa.set_title("released folds, 5 seeds", fontsize=8, fontweight="bold")
panel("a", 0.0)

# ---------------------------------------------------------------- b  24 re-partitions
axb = fig.add_axes([0.37, 0.27, 0.30, 0.60])
per = RS["resplits"]["per_partition"]
keys = ["ours", "stacked", "concat", "clinical"]
src["b"] = {k: per[k] for k in keys}
JIT = {k: np.random.default_rng(x).uniform(-0.09, 0.09, len(per[k])) for x, k in enumerate(keys)}
for j in range(len(per["ours"])):                  # each line joins one partition's own points
    axb.plot([x + JIT[k][j] for x, k in enumerate(keys)], [per[k][j] for k in keys],
             color="#D5DAE2", lw=0.6, zorder=1)
for x, (key, lab, col, mk) in enumerate(CONS):
    axb.scatter(x + JIT[key], per[key], s=11, marker=mk, color=col, zorder=3, linewidths=0)
axb.set_xticks(range(4))
axb.set_xticklabels(["ModRank", "stacked", "concat.", "clinical"], fontsize=7)
n_part = len(per["ours"])
wins = {k: sum(1 for a_, b_ in zip(per["ours"], per[k]) if a_ > b_) for k in ("stacked", "concat", "clinical")}
axb.set_title("%d re-partitions of the same patients" % n_part, fontsize=8, fontweight="bold")
src["b_wins"] = wins
panel("b", 0.305)

# ---------------------------------------------------------------- c  held-out hospitals
axc = fig.add_axes([0.745, 0.27, 0.24, 0.60])
SC = RS["site_grouped_cv"]
src["c"] = SC
for key, lab, col, mk in CONS:
    v = [f[key] for f in SC["folds"]]
    axc.plot(range(1, len(v) + 1), v, "-", color=col, marker=mk, ms=3.8, lw=0.9, label=lab)
    axc.scatter([len(v) + 1.2], [SC["pooled"][key]], s=30, marker=mk, color=col, zorder=3)
axc.set_xticks(list(range(1, len(SC["folds"]) + 1)) + [len(SC["folds"]) + 1.2])
axc.set_xticklabels([str(i) for i in range(1, len(SC["folds"]) + 1)] + ["pooled"], fontsize=7)
axc.set_xlabel("site-grouped fold")
axc.set_title("%d hospitals held out whole" % SC["sites_total"], fontsize=8, fontweight="bold")
# one legend for all three panels, under them, so no panel carries a key over its data
hs = [plt.Line2D([], [], color=col, marker=mk, lw=0.9, ms=4.5, label=lab) for _, lab, col, mk in CONS]
fig.legend(handles=hs, loc="lower center", ncol=4, fontsize=7, bbox_to_anchor=(0.5, 0.0),
           handletextpad=0.4, columnspacing=1.6)
panel("c", 0.69)

out = os.path.join(HERE, "fig6_robustness.pdf")
style.save(fig, out)
json.dump(src, open(os.path.join(HERE, "fig6_robustness.source.json"), "w"), indent=1)
print("wrote %s  %.3f x %.3f in | wins vs stacked/concat/clinical: %s" % (out, W, H, wins))
