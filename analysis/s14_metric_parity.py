#!/usr/bin/env python3
"""Pooled out-of-fold concordance against the mean of per-fold concordances, from one set of scores.

The published entries on these folds are mostly means of per-fold C, ModRank's primary is pooled over
out-of-fold predictions. This reads the reporting dump's per-case seed-mean ModRank score and writes
both, so a reader can convert. It fits nothing.

The committed metric-parity.json predates this script, whose producer was not kept. Before this file
was used to regenerate it (amendment A2), it reproduced that committed file field for field from the
committed dump (experiments/20260928-amendment-a2/legacy-checks/s14_metric_parity.txt).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import cindex                                                  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dump", required=True, help="reporting-dump.json")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cases = json.load(open(a.dump))["cases"]
    r = np.array([c["seed_mean"]["OURS"] for c in cases], float)
    t = np.array([c["months"] for c in cases], float)
    e = np.array([c["event"] for c in cases], float)
    f = np.array([c["fold"] for c in cases], int)
    pooled = cindex(r, t, e)
    per = [cindex(r[f == k], t[f == k], e[f == k]) for k in sorted(set(f.tolist()))]
    out = {"artifact_type": "s14_metric_parity", "phase_of_origin": "frozen_post_hoc",
           "why": "cold panel round 2 asked whether every row of the published table is the same "
                  "estimand. Ours is a pooled out-of-fold C-index; the convention on these folds is "
                  "generally a mean of per-fold values. Both are reported here from one set of "
                  "predictions so a reader can convert.",
           "pooled_out_of_fold": round(float(pooled), 4),
           "mean_of_per_fold": round(float(np.mean(per)), 4),
           "per_fold": [round(float(x), 4) for x in per],
           "sd_across_folds": round(float(np.std(per, ddof=1)), 4),
           "difference_mean_minus_pooled": round(float(np.mean(per) - pooled), 4)}
    json.dump(out, open(a.out, "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ("why",)}), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
