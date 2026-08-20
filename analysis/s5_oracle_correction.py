#!/usr/bin/env python3
"""Correction. The Step 0 oracle's control answered the wrong question.

Pre-freeze, exploratory. NOT reportable.

WHAT WAS RUN THIS MORNING, and what it actually measured. The headroom probe fitted a Cox model
IN-SAMPLE on the held-out fold at d principal components, scored it on that same fold, and called
the result an oracle. Its control refitted the same design on PERMUTED outcomes and scored THAT
against the REAL outcomes, landing at 0.4980 over 1,000 draws. From the pair it reported
`recoverable_above_control = 0.2686` and `estimation_loss = 0.1170` for the TITAN arm, and this
campaign then spent two candidates trying to recover that 0.1170.

THE CONTROL WAS THE WRONG ONE. Fitting to target A and scoring on target B measures whether A and
B are related -- of course a permutation destroys that, and of course the answer is 0.5. It says
NOTHING about how much a fit is inflated by being evaluated on the very outcomes it was fitted to.
Sixteen parameters against roughly 23 events will fit a great deal of any target, real or not.

THE RIGHT CONTROL, and it is a two-line change: permute the outcome, fit to the permuted outcome,
and score on THAT SAME PERMUTED OUTCOME. Zero signal by construction, identical n, identical event
count, identical censoring, identical capacity. Whatever it reaches is the pure in-sample inflation
of this fit at this capacity.

    real_in_sample  ~  null_in_sample   ->  the "oracle" was overfitting and the representation
                                            was never shown to hold more than its out-of-fold
                                            value. `estimation_loss` was an artifact, and so was
                                            the deficit two candidates were built to attack.
    real_in_sample  >>  null_in_sample  ->  the original reading survives and the headroom is real.

WHAT MADE THIS WORTH CHECKING was not suspicion of the arithmetic. It was that C2 predicted, before
running, that variance reduction would help the arms in proportion to their dimension, and the
measured result went the other way: bagging cost the 768-dimensional TITAN arm 0.056 and gained the
275-dimensional omics arm 0.024. A falsified mechanism is a reason to re-examine the measurement
the mechanism was read off.

Both controls are reported side by side so the difference between them is visible rather than
described.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, cox_fit  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--titan", required=True)
    ap.add_argument("--sp-dir", required=True)
    ap.add_argument("--perms", type=int, default=200)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    rng = np.random.default_rng(0)
    rep = {"artifact_type": "s5_oracle_correction", "reportable": False,
           "phase_of_origin": "exploratory", "perms": a.perms,
           "corrects": "s5_step0_headroom.py oracle block, whose control fitted on permuted "
                       "outcomes and scored on REAL ones, which cannot measure in-sample inflation"}

    out = {}
    for name, Xall in (("titan", co.T), ("clinical", co.CLIN)):
        bases = []
        for tri, vai in fi:
            if len(tri) < 30 or len(vai) < 20:
                bases.append(None)
                continue
            mu, sd = Xall[tri].mean(0), Xall[tri].std(0) + 1e-9
            Z = (Xall - mu) / sd
            _, _, Vt = np.linalg.svd(Z[tri], full_matrices=False)
            bases.append((Z, Vt))
        for d in (2, 4, 8, 16):
            real, null_same, null_cross, oof = [], [], [], []
            for (tri, vai), bs in zip(fi, bases):
                if bs is None:
                    continue
                Z, Vt = bs
                dd = min(d, Z.shape[1], len(tri) - 1)
                P = Z @ Vt[:dd].T
                tv, ev = co.t[vai], co.e[vai]
                ii, jj = cpairs(tv, ev)
                if ii.size == 0:
                    continue
                real.append(cidx(P[vai] @ cox_fit(P[vai], tv, ev, 1.0), ii, jj))
                oof.append(cidx(P[vai] @ cox_fit(P[tri], co.t[tri], co.e[tri], 1.0), ii, jj))
                for _ in range(a.perms):
                    q = rng.permutation(len(vai))
                    tp, ep = tv[q], ev[q]
                    bp = cox_fit(P[vai], tp, ep, 1.0)
                    # THE CORRECTED CONTROL: scored on the same permuted outcome it was fitted to
                    ip, jp = cpairs(tp, ep)
                    if ip.size:
                        null_same.append(cidx(P[vai] @ bp, ip, jp))
                    # the ORIGINAL control, kept so the two are visible together
                    null_cross.append(cidx(P[vai] @ bp, ii, jj))
            if not real:
                continue
            r, ns, nc, o = (float(np.mean(real)), float(np.mean(null_same)),
                            float(np.mean(null_cross)), float(np.mean(oof)))
            out["%s_d%d" % (name, d)] = {
                "real_in_sample": round(r, 4),
                "CORRECTED_null_fit_and_score_on_same_permutation": round(ns, 4),
                "corrected_null_sd": round(float(np.std(null_same, ddof=1)), 4),
                "original_control_fit_permuted_score_real": round(nc, 4),
                "out_of_fold": round(o, 4),
                "real_minus_corrected_null": round(r - ns, 4),
                "estimation_loss_as_originally_reported": round(r - o, 4),
            }
            print("%-12s d=%2d  real=%.4f  CORRECTED null=%.4f  old ctl=%.4f  oof=%.4f  "
                  "real-null=%+.4f" % (name, d, r, ns, nc, o, r - ns),
                  file=sys.stderr, flush=True)
    rep["blocks"] = out
    rep["reading"] = (
        "real_minus_corrected_null is the only quantity here that is evidence about the "
        "representation. If it is near zero the in-sample number was inflation, the reported "
        "estimation_loss was an artifact of comparing an inflated number to an honest one, and "
        "no headroom was ever demonstrated at that capacity.")

    text = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
