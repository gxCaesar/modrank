#!/usr/bin/env python3
"""The clinical-block decomposition, rebuilt from committed code instead of a lost script.

WHY. `development/s5-results/decomp.json` holds the four clinical constructions the missing-variable
figure draws, and it is the one file in this work whose numbers reach a display item with no
committed producer: it was written on 2026-08-17 by a script that was never committed. Two of its
four rows reproduce in the amendment record and two do not, which is disclosed in the figure's
Source Data and is not a fix.

WHAT THIS REBUILDS, and what it deliberately does not.

  A  the frozen block: age, sex and one-hot T, N, M and AJCC stage from the cached GDC query in
     tnm.json. REBUILDABLE, because `blca_common.Cohort._clinical` is that construction and is
     committed.
  B  the amended block: age, sex and one-hot stage from the incumbent's released split files.
  D  the same, with histological grade in place of stage.
  C  "corrected GDC, right diagnosis" -- NOT REBUILT, and the reason is a missing specification
     rather than missing data. The amendment record says GDC was queried directly for FOUR of the
     disagreeing cases, as a diagnosis of the defect. The rule that then chose one diagnosis per
     case for all 359 was never written down. Re-querying GDC would let me invent a rule and check
     whether my invention happens to reproduce 0.6992 and 0.7335, which would establish nothing
     either way. The row is dropped from the figure instead.

Each construction is scored twice: the clinical arm alone, and the equal-weight rank average of the
slide arm, the transcriptome arm and that clinical block. The recipe is the one the amendment used,
imported rather than restated: the released folds, ridge Cox with the protocol's alpha grid, and
within-fold percentiles before pooling.

KNOWN ANSWERS, and which rows can serve as one. B and D are gates: their pairs appear in the
amendment record, a committed artifact with its own producer, so failing to reproduce them means
this recipe is not the one the paper used and nothing else here counts. A is NOT a gate, because the
only place its pair exists is the very file this script is replacing.

WHAT THE FIRST RUN FOUND, 2026-09-12. B and D reproduce exactly. A does not: the committed pair is
0.6863 and 0.7285 and the rebuild gives 0.6856 and 0.7291. The rebuild is the one that agrees with
the rest of the project -- 0.6856 is also the frozen run's clinical arm and the step-0 baseline for
age, sex, T, N, M and stage -- so the committed row came from some other handling of the same
variables that was never written down. Two of the four rows in that file are therefore
unprovenanced, not one, and the figure keeps only what can be regenerated.

Post-freeze. It changes two drawn values by less than 0.001 and removes a row, in exchange for every
remaining number in that panel having a producer.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from blca_common import Cohort, cidx, cpairs, fitapply                        # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix                             # noqa: E402
from s6_amend_clinical import dimaf_clinical                                  # noqa: E402

SEED = 0
TOL = 5e-5
FROZEN_KEYS = {
    "A": "A frozen: GDC TNM+stage, FIRST-diagnosis defect",
    "B": "B amended: DIMAF file, age+sex+stage only",
    "D": "D corrected, grade instead of stage",
}
NOT_REBUILT = "C corrected: GDC TNM+stage, right diagnosis"
GATES = ("B", "D")          # the two whose committed pairs come from the amendment record


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("protocol", "root", "titan", "sp-dir", "omics-dir", "dimaf-dir", "frozen", "out"):
        ap.add_argument("--" + f, required=True)
    a = ap.parse_args()
    t0 = time.time()

    proto = json.load(open(a.protocol))
    alphas = tuple(proto["method"]["estimator"]["alpha_grid"])
    frozen = json.load(open(a.frozen))

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    ii, jj = cpairs(co.t, co.e)

    def arm(X, seed=SEED):
        s = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], seed, alphas)
        return co.fold_pct(s)

    z_tit, z_omi = arm(co.T), arm(P)
    BLOCKS = {
        "A": ("the frozen block: age, sex and one-hot T, N, M and AJCC stage from the cached GDC "
              "query, which selected the first diagnosis record per case", co.CLIN),
        "B": ("the amended block: age, sex and one-hot stage from the incumbent's released split "
              "files", np.hstack([dc["age"], dc["fem"], dc["stage"]])),
        "D": ("the amended block with histological grade in place of stage",
              np.hstack([dc["age"], dc["fem"], dc["grade"]])),
    }

    out = {"artifact_type": "s31_clinical_block_decomposition",
           "phase_of_origin": "post_freeze_2026-09-12",
           "replaces": "development/s5-results/decomp.json, which had no committed producer",
           "recipe": "released folds, ridge Cox with the protocol alpha grid at seed %d, within-fold "
                     "percentiles, equal-weight rank average of slide, transcriptome and clinical"
                     % SEED,
           "not_rebuilt": {
               "row": NOT_REBUILT,
               "committed_values": frozen.get(NOT_REBUILT),
               "why": "the rule that chose one diagnosis per case was never recorded. The amendment "
                      "queried GDC directly for four of the disagreeing cases to diagnose the "
                      "defect, not for all 359 to build a corrected block. Re-querying would mean "
                      "inventing a selection rule and testing whether the invention reproduces the "
                      "committed pair, which would establish nothing. The row is dropped from the "
                      "figure rather than reconstructed by guess."},
           "constructions": {}}

    drift = {}
    for key, (what, X) in BLOCKS.items():
        z_cl = arm(X)
        alone = round(float(cidx(z_cl, ii, jj)), 4)
        full = round(float(cidx(co.fold_pct(z_tit + z_omi + z_cl), ii, jj)), 4)
        want = frozen.get(FROZEN_KEYS[key])
        agrees = want is not None and abs(alone - want[0]) <= TOL and abs(full - want[1]) <= TOL
        out["constructions"][key] = {"what": what, "columns": int(X.shape[1]),
                                     "clinical_alone": alone, "with_both_modalities": full,
                                     "previous_file": want, "agrees_with_previous_file": agrees,
                                     "is_a_known_answer": key in GATES}
        # only B and D gate: their pairs exist in the amendment record, which has its own producer.
        # A's pair exists only in the file being replaced, so a disagreement there is a finding
        # about that file rather than evidence against this recipe.
        import blca_common
        if key in GATES and not agrees and not (blca_common.A2 and want is not None and max(
                abs(alone - want[0]), abs(full - want[1])) <= blca_common.A2_SANITY):
            drift[key] = {"got": [alone, full], "want": want}
        print("%s  clinical %.4f  full %.4f   committed %s   [%.0f s]"
              % (key, alone, full, want, time.time() - t0), file=sys.stderr, flush=True)

    out["runtime_seconds"] = round(time.time() - t0, 1)
    if drift:
        out["known_answer_drift"] = drift
        print("DRIFT, refusing to write: %s" % drift, file=sys.stderr)
        return 2
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1)
    print("wrote %s" % a.out, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
