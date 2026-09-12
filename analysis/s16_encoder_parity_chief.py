#!/usr/bin/env python3
"""Run ModRank on the encoder the comparators themselves consumed, and compare on that ground.

WHAT THIS ANSWERS. Cold panel round 2, finding G8: ModRank reads a TITAN slide embedding while
SurvPath and PIBD read their own tile features, so "method" and "2025 encoder" are confounded. The
existing seven-encoder sweep swaps the slide block through seven public encoders, and none of them
is the one the comparators actually used.

The comparator rerun in this study consumed CHIEF tile features. Pooling those same features to the
case and putting them where ModRank's slide arm goes makes the comparison happen on the comparator's
own ground: if the full method still clears the published entries with CHIEF, the lead is not an
artefact of having a better encoder than the methods it is compared against.

THE OTHER DIRECTION IS NOT MERELY EXPENSIVE, IT IS ILL-POSED, and that is worth writing down rather
than leaving as an unexplained omission. SurvPath and PIBD consume a bag of 4,096 patch tokens and
attend over them; the public TITAN release is one vector per slide. Feeding a bag of one destroys
the premise of the architecture. Doing it properly would mean obtaining the tile-level encoder,
re-extracting features for every slide from the whole-slide images, and retraining both models.

PHASE. Post-freeze and descriptive, like the sweep it extends. It changes no reported number and can
only make the paper's claim weaker or better supported.

WHY THE KNOWN-ANSWER CHECK IS NOT OPTIONAL. This file introduces a number into the paper, so a new
code path computing it is not enough: the path has to be shown to reproduce a number that already
exists. TITAN goes through the identical construction first and must return the sweep's own shipped
values. If it does not, the CHIEF number is not reported and this script exits 2.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs                              # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix                         # noqa: E402
from s6_amend_clinical import dimaf_clinical                              # noqa: E402
from s13_selection_null_and_encoder import oof                            # noqa: E402


def chief_percase(npz_path, keep, stem2case):
    """The pooled CHIEF matrix, aligned to the cohort's case order.

    The sweep's own loader expects `names` and `X`; this archive stores the mean-pooled block under
    `mean` alongside dispersion blocks that are a different arm and are not used here. Only the mean
    is taken, because the slide arm being swapped is a mean-pooled slide embedding and substituting
    a richer representation would change two things at once.
    """
    z = np.load(npz_path, allow_pickle=True)
    names = [str(x) for x in z["names"]]
    X = np.asarray(z["mean"], dtype=float)
    bag = {}
    for i, s in enumerate(names):
        c = stem2case.get(s)
        if c:
            bag.setdefault(c, []).append(i)
    have = [c for c in keep if c in bag]
    mu = X.mean(0)
    M = np.vstack([X[bag[c]].mean(0) if c in bag else mu for c in keep])
    return M, len(have)


def main():
    ap = argparse.ArgumentParser()
    for f in ("protocol", "root", "titan", "sp-dir", "omics-dir", "dimaf-dir", "chief", "sweep",
              "out"):
        ap.add_argument("--" + f, required=True)
    a = ap.parse_args()

    proto = json.load(open(a.protocol))
    assert proto["protocol_status"] == "frozen"

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    n = len(co.keep)
    ii, jj = cpairs(co.t, co.e)
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    # age, sex and stage, assembled exactly as the sweep assembles them. Grade is deliberately
    # absent, which is the paper's whole point about this benchmark's clinical baseline.
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    CLIN = np.hstack([dc["age"], dc["fem"], dc["stage"]])

    import csv as _csv
    stem2case = {r["slide_id"].replace(".svs", ""): r["case_id"]
                 for r in _csv.DictReader(open(a.root + "/tcga_blca_meta.csv"))}

    om = oof(P, co.t, co.e, fi, n, 0)
    cl = oof(CLIN, co.t, co.e, fi, n, 0)

    def arm(X):
        sl = oof(X, co.t, co.e, fi, n, 0)
        full = co.fold_pct(co.fold_pct(sl) + co.fold_pct(om) + co.fold_pct(cl))
        return (round(cidx(co.fold_pct(sl), ii, jj), 4), round(cidx(full, ii, jj), 4))

    # ---- known-answer check, before anything new is computed
    want = json.load(open(a.sweep))["B_encoder_swap"]["encoders"]["titan"]
    got_alone, got_full = arm(co.T)
    ok = (got_alone == want["slide_alone"] and got_full == want["full_method"])
    print("known-answer  titan slide %.4f (want %.4f)  full %.4f (want %.4f)  %s"
          % (got_alone, want["slide_alone"], got_full, want["full_method"],
             "REPRODUCED" if ok else "DOES NOT REPRODUCE"), file=sys.stderr)
    if not ok:
        print(json.dumps({"status": "error",
                          "error_code": "known_answer_check_failed_number_not_reported"}),
              file=sys.stderr)
        return 2

    # ---- the arm this file exists for
    M, cov = chief_percase(a.chief, co.keep, stem2case)
    if cov < 0.9 * n:
        print(json.dumps({"status": "error", "error_code": "chief_covers_too_few_cases"}),
              file=sys.stderr)
        return 2
    chief_alone, chief_full = arm(M)
    print("chief         slide %.4f  ModRank %.4f  (%d of %d cases covered)"
          % (chief_alone, chief_full, cov, n), file=sys.stderr)

    rep = {
        "artifact_type": "s16_encoder_parity_on_the_comparators_own_features",
        "phase_of_origin": "frozen_post_hoc",
        "answers": "cold panel round 2, finding G8",
        "n": n, "events": int(co.e.sum()), "comparable_pairs": int(len(ii)),
        "known_answer_check": {"encoder": "titan", "recomputed": [got_alone, got_full],
                               "shipped": [want["slide_alone"], want["full_method"]],
                               "reproduced": ok},
        "chief": {"slide_alone": chief_alone, "full_method": chief_full,
                  "cases_covered": cov,
                  "what_it_is": "the tile features the SurvPath rerun in this study consumed, "
                                "mean-pooled to the case"},
        "titan_for_reference": {"slide_alone": got_alone, "full_method": got_full},
        "reading": ("the comparison on the comparator's own encoder. A full-method value that "
                    "still clears the published entries means the lead is not an artefact of "
                    "reading a better slide encoder than the methods it is compared against."),
        "the_other_direction_was_not_run": (
            "comparators on TITAN. Not merely expensive: SurvPath and PIBD consume a bag of 4,096 "
            "patch tokens while the public TITAN release is one vector per slide, so the "
            "substitution is undefined without re-extracting tile-level features from the "
            "whole-slide images and retraining both models."),
    }
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(rep, open(a.out, "w"), indent=1)
    print("wrote %s" % a.out, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
