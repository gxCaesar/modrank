#!/usr/bin/env python3
"""The modality figure's redundancy and tied-pair panels, recomputed with the AMENDED clinical block.

WHY. Panels c and d of the modality figure, and the Results sentences that quote them (Spearman
0.213 and 0.166, the clinical arm at 0.5012 and the slide arm at 0.6120 on 1,281 tightly tied
pairs), came from the development atlas (analysis/s5_step1_atlas.py), which used the PRE-amendment
23-column clinical block with the first-diagnosis defect. The manuscript's clinical arm is the
amended one. Found 2026-09-11.

HOW, so that only one thing changes. The atlas's own functions are imported rather than rewritten,
and the atlas is first re-run here with its OWN clinical block as a known-answer check: if 0.213,
0.166, 0.5012 and 0.6120 do not come back, nothing is reported. Then the clinical block, and only
the clinical block, is swapped for the amended one (age, sex, AJCC pathologic stage from DIMAF's
split files), seed 0, identical folds, estimator and SurvPath risk.

Post-freeze and descriptive. No reported primary changes.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blca_common                                           # noqa: E402
from blca_common import Cohort, fitapply                     # noqa: E402
from s5_step1_atlas import conditional_probe, spearman       # noqa: E402
from s6_amend_clinical import dimaf_clinical                 # noqa: E402

KNOWN = {"spearman_titan_vs_clinical": 0.213, "spearman_survpath_vs_clinical": 0.166,
         "tied_q05_clinical": 0.5012, "tied_q05_titan": 0.6120}


def probe(co, CL):
    fi = co.fold_indices()
    p_ti, p_cl = np.full(len(co.keep), np.nan), np.full(len(co.keep), np.nan)
    for tri, vai in fi:
        p_ti[vai] = fitapply(co.T[tri], co.t[tri], co.e[tri], co.T[vai], 0)
        p_cl[vai] = fitapply(CL[tri], co.t[tri], co.e[tri], CL[vai], 0)
    z_ti, z_cl, z_sp = co.fold_pct(p_ti), co.fold_pct(p_cl), co.fold_pct(co.spr)
    arms = {"clinical": z_cl, "titan": z_ti, "survpath": z_sp,
            "titan+clinical": co.fold_pct(z_ti + z_cl)}
    return {"redundancy": {"spearman_titan_vs_clinical": round(spearman(z_ti, z_cl), 4),
                           "spearman_survpath_vs_clinical": round(spearman(z_sp, z_cl), 4),
                           "spearman_titan_vs_survpath": round(spearman(z_ti, z_sp), 4)},
            "conditional_probe": conditional_probe(co, arms, z_cl)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("root", "titan", "sp-dir", "dimaf-dir", "out"):
        ap.add_argument("--" + f, required=True)
    a = ap.parse_args()
    co = Cohort(a.root, a.titan, a.sp_dir)

    old = probe(co, co.CLIN)
    got = {"spearman_titan_vs_clinical": old["redundancy"]["spearman_titan_vs_clinical"],
           "spearman_survpath_vs_clinical": old["redundancy"]["spearman_survpath_vs_clinical"],
           "tied_q05_clinical": old["conditional_probe"]["clinically_tied_q05"]["arms"]["clinical"],
           "tied_q05_titan": old["conditional_probe"]["clinically_tied_q05"]["arms"]["titan"]}
    # amendment A2: KNOWN holds pre-A2 values, checked to 6e-4 under BLCA_A2=0 and only as a
    # gross-failure bound under A2 (the pre-A2 run itself is reproduced separately)
    tol = blca_common.A2_SANITY if blca_common.A2 else 6e-4
    ok = all(abs(round(got[k], 4) - KNOWN[k]) < tol for k in KNOWN)
    for k in KNOWN:
        print("  known answer %-32s %.4f (want %.4f)" % (k, got[k], KNOWN[k]), file=sys.stderr)
    if not ok:
        print(json.dumps({"status": "error", "error_code": "atlas_known_answer_not_reproduced"}),
              file=sys.stderr)
        return 2

    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    new = probe(co, np.hstack([dc["age"], dc["fem"], dc["stage"]]))
    rep = {"artifact_type": "s19b_modality_atlas_amended", "phase_of_origin": "post_freeze_2026-09-11",
           "reportable": True,
           "known_answer": ({"reproduced": True, "pre_amendment_block": got} if not blca_common.A2 else
                            {"reproduced": False, "pre_amendment_block": got, "amendment": "A2",
                             "checked_within": blca_common.A2_SANITY}),
           "clinical_block": "amended: age, sex, AJCC pathologic tumour stage from DIMAF's split "
                             "files; seed 0; atlas functions imported unchanged",
           **new}
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(rep, open(a.out, "w"), indent=1)
    q = new["conditional_probe"]["clinically_tied_q05"]
    print("amended: rho titan-clin %.4f, survpath-clin %.4f; q05 %d pairs clinical %.4f titan %.4f"
          % (new["redundancy"]["spearman_titan_vs_clinical"],
             new["redundancy"]["spearman_survpath_vs_clinical"], q["pairs"], q["arms"]["clinical"],
             q["arms"]["titan"]), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
