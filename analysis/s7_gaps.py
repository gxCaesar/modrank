#!/usr/bin/env python3
"""Six gaps a reviewer would find, closed in one pass.

d  MULTIPLICITY. The frozen protocol declares `multiplicity_policy: Holm within each prespecified
   comparison family` and it was never applied. Applying it after the fact is legitimate only
   because the FAMILY was named in advance; the p-values are what they are.

m  DETECTABLE MARGIN. The paper repeats that this cohort cannot resolve margins of 0.03-0.05. That
   is rhetoric until it is a number. What margin would a two-sided test at 80% power detect, given
   the case-level bootstrap SE measured here? And inverted: what power does the observed margin
   have? This makes the central limitation quantitative and it applies to every published value on
   this benchmark, not only ours.

h  MISSING STAGE. 32 of 413 rows in the benchmark's clinical file carry no stage. One-hot encoding
   gives those cases an all-zero stage block, which is an implicit "no stage information" and is a
   modelling choice rather than a neutral default. Sensitivity: drop them and refit.

k  THE CLINICIAN'S MODEL. Our clinical arm is age + sex + stage. A pathology report also carries
   grade, and a clinician building a model would use both. If stage + grade beats stage alone the
   paper should say so; if it does not, that is the cleaner statement of the finding.

g  SLIDE POOLING. A case with several slides is mean-pooled. Max-pooling is the obvious
   alternative and encodes "the worst-looking slide decides", which is closer to how a pathologist
   reads. Reported whichever way it lands.

l  SEED-MATCHED COMPETITOR METRICS. The additional survival metrics were computed for SurvPath at
   seed 1 while the headline is now seed-matched. Recomputed on the seed-averaged competitor so no
   table mixes the two.
"""

from __future__ import annotations

import argparse
import collections
import math
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import Cohort, cidx, cpairs, fitapply  # noqa: E402
from s5_omics_arm import load_rna, pathway_matrix  # noqa: E402
from s6_amend_clinical import dimaf_clinical, boot  # noqa: E402
from s7_survpath_multiseed import load_seed  # noqa: E402
from s7_survival_metrics import km_censoring, g_at, td_auc, ibs, logrank, chi2_sf  # noqa: E402


def holm(pvals: dict):
    """Holm-Bonferroni step-down. Returns adjusted p-values, order preserved."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m = len(items)
    out, running = {}, 0.0
    for i, (k, p) in enumerate(items):
        adj = min(1.0, (m - i) * p)
        running = max(running, adj)          # enforce monotonicity
        out[k] = round(running, 4)
    return {k: out[k] for k in pvals}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("root", "titan", "sp-dir", "seed-dir", "omics-dir", "dimaf-dir", "clin-dir", "out"):
        ap.add_argument("--" + f, required=True)
    a = ap.parse_args()

    co = Cohort(a.root, a.titan, a.sp_dir)
    fi = co.fold_indices()
    ii, jj = cpairs(co.t, co.e)
    G, genes, _ = load_rna(os.path.join(a.omics_dir, "rna_combine.csv"), co.keep)
    P, _ = pathway_matrix(G, genes, os.path.join(a.omics_dir, "combine_signatures.csv"))
    dc = dimaf_clinical(a.dimaf_dir, co.keep)
    AGE, FEM, ST, GR = dc["age"], dc["fem"], dc["stage"], dc["grade"]

    def arm(X, seed=0, sub=None):
        idx = np.arange(len(co.keep)) if sub is None else np.flatnonzero(sub)
        s = np.full(len(co.keep), np.nan)
        for tri, vai in fi:
            tri = np.intersect1d(tri, idx); vai = np.intersect1d(vai, idx)
            if len(tri) < 30 or not len(vai):
                continue
            s[vai] = fitapply(X[tri], co.t[tri], co.e[tri], X[vai], seed)
        return co.fold_pct(s)

    CLIN = np.hstack([AGE, FEM, ST])
    zt, zo, zc = arm(co.T), arm(P), arm(CLIN)
    OURS = co.fold_pct(zt + zo + zc)
    rep = {"artifact_type": "s7_reviewer_gaps", "phase_of_origin": "frozen_post_hoc"}

    # ---------------------------------------------------------------- l: seed-matched competitor
    zs = {}
    r, f, _, _ = load_seed(a.sp_dir, co.keep, co.idx); zs[1] = co.fold_pct(r, f)
    for s in (2, 3, 4, 5):
        d = os.path.join(a.seed_dir, "seed%d" % s)
        if os.path.isdir(d):
            r, f, _, _ = load_seed(d, co.keep, co.idx); zs[s] = co.fold_pct(r, f)
    SP_avg = co.fold_pct(np.vstack([zs[s] for s in sorted(zs)]).mean(0))
    SPC_avg = co.fold_pct(SP_avg + zc)

    uniq, g = km_censoring(co.t, co.e)
    ev_t = np.sort(co.t[co.e == 1])
    hs = [float(np.quantile(ev_t, q)) for q in (0.25, 0.5, 0.75)]
    grid = [float(x) for x in np.linspace(ev_t[0], np.quantile(ev_t, 0.9), 12)]
    kt, ks = km_censoring(co.t, 1.0 - co.e)
    base = g_at(kt, ks, np.array(grid))

    def metrics(z):
        curves = np.clip(base[None, :] ** np.exp(2.0 * (z[:, None] - 0.5)), 1e-6, 1 - 1e-6)
        q = np.quantile(z, [1 / 3, 2 / 3]); grp = np.digitize(z, q)
        c2, df, _ = logrank(co.t, co.e, grp)
        return {"cindex": round(cidx(z, ii, jj), 4),
                **{"td_auc_%dm" % round(h): round(td_auc(z, co.t, co.e, h, uniq, g), 4) for h in hs},
                "ibs": round(ibs(curves, grid, co.t, co.e, uniq, g), 4),
                "logrank_p": float("%.3g" % chi2_sf(c2, df))}
    rep["l_seed_matched_competitor_metrics"] = {
        "survpath_alone_seed_averaged": metrics(SP_avg),
        "survpath_plus_clinical_seed_averaged": metrics(SPC_avg),
        "OURS": metrics(OURS),
        "note": "recomputed on the seed-AVERAGED competitor so no table mixes a five-seed arm with "
                "a one-seed one"}
    print("l: seed-matched competitor metrics done", file=sys.stderr, flush=True)

    # ---------------------------------------------------------------- k: the clinician's model
    variants = {
        "age+sex+stage (the frozen arm)": CLIN,
        "age+sex+grade": np.hstack([AGE, FEM, GR]),
        "age+sex+stage+grade (a pathology report)": np.hstack([AGE, FEM, ST, GR]),
        "stage alone": ST,
    }
    kk = {}
    for nm, X in variants.items():
        if X.shape[1] == 0:
            continue
        z = arm(X)
        kk[nm] = {"clinical_alone": round(cidx(z, ii, jj), 4),
                  "full_arm_with_it": round(cidx(co.fold_pct(zt + zo + z), ii, jj), 4)}
    rep["k_clinicians_model"] = {
        "variants": kk,
        "reading": "if adding grade to stage does not help, the finding is cleaner: it is stage "
                   "that was missing, not 'more clinical variables'."}
    print("k: clinician variants done", file=sys.stderr, flush=True)

    # ---------------------------------------------------------------- h: missing stage
    has_stage = ST.sum(1) > 0
    rep["h_missing_stage"] = {
        "cases_without_a_stage_value": int((~has_stage).sum()),
        "of_total": len(co.keep),
        "encoding_in_the_frozen_arm": "all-zero stage block, i.e. an implicit 'no stage information'",
        "full_cohort": round(cidx(OURS, ii, jj), 4)}
    if (~has_stage).sum() > 0:
        sub = has_stage
        i2, j2 = cpairs(co.t[sub], co.e[sub])
        z2 = co.fold_pct(arm(co.T, 0, sub) + arm(P, 0, sub) + arm(CLIN, 0, sub))
        rep["h_missing_stage"]["complete_cases_only"] = {
            "n": int(sub.sum()), "events": int(co.e[sub].sum()),
            "OURS": round(cidx(z2[sub], i2, j2), 4),
            "clinical": round(cidx(arm(CLIN, 0, sub)[sub], i2, j2), 4)}
    print("h: missing-stage sensitivity done", file=sys.stderr, flush=True)

    # ---------------------------------------------------------------- g: slide pooling
    import pickle
    d0 = pickle.load(open(a.titan, "rb"))
    E = np.asarray(d0["embeddings"])
    stem = {}
    for i, fn in enumerate(d0["filenames"]):
        s = str(fn)
        for suf in (".svs", ".h5", ".pt"):
            if s.endswith(suf):
                s = s[:-len(suf)]
        stem[s] = i
    s2c = {r["slide_id"].replace(".svs", ""): r["case_id"]
           for r in csv.DictReader(open(a.root + "/tcga_blca_meta.csv"))}
    per = collections.defaultdict(list)
    for s, i in stem.items():
        c = s2c.get(s)
        if c:
            per[c].append(i)
    n_multi = sum(1 for c in co.keep if len(per.get(c, [])) > 1)
    Tmax = np.vstack([E[per[c]].max(0) if c in per else np.zeros(E.shape[1]) for c in co.keep])
    zt_max = arm(Tmax.astype(float))
    rep["g_slide_pooling"] = {
        "cases_with_more_than_one_slide": n_multi, "of_total": len(co.keep),
        "mean_pool_image_arm": round(cidx(zt, ii, jj), 4),
        "max_pool_image_arm": round(cidx(zt_max, ii, jj), 4),
        "mean_pool_full_arm": round(cidx(OURS, ii, jj), 4),
        "max_pool_full_arm": round(cidx(co.fold_pct(zt_max + zo + zc), ii, jj), 4)}
    print("g: pooling sensitivity done", file=sys.stderr, flush=True)

    # ---------------------------------------------------------------- d + m: multiplicity, power
    fam = {
        "ours vs clinical alone": boot(OURS, zc, co.t, co.e, 6000),
        "ours vs image alone": boot(OURS, zt, co.t, co.e, 6000),
        "ours vs omics alone": boot(OURS, zo, co.t, co.e, 6000),
        "ours vs SurvPath+clinical (seed-matched)": boot(OURS, SPC_avg, co.t, co.e, 6000),
        "ours vs ours-without-clinical": boot(OURS, co.fold_pct(zt + zo), co.t, co.e, 6000),
    }
    raw = {k: v["p_two_sided"] for k, v in fam.items()}
    adj = holm(raw)
    rep["d_multiplicity"] = {
        "family": "the five prespecified comparisons in the frozen protocol",
        "method": "Holm-Bonferroni, step-down, as the protocol declares",
        "comparisons": {k: {"delta": fam[k]["mean"], "ci95": fam[k]["ci95"],
                            "p_raw": raw[k], "p_holm": adj[k],
                            "survives_holm_at_0.05": bool(adj[k] < 0.05)} for k in fam}}
    print("d: Holm applied", file=sys.stderr, flush=True)

    se = float(np.mean([v["se"] for v in fam.values()]))
    z80, z975 = 0.8416, 1.9600
    rep["m_detectable_margin"] = {
        "case_level_bootstrap_se": round(se, 4),
        "detectable_at_80pc_power_two_sided_0.05": round((z80 + z975) * se, 4),
        "power_of_the_observed_input_parity_margin": round(
            float(1 - 0.5 * (1 + math.erf((z975 - abs(fam["ours vs SurvPath+clinical (seed-matched)"]["mean"]) / se) / np.sqrt(2)))), 3),
        "n": len(co.keep), "events": int(co.e.sum()),
        "reading": "a two-sided test at this cohort's SE detects a margin of about %.3f at 80%% "
                   "power. The nine published entries on this benchmark span 0.596 to 0.691 and "
                   "consecutive entries differ by 0.002 to 0.012, so almost none of the published "
                   "ORDERING is resolvable at this n -- a statement about the benchmark, not about "
                   "any one method." % ((z80 + z975) * se)}
    print("m: power done", file=sys.stderr, flush=True)

    open(a.out, "w").write(json.dumps(rep, indent=1) + "\n")
    print("wrote", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
