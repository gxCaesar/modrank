#!/usr/bin/env python3
"""The added-value inflation D against a model with no signal, on TCGA-BLCA and in the GEO cohorts.

Post-freeze, specified 2026-09-27 after analyses A, B and D showed that D is positive for a model with
no signal (experiments/20260927-field-inflation/README.md, finding 1). Nothing here changes a
committed value.

Under the paper's rank rule a model M is combined with a reference R as the within-fold percentile of
M + R. Equal weighting dilutes a stronger reference more, so a score carrying no information gains
less over the stage-based reference than over the grade-based one and has D > 0. For each
construction the no-signal benchmark keeps the reference and the combination rule and replaces the
model by random within-fold ranks: one random arm for a two-arm construction (SurvPath, PIBD and the
further architectures, and ModRank's transcriptome-plus-clinical form in GEO), two independent random
arms for ModRank's three-arm sum on TCGA-BLCA.

Reported per construction:
  D, the committed value, recomputed first as the known answer;
  the benchmark's D on the full data, mean and 2.5/97.5 percentiles over 2,000 random draws;
  the excess, D minus the benchmark's expected D, with a case-bootstrap interval: on each of 6,000
  resamples (seed 20260911, the draws of every committed D) the construction's D minus the mean D of
  20 fixed random draws evaluated on the same resample.

  s35_inflation_null.py blca --percase ... --committed ... --runs ... --models ... --out ...   (Mac)
  s35_inflation_null.py geo  --protocol ... --geo-dir ... --signatures ... --percase ... --out ... (sysu)
  s35_inflation_null.py five --titan ... --meta-dir ... --clin-dir ... --pan-dir ... --sig ...
                             --committed-b ... --out ...                                         (Mac)
The five-study form (protocol B's bladder, head and neck and stomach) uses ModRank's three-arm sum,
so its benchmark has two random arms.
Each exits 2 without writing if a known answer is not reproduced.
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
from blca_common import cidx, cpairs, pct                                      # noqa: E402

BOOT_SEED, REPS, NULL_DRAWS, NULL_PER_RESAMPLE = 20260911, 6000, 2000, 20


def summarise(x):
    x = np.asarray(x)
    return {"mean": round(float(x.mean()), 4),
            "ci95": [round(float(np.quantile(x, .025)), 4), round(float(np.quantile(x, .975)), 4)],
            "p_two_sided": round(float(2 * min((x <= 0).mean(), (x >= 0).mean())), 4)}


def d_of(w, s, g, c, ii, jj):
    return (cidx(w, ii, jj) - cidx(g, ii, jj)) - (cidx(s, ii, jj) - cidx(c, ii, jj))


def against_null(w, s, g, c, t, e, make_null, rng_null):
    """w, s: the construction with grade and with stage; g, c: the references; make_null(rng) returns
    one random construction's (with grade, with stage) pair."""
    n = len(t)
    ii, jj = cpairs(t, e)
    D = d_of(w, s, g, c, ii, jj)
    full = [d_of(*make_null(rng_null), g, c, ii, jj) for _ in range(NULL_DRAWS)]
    fixed = [make_null(rng_null) for _ in range(NULL_PER_RESAMPLE)]
    rng = np.random.default_rng(BOOT_SEED)
    ex = []
    for _ in range(REPS):
        b = rng.choice(n, size=n, replace=True)
        i2, j2 = cpairs(t[b], e[b])
        if not i2.size:
            continue
        dn = np.mean([d_of(nw[b], ns[b], g[b], c[b], i2, j2) for nw, ns in fixed])
        ex.append(d_of(w[b], s[b], g[b], c[b], i2, j2) - dn)
    return {"D": round(float(D), 4),
            "no_signal_D": {"mean": round(float(np.mean(full)), 4),
                            "q025_q975": [round(float(np.quantile(full, .025)), 4),
                                          round(float(np.quantile(full, .975)), 4)]},
            "excess_over_no_signal": {"point": round(float(D - np.mean(full)), 4), **summarise(ex)}}


# ------------------------------------------------------------------------------------------ TCGA-BLCA

def run_blca(a):
    from s33_field_inflation import exact, fold_pct, load_model
    rows = [json.loads(x) for x in open(a.percase) if x.strip()]
    n = len(rows)
    idx = {r["case_id"]: i for i, r in enumerate(rows)}
    t = np.array([r["months"] for r in rows], float)
    e = np.array([r["event"] for r in rows], float)
    fold = np.array([r["fold"] for r in rows], int)
    col = lambda k: exact(np.array([r[k] for r in rows], float), fold)
    zg, zc = col("clinical_grade"), col("clinical_stage")
    ii, jj = cpairs(t, e)
    committed = json.load(open(a.committed))["added_value_inflation"]["constructions"]
    committed_a = json.load(open(a.committed_a))["models"]

    def null2(rng):
        u = fold_pct(rng.random(n), fold)
        return fold_pct(u + zg, fold), fold_pct(u + zc, fold)

    def null3(rng):
        u1, u2 = fold_pct(rng.random(n), fold), fold_pct(rng.random(n), fold)
        return fold_pct(u1 + u2 + zg, fold), fold_pct(u1 + u2 + zc, fold)

    cons = {"ModRank": (col("ours_grade"), col("ours"), null3, committed["ours"]["D"]),
            "SurvPath": (col("survpath_plus_grade"), col("survpath_plus_stage"), null2,
                         committed["survpath"]["D"]),
            "PIBD": (col("pibd_best_val_plus_grade"), col("pibd_best_val_plus_stage"), null2,
                     committed["pibd_best_val"]["D"])}
    for m in a.models:
        seeds = [r for _, r, _ in load_model(os.path.join(a.runs, m), idx, n, t, e)]
        M = fold_pct(np.vstack([fold_pct(r, fold) for r in seeds]).mean(0), fold)
        cons[m] = (fold_pct(M + zg, fold), fold_pct(M + zc, fold), null2, committed_a[m]["D"])

    drift = {k: (round(d_of(w, s, zg, zc, ii, jj), 4), want) for k, (w, s, _, want) in cons.items()
             if round(d_of(w, s, zg, zc, ii, jj), 4) != want}
    if drift:
        print("s35 blca: known answer not reproduced: %s" % drift, file=sys.stderr)
        return 2
    rng_null = np.random.default_rng(35)
    out = {}
    for k, (w, s, mk, _) in cons.items():
        out[k] = against_null(w, s, zg, zc, t, e, mk, rng_null)
        out[k]["no_signal_arms"] = 2 if mk is null3 else 1
        print("%-22s D %+.4f | no-signal D %+.4f %s | excess %+.4f %s p %.4f"
              % (k, out[k]["D"], out[k]["no_signal_D"]["mean"], out[k]["no_signal_D"]["q025_q975"],
                 out[k]["excess_over_no_signal"]["point"], out[k]["excess_over_no_signal"]["ci95"],
                 out[k]["excess_over_no_signal"]["p_two_sided"]), file=sys.stderr, flush=True)
    return out


# ------------------------------------------------------------------------------------------ GEO

def run_geo(a):
    import s24_geo_external as s24
    from s5_omics_arm import pathway_matrix
    ann = {"GPL570": s24.read_annot(os.path.join(a.geo_dir, "GPL570.annot.gz")),
           "GPL6947": s24.read_annot(os.path.join(a.geo_dir, "GPL6947.annot.gz"))}

    def prep(fname, gpl):
        s, ch, pr, X = s24.read_series(os.path.join(a.geo_dir, fname))
        G, genes = s24.gene_matrix(pr, X, ann[gpl])
        P, _ = pathway_matrix(G, genes, a.signatures)
        return s, ch, P

    def num(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return np.nan

    cohorts = {}
    # GSE31684, parsed exactly as s24 parses it
    s, ch, P = prep("GSE31684_series_matrix.txt.gz", "GPL570")
    t = np.array([num(ch["survival.months"].get(x)) for x in s])
    e = np.array([1.0 if ch["last known status"].get(x, "") == "DOD" else 0.0 for x in s])
    age = np.array([num(ch["age at rc"].get(x)) for x in s])[:, None]
    fem = np.array([1.0 if ch["gender"].get(x, "").lower() == "female" else 0.0 for x in s])[:, None]
    grade, _ = s24.onehot([ch["rc grade"].get(x, "").capitalize() for x in s])
    stage, _ = s24.onehot([ch["rc_stage"].get(x, "") for x in s])
    ok = np.isfinite(t) & (t > 0) & np.isfinite(age[:, 0])
    cohorts["GSE31684"] = (P[ok], np.hstack([age, fem, stage])[ok], np.hstack([age, fem, grade])[ok],
                           t[ok], e[ok], np.array(s)[ok])
    # GSE32894
    s, ch, P = prep("GSE32894_series_matrix.txt.gz", "GPL6947")
    t = np.array([num(ch["time_to_dod_(months)"].get(x)) for x in s])
    ev = [ch["dod_event_(yes/no)"].get(x, "") for x in s]
    exc = np.array([ch["dod_excluded_(yes/no)"].get(x, "") == "yes" for x in s])
    e = np.array([1.0 if v == "yes" else 0.0 for v in ev])
    ok = np.isfinite(t) & (t > 0) & ~exc & np.array([v in ("yes", "no") for v in ev])
    age = np.array([num(ch["age"].get(x)) for x in s])[:, None]
    fem = np.array([1.0 if ch["gender"].get(x, "") == "F" else 0.0 for x in s])[:, None]
    grade, _ = s24.onehot([ch["tumor_grade"].get(x, "") for x in s])
    stage, _ = s24.onehot([ch["tumor_stage"].get(x, "") for x in s])
    ok &= np.isfinite(age[:, 0])
    cohorts["GSE32894"] = (P[ok], np.hstack([age, fem, stage])[ok], np.hstack([age, fem, grade])[ok],
                           t[ok], e[ok], np.array(s)[ok])

    rows = [json.loads(x) for x in open(a.percase) if x.strip()]
    out, drift = {}, {}
    rng_null = np.random.default_rng(35)
    for name, (P, cs, cw, t, e, ids) in cohorts.items():
        n = len(t)
        # s24.evaluate's splits and arms, repeat by repeat
        rng = np.random.default_rng(0)
        splits = [s24.strat_folds(e, s24.K, rng) for _ in range(s24.REPEATS)]
        reps = []
        for r, folds in enumerate(splits):
            zo = s24.fpct(s24.oof(P, t, e, folds, r), folds)
            zs = s24.fpct(s24.oof(cs, t, e, folds, r), folds)
            zw = s24.fpct(s24.oof(cw, t, e, folds, r), folds)
            reps.append((folds, zo, zs, zw))
        avg = lambda f: np.mean(np.vstack([f(*x) for x in reps]), 0)
        w = avg(lambda folds, zo, zs, zw: s24.fpct(zo + zw, folds))
        st = avg(lambda folds, zo, zs, zw: s24.fpct(zo + zs, folds))
        g = avg(lambda folds, zo, zs, zw: zw)
        c = avg(lambda folds, zo, zs, zw: zs)
        # known answer: the committed per-case rows of this analysis
        mine = {k: v for k, v in zip(ids, zip(w, st, g, c))}
        for row in rows:
            if row["analysis"] != name:
                continue
            got = mine[row["sample"]]
            want = (row["modrank_with_weak"], row["modrank"], row["clinical_weak"], row["clinical_stage"])
            if max(abs(x - y) for x, y in zip(got, want)) > 1e-9:
                drift[name] = "per-case vectors differ from the committed rows"
                break

        def make_null(rng_):
            us = [rng_.random(n) for _ in reps]
            nw = np.mean(np.vstack([s24.fpct(s24.fpct(u, f) + zw_, f) for u, (f, _, _, zw_) in zip(us, reps)]), 0)
            ns = np.mean(np.vstack([s24.fpct(s24.fpct(u, f) + zs_, f) for u, (f, _, zs_, _) in zip(us, reps)]), 0)
            return nw, ns
        out[name] = {"n": int(n), "events": int(e.sum()),
                     **against_null(w, st, g, c, t, e, make_null, rng_null), "no_signal_arms": 1}
        print("%-9s n %d ev %d | D %+.4f | no-signal D %+.4f %s | excess %+.4f %s p %.4f"
              % (name, n, int(e.sum()), out[name]["D"], out[name]["no_signal_D"]["mean"],
                 out[name]["no_signal_D"]["q025_q975"], out[name]["excess_over_no_signal"]["point"],
                 out[name]["excess_over_no_signal"]["ci95"], out[name]["excess_over_no_signal"]["p_two_sided"]),
              file=sys.stderr, flush=True)
    if drift:
        print("s35 geo: known answer not reproduced: %s" % drift, file=sys.stderr)
        return 2
    return out


# ------------------------------------------------------------------------------------------ five studies

def run_five(a):
    import pickle
    from s5_stage_five_cohorts import fold_pct as fpct5, load_cohort
    from s32_five_cohort_fusion import arm
    committed = json.load(open(a.committed_b))["inflation_B"]
    d0 = pickle.load(open(a.titan, "rb"))
    E = np.asarray(d0["embeddings"])
    stem2emb = {}
    for i, f in enumerate(d0["filenames"]):
        s_ = str(f)
        for suf in (".svs", ".h5", ".pt"):
            if s_.endswith(suf):
                s_ = s_[:-len(suf)]
        stem2emb[s_] = E[i]
    out, drift = {}, {}
    rng_null = np.random.default_rng(35)
    for c in ("blca", "hnsc", "stad"):
        co = load_cohort(os.path.join(a.meta_dir, "tcga_%s.csv" % c),
                         os.path.join(a.clin_dir, "tcga_%s_clinical.csv" % c),
                         os.path.join(a.pan_dir, "splits", c), os.path.join(a.pan_dir, "rna_%s.csv" % c),
                         a.sig, stem2emb)
        AS = np.hstack([co["age"], co["fem"]])
        zt, zo = arm(co["T"], co), arm(co["P"], co)
        ass, asg = arm(np.hstack([AS, co["stage"]]), co), arm(np.hstack([AS, co["grade"]]), co)
        w, st = fpct5(zt + zo + asg, co["assign"]), fpct5(zt + zo + ass, co["assign"])
        ii, jj = cpairs(co["t"], co["e"])
        got = round(d_of(w, st, asg, ass, ii, jj), 4)
        if got != committed[c]["D"]:
            drift[c] = (got, committed[c]["D"])
            continue
        n = len(co["keep"])

        def make_null(rng_):
            u1, u2 = fpct5(rng_.random(n), co["assign"]), fpct5(rng_.random(n), co["assign"])
            return fpct5(u1 + u2 + asg, co["assign"]), fpct5(u1 + u2 + ass, co["assign"])
        out[c] = {"n": int(n), "events": int(co["e"].sum()),
                  **against_null(w, st, asg, ass, co["t"], co["e"], make_null, rng_null), "no_signal_arms": 2}
        print("%-5s D %+.4f | no-signal D %+.4f %s | excess %+.4f %s p %.4f"
              % (c, out[c]["D"], out[c]["no_signal_D"]["mean"], out[c]["no_signal_D"]["q025_q975"],
                 out[c]["excess_over_no_signal"]["point"], out[c]["excess_over_no_signal"]["ci95"],
                 out[c]["excess_over_no_signal"]["p_two_sided"]), file=sys.stderr, flush=True)
    if drift:
        print("s35 five: known answer not reproduced: %s" % drift, file=sys.stderr)
        return 2
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("blca")
    for f in ("percase", "committed", "committed-a", "runs", "out"):
        b.add_argument("--" + f, required=True)
    b.add_argument("--models", nargs="*", default=[])
    g = sub.add_parser("geo")
    for f in ("protocol", "geo-dir", "signatures", "percase", "out"):
        g.add_argument("--" + f, required=True)
    f5 = sub.add_parser("five")
    for f in ("titan", "meta-dir", "clin-dir", "pan-dir", "sig", "committed-b", "out"):
        f5.add_argument("--" + f, required=True)
    a = ap.parse_args()
    t0 = time.time()
    res = {"blca": run_blca, "geo": run_geo, "five": run_five}[a.cmd](a)
    if res == 2:
        return 2
    json.dump({"artifact_type": "s35_inflation_null_%s" % a.cmd,
               "phase_of_origin": "post_freeze_2026-09-27",
               "basis": "experiments/20260927-field-inflation/README.md, finding 1",
               "bootstrap": {"replicates": REPS, "seed": BOOT_SEED, "unit": "cases"},
               "no_signal": {"draws_full_data": NULL_DRAWS, "draws_per_resample": NULL_PER_RESAMPLE,
                             "rng_seed": 35},
               "constructions": res, "runtime_seconds": round(time.time() - t0, 1)},
              open(a.out, "w"), indent=1)
    print("wrote %s (%.0f s)" % (a.out, time.time() - t0), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
