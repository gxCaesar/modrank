#!/usr/bin/env python3
"""The combination rule and the added-value inflation, in independent public GEO bladder cohorts.

Implements experiments/20260911-geo-external/protocol.json, frozen (f0fc6bf) before any GEO outcome
was analysed. Three cohorts with expression and clinical data and no histology: GSE31684 (primary,
cystectomy), GSE32894 and GSE48075 (secondary). The recipe is re-fitted inside each cohort by
repeated cross-validation; nothing fitted on TCGA is transferred across platforms.

For each cohort: the two-arm ModRank (transcriptome pathways and the clinical stage block, equal
weight rank average) against its components and the two fusion alternatives, the added-value
inflation D where a grade field exists, and a split-reseed floor.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import hashlib
import json
import os
import re
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blca_common import ALPHAS, cidx, cox_fit, cpairs, fitapply, pct       # noqa: E402
from s5_omics_arm import pathway_matrix                                   # noqa: E402

EXT = ALPHAS + (32768.0,)
REPEATS, K, REPS, BOOT_SEED = 5, 5, 6000, 20260911
UNKNOWN = {"", "na", "nan", "unknown", "gx", "tx", "nx", "none", "not available"}


def strat_folds(e, k, rng):
    """Stratified by event: the same construction analysis/s18_chimera_two_arm.py uses."""
    idx = {c: rng.permutation(np.flatnonzero(e == c)) for c in (0, 1)}
    folds = [[] for _ in range(k)]
    for c in (0, 1):
        for j, i in enumerate(idx[c]):
            folds[j % k].append(i)
    return [np.array(sorted(f)) for f in folds]


def onehot(vals):
    lv = sorted({v for v in vals if v.lower() not in UNKNOWN})
    if len(lv) < 2:
        return np.zeros((len(vals), 0)), lv
    return np.array([[1.0 if v == l else 0.0 for l in lv[1:]] for v in vals]), lv


def read_series(path):
    samples, chars, probes, rows, intab = [], collections.defaultdict(dict), [], [], False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                samples = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                for s, v in zip(samples, line.rstrip("\n").split("\t")[1:]):
                    v = v.strip().strip('"')
                    if ":" in v:
                        k_, val = v.split(":", 1)
                        chars[k_.strip().lower()][s] = val.strip()
            elif line.startswith("!series_matrix_table_begin"):
                intab = True
                next(fh)
            elif line.startswith("!series_matrix_table_end"):
                break
            elif intab:
                p = line.rstrip("\n").split("\t")
                probes.append(p[0].strip('"'))
                rows.append([float(x) if x.strip('"') not in ("", "null", "NA") else np.nan for x in p[1:]])
    return samples, chars, probes, np.asarray(rows, dtype=float)


def read_annot(path):
    """probe -> gene symbol from a GEO .annot table; ambiguous or empty symbols are dropped."""
    out, intab, col = {}, False, None
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                intab = True
                hdr = next(fh).rstrip("\n").split("\t")
                col = hdr.index("Gene symbol")
                continue
            if line.startswith("!platform_table_end"):
                break
            if intab:
                p = line.rstrip("\n").split("\t")
                if len(p) > col and p[col] and "///" not in p[col]:
                    out[p[0]] = p[col]
    return out


def gene_matrix(probes, X, annot):
    by = collections.defaultdict(list)
    for i, pr in enumerate(probes):
        g = annot.get(pr)
        if g:
            by[g].append(i)
    genes = sorted(by)
    G = np.vstack([np.nanmean(X[by[g]], axis=0) for g in genes]).T     # samples x genes
    G = np.where(np.isfinite(G), G, np.nanmean(G, axis=0))
    return G, genes


def fpct(v, folds):
    out = np.full(len(v), np.nan)
    for va in folds:
        out[va] = pct(v[va])
    return out


def oof(X, t, e, folds, seed, alphas=ALPHAS):
    s = np.full(len(t), np.nan)
    for va in folds:
        tr = np.setdiff1d(np.arange(len(t)), va)
        s[va] = fitapply(X[tr], t[tr], e[tr], X[va], seed, alphas)
    return s


def stacked(blocks, t, e, folds, seed, zval):
    s = np.full(len(t), np.nan)
    for va in folds:
        tr = np.setdiff1d(np.arange(len(t)), va)
        m = len(tr)
        cols = []
        for X in blocks:
            io = np.full(m, np.nan)
            for ite in np.array_split(np.random.default_rng(seed).permutation(m), 3):
                itr = np.setdiff1d(np.arange(m), ite)
                io[ite] = pct(fitapply(X[tr][itr], t[tr][itr], e[tr][itr], X[tr][ite], seed))
            cols.append(io)
        Zi = np.column_stack(cols)
        w = cox_fit(Zi - Zi.mean(0), t[tr], e[tr], 1.0)
        s[va] = np.column_stack([z[va] for z in zval]) @ w
    return fpct(s, folds)


def boot_pair(x, y, t, e):
    rng = np.random.default_rng(BOOT_SEED)
    v = []
    for _ in range(REPS):
        bs = rng.choice(len(t), size=len(t), replace=True)
        i2, j2 = cpairs(t[bs], e[bs])
        v.append(cidx(x[bs], i2, j2) - cidx(y[bs], i2, j2))
    v = np.asarray(v)
    return {"mean": round(float(v.mean()), 4),
            "ci95": [round(float(np.quantile(v, .025)), 4), round(float(np.quantile(v, .975)), 4)],
            "p_two_sided": round(float(2 * min((v <= 0).mean(), (v >= 0).mean())), 4)}


def holm(p):
    keys = sorted(p, key=lambda k_: p[k_])
    run, adj = 0.0, {}
    for i, k_ in enumerate(keys):
        run = max(run, min(1.0, (len(keys) - i) * p[k_]))
        adj[k_] = round(run, 4)
    return adj


def evaluate(name, P, clin_stage, clin_weak, t, e):
    ii, jj = cpairs(t, e)
    rng = np.random.default_rng(0)
    splits = [strat_folds(e, K, rng) for _ in range(REPEATS)]
    per, vec = collections.defaultdict(list), collections.defaultdict(list)
    for r, folds in enumerate(splits):
        zo = fpct(oof(P, t, e, folds, r), folds)
        zs = fpct(oof(clin_stage, t, e, folds, r), folds)
        arms = {"transcriptome": zo, "clinical_stage": zs,
                "modrank": fpct(zo + zs, folds),
                "concatenated": fpct(oof(np.hstack([P, clin_stage]), t, e, folds, r, EXT), folds),
                "stacked": stacked([P, clin_stage], t, e, folds, r, [zo, zs])}
        if clin_weak is not None:
            zw = fpct(oof(clin_weak, t, e, folds, r), folds)
            arms["clinical_weak"] = zw
            arms["modrank_with_weak"] = fpct(zo + zw, folds)
        for k_, v in arms.items():
            per[k_].append(cidx(v, ii, jj))
            vec[k_].append(v)
    avg = {k_: np.mean(np.vstack(v), 0) for k_, v in vec.items()}
    res = {"n": int(len(t)), "events": int(e.sum()), "comparable_pairs": int(len(ii)),
           "cv": {k_: {"mean": round(float(np.mean(v)), 4), "sd": round(float(np.std(v)), 4),
                       "repeat_averaged_vector_c": round(cidx(avg[k_], ii, jj), 4)}
                  for k_, v in per.items()}}
    comps = [c for c in ("clinical_stage", "transcriptome", "concatenated", "stacked") if c in avg]
    pb = {c: boot_pair(avg["modrank"], avg[c], t, e) for c in comps}
    adj = holm({c: v["p_two_sided"] for c, v in pb.items()})
    res["paired_vs_modrank"] = {c: {**v, "p_holm": adj[c]} for c, v in pb.items()}
    deltas = []
    for k_ in range(24):
        f2 = strat_folds(e, K, np.random.default_rng(3000 + k_))
        zo = fpct(oof(P, t, e, f2, 0), f2)
        zs = fpct(oof(clin_stage, t, e, f2, 0), f2)
        deltas.append(cidx(fpct(zo + zs, f2), ii, jj) - cidx(zs, ii, jj))
    bar = 2 * float(np.std(deltas, ddof=1))
    res["floor"] = {"sd": round(bar / 2, 4), "bar": round(bar, 4)}
    margins = {c: round(res["cv"]["modrank"]["mean"] - res["cv"][c]["mean"], 4) for c in comps}
    res["decision"] = {"margins": margins, "clears_every_comparator": all(m > bar for m in margins.values())}
    if clin_weak is not None:
        rng_d = np.random.default_rng(BOOT_SEED)
        dv = []
        for _ in range(REPS):
            bs = rng_d.choice(len(t), size=len(t), replace=True)
            i2, j2 = cpairs(t[bs], e[bs])
            dv.append((cidx(avg["modrank_with_weak"][bs], i2, j2) - cidx(avg["clinical_weak"][bs], i2, j2))
                      - (cidx(avg["modrank"][bs], i2, j2) - cidx(avg["clinical_stage"][bs], i2, j2)))
        res["D"] = {"point": round((cidx(avg["modrank_with_weak"], ii, jj) - cidx(avg["clinical_weak"], ii, jj))
                                   - (cidx(avg["modrank"], ii, jj) - cidx(avg["clinical_stage"], ii, jj)), 4),
                    "ci95": [round(float(np.quantile(dv, .025)), 4), round(float(np.quantile(dv, .975)), 4)],
                    "p_two_sided": round(float(2 * min((np.asarray(dv) <= 0).mean(), (np.asarray(dv) >= 0).mean())), 4)}
    print("  %s: n %d events %d | modrank %.4f stage %.4f omics %.4f concat %.4f stacked %.4f%s"
          % (name, len(t), int(e.sum()), res["cv"]["modrank"]["mean"], res["cv"]["clinical_stage"]["mean"],
             res["cv"]["transcriptome"]["mean"], res["cv"]["concatenated"]["mean"], res["cv"]["stacked"]["mean"],
             "" if clin_weak is None else " | D %+.4f" % res["D"]["point"]), file=sys.stderr, flush=True)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for f in ("protocol", "geo-dir", "signatures", "out"):
        ap.add_argument("--" + f, required=True)
    a = ap.parse_args()
    t0 = time.time()
    rep = {"artifact_type": "s24_geo_external", "phase_of_origin": "post_freeze_2026-09-11",
           "protocol_sha256": hashlib.sha256(open(a.protocol, "rb").read()).hexdigest(),
           "reportable": True, "cohorts": {}}
    ann = {"GPL570": read_annot(os.path.join(a.geo_dir, "GPL570.annot.gz")),
           "GPL6947": read_annot(os.path.join(a.geo_dir, "GPL6947.annot.gz"))}

    def prep(fname, gpl):
        s, ch, pr, X = read_series(os.path.join(a.geo_dir, fname))
        G, genes = gene_matrix(pr, X, ann[gpl])
        P, kept = pathway_matrix(G, genes, a.signatures)
        return s, ch, P, kept, len(genes)

    def num(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return np.nan

    # ---- GSE31684, primary
    s, ch, P, kept, ng = prep("GSE31684_series_matrix.txt.gz", "GPL570")
    t = np.array([num(ch["survival.months"].get(x)) for x in s])
    e = np.array([1.0 if ch["last known status"].get(x, "") == "DOD" else 0.0 for x in s])
    age = np.array([num(ch["age at rc"].get(x)) for x in s])[:, None]
    fem = np.array([1.0 if ch["gender"].get(x, "").lower() == "female" else 0.0 for x in s])[:, None]
    grade, glv = onehot([ch["rc grade"].get(x, "").capitalize() for x in s])
    stage, slv = onehot([ch["rc_stage"].get(x, "") for x in s])
    ok = np.isfinite(t) & (t > 0) & np.isfinite(age[:, 0])
    rep["cohorts"]["GSE31684"] = {"genes_mapped": ng, "pathways_kept": len(kept),
                                  "grade_levels": glv, "stage_levels": slv, "excluded_incomplete": int((~ok).sum()),
                                  **evaluate("GSE31684", P[ok], np.hstack([age, fem, stage])[ok],
                                             np.hstack([age, fem, grade])[ok], t[ok], e[ok])}
    chemo = np.array([ch["prerc_chemo"].get(x, "").lower() == "yes" for x in s])
    ok2 = ok & ~chemo
    rep["cohorts"]["GSE31684"]["sensitivity_without_pre_cystectomy_chemotherapy"] = evaluate(
        "GSE31684 no pre-RC chemo", P[ok2], np.hstack([age, fem, stage])[ok2],
        np.hstack([age, fem, grade])[ok2], t[ok2], e[ok2])

    # ---- GSE32894, secondary
    s, ch, P, kept, ng = prep("GSE32894_series_matrix.txt.gz", "GPL6947")
    t = np.array([num(ch["time_to_dod_(months)"].get(x)) for x in s])
    ev = [ch["dod_event_(yes/no)"].get(x, "") for x in s]
    exc = np.array([ch["dod_excluded_(yes/no)"].get(x, "") == "yes" for x in s])
    e = np.array([1.0 if v == "yes" else 0.0 for v in ev])
    ok = np.isfinite(t) & (t > 0) & ~exc & np.array([v in ("yes", "no") for v in ev])
    age = np.array([num(ch["age"].get(x)) for x in s])[:, None]
    fem = np.array([1.0 if ch["gender"].get(x, "") == "F" else 0.0 for x in s])[:, None]
    grade, glv = onehot([ch["tumor_grade"].get(x, "") for x in s])
    stage, slv = onehot([ch["tumor_stage"].get(x, "") for x in s])
    ok &= np.isfinite(age[:, 0])
    res = evaluate("GSE32894", P[ok], np.hstack([age, fem, stage])[ok], np.hstack([age, fem, grade])[ok],
                   t[ok], e[ok])
    rule = evaluate("GSE32894 rule block", P[ok], np.hstack([age, fem, stage, grade])[ok], None, t[ok], e[ok])
    rep["cohorts"]["GSE32894"] = {"genes_mapped": ng, "pathways_kept": len(kept), "grade_levels": glv,
                                  "stage_levels": slv, **res, "with_the_rule_block": rule["cv"]}

    # ---- GSE48075, secondary, method only
    s, ch, P, kept, ng = prep("GSE48075_series_matrix.txt.gz", "GPL6947")
    t = np.array([num(ch["survival (mo)"].get(x)) for x in s])
    cen = [ch["dss censor"].get(x, "") for x in s]
    e = np.array([1.0 if c == "uncensored" else 0.0 for c in cen])
    age = np.array([num(ch["age (at specimen collection)"].get(x)) for x in s])[:, None]
    fem = np.array([1.0 if ch["gender"].get(x, "").lower() == "female" else 0.0 for x in s])[:, None]

    def tcat(v):
        m = re.match(r"^c?T(is|a|\d)", v.strip())
        return ("T" + m.group(1)) if m else ""
    tst, tlv = onehot([tcat(ch["cstage"].get(x, "")) for x in s])
    ok = np.isfinite(t) & (t > 0) & np.array([c in ("censored", "uncensored") for c in cen]) & np.isfinite(age[:, 0])
    rep["cohorts"]["GSE48075"] = {"genes_mapped": ng, "pathways_kept": len(kept), "t_levels": tlv,
                                  **evaluate("GSE48075", P[ok], np.hstack([age, fem, tst])[ok], None, t[ok], e[ok])}
    rep["runtime_seconds"] = round(time.time() - t0, 1)
    json.dump(rep, open(a.out, "w"), indent=1)
    print("wrote %s (%.0f s)" % (a.out, time.time() - t0), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
