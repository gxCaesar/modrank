#!/usr/bin/env python3
"""The stage-versus-grade reference gap in external bladder cohorts (analysis C).

Implements experiments/20260927-field-inflation/protocol-c.md, written after the search and before any
download. Post-freeze and exploratory. Nothing here changes a committed value.

  s36_reference_gap_external.py gate0   --geo-c DIR --retained DIR --out gate0.json
  s36_reference_gap_external.py analyze --geo-c DIR --retained DIR --gate0 gate0.json \\
        --geo-dir DIR --committed-geo geo-external.json --out gap.json

gate0 parses each cohort's clinical fields into one row per patient and counts, before any outcome is
analysed, the patients with age, sex, grade, stage, a time and an event indicator, the events among
them, and the levels of grade and stage. A cohort is admitted with at least 50 such patients, 20
events, and two or more levels of both grade and stage.

analyze fits, in each admitted cohort, two clinical models by the GEO protocol's recipe (ridge Cox,
penalty by inner three-fold cross-validation, five repeats of stratified five-fold cross-validation,
analysis/s24_geo_external.py's own functions): age, sex and grade, and age, sex and stage. The gap is
the mean over repeats of C(stage model) minus that of C(grade model), and its interval is the paired
case bootstrap of the repeat-averaged scores (6,000 resamples, seed 20260911). Before any new cohort
is scored, the same code must reproduce the committed stage-block and grade-block concordances of the
two GEO cohorts analysed before (GSE31684, GSE32894).
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s24_geo_external as s24                                                 # noqa: E402
from blca_common import cidx, cpairs                                           # noqa: E402

MIN_N, MIN_EVENTS = 50, 20


def num(x):
    try:
        return float(re.sub(r"[^0-9.eE+-]", "", str(x)))
    except ValueError:
        return float("nan")


def tcat(v):
    m = re.search(r"T(is|a|\d)", str(v))
    return ("T" + m.group(1)) if m else ""


def series_chars(path):
    """One dict per sample of 'key: value' characteristics, which GEO does not align by row."""
    samples, chars = [], collections.defaultdict(dict)
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                samples = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith(("!Sample_characteristics_ch1", "!Sample_source_name_ch1",
                                  "!Sample_title")):
                tag = line.split("\t", 1)[0]
                for s_, v in zip(samples, line.rstrip("\n").split("\t")[1:]):
                    v = v.strip('"')
                    if tag != "!Sample_characteristics_ch1":
                        chars[s_][tag] = v
                    elif ":" in v:
                        k, val = v.split(":", 1)
                        chars[s_][k.strip().lower()] = val.strip()
            elif line.startswith("!series_matrix_table_begin"):
                break
    return samples, chars


def sdrf_rows(path, key):
    """One row per individual from an SDRF, which repeats an individual once per assay."""
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    seen = collections.OrderedDict()
    for r in rows:
        seen.setdefault(r[key], r)
    return list(seen.values())


# ------------------------------------------------------------------------------------------ cohorts

def gse19915(d):
    s, ch = series_chars(os.path.join(d, "GSE19915-GPL5186_series_matrix.txt.gz"))
    out = []
    for x in s:
        c = ch[x]
        if "tumor stage" not in c and "tumor grade" not in c:
            continue                                           # normal bladder
        dod = c.get("dead of disease", "").upper()
        out.append({"id": c.get("!Sample_title", x), "age": float("nan"), "fem": float("nan"),
                    # "G3^Nested" is grade 3 with the nested variant appended (protocol-c.md, gate 0)
                    "grade": re.sub(r"^(G\d)\^.*$", r"\1", c.get("tumor grade", "")),
                    "stage": tcat(c.get("tumor stage", "")),
                    "time": num(c.get("follow up time (months)", "")),
                    "event": 1.0 if dod == "YES" else 0.0 if dod == "NO" else float("nan")})
    return out, "disease-specific survival", False


def emtab1803(d):
    out = []
    for r in sdrf_rows(os.path.join(d, "E-MTAB-1803.sdrf.txt"), "Characteristics[individual]"):
        v = r.get("Characteristics[viability]", "").strip().lower()
        sex = r.get("Characteristics[sex]", "").strip().lower()
        out.append({"id": r["Characteristics[individual]"], "age": num(r.get("Characteristics[age]", "")),
                    "fem": 1.0 if sex == "female" else 0.0 if sex == "male" else float("nan"),
                    "grade": r.get("Comment[Tumor grade73]", "").strip(),
                    "stage": tcat(r.get("Characteristics[tumor stage]", "")),
                    "time": num(r.get("Comment[Follow up time]", "")),
                    "event": 1.0 if v == "dead" else 0.0 if v == "alive" else float("nan")})
    return out, "overall survival", True


def gse13507(d):
    s, ch = series_chars(os.path.join(d, "GSE13507_series_matrix.txt.gz"))
    out = []
    for x in s:
        c = ch[x]
        if c.get("biological source", "").lower() != "primary bladder cancer":
            continue
        css = c.get("cancer specific survival", "").lower()
        sex = c.get("sex", "").upper()
        out.append({"id": x, "age": num(c.get("age", "")),
                    "fem": 1.0 if sex == "F" else 0.0 if sex == "M" else float("nan"),
                    "grade": c.get("grade", ""), "stage": tcat(c.get("stage", "")),
                    "time": num(c.get("survival month", "")),
                    "event": 1.0 if css.startswith("death") else 0.0 if css == "survival" else float("nan")})
    return out, "cancer-specific survival", True


def emtab4321(d):
    out = []
    for r in sdrf_rows(os.path.join(d, "E-MTAB-4321.sdrf.txt"), "Source Name"):
        p = r.get("Characteristics[progression to T2+]", "").strip().lower()
        sex = r.get("Characteristics[sex]", "").strip().lower()
        st = r.get("Characteristics[disease staging]", "").strip()
        out.append({"id": r["Source Name"], "age": num(r.get("Characteristics[age]", "")),
                    "fem": 1.0 if sex == "female" else 0.0 if sex == "male" else float("nan"),
                    "grade": r.get("Characteristics[tumor grading]", "").strip(),
                    "stage": st if st.upper() == "CIS" else tcat(st),
                    "time": num(r.get("Characteristics[progression free survival]", "")),
                    "event": 1.0 if p in ("yes", "1", "true") else 0.0 if p in ("no", "0", "false")
                    else float("nan")})
    return out, "progression-free survival (progression to T2+)", True


COHORTS = {"GSE19915": ("geo-c", gse19915), "E-MTAB-1803": ("geo-c", emtab1803),
           "GSE13507": ("retained", gse13507), "E-MTAB-4321": ("retained", emtab4321)}
NO_SURVIVAL = {"GSE5479": "GSE5479_clinical_information.txt carries stage, grade, sex and age and no "
                          "time or event field; the series has no matrix",
               "GSE1827": "the series matrix carries stage and some grade in free-text descriptions "
                          "and no characteristics, time or event field"}


def entropy(vals):
    c = collections.Counter(vals)
    n = sum(c.values())
    if len(c) < 2:
        return 0.0
    return -sum(v / n * math.log(v / n) for v in c.values()) / math.log(len(c))


def complete(rows, with_age_sex):
    keys = ("grade", "stage", "time", "event") + (("age", "fem") if with_age_sex else ())
    ok = []
    for r in rows:
        if any((isinstance(r[k], float) and not np.isfinite(r[k])) for k in keys if k != "grade" and k != "stage"):
            continue
        if r["grade"].lower() in s24.UNKNOWN or r["stage"].lower() in s24.UNKNOWN:
            continue
        if not r["time"] > 0:
            continue
        ok.append(r)
    return ok


def load(a):
    out = {}
    for name, (where, fn) in COHORTS.items():
        base = a.geo_c if where == "geo-c" else os.path.join(a.retained, name)
        rows, endpoint, with_age_sex = fn(base)
        out[name] = (complete(rows, with_age_sex), rows, endpoint, with_age_sex)
    return out


# ------------------------------------------------------------------------------------------ gate 0

def run_gate0(a):
    res = {"artifact_type": "s36_gate0", "phase_of_origin": "post_freeze_2026-09-27",
           "rule": "admit with >= %d complete patients, >= %d events, >= 2 levels of grade and of stage"
                   % (MIN_N, MIN_EVENTS), "cohorts": {}, "not_admissible_from_deposited_files": NO_SURVIVAL}
    for name, (ok, rows, endpoint, was) in load(a).items():
        g = collections.Counter(r["grade"] for r in ok)
        st = collections.Counter(r["stage"] for r in ok)
        ev = int(sum(r["event"] for r in ok))
        admit = len(ok) >= MIN_N and ev >= MIN_EVENTS and len(g) >= 2 and len(st) >= 2
        res["cohorts"][name] = {"endpoint": endpoint, "age_and_sex_in_models": was,
                                "patients_parsed": len(rows), "patients_complete": len(ok), "events": ev,
                                "grade_levels": dict(sorted(g.items())), "stage_levels": dict(sorted(st.items())),
                                "grade_entropy": round(entropy([r["grade"] for r in ok]), 3),
                                "stage_entropy": round(entropy([r["stage"] for r in ok]), 3),
                                "admitted": admit}
        print("%-12s parsed %4d complete %4d events %3d | grade %s | stage %s | admitted %s"
              % (name, len(rows), len(ok), ev, dict(g), dict(st), admit), file=sys.stderr)
    json.dump(res, open(a.out, "w"), indent=1)
    print("wrote %s" % a.out, file=sys.stderr)
    return 0


# ------------------------------------------------------------------------------------------ analysis

def fit_two(cs, cw, t, e):
    """s24.evaluate's splits and clinical arms, repeat by repeat."""
    rng = np.random.default_rng(0)
    splits = [s24.strat_folds(e, s24.K, rng) for _ in range(s24.REPEATS)]
    ii, jj = cpairs(t, e)
    per = {"stage": [], "grade": []}
    vec = {"stage": [], "grade": []}
    for r, folds in enumerate(splits):
        for k, X in (("stage", cs), ("grade", cw)):
            z = s24.fpct(s24.oof(X, t, e, folds, r), folds)
            per[k].append(cidx(z, ii, jj))
            vec[k].append(z)
    avg = {k: np.mean(np.vstack(v), 0) for k, v in vec.items()}
    return {k: float(np.mean(v)) for k, v in per.items()}, avg


def run_analyze(a):
    t0 = time.time()
    # known answers: the two GEO cohorts' stage and grade blocks, through this code path
    committed = json.load(open(a.committed_geo))["cohorts"]
    ann = {"GPL570": s24.read_annot(os.path.join(a.geo_dir, "GPL570.annot.gz")),
           "GPL6947": s24.read_annot(os.path.join(a.geo_dir, "GPL6947.annot.gz"))}
    known, drift = {}, {}
    for fname, gpl, name in (("GSE31684_series_matrix.txt.gz", "GPL570", "GSE31684"),
                             ("GSE32894_series_matrix.txt.gz", "GPL6947", "GSE32894")):
        s, ch, pr, X = s24.read_series(os.path.join(a.geo_dir, fname))
        if name == "GSE31684":
            t = np.array([num(ch["survival.months"].get(x)) for x in s])
            e = np.array([1.0 if ch["last known status"].get(x, "") == "DOD" else 0.0 for x in s])
            age = np.array([num(ch["age at rc"].get(x)) for x in s])[:, None]
            fem = np.array([1.0 if ch["gender"].get(x, "").lower() == "female" else 0.0 for x in s])[:, None]
            glab = [ch["rc grade"].get(x, "").capitalize() for x in s]
            grade, _ = s24.onehot(glab)
            stage, _ = s24.onehot([ch["rc_stage"].get(x, "") for x in s])
            ok = np.isfinite(t) & (t > 0) & np.isfinite(age[:, 0])
        else:
            t = np.array([num(ch["time_to_dod_(months)"].get(x)) for x in s])
            ev = [ch["dod_event_(yes/no)"].get(x, "") for x in s]
            exc = np.array([ch["dod_excluded_(yes/no)"].get(x, "") == "yes" for x in s])
            e = np.array([1.0 if v == "yes" else 0.0 for v in ev])
            age = np.array([num(ch["age"].get(x)) for x in s])[:, None]
            fem = np.array([1.0 if ch["gender"].get(x, "") == "F" else 0.0 for x in s])[:, None]
            glab = [ch["tumor_grade"].get(x, "") for x in s]
            grade, _ = s24.onehot(glab)
            stage, _ = s24.onehot([ch["tumor_stage"].get(x, "") for x in s])
            ok = np.isfinite(t) & (t > 0) & ~exc & np.array([v in ("yes", "no") for v in ev]) & np.isfinite(age[:, 0])
        cv, avg = fit_two(np.hstack([age, fem, stage])[ok], np.hstack([age, fem, grade])[ok], t[ok], e[ok])
        want = {"stage": committed[name]["cv"]["clinical_stage"]["mean"],
                "grade": committed[name]["cv"]["clinical_weak"]["mean"]}
        for k in ("stage", "grade"):
            if round(cv[k], 4) != want[k]:
                drift["%s/%s" % (name, k)] = (round(cv[k], 4), want[k])
        gap = s24.boot_pair(avg["stage"], avg["grade"], t[ok], e[ok])
        known[name] = {"n": int(ok.sum()), "events": int(e[ok].sum()),
                       "grade_entropy": round(entropy([g for g, k in zip(glab, ok)
                                                       if k and g.lower() not in s24.UNKNOWN]), 3),
                       "C_stage_model": round(cv["stage"], 4), "C_grade_model": round(cv["grade"], 4),
                       "gap": round(cv["stage"] - cv["grade"], 4), "gap_boot": gap}
        print("known %-9s stage %.4f grade %.4f (committed %.4f %.4f)"
              % (name, cv["stage"], cv["grade"], want["stage"], want["grade"]), file=sys.stderr, flush=True)
    if drift:
        print("s36: known answer not reproduced: %s" % drift, file=sys.stderr)
        return 2

    gate = json.load(open(a.gate0))["cohorts"]
    out = {}
    for name, (ok, rows, endpoint, was) in load(a).items():
        if not gate[name]["admitted"]:
            out[name] = {"admitted": False}
            continue
        t = np.array([r["time"] for r in ok])
        e = np.array([r["event"] for r in ok])
        glv, slv = sorted({r["grade"] for r in ok}), sorted({r["stage"] for r in ok})
        G = np.array([[1.0 if r["grade"] == l else 0.0 for l in glv[1:]] for r in ok])
        S = np.array([[1.0 if r["stage"] == l else 0.0 for l in slv[1:]] for r in ok])
        base = np.array([[r["age"], r["fem"]] for r in ok]) if was else np.zeros((len(ok), 0))
        cv, avg = fit_two(np.hstack([base, S]), np.hstack([base, G]), t, e)
        gap = s24.boot_pair(avg["stage"], avg["grade"], t, e)
        out[name] = {"admitted": True, "endpoint": endpoint, "n": len(ok), "events": int(e.sum()),
                     "age_and_sex_in_models": was, "grade_entropy": gate[name]["grade_entropy"],
                     "stage_entropy": gate[name]["stage_entropy"],
                     "C_stage_model": round(cv["stage"], 4), "C_grade_model": round(cv["grade"], 4),
                     "gap": round(cv["stage"] - cv["grade"], 4), "gap_boot": gap}
        print("%-12s n %d ev %d | stage %.4f grade %.4f | gap %+.4f %s | grade entropy %.3f"
              % (name, len(ok), int(e.sum()), cv["stage"], cv["grade"], cv["stage"] - cv["grade"],
                 gap["ci95"], gate[name]["grade_entropy"]), file=sys.stderr, flush=True)
    json.dump({"artifact_type": "s36_reference_gap_external", "phase_of_origin": "post_freeze_2026-09-27",
               "protocol": "experiments/20260927-field-inflation/protocol-c.md",
               "bootstrap": {"replicates": s24.REPS, "seed": s24.BOOT_SEED, "unit": "cases"},
               "known_answers_geo": known, "cohorts": out, "runtime_seconds": round(time.time() - t0, 1)},
              open(a.out, "w"), indent=1)
    print("wrote %s (%.0f s)" % (a.out, time.time() - t0), file=sys.stderr)
    return 0


def run_summary(a):
    """Every bladder-or-benchmark cohort with a usable grade, one row each: the gap and grade's entropy.

    TCGA rows are the committed five-study values (released folds, the benchmark's clinical file, one
    penalty selection): the point from stage5.json and the interval from five-cohort-intervals.json.
    Breast and colorectal carry no usable grade and are left out, as in the paper's Figure 2d.
    """
    from scipy.stats import spearmanr
    ex = json.load(open(a.gap))
    s5 = json.load(open(a.stage5))["cohorts"]
    iv = json.load(open(a.intervals))["cohorts"]
    why = json.load(open(a.why))["cohorts"]
    rows = []
    for c in ("blca", "hnsc", "stad"):
        d = iv[c]["differences"]["stage_minus_grade"]
        rows.append({"cohort": "TCGA-" + c.upper(), "source": "benchmark, five-study protocol",
                     "gap": s5[c]["stage_minus_grade_same_pipeline"], "ci95": d["ci95"],
                     "grade_entropy": why[c]["grade"]["normalised_entropy"]})
    for k, v in ex["known_answers_geo"].items():
        rows.append({"cohort": k, "source": "external, analysed before", "gap": v["gap"],
                     "ci95": v["gap_boot"]["ci95"], "grade_entropy": v["grade_entropy"]})
    for k, v in ex["cohorts"].items():
        if v.get("admitted"):
            rows.append({"cohort": k, "source": "external, analysis C", "gap": v["gap"],
                         "ci95": v["gap_boot"]["ci95"], "grade_entropy": v["grade_entropy"]})
    rho, p = spearmanr([r["grade_entropy"] for r in rows], [r["gap"] for r in rows])
    out = {"artifact_type": "s36_reference_gap_summary", "rows": rows, "n": len(rows),
           "gap_positive": sum(r["gap"] > 0 for r in rows),
           "interval_excludes_zero": sum(r["ci95"][0] > 0 for r in rows),
           "spearman_gap_vs_grade_entropy": {"rho": round(float(rho), 3), "p": round(float(p), 3)},
           "external_only": {"n": sum(r["source"].startswith("external") for r in rows),
                             "gap_positive": sum(r["gap"] > 0 for r in rows if r["source"].startswith("external")),
                             "interval_excludes_zero": sum(r["ci95"][0] > 0 for r in rows
                                                           if r["source"].startswith("external"))}}
    for r in rows:
        print("%-14s %-32s gap %+.4f %s grade entropy %.3f" % (r["cohort"], r["source"], r["gap"],
                                                               r["ci95"], r["grade_entropy"]), file=sys.stderr)
    print("positive %d of %d, interval excluding zero %d; Spearman rho %.3f (p %.3f)"
          % (out["gap_positive"], out["n"], out["interval_excludes_zero"], rho, p), file=sys.stderr)
    json.dump(out, open(a.out, "w"), indent=1)
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gate0")
    for f in ("geo-c", "retained", "out"):
        g.add_argument("--" + f, required=True)
    z = sub.add_parser("analyze")
    for f in ("geo-c", "retained", "gate0", "geo-dir", "committed-geo", "out"):
        z.add_argument("--" + f, required=True)
    m = sub.add_parser("summary")
    for f in ("gap", "stage5", "intervals", "why", "out"):
        m.add_argument("--" + f, required=True)
    a = ap.parse_args()
    return {"gate0": run_gate0, "analyze": run_analyze, "summary": run_summary}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
