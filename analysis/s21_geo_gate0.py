#!/usr/bin/env python3
"""GATE 0 for public GEO bladder cohorts: do expression, stage, grade, age, sex and a timed outcome
CO-EXIST in the same downloadable samples?

RUNS BEFORE ANY MODELLING CODE, and reads no outcome VALUE beyond whether it is present. The output
is a count table per series, not a verdict: a field is only counted as present when its values are
non-missing, and the mapping from GEO's free-text characteristics keys to the variables the design
needs is printed in full, so a person adjudicates it rather than trusting a keyword match.

Every characteristic key is reported with its level counts (top ten), because GEO keys are free text
and a heuristic that matched 'stage' to a key holding 'stage of treatment' would pass silently.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
import os
import re
import sys

MISSING = {"", "na", "n/a", "nan", "none", "unknown", "not available", "--", "-", "x", "nx", "tx",
           "[not available]", "[unknown]", "not evaluated"}
PATTERNS = {
    "stage": r"(^|[\s_])(stage|pt|t[\s_-]?stage|tumou?r stage|clinical t|pathologic.*t|tnm)",
    "grade": r"grade",
    "age": r"(^|[\s_])age",
    "sex": r"(sex|gender)",
    "time": r"(time|months|follow|survival.*(time|month)|days)",
    "event": r"(status|dead|death|event|died|vital|dss|progression|recurrence|censor)",
}


def parse(path):
    samples, chars, n_probes, platform = [], collections.defaultdict(dict), 0, None
    in_table = False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Series_platform_id"):
                platform = line.split("\t")[1].strip().strip('"')
            elif line.startswith("!Sample_geo_accession"):
                samples = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                for s, v in zip(samples, vals):
                    if ":" in v:
                        k, val = v.split(":", 1)
                        chars[k.strip().lower()][s] = val.strip()
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
                next(fh)                                   # the header row of sample ids
            elif line.startswith("!series_matrix_table_end"):
                in_table = False
            elif in_table:
                n_probes += 1
    return samples, chars, n_probes, platform


def present(v):
    return v is not None and v.strip().lower() not in MISSING


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rep = {"artifact_type": "geo_bladder_gate0_joint_count", "phase_of_origin": "pre_analysis",
           "reportable": False, "series": {}}
    for f in sorted(os.listdir(a.dir)):
        if not f.endswith("_series_matrix.txt.gz"):
            continue
        samples, chars, n_probes, platform = parse(os.path.join(a.dir, f))
        keys = {}
        for k, d in chars.items():
            c = collections.Counter(d.get(s, "") for s in samples)
            keys[k] = {"present": sum(1 for s in samples if present(d.get(s))),
                       "n_levels": len(c), "top_levels": c.most_common(10)}
        cand = {var: sorted(k for k in chars if re.search(rx, k)) for var, rx in PATTERNS.items()}
        rep["series"][f] = {"platform": platform, "samples": len(samples), "probes": n_probes,
                            "characteristics": keys, "candidate_keys_by_variable": cand}
        print("\n=== %s  platform %s  %d samples  %d probes" % (f, platform, len(samples), n_probes))
        for k, v in keys.items():
            print("  %-40s present %4d  levels %4d  top %s"
                  % (k[:40], v["present"], v["n_levels"], v["top_levels"][:4]))
        print("  candidates:", json.dumps(cand))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(rep, open(a.out, "w"), indent=1)
    print("\nwrote", a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
