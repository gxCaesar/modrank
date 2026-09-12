#!/usr/bin/env python3
"""Reduce six public pan-TCGA feature releases to per-slide vectors for the five benchmark cohorts.

Only the slides those cohorts' metadata names are kept, so what comes back to the local machine is
a few megabytes per encoder rather than gigabytes. Nothing here fits, selects or looks at an
outcome; it is a format conversion.
"""
import glob, io, json, os, sys
import numpy as np

import os
# default is the path this ran at on the compute host; override with $BLCA_FM_DIR
SRC = os.environ.get("BLCA_FM_DIR", "/data1/guanxing/bladder_cancer/data/fm-slide-features")
META = os.environ.get("BLCA_META_DIR",
                      "/data1/guanxing/bladder_cancer/SurvPath/datasets_csv/metadata")
OUT = os.path.join(SRC, "extracted")
os.makedirs(OUT, exist_ok=True)

import csv
wanted, case_of = set(), {}
for f in glob.glob(os.path.join(META, "tcga_*.csv")):
    coh = os.path.basename(f)[5:-4]
    for r in csv.DictReader(open(f)):
        stem = r["slide_id"].replace(".svs", "")
        wanted.add(stem)
        case_of[stem] = (coh, r["case_id"])
print("slides wanted across five cohorts:", len(wanted), flush=True)

DROPPED = {}

def save(tag, names, M, extra=None):
    # Some source entries are malformed -- a handful of GigaSSL values come back with dimension 1
    # instead of 512. Keeping only the modal dimension is right; doing it SILENTLY is not, because
    # a loader that skips what it cannot parse and saves what is left looks identical to one that
    # found everything. The count is recorded and printed.
    from collections import Counter
    dims = Counter(len(v) for v in M)
    keep_d = dims.most_common(1)[0][0]
    ok = [i for i, v in enumerate(M) if len(v) == keep_d]
    n_drop = len(M) - len(ok)
    DROPPED[tag] = {"dropped": n_drop, "kept": len(ok), "modal_dim": int(keep_d),
                    "dim_histogram": {int(k): int(v) for k, v in dims.items()}}
    if n_drop:
        print("  !! %s: dropping %d of %d entries whose dimension is not %d -- histogram %s"
              % (tag, n_drop, len(M), keep_d, dict(dims)), flush=True)
    names = [names[i] for i in ok]
    M = [M[i] for i in ok]
    if extra:
        extra = {k: np.asarray(v)[ok] if len(np.asarray(v)) == len(ok) + n_drop else v
                 for k, v in extra.items()}
    d = {"names": np.array(names), "X": np.vstack(M).astype(np.float32)}
    if extra: d.update(extra)
    np.savez_compressed(os.path.join(OUT, tag + ".npz"), **d)
    print("  %-14s %d slides x %d dims" % (tag, len(names), d["X"].shape[1]), flush=True)

for tag in ("phikon", "gigapath", "optimus", "ctranspath"):
    p = os.path.join(SRC, "gigassl_%s.npy" % tag)
    if not os.path.exists(p): continue
    d = np.load(p, allow_pickle=True).item()
    hit = [(k, np.asarray(v, dtype=np.float32).ravel()) for k, v in d.items() if k in wanted]
    print("gigassl_%s: %d of %d source keys matched" % (tag, len(hit), len(d)), flush=True)
    if hit: save("gigassl_" + tag, [k for k, _ in hit], [v for _, v in hit])

import pyarrow.parquet as pq
p = os.path.join(SRC, "provgigapath.parquet")
if os.path.exists(p):
    t = pq.read_table(p)
    cols = t.column_names
    fn = [str(x) for x in t.column("filename").to_pylist()]
    emb = t.column("embedding").to_pylist()
    stage = ([str(x) for x in t.column("ajcc_pathologic_tumor_stage").to_pylist()]
             if "ajcc_pathologic_tumor_stage" in cols else None)
    grade = ([str(x) for x in t.column("histological_grade").to_pylist()]
             if "histological_grade" in cols else None)
    names, M, st, gr = [], [], [], []
    for i, f in enumerate(fn):
        stem = f.replace(".svs", "").replace(".h5", "")
        if stem in wanted:
            names.append(stem); M.append(np.asarray(emb[i], dtype=np.float32).ravel())
            st.append(stage[i] if stage else ""); gr.append(grade[i] if grade else "")
    print("provgigapath: %d of %d matched" % (len(names), len(fn)), flush=True)
    if names:
        save("provgigapath", names, M,
             {"ajcc_stage": np.array(st), "histological_grade": np.array(gr)})

p = os.path.join(SRC, "uni_slide.parquet")
if os.path.exists(p):
    pf = pq.ParquetFile(p)
    names, M = [], []
    for rg in range(pf.num_row_groups):
        t = pf.read_row_group(rg, columns=["file_name", "embedding", "embedding_shape"])
        d = t.to_pydict()
        for f, b, sh in zip(d["file_name"], d["embedding"], d["embedding_shape"]):
            stem = str(f).replace(".svs", "")
            if stem not in wanted: continue
            a = np.frombuffer(b, dtype=np.float32).reshape(tuple(int(x) for x in sh))
            names.append(stem); M.append(a.mean(0))
        if rg % 5 == 0: print("  uni row-group %d/%d, kept %d" % (rg, pf.num_row_groups, len(names)), flush=True)
    print("uni: %d matched" % len(names), flush=True)
    if names: save("uni_meanpool", names, M)
print("DROP REPORT:", json.dumps(DROPPED, indent=1))
open(os.path.join(OUT, "drop_report.json"), "w").write(json.dumps(DROPPED, indent=1))
print("DONE")
