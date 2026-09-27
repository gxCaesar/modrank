#!/usr/bin/env python3
"""Run amendment A2's producers in one of two modes.

legacy: BLCA_A2=0, upstream inputs are the committed repository files, outputs go to legacy2/<job>/,
        and each output is deep-diffed against its committed file.
a2:     BLCA_A2=1, upstream inputs are taken from the staging tree when a staged copy exists,
        outputs go to the staging tree at the committed file's repository path.

usage: run_jobs.py legacy|a2 job [job ...]      (or "all")
"""
import concurrent.futures as cf
import json
import os
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
S = os.environ.get("A2_SCRATCH", os.path.dirname(os.path.abspath(__file__)))
STAGE = os.path.join(S, "stage")
# the interpreter of results committed on or before 2026-09-11 (Python 3.8.20, NumPy 1.24.4, SciPy 1.10.1)
ML = os.environ.get("A2_LEGACY_PYTHON", "python3.8")
# every later result and every A2 run (Python 3.9.6, NumPy 2.0.2, SciPy 1.13.1)
SYS = os.environ.get("A2_PYTHON", "python3")
D = "data/frozen-inputs"
AUX = os.path.join(S, "legacy", "aux")
FI = "experiments/20260927-field-inflation"
CONF = "experiments/20260817-blca-confirm/results"
POST = "experiments/20260911-blca-posthoc/results"
BLCA = ["--root", D, "--titan", D + "/TCGA_TITAN_features.pkl", "--sp-dir", D + "/survpath_preds",
        "--omics-dir", D + "/omics", "--dimaf-dir", D + "/dimaf_splits"]
FIVE = ["--titan", D + "/TCGA_TITAN_features.pkl", "--meta-dir", D + "/meta", "--clin-dir", D + "/clin",
        "--pan-dir", D + "/pan", "--sig", D + "/omics/combine_signatures.csv"]
ARCH = ["coattn", "abmil_wsi_pathways", "transmil_wsi_pathways", "deepmisl_wsi_pathways"]

# name: (legacy interpreter, script, argv with {UP:path} upstream inputs, committed output path)
JOBS = {
    "s8": (ML, "analysis/s8_biology.py", BLCA, CONF + "/biology.json"),
    "s9": (ML, "analysis/s9_pibd_parity.py", BLCA + ["--seed-dir", AUX + "/seed-dir", "--pibd-dir", AUX + "/pibd-dir"],
           CONF + "/pibd-parity.json"),
    "s11": (ML, "analysis/s11_cold_panel_round1.py", BLCA, CONF + "/cold-panel-round1-measurements.json"),
    "s7interp": (ML, "analysis/s7_interpret_calibrate_cases.py", BLCA, CONF + "/interpret-calibrate-cases.json"),
    "s7gaps": (ML, "analysis/s7_gaps.py", BLCA + ["--seed-dir", AUX + "/seed-dir", "--clin-dir", D + "/clin"],
               CONF + "/reviewer-gaps.json"),
    "s7surv": (ML, "analysis/s7_survival_metrics.py", BLCA, CONF + "/survival-metrics.json"),
    "s7multiseed": (ML, "analysis/s7_survpath_multiseed.py", BLCA + ["--seed-dir", AUX + "/seed-dir"],
                    CONF + "/survpath-multiseed.json"),
    "s10": (ML, "analysis/s10_figure_data.py", BLCA + ["--seed-dir", AUX + "/seed-dir", "--pibd-dir", AUX + "/pibd-dir",
            "--parity", "{UP:" + CONF + "/pibd-parity.json}", "--a1", "{UP:" + CONF + "/amendment-A1-clinical-provenance.json}",
            "--published", CONF + "/published-benchmark-table.json"],
            "paper/figures/figure-source-data.json"),
    "s16": (ML, "analysis/s16_encoder_parity_chief.py",
            ["--protocol", "development/benchmark-protocol.json"] + BLCA +
            ["--chief", D + "/chief_pooled.npz", "--sweep", "{UP:experiments/20260818-selection-null/results/selection-null.json}"],
            "experiments/20260820-encoder-parity/results/encoder-parity.json"),
    "s31": (SYS, "analysis/s31_decomp_emitter.py",
            ["--protocol", "development/benchmark-protocol.json"] + BLCA + ["--frozen", "development/s5-results/decomp.json"],
            "experiments/20260912-decomp-emitter/results/clinical-block-decomposition.json"),
    "stage5": (ML, "analysis/s5_stage_five_cohorts.py", FIVE,
               "experiments/20260911-five-study-repro/results/stage5.json"),
    "s27": (SYS, "analysis/s27_five_cohort_intervals.py",
            FIVE + ["--frozen", "{UP:experiments/20260911-five-study-repro/results/stage5.json}"],
            "experiments/20260912-five-cohort-intervals/results/five-cohort-intervals.json"),
    "s27g": (SYS, "analysis/s27_five_cohort_intervals.py",
             FIVE + ["--slide-npz", D + "/fmfeat/gigassl_gigapath.npz",
                     "--frozen", "{UP:experiments/20260911-five-study-repro/results/stage5.json}"],
             "experiments/20260913-encoder-sensitivity/results/five-cohort-gigassl.json"),
    "s32": (SYS, "analysis/s32_five_cohort_fusion.py",
            FIVE + ["--frozen", "{UP:experiments/20260911-five-study-repro/results/stage5.json}"],
            "experiments/20260913-five-cohort-fusion/results/five-cohort-fusion.json"),
    "s28": (SYS, "analysis/s28_subgroup_fairness.py",
            ["--vectors", "{UP:" + POST + "/percase-canonical-vectors.json}", "--meta", D + "/meta/tcga_blca.csv",
             "--race", D + "/gdc_demographics.json"],
            "experiments/20260912-subgroup-fairness/results/subgroup-fairness.json"),
    "s34": (SYS, "analysis/s34_five_study_extensions.py",
            FIVE + ["--frozen", "{UP:experiments/20260911-five-study-repro/results/stage5.json}",
                    "--frozen-fusion", "{UP:experiments/20260913-five-cohort-fusion/results/five-cohort-fusion.json}"],
            FI + "/analysis-results/five-study-extensions-bd.json"),
    "s33": (SYS, "analysis/s33_field_inflation.py",
            ["--percase", "{UP:" + POST + "/percase-canonical-vectors.jsonl}",
             "--committed", "{UP:" + POST + "/unified-fusion-and-added-value.json}",
             "--runs", FI + "/results", "--models", "coattn", "coattn_motcat"] + ARCH[1:],
            FI + "/analysis-results/field-inflation-a.json"),
    "s35b": (SYS, "analysis/s35_inflation_null.py",
             ["blca", "--percase", "{UP:" + POST + "/percase-canonical-vectors.jsonl}",
              "--committed", "{UP:" + POST + "/unified-fusion-and-added-value.json}",
              "--committed-a", "{UP:" + FI + "/analysis-results/field-inflation-a.json}",
              "--runs", FI + "/results", "--models"] + ARCH,
             FI + "/analysis-results/inflation-null-blca.json"),
    "s35f": (SYS, "analysis/s35_inflation_null.py",
             ["five"] + FIVE[:2] + FIVE[2:] + ["--committed-b", "{UP:" + FI + "/analysis-results/five-study-extensions-bd.json}"],
             FI + "/analysis-results/inflation-null-five.json"),
}
# s35's subcommand comes first and its option names differ from FIVE's for the signature file
JOBS["s35f"] = (SYS, "analysis/s35_inflation_null.py",
                ["five", "--titan", D + "/TCGA_TITAN_features.pkl", "--meta-dir", D + "/meta", "--clin-dir", D + "/clin",
                 "--pan-dir", D + "/pan", "--sig", D + "/omics/combine_signatures.csv",
                 "--committed-b", "{UP:" + FI + "/analysis-results/five-study-extensions-bd.json}"],
                FI + "/analysis-results/inflation-null-five.json")
JOBS.update({
    "s13": (ML, "analysis/s13_selection_null_and_encoder.py",
            ["--protocol", "development/benchmark-protocol.json"] + BLCA + ["--fm-dir", D + "/fmfeat", "--perms", "200"],
            "experiments/20260818-selection-null/results/selection-null.json"),
    # s22 reads unified-fusion-and-added-value.json as a sibling of --out: under a2 that is the staged s19 file
    "s22": (ML, "analysis/s22_resplit_and_site_cv.py", BLCA, POST + "/resplit-and-site-cv.json"),
    "s19b": (ML, "analysis/s19b_modality_atlas_amended.py",
             ["--root", D, "--titan", D + "/TCGA_TITAN_features.pkl", "--sp-dir", D + "/survpath_preds",
              "--dimaf-dir", D + "/dimaf_splits"], POST + "/modality-atlas-amended.json"),
    "s7abl": (ML, "analysis/s7_ablation_and_generalisation.py",
              BLCA + ["--meta-dir", D + "/meta", "--pan-dir", D + "/pan", "--clin-dir", D + "/clin"],
              CONF + "/ablation-generalisation.json"),
    "s29": (SYS, "analysis/s29_five_cohort_floors.py",
            FIVE + ["--frozen", "{UP:experiments/20260911-five-study-repro/results/stage5.json}"],
            "experiments/20260912-five-cohort-floors/results/five-cohort-floors.json"),
})
JOBS["s36s"] = (SYS, "analysis/s36_reference_gap_external.py",
                ["summary", "--gap", "{UP:" + FI + "/analysis-results/reference-gap-external.json}",
                 "--stage5", "{UP:experiments/20260911-five-study-repro/results/stage5.json}",
                 "--intervals", "{UP:experiments/20260912-five-cohort-intervals/results/five-cohort-intervals.json}",
                 "--why", "{UP:" + CONF + "/why-grade-fails.json}"],
                FI + "/analysis-results/reference-gap-summary.json")
# the local jsonl the committed s33/s35 read is the committed table (legacy); a2 uses the staged export
IGNORE = {"runtime_seconds", "runtime", "executed_at_utc", "generated_at", "git_head", "elapsed", "host"}


def resolve(arg, mode):
    if not (arg.startswith("{UP:") and arg.endswith("}")):
        return arg
    rel = arg[4:-1]
    if mode == "a2" and os.path.exists(os.path.join(STAGE, rel)):
        return os.path.join(STAGE, rel)
    if mode == "a2":
        raise SystemExit("A2 run needs the staged upstream %s first" % rel)
    return rel


def diff(a, b, p=""):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k in IGNORE:
                continue
            if k not in a or k not in b:
                out.append("%s/%s missing" % (p, k))
            else:
                out += diff(a[k], b[k], p + "/" + k)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [p + " length %d vs %d" % (len(a), len(b))]
        for i, (u, v) in enumerate(zip(a, b)):
            out += diff(u, v, "%s[%d]" % (p, i))
    elif a != b:
        out.append("%s: %r vs %r" % (p, a, b))
    return out


def run(name, mode):
    py, script, argv, committed = JOBS[name]
    if mode == "legacy":
        out = os.path.join(S, "legacy2", name, os.path.basename(committed))
        interp = py
    else:
        out = os.path.join(STAGE, committed)
        interp = SYS
    os.makedirs(os.path.dirname(out), exist_ok=True)
    args = [resolve(x, mode) for x in argv]
    cmd = [interp, "-W", "ignore", script] + args + ["--out", out]
    env = dict(os.environ, BLCA_A2="0" if mode == "legacy" else "1",
               OMP_NUM_THREADS="3", MKL_NUM_THREADS="3", OPENBLAS_NUM_THREADS="3")
    logd = os.path.join(S, mode + "-logs")
    os.makedirs(logd, exist_ok=True)
    t0 = time.time()
    with open(os.path.join(logd, name + ".log"), "w") as fh:
        fh.write(" ".join(cmd) + "\n")
        fh.flush()
        rc = subprocess.run(cmd, cwd=REPO, env=env, stdout=fh, stderr=subprocess.STDOUT).returncode
    res = {"job": name, "mode": mode, "exit": rc, "seconds": round(time.time() - t0, 1),
           "interpreter": interp, "out": out}
    if mode == "legacy" and rc == 0 and os.path.exists(out):
        d = diff(json.load(open(os.path.join(REPO, committed))), json.load(open(out)))
        res["differences"] = len(d)
        res["first"] = d[:5]
    return res


def main():
    mode, names = sys.argv[1], sys.argv[2:]
    if names == ["all"]:
        names = list(JOBS)
    with cf.ThreadPoolExecutor(max_workers=int(os.environ.get("JOBS", "4"))) as ex:
        for r in ex.map(lambda n: run(n, mode), names):
            print(json.dumps(r), flush=True)
            with open(os.path.join(S, mode + "-results.jsonl"), "a") as fh:
                fh.write(json.dumps(r) + "\n")


if __name__ == "__main__":
    main()
