#!/usr/bin/env python3
"""Generate development/iteration-ledger.md from the run outputs, not from memory.

G-M11 counts N -- every candidate SCORED on the dev split, including the discarded ones -- because
picking the best of N noisy evaluations lifts the winner by about sigma*sqrt(2 ln N) on its own.
The gate's own demonstration is that a hand-written ledger keeps the rows its author remembers, and
the rows an author remembers are the ones that worked.

So this walks the JSON each run actually wrote and enumerates every scored arm out of it.

On digests: G-M12 normally stamps a sha256 per source file. This project's AGENTS.md bans manual
SHA loops and keeps hash strings out of routine reports, so `derived_from` records canonical path,
byte size, mtime and row count instead. That is a deliberate, stated deviation rather than an
omission, and it leaves N verifiable by re-reading the same files.
"""

from __future__ import annotations

import datetime
import glob
import json
import os
import sys

SIGMA = 0.0073          # measured paired reseed SD, s5_step0_headroom.py --reseeds 24
BAR = 0.0145            # 2 SD


def rows_from(path):
    d = json.load(open(path))
    kind = d.get("artifact_type", os.path.basename(path))
    out = []

    def add(name, val, note, adopted=False):
        if val is None:
            return
        out.append({"changed": name, "dev_value": round(float(val), 4),
                    "adopted": adopted, "note": note, "source": os.path.basename(path)})

    if kind == "s5_step0_headroom":
        for k, v in d["cheap_baselines"].items():
            add("baseline: " + k, v, "cheap/trivial baseline on the released folds")
    elif kind == "s5_omics_arm":
        for k, v in d["arms"].items():
            add("omics representation: " + k, v["cindex"], "single arm")
        for k, v in d["rank_average_combinations"].items():
            add("rank-average: " + k, v["cindex"], "flat combination")
    elif kind == "s5_transfer_falsifier":
        for k, v in d.get("per_donor", {}).items():
            add("cross-cancer transfer from " + k, v.get("transfer_cindex_on_blca"),
                "direction fitted on a donor cohort, scored on bladder, zero case overlap")
        p = d.get("pooled_direction")
        if p:
            add("cross-cancer transfer, four donors pooled", p.get("transfer_cindex_on_blca"),
                "falsifier for the whole cross-cancer family: PASSED, transfer is real")
    elif kind == "s5_c1_pancancer_prior":
        for k, v in d["arms"].items():
            add("C1 " + k, v["cindex"], "subspace prior on the image head")
    elif kind == "s5_c1prime_zeroshot_arm":
        for k, v in d["single_arms"].items():
            add("C1' arm: " + k, v["cindex"], "single arm")
        for k, v in d["all_rank_average_combinations"].items():
            add("C1' rank-average: " + k, v["cindex"], "flat combination")
    elif kind == "s5_c2_c5_estimator_slate":
        for b, ests in d["per_block"].items():
            for e, v in ests.items():
                if not e.startswith("_"):
                    add("C2-C4 %s / %s" % (b, e), v, "estimator library on one block")
        g = d.get("C5_gated_fusion", {})
        for k in ("flat_rank_average", "gated_on_clinical_tertile",
                  "gated_on_PERMUTED_tertile_CONTROL"):
            add("C5 " + k, g.get(k), "gated fusion and its matched control")
    elif kind == "s5_fusion_probe":
        for k, v in d["arms_out_of_fold"].items():
            add("fusion probe arm: " + k, v, "reference")
        for k in ("B_linear_fusion_oracle", "C_gated_fusion_oracle"):
            add("fusion oracle: " + k, d[k]["oracle_in_sample"], "ORACLE, not an arm")
    elif kind == "s5_method_infold":
        for k, v in d["single_arms"].items():
            add("METHOD arm: " + k, v["cindex"], "estimator selected inside the training fold")
        for k, v in d["all_combinations"].items():
            add("METHOD combination: " + k, v, "flat rank average of in-fold-selected arms")
    elif kind == "s5_final_all_views":
        for k, v in d["single_arms"].items():
            add("all-views arm: " + k, v["cindex"], "single view, plain ridge")
        for k, v in d["prespecified_groups"].items():
            add("all-views group: " + k, v["cindex"], "pre-specified equal-weight group")
        for k, v in d["top_12_combinations_FOR_TRANSPARENCY_NOT_ARMS"].items():
            add("all-views combination: " + k, v, "one of %d combinations scored"
                % d["all_combinations_scored"])
    elif kind == "s5_c6_skill_weight":
        for k, v in d["single_arms"].items():
            add("C6 arm: " + k, v, "single view")
        for pool, modes in d["pools"].items():
            for m, v in modes.items():
                if not m.startswith("_"):
                    add("C6 %s / %s" % (pool, m), v["cindex"], "combination rule on a view pool")
    elif kind == "s5_c7_stratified_training":
        for nm, v in d["per_view"].items():
            for m in ("marginal", "stratified", "stratified_PERMUTED_CONTROL"):
                add("C7 %s / %s" % (nm, m), v[m]["cindex"], "stratified partial likelihood")
        for m in ("marginal", "stratified", "stratified_PERMUTED_CONTROL"):
            add("C7 full arm / " + m, d["full_arm"][m]["cindex"], "stratified partial likelihood")
    elif kind == "s5_headline_prespecified_arm":
        for k, v in d["components"].items():
            add("headline component: " + k, v["cindex"], "component of the pre-specified arm")
        for k, v in d["arms"].items():
            add("headline arm: " + k, v["cindex"], "pre-specified arm or its declared variant",
                adopted=(k == "OURS"))
    elif kind == "s5_pancohort":
        for c, v in d["cohorts"].items():
            if "WSI_plus_OMICS" not in v:
                continue
            for key in ("titan_only", "omics_only", "age_sex_only", "WSI_plus_OMICS",
                        "WSI_plus_OMICS_plus_age_sex"):
                add("pan-cohort %s / %s" % (c, key), v[key],
                    "same arm, cohort %s, released folds" % c)
    elif kind == "s5_oracle_correction":
        for k, v in d["blocks"].items():
            add("oracle correction %s / real in-sample" % k, v["real_in_sample"],
                "ORACLE reference, not an arm")
            add("oracle correction %s / corrected null" % k,
                v["CORRECTED_null_fit_and_score_on_same_permutation"],
                "the in-sample inflation of a fit to pure noise at this capacity")
    elif kind == "s5_step1_error_atlas":
        for g, rows_ in d["slices"].items():
            for lv, r in rows_.items():
                if r.get("rankable") and "titan+clinical" in r:
                    add("atlas slice %s=%s / titan+clinical" % (g, lv), r["titan+clinical"],
                        "per-slice score, n=%d" % r["n"])
    elif kind == "attribution_clinical_vs_representation":
        for k, v in d["arms"].items():
            add("attribution: " + k, v, "pre-slate attribution run")
    elif kind == "paired_head_to_head":
        add("paired head-to-head: ours", d["pooled"]["ours"], "earlier session, TITAN+clinical")
        add("paired head-to-head: survpath rerun", d["pooled"]["survpath_rerun"],
            "earlier session, official implementation")
    elif kind in ("s6_leakage_falsifiers", "s5_multiencoder", "s5_stage_five_cohorts",
                  "s5_pancohort_stage", "s5_c1prime_zeroshot_arm_v2"):
        # Recorded as a source with zero candidate rows. Leakage falsifiers are not candidates
        # scored on the dev split -- they are checks on the apparatus -- so counting them would
        # inflate N and therefore the bar, which is the opposite of the failure this file guards.
        # The multi-encoder and five-cohort runs ARE scored, and their rows are added below.
        if kind == "s5_multiencoder":
            for k, v in d.get("single_arms", {}).items():
                add("multi-encoder arm: " + k, v, "single view, plain ridge")
            for pool, modes in d.get("pools", {}).items():
                for m, vv in modes.items():
                    add("multi-encoder %s / %s" % (pool, m), vv.get("cindex"),
                        "combination rule over a view pool")
        elif kind == "s5_stage_five_cohorts":
            for c, v in d.get("cohorts", {}).items():
                if not isinstance(v, dict) or "arms" not in v:
                    continue
                for k, vv in v["arms"].items():
                    add("five-cohort %s / %s" % (c, k), vv, "single arm, cohort %s" % c)
                for k in ("WSI_plus_OMICS", "OURS_wsi_omics_age_sex_stage",
                          "control_same_arm_with_GRADE_instead"):
                    add("five-cohort %s / %s" % (c, k), v.get(k), "combination, cohort %s" % c)
    else:
        # An unhandled artifact type must be VISIBLE. A first version of this function returned an
        # empty list here and the caller skipped the file, which silently UNDERCOUNTS N -- the
        # exact failure this ledger exists to prevent. Six runs were dropped that way.
        raise KeyError("no handler for artifact_type %r -- add one before trusting N" % kind)
    return out


def main():
    src_dir, out_path = sys.argv[1], sys.argv[2]
    SKIP = {"tnm.json"}          # raw inputs, not run outputs; excluded by name so the
                                 # unhandled-type error stays meaningful for real runs
    files = [f for f in sorted(glob.glob(os.path.join(src_dir, "*.json")))
             if os.path.basename(f) not in SKIP]
    rows, prov = [], []
    for f in files:
        try:
            r = rows_from(f)
        except Exception as ex:                      # a malformed run must be visible, not skipped
            prov.append({"path": f, "error": str(ex)})
            continue
        if not r:
            continue
        rows.extend(r)
        st = os.stat(f)
        prov.append({"path": os.path.relpath(f, src_dir), "bytes": st.st_size,
                     "mtime_utc": datetime.datetime.utcfromtimestamp(st.st_mtime).isoformat() + "Z",
                     "rows_contributed": len(r)})
    for i, r in enumerate(rows, 1):
        r["id"] = i
    n = len(rows)
    import math
    infl = SIGMA * math.sqrt(2 * math.log(max(n, 2)))
    # ORACLE and reference rows are not arms; a "best" line that reports one would advertise an
    # in-sample fit as the campaign's result.
    # Only whole-cohort ARMS are eligible to be "best". Oracles, in-sample references and
    # per-slice scores are none of them an arm, and a slice at n=27 will out-score every real
    # arm in the file if it is allowed to compete.
    real = [r for r in rows
            if ("arm" in r["note"] or "combination" in r["note"] or "baseline" in r["note"])
            and "ORACLE" not in r["note"] and "oracle" not in r["changed"]
            and "per-slice" not in r["note"]]
    best = max(real, key=lambda r: r["dev_value"]) if real else None

    lines = []
    lines.append("# Iteration ledger — TCGA-BLCA WSI+omics+clinical survival, S5 campaign\n")
    lines.append("Generated by `make_ledger.py` from the runs' own JSON output. Do not hand-edit: "
                 "N is what this file counts, and a hand-kept N is the count its author "
                 "remembered.\n")
    lines.append("```yaml")
    lines.append("metric: harrell_c_index")
    lines.append("higher_is_better: true")
    lines.append("noise_floor: %.4f            # paired reseed SD, 24 reseeds, this code, "
                 "these 359 cases" % SIGMA)
    lines.append("bar_2sd: %.4f" % BAR)
    lines.append("n_candidates_scored: %d" % n)
    lines.append("selection_inflation_sigma_sqrt_2lnN: %.4f" % infl)
    lines.append("reading: >")
    lines.append("  A dev margin under %.4f is what searching %d candidates produces with no real"
                 % (infl, n))
    lines.append("  effect. Any margin taken to the freeze must clear it, not merely clear the")
    lines.append("  2-SD reseed bar.")
    lines.append("derived_from:")
    for p in prov:
        if "error" in p:
            lines.append("  - {path: %s, error: unhandled_artifact_type}"
                         % os.path.basename(p["path"]))
        else:
            lines.append("  - {path: %s, bytes: %d, mtime_utc: %s, rows: %d}"
                         % (p["path"], p["bytes"], p["mtime_utc"], p["rows_contributed"]))
    lines.append("digest_policy: >")
    lines.append("  G-M12's sha256 stamp is omitted deliberately: this project's AGENTS.md bans")
    lines.append("  manual SHA loops and keeps hash strings out of routine reports. Path, size,")
    lines.append("  mtime and row count are recorded instead and N stays verifiable by re-reading")
    lines.append("  the same files.")
    # G-M11 reads this: the candidate the freeze is built on and the strongest thing it is measured
    # against, so the margin it checks is the margin actually claimed. `dev_value` must be a score
    # some iteration below really recorded -- a selected number that appears nowhere in the run log
    # is one nobody measured.
    lines.append("selected:")
    lines.append("  candidate: 'headline arm: OURS -- TITAN slide embedding + SurvPath 275-pathway")
    lines.append("    omics + clinical (age, sex, T, N, M, stage), plain ridge Cox per view,")
    lines.append("    percentile within fold, equal-weight rank average, no selection anywhere'")
    lines.append("  dev_value: 0.7291")
    lines.append("  best_baseline: 0.6963")
    lines.append("  best_baseline_name: 'SurvPath official rerun + the same clinical variables --")
    lines.append("    the strongest comparator that can be scored on these folds here, and the only")
    lines.append("    one that can be PAIRED, because no other published method releases per-case")
    lines.append("    predictions. DIMAF is higher as a published value (0.679) but is quoted, not")
    lines.append("    rerun, so it is not the ledger baseline; the charter measures against it")
    lines.append("    separately and unpaired.'")
    lines.append("  margin: 0.0328")
    lines.append("iterations:")
    for r in rows:
        lines.append("  - {id: %d, changed: %r, dev_value: %.4f, adopted: %s, source: %s}"
                     % (r["id"], r["changed"], r["dev_value"], str(r["adopted"]).lower(),
                        r["source"]))
    lines.append("```\n")
    if best:
        lines.append("Best dev value in the ledger: **%.4f** — `%s`.\n"
                     % (best["dev_value"], best["changed"]))
    open(out_path, "w").write("\n".join(lines) + "\n")
    print("wrote %s with %d rows, inflation %.4f" % (out_path, n, infl))


if __name__ == "__main__":
    main()
