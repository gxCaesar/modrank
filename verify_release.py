#!/usr/bin/env python3
"""Recompute the manuscript's load-bearing numbers from the result files in this bundle.

This is the first thing to run, and it needs nothing installed beyond Python 3.9. It does not read
the manuscript: every check below RECOMPUTES a quantity from the file that is supposed to contain
it, and fails if the file disagrees with itself. A results file that merely records a number is a
record; one whose derived fields can be regenerated from its own primitives is evidence.

Every one of these caught a real error while the paper was being written. The published-gap check
found the manuscript describing the range as 0.002-0.012 when two of the eight gaps are 0.021 and
0.029; the component check found a table listing eight rows against a text that said seven; the
fusion check found a gain quoted against the wrong baseline.
"""

from __future__ import annotations

import json
import math
import os
import re
import sys

R = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "experiments", "20260817-blca-confirm", "results")
F = os.path.join(os.path.dirname(os.path.abspath(__file__)), "paper", "figures")
load = lambda d, n: json.load(open(os.path.join(d, n)))

ok, fail = [], []


def check(name, cond, detail=""):
    (ok if cond else fail).append((name, detail))


# --- 1. the published table's gaps must regenerate from its own entries
PUB = load(R, "published-benchmark-table.json")
v = [e["cindex"] for e in PUB["entries"]]
gaps = [round(a - b, 4) for a, b in zip(v[:-1], v[1:])]
check("published gaps regenerate", gaps == PUB["consecutive_gaps"],
      "%s vs %s" % (gaps, PUB["consecutive_gaps"]))
check("published gap range regenerates",
      [min(gaps), max(gaps)] == PUB["consecutive_gap_range"],
      "range %.3f-%.3f" % (min(gaps), max(gaps)))
check("nine published entries", len(PUB["entries"]) == 9, "%d" % len(PUB["entries"]))

# --- 2. the capacity sweep's excess is in_sample minus the matched null, every row
REG = load(R, "regime-evidence.json")
bad = [a for a in REG["capacity_sweep"]["arms"]
       if abs((a["in_sample"] - a["corrected_null"]) - a["excess"]) > 5e-4]
check("capacity sweep excess regenerates", not bad, "%d rows disagree" % len(bad))
t16 = [a for a in REG["capacity_sweep"]["arms"] if a["arm"] == "TITAN" and a["d"] == 16][0]
t8 = [a for a in REG["capacity_sweep"]["arms"] if a["arm"] == "TITAN" and a["d"] == 8][0]
check("the excess SHRINKS from d=8 to d=16", t16["excess"] < t8["excess"],
      "%+.4f at d=16 against %+.4f at d=8" % (t16["excess"], t8["excess"]))

# --- 3. the fusion gains, against BOTH baselines, because quoting one without saying which was
#        an ambiguity a self-consistency check caught
FU = REG["fusion"]
check("gated oracle gain over the linear oracle",
      abs(FU["gated_fusion_oracle"]["value"] - FU["linear_fusion_oracle"]["value"]
          - FU["gated_fusion_oracle"]["gain_over_the_LINEAR_oracle"]) < 1e-3)
check("gated oracle gain over the deployed rule",
      abs(FU["gated_fusion_oracle"]["value"] - FU["deployed_equal_weight"]
          - FU["gated_fusion_oracle"]["gain_over_the_deployed_equal_weight"]) < 1e-3)
check("the weight sweep peaks at equal weight",
      abs(float(FU["best_fixed_weight_on_a_21_point_sweep"]["at_w"]) - 0.50) < 1e-9)

# --- 4. seven components, and not one of them cleared the bar
CS = load(R, "component-slate.json")
check("seven pre-registered components", CS["n_formal_components"] == 7,
      "%d" % CS["n_formal_components"])
cleared = [c["id"] for c in CS["components"] if c["observed"] > CS["bar"]]
check("none cleared the bar", not cleared, "cleared: %s" % cleared)
voided = [c for c in CS["components"]
          if c["status"] == "void" and c["control"] is not None
          and c["observed"] > 0 and c["control"] > c["observed"]]
check("at least one component is void on its own control", bool(voided),
      "%s" % [c["id"] for c in voided])

# --- 5. the six-comparison family's deltas regenerate from the canonical arm table
PAR = load(R, "pibd-parity.json")
A = PAR["canonical_arms_seed_averaged_scores"]
FAM = PAR["holm_family_of_six"]["comparisons"]
MAP = {"ours vs slide alone": "slide_alone", "ours vs omics alone": "omics_alone",
       "ours vs clinical alone": "clinical_alone",
       "ours vs ours-without-clinical": "ours_without_clinical",
       "ours vs SurvPath+clinical (seed-matched)": "survpath_plus_clinical",
       "ours vs PIBD+clinical (best-val checkpoint)": "pibd_plus_clinical_best_val"}
check("the family is six comparisons", len(FAM) == 6, "%d" % len(FAM))
off = [(k, round(A["ours"] - A[MAP[k]], 4), FAM[k]["delta"]) for k in FAM
       if abs((A["ours"] - A[MAP[k]]) - FAM[k]["delta"]) > 1.5e-3]
check("every delta regenerates from the arm table", not off, "%s" % off)
surv = [k for k, x in FAM.items() if x["survives_at_0.05"]]
check("exactly two survive Holm", len(surv) == 2, "%s" % surv)
check("NEITHER input-parity comparison survives",
      not any("SurvPath" in k or "PIBD" in k for k in surv), "%s" % surv)

# --- 6. checkpoint selection is the difference between the two checkpoints, not an assertion
P = PAR["pibd"]
check("checkpoint effect regenerates",
      abs((P["best_validation_epoch"]["fold_mean"] - P["final_epoch"]["fold_mean"])
          - PAR["published"]["checkpoint_selection_is_worth"]) < 1e-3,
      "%.4f" % PAR["published"]["checkpoint_selection_is_worth"])

# --- 7. the power curve passes through its own reported detectable margin
FSD = json.load(open(os.path.join(F, "figure-source-data.json")))
PW = FSD["F3_power"]
z = 1.959963984540054
power = lambda d: (0.5 * math.erfc(-(d / PW["se"] - z) / math.sqrt(2))
                   + 0.5 * math.erfc(-(-d / PW["se"] - z) / math.sqrt(2)))
check("80%% power lands at the reported margin",
      abs(power(PW["detectable_at_80pc"]) - 0.80) < 0.01,
      "power(%.4f) = %.3f" % (PW["detectable_at_80pc"], power(PW["detectable_at_80pc"])))
check("the published gaps sit far below detectability",
      power(max(PUB["consecutive_gaps"])) < 0.30,
      "largest gap %.3f -> power %.3f" % (max(PUB["consecutive_gaps"]),
                                          power(max(PUB["consecutive_gaps"]))))

# --- 8. the risk ladder survives greyscale
L = FSD["F4_risk_group_colours"]
lum = L["luminance"]
check("risk colours run dark to pale, monotonically",
      all(b > a for a, b in zip(lum, lum[1:])), "%s" % lum)
check("minimum adjacent luminance gap clears 0.05",
      min(b - a for a, b in zip(lum, lum[1:])) >= 0.05,
      "%.4f" % min(b - a for a, b in zip(lum, lum[1:])))

# --- 9. the parameter ratios
PC = load(R, "parameter-counts.json")
for k in ("survpath", "pibd"):
    check("%s parameter ratio regenerates" % k,
          abs(round(PC[k]["parameters"] / PC["ours"]["parameters"]) - PC[k]["ratio_to_ours"]) <= 1,
          "%d" % PC[k]["ratio_to_ours"])

# --- 10. THE HEADLINE, which round 2's artifact executor correctly observed this script never
# touched. "22 checks passed" was an internal-consistency result that said nothing about the primary.
# Three numbers exist and the distinction is the point: the frozen confirmatory run produced 0.7260;
# amendment A1 rebuilt the clinical block from DIMAF's released split files, which moved the primary to
# 0.7212; amendment A2 corrected the Cox risk sets and tied ranks and regenerated every later result,
# which moved it to 0.7214. The first and last are checked against the files that record them, the
# middle one against the amendment's before-after table.
CONF = load(R, "confirmatory.json")
check("the frozen confirmatory run recorded 0.7260",
      abs(CONF["primary"]["value"] - 0.7260) < 5e-5, "%.4f" % CONF["primary"]["value"])
A1 = load(R, "amendment-A1-clinical-provenance.json")
REPORTED = 0.7214
check("the per-seed values average to the reported %.4f" % REPORTED,
      abs(sum(r["OURS"] for r in A1["per_seed"]) / len(A1["per_seed"]) - REPORTED) < 5e-5,
      "%.4f over %d seeds" % (sum(r["OURS"] for r in A1["per_seed"]) / len(A1["per_seed"]),
                              len(A1["per_seed"])))
check("the manuscript's headline is the AMENDED value, not the frozen one",
      abs(A1["primary"]["value"] - REPORTED) < 5e-5 and CONF["primary"]["value"] != A1["primary"]["value"],
      "frozen %.4f -> amended %.4f" % (CONF["primary"]["value"], A1["primary"]["value"]))

# --- AMENDMENT A2: the before-after table records the pre-A2 values beside the regenerated ones
BAF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "experiments", "20260928-amendment-a2",
                   "results", "before-after.json")
if os.path.exists(BAF):
    BAR = {r["quantity"]: r for r in json.load(open(BAF))["rows"]}
    _p = BAR["ModRank, mean over five seeds"]
    check("amendment A2 moved the primary from A1's 0.7212 to the reported value",
          abs(_p["before"] - 0.7212) < 5e-5 and abs(_p["after"] - A1["primary"]["value"]) < 5e-5,
          "%.4f -> %.4f" % (_p["before"], _p["after"]))

DUMP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "experiments",
                    "20260818-reporting-dump", "results", "reporting-dump.json")
if os.path.exists(DUMP):
    D = json.load(open(DUMP))
    check("the reporting dump re-executed the arm and reproduced %.4f" % REPORTED,
          (D.get("frozen_primary_matched") or D.get("amendment") == "A2")
          and abs(D["primary"] - REPORTED) < 5e-5,
          "%.4f from %d cases" % (D["primary"], len(D["cases"])))
    check("its per-seed values equal amendment A1's",
          [round(x, 4) for x in D["per_seed_OURS"]] == [r["OURS"] for r in A1["per_seed"]],
          "%s" % D["per_seed_OURS"])

# --- ROUND 2's G1, the finding that changed the paper's verdict, and the one this file could not
# see. G4 was raised against exactly this shape: 22 checks passed while the verifier had never
# touched 0.7212. It touches 0.7212 now, and until this block it had never touched 0.7716 -- the
# bar the primary FAILS. A verifier that only recomputes the numbers its paper passes has the same
# defect wearing the other sign.
HERE = os.path.dirname(os.path.abspath(__file__))
NULL = os.path.join(HERE, "experiments", "20260818-selection-null", "results", "selection-null.json")
if os.path.exists(NULL):
    S = json.load(open(NULL))
    A = S["A_selection_null"]
    INC = S["B_encoder_swap"]["incumbent_published"]

    # The inflation term is sigma*sqrt(2 ln N) in every construction; only sigma and N differ. So
    # regenerating each from its own stated inputs shows the arithmetic was never the problem and
    # isolates the defect to the ESTIMAND sigma measures.
    #
    # The PRE-REGISTERED figures come from the frozen protocol and from nowhere else. This file
    # previously took them from the block below named "as_preregistered", which carries the
    # ledger's FINAL count of 269 rather than the 211 the protocol froze. Both are real; only one
    # was pre-registered, and reading the paper's own account of its pre-registration back out of a
    # post-hoc results file is how the two came to be confused.
    PROTO = os.path.join(HERE, "development", "benchmark-protocol.json")
    if os.path.exists(PROTO):
        pr = json.load(open(PROTO))
        unc = pr["metrics"]["uncertainty"]
        n_frozen = int(re.search(r"N\s*=\s*(\d+)", unc["selection_inflation_basis"]).group(1))
        sig, infl = unc["paired_split_reseed_sd"], unc["selection_inflation"]
        tgt = pr["amendments"][0]["effect_on_the_frozen_primary"]["target"]
        got = sig * math.sqrt(2 * math.log(n_frozen))
        check("the FROZEN selection inflation regenerates from the protocol",
              abs(got - infl) < 5e-4 and n_frozen == 211,
              "%.4f x sqrt(2 ln %d) = %.4f" % (sig, n_frozen, got))
        check("the frozen target is the incumbent plus that inflation",
              abs(round(INC + infl, 4) - tgt) < 5e-5,
              "%.3f + %.4f = %.4f" % (INC, infl, tgt))
        check("the results file's \"as_preregistered\" block is the LEDGER-FINAL count, not this",
              A["bars"]["as_preregistered"]["N"] == 269 and n_frozen == 211,
              "protocol N=%d, that block N=%d" % (n_frozen, A["bars"]["as_preregistered"]["N"]))

    b = A["bars"]["gaussian_with_the_null_sigma"]
    got = b["sigma"] * math.sqrt(2 * math.log(b["N"]))
    check("gaussian with the null sigma inflation regenerates",
          abs(got - b["inflation"]) < 5e-4 and abs(b["inflation"] - 0.1299) < 5e-5,
          "%.4f x sqrt(2 ln %d) = %.4f" % (b["sigma"], b["N"], got))

    check("the null's sigma is 5x the sigma the rule used",
          round(A["null_single_candidate"]["sd"] / A["bars"]["as_preregistered"]["sigma"]) == 5,
          "%.4f measured against %.4f assumed"
          % (A["null_single_candidate"]["sd"], A["bars"]["as_preregistered"]["sigma"]))

    # The empirical bar assumes neither normality nor independent looks: the null max is centred on
    # chance, so its excess over 0.5 is what a search of this family buys.
    emp = INC + (A["null_max_over_the_family"]["q95"] - 0.5)
    check("the empirical-max bar regenerates from the null's own q95",
          abs(emp - A["bars"]["empirical_max_q95"]["bar"]) < 5e-4 and abs(emp - 0.7708) < 5e-4,
          "%.4f + (%.4f - 0.5) = %.4f" % (INC, A["null_max_over_the_family"]["q95"], emp))

    P = A1["primary"]["value"]
    check("the reported primary FAILS both correctly-specified bars",
          P < A["bars"]["empirical_max_q95"]["bar"]
          and P < A["bars"]["gaussian_with_the_null_sigma"]["bar"],
          "%.4f < %.4f and < %.4f" % (P, A["bars"]["empirical_max_q95"]["bar"],
                                      A["bars"]["gaussian_with_the_null_sigma"]["bar"]))
    check("and CLEARS the frozen bar, which is the contradiction reported",
          P > tgt and P > A["bars"]["as_preregistered"]["bar"],
          "%.4f > %.4f frozen and > %.4f at the ledger-final count"
          % (P, tgt, A["bars"]["as_preregistered"]["bar"]))

    # G8, the encoder confound. Hold the method fixed, swap the slide block through all seven.
    E = S["B_encoder_swap"]
    slide = [e["slide_alone"] for e in E["encoders"].values()]
    full = [e["full_method"] for e in E["encoders"].values()]
    check("the encoder sweep's two spans regenerate",
          abs((max(slide) - min(slide)) - E["slide_arm_span"]) < 5e-4
          and abs((max(full) - min(full)) - E["full_method_span"]) < 5e-4,
          "slide %.4f, full method %.4f" % (max(slide) - min(slide), max(full) - min(full)))
    beat = sum(1 for x in full if x > INC)
    check("five of the seven encoders still beat the incumbent",
          beat == E["encoders_whose_full_method_still_beats_the_incumbent"] == 5,
          "%d of %d above %.3f" % (beat, len(full), INC))

# --- and the rank is RE-EXECUTED rather than read back. Its emitter reads one markdown file and
# needs no data, so the strongest available check is to run the shipped script against the shipped
# ledger in a scratch directory and require the output to equal the shipped artifact.
S15 = os.path.join(HERE, "analysis", "s15_selection_rank.py")
LEDGER = os.path.join(HERE, "development", "iteration-ledger.md")
RANK = os.path.join(HERE, "experiments", "20260818-selection-null", "results", "selection-rank.json")
if all(os.path.exists(p) for p in (S15, LEDGER, RANK)):
    import subprocess
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "selection-rank.json")
        r = subprocess.run([sys.executable, S15, "--ledger", LEDGER, "--out", out],
                           capture_output=True)
        got = json.load(open(out)) if os.path.exists(out) else {}
    want = json.load(open(RANK))
    check("the frozen arm's rank re-executes from the shipped ledger",
          r.returncode == 0 and got == want,
          "rank %s of %s comparable, %s strictly above"
          % (got.get("rank_of_the_frozen_arm"), got.get("full_cohort_out_of_fold_candidates"),
             got.get("candidates_strictly_above")))
    check("the rank equals the strictly-higher count plus its ties",
          want["candidates_strictly_above"] + want["ledger_rows_at_the_frozen_value"]
          == want["rank_of_the_frozen_arm"],
          "%d + %d = %d" % (want["candidates_strictly_above"],
                            want["ledger_rows_at_the_frozen_value"],
                            want["rank_of_the_frozen_arm"]))
    check("the five strata account for every scored candidate",
          sum(want["strata"].values()) == want["candidates_scored"],
          "%s = %d" % (" + ".join(str(v) for v in want["strata"].values()),
                       want["candidates_scored"]))


# --- 12. every artifact the manuscript's Availability paragraph PROMISES is actually here.
#         G50 checks that identifiers resolve and this project has no gate for the other half: a
#         promise can name a real thing and the bundle can still not carry it. It shipped that way
#         once. The paragraph promises six artifacts and an environment record; the run LOGS were
#         the missing one, which is also the only place the frozen decision rule appears as it
#         executed rather than as the manuscript describes it.
HERE = os.path.dirname(os.path.abspath(__file__))
PROMISED = [
    ("the frozen protocol", "development/benchmark-protocol.json"),
    ("the split manifest with its hash", "development/split-manifest.json"),
    ("the iteration ledger covering all 269 scored candidates", "development/iteration-ledger.md"),
    ("the error atlas", "development/error-atlas.json"),
    ("the component pre-registrations with their controls", "development/contribution-design.md"),
    ("the confirmatory run's environment record", "experiments/20260817-blca-confirm/env.txt"),
    ("the confirmatory run's logs", "experiments/20260817-blca-confirm/logs/confirm.log"),
]
missing = [t for t, rel in PROMISED if not os.path.isfile(os.path.join(HERE, rel))]
check("the Availability paragraph's promises are all present", not missing,
      "%d of %d present" % (len(PROMISED) - len(missing), len(PROMISED))
      + ("" if not missing else "; MISSING: " + "; ".join(missing)))

# The frozen rule is quoted in the manuscript. The log is where it was printed by the run itself,
# so the two are checked against each other rather than the log being shipped as decoration.
_logp = os.path.join(HERE, "experiments/20260817-blca-confirm/logs/confirm.log")
_log = open(_logp).read() if os.path.isfile(_logp) else ""
check("the run log records the frozen rule at the protocol's own selection term",
      "0.679 + 0.0239 = 0.7029" in _log,
      "the executed rule, printed by the confirmatory run before any of this was written up")

print("recomputed from the bundle's own result files:\n")
for n, d in ok:
    print("  ok    %-48s %s" % (n, d))
for n, d in fail:
    print("  FAIL  %-48s %s" % (n, d))
print("\n%d checks, %d passed, %d failed" % (len(ok) + len(fail), len(ok), len(fail)))
sys.exit(1 if fail else 0)
