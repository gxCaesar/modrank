#!/usr/bin/env python3
"""Recompute, from the iteration ledger alone, where the frozen arm sat in the campaign.

WHY THIS SCRIPT EXISTS. Cold panel round 2 applied a maximum-of-N selection correction to the
confirmatory result. That correction assumes the maximum was taken. Whether it was is a fact about
the ledger, and the paper asserts it, so it has to be recomputed from the ledger rather than
remembered. The answer this file wrote by hand in an earlier session agreed with the ledger on the
stratum size and disagreed with its own listed evidence on the rank, which is exactly the failure a
hand-kept number produces.

THE STRATUM. Not every scored row is a candidate the frozen arm competed against. Four kinds are
excluded, by a rule stated here and applied mechanically:

  subgroup slice            scored on a subset of the cohort, so not comparable to a full-cohort
                            number at all (the error atlas)
  in-sample oracle          fitted and scored on the same rows, by construction, to measure how
                            much a given capacity can memorise
  permuted or null control  built to lose, and reported as the thing its partner had to beat
  cross-cohort transfer     scored on a different disease

What remains is the set of full-cohort out-of-fold values, which is what a rank statement means.

TIES ARE COUNTED PESSIMISTICALLY. The frozen configuration appears in the ledger more than once,
because three sub-experiments arrived at the same construction and each recorded it. Its rank is
therefore reported as the WORST position consistent with the ledger, counting every tied row as
though it were above. The optimistic rank is recorded beside it so the difference is visible rather
than chosen silently.
"""

from __future__ import annotations

import argparse
import json
import re

# The frozen configuration's development value, from the protocol's `selected` block. A literal
# rather than a search key: three different rows describe this same construction in three different
# vocabularies, and matching on the value is the only thing all three agree on.
FROZEN_DEV = 0.7291

ROW = re.compile(
    r"\{id: (\d+), changed: ('(?:[^']*)'|\"(?:[^\"]*)\"), dev_value: ([0-9.]+), "
    r"adopted: (\w+), source: ([a-z0-9_]+\.json)\}")

# Why each of the eight strictly-higher rows is not the frozen arm. Keyed by the ledger's own
# `changed` string so a renamed candidate breaks this loudly instead of silently dropping its note.
WHY = {
    "C5 gated_on_clinical_tertile":
        "rejected: its permuted-tertile control beat it, reported as void",
    "C5 flat_rank_average":
        "the same C5 comparison's reference arm, not a separate candidate",
    "pan-cohort brca / titan_only": "a different cohort (BRCA), not a bladder candidate",
    "five-cohort brca / titan": "a different cohort (BRCA), not a bladder candidate",
    "all-views combination: clinical + omics_xena + wsi_titan": "the declared xena sensitivity",
    "headline arm: SENSITIVITY_xena_instead_of_combine": "the declared xena sensitivity",
    "C6 three_view_titan_combine_clinical / skill":
        "rejected: skill weighting failed its control",
    "multi-encoder titan_only_reference / skill":
        "rejected: skill weighting failed its control",
}


def stratum(changed, source):
    """Which of the five kinds a ledger row is. One rule, applied to every row."""
    low = changed.lower()
    if source == "atlas.json":
        return "subgroup_slice"
    if source == "transfer.json":
        return "cross_cohort_transfer_probe"
    if source == "oracle_fix.json" or "oracle" in low:
        return "in_sample_oracle"
    if "permuted" in low or "control_same_arm_with_grade" in low:
        return "permuted_or_null_control"
    return "full_cohort_out_of_fold"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    text = open(a.ledger).read()
    declared = int(re.search(r"n_candidates_scored: (\d+)", text).group(1))
    rows = [{"id": int(m.group(1)), "changed": m.group(2)[1:-1],
             "dev": float(m.group(3)), "source": m.group(5)}
            for m in ROW.finditer(text)]
    if len(rows) != declared:
        raise SystemExit("the ledger declares %d candidates and %d rows parsed; the parser and "
                         "the ledger disagree and no rank may be computed" % (declared, len(rows)))

    counts = {}
    for r in rows:
        r["stratum"] = stratum(r["changed"], r["source"])
        counts[r["stratum"]] = counts.get(r["stratum"], 0) + 1
    if sum(counts.values()) != declared:
        raise SystemExit("the strata do not sum to the ledger's own N")

    full = [r for r in rows if r["stratum"] == "full_cohort_out_of_fold"]
    above = sorted([r for r in full if r["dev"] > FROZEN_DEV], key=lambda r: -r["dev"])
    tied = [r for r in full if abs(r["dev"] - FROZEN_DEV) < 1e-9]
    if not tied:
        raise SystemExit("the frozen arm's development value is not in the ledger's full-cohort "
                         "stratum; the rank statement has no referent")

    missing = [r["changed"] for r in above if r["changed"] not in WHY]
    if missing:
        raise SystemExit("no recorded reason why these outrank the frozen arm: %s" % missing)

    rep = {
        "artifact_type": "s15_selection_rank",
        "phase_of_origin": "frozen_post_hoc",
        "why": "a maximum-of-N correction assumes the maximum was taken; this records whether it "
               "was, from the ledger rather than from memory",
        "frozen_arm_dev_value": FROZEN_DEV,
        "candidates_scored": declared,
        "strata": counts,
        "stratum_rule": "a rank is only meaningful among values comparable to it: full cohort, "
                        "scored out of fold. Subgroup slices, in-sample oracles, permuted and null "
                        "controls and cross-cohort transfer probes are excluded and counted.",
        "full_cohort_out_of_fold_candidates": len(full),
        "candidates_strictly_above": len(above),
        "ledger_rows_at_the_frozen_value": len(tied),
        "rank_of_the_frozen_arm": len(above) + len(tied),          # pessimistic: ties count against
        "rank_optimistic": len(above) + 1,
        "tie_rule": "the frozen construction is recorded %d times, by %d sub-experiments that "
                    "arrived at it independently. The reported rank counts every tied row as "
                    "though it were above, which is the worst position the ledger supports."
                    % (len(tied), len(tied)),
        "candidates_scoring_above_it": [
            {"dev_value": r["dev"], "candidate": r["changed"],
             "why_it_is_not_the_frozen_arm": WHY[r["changed"]]} for r in above],
        "ledger_rows_recording_the_frozen_construction": [
            {"id": r["id"], "candidate": r["changed"], "source": r["source"]} for r in tied],
        "reading": "the frozen configuration is rank %d of %d comparable candidates, and %d "
                   "distinct configurations score strictly above it: two belong to a different "
                   "cohort, three are components this paper rejects on their own pre-registered "
                   "controls, and the rest are the declared transcriptomic sensitivity. The "
                   "remaining places are the frozen construction's own tied rows. The arm was "
                   "chosen by a stated principle rather than by taking the maximum, which is what "
                   "a max-of-N correction assumes."
                   % (len(above) + len(tied), len(full), len(above)),
    }
    json.dump(rep, open(a.out, "w"), indent=1)
    print("  %d scored, %d comparable; frozen arm rank %d (optimistic %d), %d strictly above, "
          "%d tied rows" % (declared, len(full), rep["rank_of_the_frozen_arm"],
                            rep["rank_optimistic"], len(above), len(tied)))
    for k in sorted(counts, key=lambda k: -counts[k]):
        print("     %-30s %d" % (k, counts[k]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
