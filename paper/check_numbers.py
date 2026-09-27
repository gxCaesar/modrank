#!/usr/bin/env python3
"""Every load-bearing number in main.tex, checked against the file it came from.

This exists because two defects got into the manuscript that no gate could see. The first was a
SIGN: amendment A1 was described as moving the headline +0.0050 when the reported primary fell
-0.0048, because the sentence kept the flattering half of a two-part change. The second was a
CONTRADICTION: +0.0318 in the results text against +0.0334 in the multiplicity table, because our
arm was built at seed 0 in one place and seed-averaged in the other.

Neither is caught by "does the number appear in a results file" -- both numbers did. What catches
them is asserting the number against the SPECIFIC field it claims to be, and asserting that the
manuscript does not contain two values for one quantity. So each row below names a path into a
results JSON, and a missing path is a failure rather than a skip.

Run:  /usr/bin/python3 paper/check_numbers.py
"""

from __future__ import annotations

import ast
import collections
import json
import math
import os
import re
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "experiments", "20260817-blca-confirm", "results")
TEX = os.path.join(ROOT, "paper", "main.tex")

# The Briefings in Bioinformatics (BiB) build was retired on 2026-09-27 when the paper moved to
# Nature Communications (paper/nc-submission/). main() never opens or checks these five documents
# again; the live manuscript is checked by check_nc() and check_nc_si() instead. Two derived files
# of the same retired build (paper/bib-submission/main.aux and supplementary.aux) are likewise no
# longer read; they are not manuscripts themselves so they are not added to this set.
RETIRED_DOCS = frozenset({
    TEX,                                                            # paper/main.tex
    os.path.join(ROOT, "paper", "methods.tex"),
    os.path.join(ROOT, "paper", "bib-submission", "main.tex"),
    os.path.join(ROOT, "paper", "bib-submission", "supplementary.tex"),
    os.path.join(ROOT, "paper", "cover_letter.tex"),                # the BiB cover letter
})


def load(name):
    with open(os.path.join(RES, name)) as fh:
        return json.load(fh)


S5 = os.path.join(ROOT, "development", "s5-results")


def load5(name):
    """The modality figure reads the S5 evidence directory, which the confirmatory checks do not."""
    with open(os.path.join(S5, name)) as fh:
        return json.load(fh)


def dig(obj, path):
    """Walk a slash-separated path; raise rather than return a default."""
    cur = obj
    for part in path.split("/"):
        if isinstance(cur, list):
            cur = cur[int(part)]
        else:
            if part not in cur:
                raise KeyError("no such field: %s (stopped at %r)" % (path, part))
            cur = cur[part]
    return cur


# (label, file, path, expected-as-written-in-the-paper, tolerance)
CHECKS = [
    ("primary metric", "amendment-A1-clinical-provenance.json", "primary/value", 0.7212, 5e-5),
    ("primary SD", "amendment-A1-clinical-provenance.json", "primary/sd_over_seeds", 0.0048, 5e-5),
    ("amendment moved the reported primary", "amendment-A1-clinical-provenance.json",
     "decision/moved_from_frozen_run", -0.0048, 5e-5),
    ("frozen-run primary", "amendment-A1-clinical-provenance.json",
     "frozen_run_for_comparison/primary", 0.7260, 5e-5),
    ("clinical with stage", "amendment-A1-clinical-provenance.json",
     "arms_seed0/clinical_age_sex_stage_from_DIMAF_file", 0.6638, 5e-5),
    ("clinical with grade", "amendment-A1-clinical-provenance.json",
     "arms_seed0/clinical_age_sex_GRADE_from_DIMAF_file", 0.5666, 5e-5),
    ("stage minus grade", "amendment-A1-clinical-provenance.json",
     "arms_seed0/stage_minus_grade", 0.0972, 5e-5),
    ("ours at seed 0", "amendment-A1-clinical-provenance.json", "arms_seed0/OURS", 0.7225, 5e-5),
    ("ours, seed-averaged score", "pibd-parity.json",
     "canonical_arms_seed_averaged_scores/ours", 0.7237, 5e-5),
    ("SurvPath + clinical, seed-averaged", "pibd-parity.json",
     "canonical_arms_seed_averaged_scores/survpath_plus_clinical", 0.6897, 5e-5),
    ("PIBD + clinical, best-val", "pibd-parity.json",
     "canonical_arms_seed_averaged_scores/pibd_plus_clinical_best_val", 0.6989, 5e-5),
    ("PIBD + clinical, final epoch", "pibd-parity.json",
     "canonical_arms_seed_averaged_scores/pibd_plus_clinical_final_epoch", 0.6639, 5e-5),
    ("margin vs SurvPath", "pibd-parity.json",
     "holm_family_of_six/comparisons/ours vs SurvPath+clinical (seed-matched)/delta", 0.0340, 5e-5),
    ("p vs SurvPath", "pibd-parity.json",
     "holm_family_of_six/comparisons/ours vs SurvPath+clinical (seed-matched)/p_raw", 0.0843, 5e-5),
    ("margin vs PIBD", "pibd-parity.json",
     "comparison_ours_vs_pibd_plus_clinical/paired_vs_best_val/mean", 0.0249, 5e-5),
    ("p vs PIBD", "pibd-parity.json",
     "comparison_ours_vs_pibd_plus_clinical/paired_vs_best_val/p_two_sided", 0.2413, 5e-5),
    ("margin vs PIBD final epoch", "pibd-parity.json",
     "comparison_ours_vs_pibd_plus_clinical/paired_vs_final/mean", 0.0595, 5e-5),
    ("PIBD reproduction, fold mean", "pibd-parity.json",
     "pibd/best_validation_epoch/fold_mean", 0.6609, 5e-5),
    ("checkpoint selection is worth", "pibd-parity.json",
     "published/checkpoint_selection_is_worth", 0.0709, 5e-5),
    ("our parameter count", "parameter-counts.json", "ours/parameters", 1049, 0),
    ("SurvPath parameter count", "parameter-counts.json", "survpath/parameters", 24702532, 0),
    ("PIBD parameter count", "parameter-counts.json", "pibd/parameters", 26930702, 0),
    ("detectable margin at 80% power", "reviewer-gaps.json",
     "m_detectable_margin/detectable_at_80pc_power_two_sided_0.05", 0.0589, 2e-4),
    ("case-level bootstrap SE", "reviewer-gaps.json",
     "m_detectable_margin/case_level_bootstrap_se", 0.0210, 5e-5),
    # The Results described the ten-case series with the stage composition and event times of a
    # DIFFERENT analysis -- "all five down-revised are stage IV" when none of them is. The prose is
    # generated from these fields now; these pin them.
    ("case series A, band size", "interpret-calibrate-cases.json",
     "I_case_studies/band_size", 119, 0),
    ("case series A, up-revised events", "interpret-calibrate-cases.json",
     "I_case_studies/summary/high_group/events", 2, 0),
    ("case series A, down-revised events", "interpret-calibrate-cases.json",
     "I_case_studies/summary/low_group/events", 1, 0),
    # The published table now lives in a results file whose gap list is RECOMPUTED from its own
    # entries. That recomputation is what caught the manuscript describing the gaps as 0.002-0.012
    # when two of the eight are 0.021 and 0.029.
    ("largest published consecutive gap", "published-benchmark-table.json",
     "consecutive_gap_range/1", 0.029, 5e-5),
    ("smallest published consecutive gap", "published-benchmark-table.json",
     "consecutive_gap_range/0", 0.002, 5e-5),
]


def main():
    fails = []

    # CHECKS binds each result-file value to the literal exactly as the BiB paper wrote it. The
    # BiB paper is retired (RETIRED_DOCS, 2026-09-27), so this table is no longer walked; the list
    # itself is left in place rather than deleted. The NC build binds its own numbers in check_nc().
    print("  --  %-38s 29 literal bindings retired with the BiB build" % "CHECKS table")

    # --- the contradiction check: a superseded value must not survive anywhere in the manuscript
    SUPERSEDED = {
        "0.0318": "margin vs SurvPath under the OLD mixed construction",
        "0.0334": "margin vs SurvPath under the OLD seed-0 numerator",
        "0.116":  "p vs SurvPath, superseded",
        "0.263":  "p vs PIBD, superseded",
        "0.6982": "PIBD+clinical under the OLD seed-0 clinical block",
        "$p=0.063$": "a p for the clinical-alone comparison that matches no source (0.0163 raw, 0.0652 Holm)",
        "0.5012": "tied-pair clinical C under the PRE-amendment clinical block (now 0.4959)",
        "0.6120": "tied-pair slide C under the PRE-amendment clinical block (now 0.6220)",
        "1,281": "tied-pair count under the PRE-amendment clinical block (now 1,213)",
        "is 0.213": "slide-clinical Spearman under the PRE-amendment block (now 0.2445)",
        "arm 0.166": "incumbent-clinical Spearman under the PRE-amendment block (now 0.1924)",
    }

    # --- a CONTEXT check, not a ban. "0.002 to 0.012" is a true statement about six of the eight
    # published gaps and a false one about the range, and a blunt substring ban cannot tell those
    # apart -- it flagged the corrected sentence. So: wherever the narrow band appears, the largest
    # gap must appear near it, which is exactly the thing whose absence made the claim wrong.
    # --- a DERIVED count, recomputed rather than trusted. "the corrected clinical baseline clears
    # N of the M published entries" is the kind of claim that is written once and then survives
    # every later change to either number. It was wrong twice over -- six of ten, against a true
    # five of nine -- until this check was written.
    pub = load("published-benchmark-table.json")
    a1 = load("amendment-A1-clinical-provenance.json")
    stage_c = a1["arms_seed0"]["clinical_age_sex_stage_from_DIMAF_file"]
    n_below = sum(1 for e in pub["entries"] if e["cindex"] < stage_c)
    n_pub = len(pub["entries"])
    WORD = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven",
            8: "eight", 9: "nine", 10: "ten", 11: "eleven"}
    # TEX, methods.tex and the bib-submission pair were retired from this scan on 2026-09-27
    # (RETIRED_DOCS); the live manuscript is the NC build below.
    SCAN = (os.path.join(ROOT, "paper", "nc-submission", "main.tex"),
            os.path.join(ROOT, "paper", "nc-submission", "supplementary.tex"),
            os.path.join(ROOT, "paper", "nc-submission", "cover_letter.tex"))
    # The released archive omits the cover letter (editor correspondence, not code), and the gate has
    # to run to completion there too. A release tree has no research_state.yaml; in the project the
    # letter must exist and a missing one still fails below.
    if not os.path.isfile(os.path.join(ROOT, "research_state.yaml")):
        _absent = [p_ for p_ in SCAN if not os.path.isfile(p_)]
        for p_ in _absent:
            print("  --  %-38s not in this archive, skipped" % os.path.relpath(p_, ROOT))
        SCAN = tuple(p_ for p_ in SCAN if p_ not in _absent)
    # --- STALE CLAIMS, banned as phrases. Found by an external review on 2026-09-11: after the
    # selection correction was shown to be mis-specified (0.0397, not 0.0073), four sentences still
    # said the result "clears the selection-inflation term", and the Limitations quoted a p that
    # matched no source. The derived +0.0178 margin that used to be checked here was one of those
    # sentences, so the check that guarded its arithmetic now guards its absence. Each phrase below
    # is a claim the evidence no longer supports; none of them may come back by a later edit.
    STALE = {
        "played no part in developing": "the other four studies were consulted during development (cold panel round 2)",
        "never developed on": "the same",
        "correctly specified correction": "the selection bar is a sensitivity bound measured under a permuted null",
        "correctly measured bar": "the same",
        "gives calibrated two-year risks": "calibrated in the large, slope 0.770",
        "The only development cohort with slides": "the other four studies' slide embeddings were used in a rejected candidate",
        # cold panel round 3: prose that survived amendment A2 unchanged while its result reversed
        "moved toward the stage": "after A2 the screened variant fell below the stage block in GSE48075",
        "moved the average toward the stage block": "the same",
        "exceeded the stage block by less than the bar in both": "the same",
        "only published number computed on this partition": "PIBD's released split files are the same, and it reports higher values outside bladder",
        "replicated in an independent cohort": "D against zero is passed by a score without information; see SI Note 8",
        "The inflation replicates": "the same",
        "The inflation is significant for all three": "D was tested against zero; the benchmark is a score without information",
        "a grade-based reference inflated the": "the same, in the external-cohort paragraph",
        "neither concatenated nor stacked": "five fitted fusions were scored, not two",
        "neither fitted fusion": "the same",
        "no fitted fusion separates": "the same",
        "clears the selection-inflation": "the frozen selection term is mis-specified; nothing clears it as a test",
        "and the selection-inflation term": "the same, in a list of bars the result clears",
        "clears it by $+0.0178$": "the derived margin over the mis-specified term",
        "verdict is unchanged": "the corrected bar reverses the pre-registered verdict",
        "still clears the pre-registered": "the same",
        "cannot disadvantage": "late fusion weights the clinical arm 1/2 for a competitor and 1/3 for ours",
        "present in neither": "a pair restriction is not a conditional-information result",
        "what all 11 method papers": "three of the eleven were never read",
        "none of the eleven uses": "the same",
        "Each paper was read in full": "three could not be read",
        "Every one of those papers is measured against a clinical baseline": "four report none",
        "whether neoadjuvant chemotherapy is given": "neoadjuvant chemotherapy precedes the specimen that is staged",
        "roughly six times the largest gap": "0.0709 is 2.4 times the largest gap and 7 times the median",
        "cannot resolve concordance differences": "a power statement, not an impossibility",
        "two-centroid rule": "the sign-split table had no committed emitter; the median split does",
        # --- the 2026-09-12 desk-reject review. Each of these overstated what was measured, and the
        # replacement is in the same sentence. A later edit may not restore any of them.
        "the pathology report holds nothing that separates them":
            "the staging variables do not separate them; the report holds morphology the slide arm reads",
        "and neither has been examined":
            "we examined the method papers reported on the benchmark, which is a narrower claim",
        "ModRank leads the benchmark":
            "it has the highest reported concordance and cannot be distinguished from two competitors",
        "the added value of every multimodal model":
            "three constructions were measured, not every model",
        "both input-parity comparisons":
            "parity is in the clinical variables only; ModRank reads a slide encoder they do not",
        "The three arms are close to independent":
            "0.2445 and 0.1924 are weak correlations, not independence",
        "correlated in exactly the way the real search was":
            "they share folds, features and estimator; 'exactly' was not measured",
        "It cannot make a combination exceed":
            "in the Lund cohort the combination reaches 0.8779 against its best arm's 0.8672",
        # --- Phase A4 claim scoping, 2026-09-27. A sign count over 24 re-partitions measures
        # stability, not superiority; the margin over stacking (+0.0111) is under the 0.0145 bar; and
        # across the five studies no fitted fusion separates in either direction (breast leans the
        # other way, p=0.061). None of these may return.
        "beat concatenated and learned fusion":
            "the margin over stacking does not clear 0.0145, and five studies show no separation",
        "the strongest construction on identical inputs":
            "a point-estimate ranking by less than the benchmark resolves against stacking",
        "exceeds concatenated and learned fusion on identical inputs in every one of 24":
            "a sign count presented as superiority",
        "The simplest combination wins here":
            "it is not separable from stacking on bladder and from either fusion in five studies",
        "fitting the combination buys nothing":
            "breast: concatenation ahead by 0.0679, p=0.061; the title that said so was retired",
        "method beats a correctly specified clinical reference outside the development cohort":
            "two of the four other studies separate, two do not",
        "a property of the recipe rather than of the encoder":
            "two encoders agree on which studies separate; that is all it shows",
        "three published multimodal constructions":
            "ModRank is one of the three and is not published",
        "the encoder's contribution near 0.022":
            "0.7225 and 0.7009 are the two values on the shared seed, a difference of 0.0216",
    }
    for _doc in SCAN:
        _raw = open(_doc).read()
        _flat = " ".join(_raw.split())
        for _p, _why in STALE.items():
            _pp = " ".join(_p.split())
            if _pp in _flat:
                fails.append("STALE CLAIM %r in %s -- %s" % (_p, os.path.relpath(_doc, ROOT), _why))
    print("  ok  %-38s %d phrases, %d documents" % ("stale claims absent", len(STALE), len(SCAN)))

    # --- TABLE STRUCTURE. Found 2026-09-11: an automated em-dash rewrite (commit 32e121d) turned a
    # "---" cell pair into "(" and ")" and merged two rows of the Holm table, and it shipped in the
    # PDF and in the archived v1.0.0 for three weeks, because every gate here read numbers and none
    # read the shape of a table. Every row of every tabular must carry its declared column count, and
    # no cell may be a bare parenthesis.
    _ntab = 0
    for _doc in SCAN:
        _src = open(_doc).read()
        for _m in re.finditer(r"\\begin\{tabular\}\{([^}]*)\}(.*?)\\end\{tabular\}", _src, re.S):
            _ntab += 1
            _ncol = len(re.findall(r"[lcrp]", re.sub(r"\{[^}]*\}", "", _m.group(1))))
            _body = re.sub(r"\\(toprule|midrule|bottomrule|hline)", "", _m.group(2))
            for _row in _body.split("\\\\"):
                _row = _row.strip()
                if not _row or _row.startswith("%") or "\\multicolumn" in _row:
                    continue
                _cells = [c.strip() for c in _row.split("&")]
                _line = _src[:_m.start()].count("\n") + 1
                if len(_cells) != _ncol:
                    fails.append("%s: a row of the table at line %d has %d cells where the table "
                                 "declares %d: %r" % (os.path.relpath(_doc, ROOT), _line,
                                                      len(_cells), _ncol, _row[:90]))
                for _c in _cells:
                    if _c in ("(", ")"):
                        fails.append("%s: a bare %r cell in the table at line %d"
                                     % (os.path.relpath(_doc, ROOT), _c, _line))
    print("  ok  %-38s %d tables, every row the declared width" % ("table structure", _ntab))

    # --- 0.6856 IS TWO QUANTITIES. It is ModRank with the grade arm (amendment-A1) and, by
    # coincidence, the development-era 23-column clinical block. The prose once called the second
    # "the corrected clinical block". Every occurrence must now sit beside what disambiguates it:
    # 0.7225 (the grade-for-stage swap inside ModRank) or the words naming the development-era block.
    # Amendment A2 moved ModRank with the grade arm to 0.6870, so the coincidence is gone; the scan
    # still requires any 0.6856 left in the documents to name the development-era block.
    for _doc in SCAN:
        _b = " ".join(open(_doc).read().split())
        for _m in re.finditer(r"0\.6856", _b):
            _w = _b[max(0, _m.start() - 260):_m.end() + 260]
            if "0.7225" not in _w and "Amendment~A1 replaced" not in _w:
                fails.append("%s quotes 0.6856 without saying which of its two quantities it is: %r"
                             % (os.path.relpath(_doc, ROOT), _w[200:330]))
    print("  ok  %-38s every 0.6856 disambiguated" % "the two 0.6856s")

    # --- THE SUBTYPE TABLE, bound to the one committed emitter. The previous table (150/59 vs
    # 209/54, a sign split) was produced by no committed script; the figure beside it used a median
    # split. Both now come from biology.json B5, written by analysis/s8_biology.py. The TABLE ROWS
    # are parsed cell by cell: a substring test passed a mutated cell because the same value also
    # appears in the prose (caught by this gate's own mutation test on 2026-09-11).
    # This checked the subtype table (tab:subtype) and its prose paragraph in bib-submission/main.tex
    # and paper/main.tex only -- both retired with the BiB build on 2026-09-27 (RETIRED_DOCS).
    print("  --  %-38s retired with the BiB build" % "within-subtype table")

    # --- THE LIMITATIONS' p FOR THE CLINICAL-ALONE COMPARISON, bound to the Holm family. This
    # checked bib-submission/main.tex and paper/main.tex only -- both retired with the BiB build
    # on 2026-09-27 (RETIRED_DOCS).
    print("  --  %-38s retired with the BiB build" % "clinical-alone p in Limitations")

    # --- THE PRE-REGISTERED SELECTION TERM, bound to the FROZEN PROTOCOL rather than to the
    # ledger's header. Found by the 2026-08-20 carpet sweep. The protocol froze
    # sigma*sqrt(2 ln N) at N = 211, giving 0.0239 and a target of 0.7029; the manuscript
    # described the SAME term at N = 269 and 0.0244 and quoted 0.7034 as the pre-registered bar in
    # one sentence while five others quoted 0.7029. Both arithmetics are right. Attributing the
    # second to the pre-registration is not, and it is the one misstatement a paper arguing from
    # its own pre-registration cannot carry. The ledger FINISHED at 269 and its header recomputes
    # the term over that final count, which is where the prose picked the wrong number up.
    _proto = json.load(open(os.path.join(ROOT, "development", "benchmark-protocol.json")))
    _unc = _proto["metrics"]["uncertainty"]
    _n_frozen = int(re.search(r"N\s*=\s*(\d+)", _unc["selection_inflation_basis"]).group(1))
    _sig = _unc["paired_split_reseed_sd"]
    _infl = _unc["selection_inflation"]
    _target = _proto["amendments"][0]["effect_on_the_frozen_primary"]["target"]
    if abs(_sig * math.sqrt(2 * math.log(_n_frozen)) - _infl) > 5e-4:
        fails.append("the protocol's own selection inflation does not regenerate: %.4f x "
                     "sqrt(2 ln %d) is not %.4f" % (_sig, _n_frozen, _infl))
    if abs(round(0.679 + _infl, 4) - _target) > 5e-5:
        fails.append("the protocol's target is not the incumbent plus its inflation: "
                     "0.679 + %.4f is not %.4f" % (_infl, _target))
    # the frozen rule stated with the ledger-final term is the specific defect, so it is banned
    # as a literal. 0.0244 remains legal everywhere it is labelled as the final-count figure.
    # This scanned bib-submission/main.tex, bib-submission/supplementary.tex, paper/main.tex and
    # paper/methods.tex -- all four retired with the BiB build on 2026-09-27 (RETIRED_DOCS). The
    # protocol arithmetic above (regeneration and target) is unaffected and still runs.
    print("  --  %-38s retired with the BiB build" % "frozen decision-rule text")
    print("  ok  %-38s frozen N=%d, inflation %.4f, target %.4f"
          % ("pre-registered selection term", _n_frozen, _infl, _target))

    # --- A CROSS-DOCUMENT DISPLAY-ITEM POINTER. The body and the supplement are separate compiled
    # documents, so a \ref between them renders as ?? and the only workable pointer is a typed
    # number. A typed number is checked by nobody: the supplement sent readers to "Figure 5e" for
    # the twelve-case series, which is Figure 7e, and Figure 5 has a real panel e about something
    # else, so the reader landed somewhere plausible and wrong. Found 2026-08-20. The number is
    # therefore bound to the main document's own .aux, and the figure is identified by what its
    # caption says rather than by a label spelled twice.
    # This bound the pointer to bib-submission/main.aux, main.tex and supplementary.tex -- all
    # retired with the BiB build on 2026-09-27 (RETIRED_DOCS plus its derived .aux).
    print("  --  %-38s retired with the BiB build" % "cross-document figure pointer")

    # --- THE COVER LETTER gets the two prose rules the loop above applies to the manuscript and
    # skipped for it. Its register is different (it is addressed to an editor, so the second person
    # is correct there and is not checked), but a path or a clause semicolon is wrong in a letter
    # for the same reason it is wrong in the paper.
    for _cl in (os.path.join(ROOT, "paper", "nc-submission", "cover_letter.tex"),):
      if os.path.exists(_cl):
        _cb = open(_cl).read()
        for _m in re.finditer(r"\\texttt\{([^}]*)\}", _cb):
            if "@" not in _m.group(1) and ("/" in _m.group(1) or _m.group(1).endswith(
                    (".json", ".csv", ".py", ".md"))):
                fails.append("the cover letter puts a path in the prose: %r" % _m.group(1))
        for _m in re.finditer(";", _cb):
            _ls = _cb.rfind("\n", 0, _m.start()) + 1
            _line = _cb[_ls:_cb.find("\n", _m.start())]
            if _line.lstrip().startswith("%") or "&" in _line:
                continue
            if re.match(r"\s*(\\item|\\end\{(itemize|enumerate)\})", _cb[_m.end():_m.end() + 60]):
                continue
            fails.append("the cover letter joins clauses with a semicolon: %r"
                         % " ".join(_cb[max(0, _m.start() - 46):_m.start() + 44].split()))
        print("  ok  %-38s %s" % ("cover letter prose", os.path.relpath(_cl, ROOT)))

    # --- THE MODALITY FIGURE'S NUMBERS. Its panels come from the S5 evidence directory rather
    # than the confirmatory results, so the table above never reaches them. Every quantity the
    # new Results subsection states in prose is recomputed here from the frozen source.
    # Both prose targets here were BiB documents (bib-submission/main.tex, paper/main.tex),
    # retired 2026-09-27 (RETIRED_DOCS); the NC build's own prose is bound in check_nc(). The
    # numeric checks below, which do not read any document, are unaffected.
    PROSE_MODALITY_DOCS = []
    _atlas, _c6 = load5("atlas.json"), load5("c6.json")
    _am = json.load(open(os.path.join(ROOT, "experiments", "20260911-blca-posthoc", "results",
                                      "modality-atlas-amended.json")))
    # under amendment A2 the pre-A2 known answer is checked only as a gross bound by s19b itself
    if not (_am["known_answer"]["reproduced"] or _am["known_answer"].get("amendment") == "A2"):
        fails.append("the amended atlas did not reproduce the development atlas first")
    _tie = _am["conditional_probe"]["clinically_tied_q05"]
    _img = [_c6["single_arms"][k] for k in ("wsi_titan", "wsi_chief_mean", "wsi_chief_dispersion")]
    _om = [_c6["single_arms"][k] for k in ("omics_xena", "omics_combine", "omics_hallmarks")]
    MOD_CHECKS = [
        # the AMENDED block. c6.json holds the development-era 23-column block (0.6856), which the
        # prose once called "the corrected clinical block"; that value now appears only where it
        # is labelled as the development-era block, and the 0.6856 context gate below enforces it.
        ("clinical alone, amended, seed 0",
         load("amendment-A1-clinical-provenance.json")["arms_seed0"]
         ["clinical_age_sex_stage_from_DIMAF_file"], "0.6637"),
        ("slide alone", _c6["single_arms"]["wsi_titan"], "0.6596"),
        ("pathway means alone", _c6["single_arms"]["omics_combine"], "0.6510"),
        ("spread within the image family", max(_img) - min(_img), "0.104"),
        ("spread within the omics family", max(_om) - min(_om), "0.069"),
        # redundancy and tied pairs from the AMENDED atlas (s19b), which first reproduces the
        # development atlas on its own pre-amendment block and then swaps only the clinical block
        ("rho slide vs clinical, amended", _am["redundancy"]["spearman_titan_vs_clinical"], "0.2443"),
        ("rho incumbent vs clinical, amended", _am["redundancy"]["spearman_survpath_vs_clinical"],
         "0.1927"),
        ("clinical on the tightest ties, amended", _tie["arms"]["clinical"], "0.4975"),
        ("slide on the tightest ties, amended", _tie["arms"]["titan"], "0.6214"),
        ("tightest-tie pair count, amended", _tie["pairs"], "1,211"),
    ]
    for label, value, written in MOD_CHECKS:
        if "," in written:                                   # a count written with a thousands comma
            if "{:,}".format(int(value)) != written:
                fails.append("modality figure: %s is %s at source, written as %s"
                             % (label, value, written))
            for doc in PROSE_MODALITY_DOCS:
                if written not in open(doc).read():
                    fails.append("%s never states %s (%s)" % (os.path.relpath(doc, ROOT), written, label))
            continue
        dp = len(written.split(".")[1])
        if "%.*f" % (dp, value) != written:
            fails.append("modality figure: %s is %.4f at source, written as %s"
                         % (label, value, written))
        for doc in PROSE_MODALITY_DOCS:
            if written not in open(doc).read():
                fails.append("%s never states %s (%s)"
                             % (os.path.relpath(doc, ROOT), written, label))
    if _tie["arms"]["clinical"] >= 0.52:
        fails.append("the tie panel's claim is that the clinical arm is a COIN FLIP there; "
                     "it reads %.4f" % _tie["arms"]["clinical"])
    print("  ok  %-38s %d quantities recomputed from source" % ("modality figure", len(MOD_CHECKS)))

    # --- THE REPORTING DUMP. Two supplementary figures and one Results paragraph are built on
    # it, and it is a SECOND execution of the frozen arm construction, so the first thing checked
    # is that it reproduced the frozen primary. Everything else here is recomputed from its own
    # per-case table rather than trusted from a summary field.
    dump = json.load(open(os.path.join(ROOT, "experiments", "20260818-reporting-dump", "results",
                                       "reporting-dump.json")))
    _a1p = load("amendment-A1-clinical-provenance.json")["primary"]["value"]
    if not (dump.get("frozen_primary_matched") or dump.get("amendment") == "A2") \
            or abs(dump["primary"] - _a1p) > 5e-5:
        fails.append("the reporting dump does not reproduce the frozen primary (%s)"
                     % dump.get("primary"))
    if len(dump["cases"]) != 359 or len(dump["pathways"]) != 275:
        fails.append("the dump should carry 359 cases and 275 pathways; it has %d and %d"
                     % (len(dump["cases"]), len(dump["pathways"])))
    if sum(1 for p in dump["pathways"] if p["q_BH"] < 0.10) != 0:
        fails.append("the supplement states that NO pathway clears BH q<0.10; some do")
    _folds = sorted({c["fold"] for c in dump["cases"]})
    _pf = []
    for fk in _folds:
        idx = [i for i, c in enumerate(dump["cases"]) if c["fold"] == fk]
        r = [dump["cases"][i]["seed_mean"]["OURS"] for i in idx]
        tt = [dump["cases"][i]["months"] for i in idx]
        ee = [dump["cases"][i]["event"] for i in idx]
        num = den = 0.0
        for a_ in range(len(idx)):
            if not ee[a_]:
                continue
            for b_ in range(len(idx)):
                if tt[b_] > tt[a_]:
                    den += 1.0
                    num += 1.0 if r[b_] < r[a_] else (0.5 if r[b_] == r[a_] else 0.0)
        _pf.append(num / den if den else float("nan"))
    for label, val, written in (("per-fold minimum", min(_pf), "0.680"),
                                ("per-fold maximum", max(_pf), "0.780"),
                                ("per-fold range", max(_pf) - min(_pf), "0.100")):
        if "%.3f" % val != written:
            fails.append("supplement: %s is %.4f recomputed, written as %s"
                         % (label, val, written))
        # the retired BiB supplement was the only prose target here (PROSE_MODALITY_DOCS is now
        # empty); the numeric check above, which reads only the reporting dump, is unaffected.
        for doc in PROSE_MODALITY_DOCS:
            if written not in open(doc).read():
                fails.append("%s never states %s (%s)" % (os.path.relpath(doc, ROOT), written,
                                                          label))
    print("  ok  %-38s primary reproduced, %d folds recomputed"
          % ("reporting dump", len(_folds)))
    # The two correlations the BiB supplement's pathway paragraph stated, recomputed from the
    # dump, are retired with the BiB build on 2026-09-27 (RETIRED_DOCS); the NC supplement's own
    # numbers are checked in check_nc_si().
    print("  --  %-38s retired with the BiB build" % "reporting-dump pathway correlations")


    # --- METRIC PARITY, answering round 2's G10. The paper now states both conventions on the
    # same predictions, so both are recomputed from the dump's own per-case table rather than
    # read back from the summary that produced them.
    _mp = os.path.join(ROOT, "experiments", "20260818-reporting-dump", "results",
                       "metric-parity.json")
    if os.path.exists(_mp):
        mp = json.load(open(_mp))
        for label, key, written in (("pooled out-of-fold", "pooled_out_of_fold", "0.7235"),
                                    ("mean of per-fold", "mean_of_per_fold", "0.7276"),
                                    ("their difference", "difference_mean_minus_pooled",
                                     "+0.0041")):
            v = mp[key]
            got = ("%+.4f" % v) if written.startswith("+") else ("%.4f" % v)
            if got != written:
                fails.append("metric parity: %s is %s at source, written as %s"
                             % (label, got, written))
            for doc in PROSE_MODALITY_DOCS:
                if written.lstrip("+") not in open(doc).read():
                    fails.append("%s never states %s (%s)"
                                 % (os.path.relpath(doc, ROOT), written, label))
        print("  ok  %-38s pooled %.4f vs fold-mean %.4f"
              % ("metric parity", mp["pooled_out_of_fold"], mp["mean_of_per_fold"]))

    # --- THE SELECTION NULL AND THE ENCODER SWAP, answering round 2's G1 and G8. G1 is the
    # paper's own decision rule failing its own correction, so every number in it is recomputed
    # from the run's output rather than trusted from the prose that reports it.
    _sn = os.path.join(ROOT, "experiments", "20260818-selection-null", "results",
                       "selection-null.json")
    if os.path.exists(_sn):
        sn = json.load(open(_sn))
        A, B = sn["A_selection_null"], sn["B_encoder_swap"]
        SEL = [("null single-candidate sd", A["null_single_candidate"]["sd"], "0.0398"),
               ("null single-candidate mean", A["null_single_candidate"]["mean"], "0.5011"),
               ("null max q95", A["null_max_over_the_family"]["q95"], "0.5918"),
               ("gaussian bar with the null sigma", A["bars"]["gaussian_with_the_null_sigma"]["bar"],
                "0.8089"),
               ("empirical max bar", A["bars"]["empirical_max_q95"]["bar"], "0.7708"),
               ("inflation of the null maximum above chance",
                A["null_max_over_the_family"]["q95"] - 0.5, "0.0918"),
               ("slide-arm span across encoders", B["slide_arm_span"], "0.1369"),
               ("full-method span across encoders", B["full_method_span"], "0.0529"),
               ("full method, weakest encoder", B["full_method_min"], "0.6695"),
               ("full method, strongest encoder", B["full_method_max"], "0.7224")]
        # the NC Supplementary Note 2 states the first six in prose (the BiB documents that also
        # stated the encoder spans are retired)
        _sel_si = " ".join(open(os.path.join(ROOT, "paper", "nc-submission",
                                             "supplementary.tex")).read().split())
        for label, value, written in SEL:
            if "%.4f" % value != written:
                fails.append("selection null: %s is %.4f at source, written as %s"
                             % (label, value, written))
            if "encoder" not in label and written not in _sel_si:
                fails.append("NC SI never states %s (%s)" % (written, label))
            for doc in PROSE_MODALITY_DOCS:
                if written not in open(doc).read():
                    fails.append("%s never states %s (%s)"
                                 % (os.path.relpath(doc, ROOT), written, label))
        # the two claims the paragraphs rest on, asserted rather than assumed
        if load("amendment-A1-clinical-provenance.json")["primary"]["value"] >= \
                A["bars"]["empirical_max_q95"]["bar"]:
            fails.append("the paper says the primary does NOT clear the corrected bar; it now does")
        if B["encoders_whose_full_method_still_beats_the_incumbent"] != 5:
            fails.append("the paper says five of seven encoders still beat 0.679; the run says %d"
                         % B["encoders_whose_full_method_still_beats_the_incumbent"])
        _sr = os.path.join(ROOT, "experiments", "20260818-selection-null", "results",
                           "selection-rank.json")
        if os.path.exists(_sr):
            sr = json.load(open(_sr))
            if sr["rank_of_the_frozen_arm"] != 11 or sr["full_cohort_out_of_fold_candidates"] != 206:
                fails.append("the paper says the frozen arm is rank 11 of 206; the ledger says "
                             "%d of %d" % (sr["rank_of_the_frozen_arm"],
                                           sr["full_cohort_out_of_fold_candidates"]))
            # The arithmetic a reader does out loud. Eight named above and a rank of 11 only
            # reconcile through the ties, so the tie count is load-bearing prose and is checked.
            if sr["candidates_strictly_above"] != 8 or sr["ledger_rows_at_the_frozen_value"] != 3:
                fails.append("the paper says eight configurations score strictly above and the "
                             "frozen construction is recorded three times; the ledger says %d "
                             "and %d" % (sr["candidates_strictly_above"],
                                         sr["ledger_rows_at_the_frozen_value"]))
            if (sr["candidates_strictly_above"] + sr["ledger_rows_at_the_frozen_value"]
                    != sr["rank_of_the_frozen_arm"]):
                fails.append("the rank does not equal the strictly-higher count plus the ties")
            # The stratification is quoted in the supplement, so its parts must sum to the N the
            # selection correction is computed over. A stratum that does not close is a stratum
            # with a silent residual, which is where a discarded candidate hides.
            if (sum(sr["strata"].values()) != sr["candidates_scored"]
                    or sr["candidates_scored"] != 269):
                fails.append("the ledger strata do not sum to 269: %r" % (sr["strata"],))
            # The strata-to-prose binding below checked only the BiB supplementary.tex; retired
            # with the BiB build on 2026-09-27 (RETIRED_DOCS). The sr-field checks above, which
            # read only the result files, are unaffected.
            print("  --  %-38s retired with the BiB build" % "strata bound to supplement prose")
        print("  ok  %-38s corrected bar %.4f, %d of 7 encoders clear 0.679"
              % ("selection null and encoder swap",
                 A["bars"]["empirical_max_q95"]["bar"],
                 B["encoders_whose_full_method_still_beats_the_incumbent"]))

    # --- DISPLAY-ITEM ORDER in the submission document. This was gated on the BiB main only
    # (bib-submission/main.tex): the internal draft is Nature-shaped, its abstract cites Figure 2
    # first, and that ordering was deliberate there. Retired with the BiB build on 2026-09-27
    # (RETIRED_DOCS).
    print("  --  %-38s retired with the BiB build" % "display-item order")

    # --- THE COVER LETTER. It is the one document an editor reads beside the paper, and a letter
    # that disagrees with its manuscript is the inconsistency they are guaranteed to see. Every
    # number in it must appear in the submission document, and the letter carries the same prose
    # rules as the paper.
    for _cl, _subdoc in ((os.path.join(ROOT, "paper", "nc-submission", "cover_letter.tex"),
                          os.path.join(ROOT, "paper", "nc-submission", "main.tex")),):
      if os.path.exists(_cl):
        letter = open(_cl).read()
        sub = open(_subdoc).read()
        body = re.sub(r"^\s*%.*$", "", letter, flags=re.M)          # skip the file's own comments
        nums = sorted(set(re.findall(r"\b\d[\d,]*\.\d+\b", body)))
        # vspace and geometry lengths are typography, not claims
        nums = [n for n in nums if n not in ("0.9", "0.95", "1.0", "1.2", "1.4", "1.6", "2.2")]
        for n in nums:
            if n not in sub:
                fails.append("the cover letter states %s, which the submission document does not"
                             % n)
        for pat, what in ((r"---", "an em-dash"), (r"(?<!-) -- (?!-)", "the dash construction")):
            if re.search(pat, body):
                fails.append("the cover letter uses %s" % what)
        print("  ok  %-38s %d numbers, all in the manuscript" % ("cover letter", len(nums)))

    # --- REGISTER, set by the author 2026-08-18. A journal body describes objectively. First
    # person plural is fine and standard here; what is not is the conversational construction --
    # "we read the papers one at a time", "we think", "nobody's contribution", "so we say which is
    # which" -- and contractions, second person, and loose indefinite nouns. NOTE that `cannot` is
    # a word, not a contraction: the first version of this check flagged 19 of them and every one
    # was a false positive.
    # The four BiB documents were removed from this list on 2026-09-27 (RETIRED_DOCS).
    PROSE_DOCS = [os.path.join(ROOT, "paper", "nc-submission", "main.tex"),
                  os.path.join(ROOT, "paper", "nc-submission", "supplementary.tex")]
    REGISTER = [(r"\b(?:can't|won't|isn't|aren't|doesn't|didn't|it's|that's|we're|don't|"
                 r"there's|let's)\b", "a contraction"),
                (r"\b(?:you|your|yours)\b", "the second person"),
                (r"\b(?:somebody|nobody|anybody|everybody)\b", "an indefinite person"),
                (r"\b[Ww]e (?:think|believe|feel|wanted|hoped|liked|would rather)\b",
                 "a first-person opinion"),
                (r"\b(?:a lot of|lots of|kind of|sort of|basically|obviously)\b",
                 "an informal intensifier")]
    n_reg = 0
    for doc in PROSE_DOCS:
        text = re.sub(r"^\s*%.*$", "", open(doc).read(), flags=re.M)
        rel = os.path.relpath(doc, ROOT)
        for pat, what in REGISTER:
            for m in re.finditer(pat, text):
                n_reg += 1
                fails.append("%s uses %s: %r" % (rel, what,
                             " ".join(text[max(0, m.start()-40):m.end()+34].split())))
    print("  ok  %-38s %d documents, no conversational construction" % ("prose register",
                                                                       len(PROSE_DOCS)))

    # --- MANUSCRIPT PROSE, two standing rules set by the author 2026-08-18.
    #   No file paths, repository paths or code identifiers in the body: a reader finds the artefact
    #   through Data availability, and the prose should read like a paper.
    #   No em-dash sentence, and no clause-joining semicolon. `---` joining two clauses, the
    #   paired `--- ... ---` parenthetical, and `;` linking two independent clauses all read as
    #   machine-written. This is a STYLE ruling, not a LaTeX one, so only the .tex bodies are
    #   scanned and the glyph stays legal inside a real paper title in the bibliography. A
    #   semicolon doing a DIFFERENT job stays: separating keywords, separating the items of an
    #   enumerate, separating several citations in one parenthesis, or inside math.
    #   Also scanned: a space before a comma or full stop. That is not a style opinion, it is the
    #   residue a bulk rewriter leaves behind, and 20 of them survived one such pass unnoticed.
    # The four BiB documents were removed from this list on 2026-09-27 (RETIRED_DOCS).
    PROSE_DOCS = [os.path.join(ROOT, "paper", "nc-submission", "main.tex"),
                  os.path.join(ROOT, "paper", "nc-submission", "supplementary.tex")]
    n_dash = n_path = n_semi = n_space = n_dup = 0
    for doc in PROSE_DOCS:
        body = open(doc).read()
        rel = os.path.relpath(doc, ROOT)
        # \TODO{...} holds text that is not yet prose and is gated separately, so its spans are
        # computed once with real brace matching rather than by counting braces in a window.
        todo_spans = []
        for t in re.finditer(r"\\newcommand\{|\\TODO(?:author)?\{|\\newcommand\{\\TODO", body):
            depth, k = 0, t.end() - 1
            while k < len(body):
                if body[k] == "{":
                    depth += 1
                elif body[k] == "}":
                    depth -= 1
                    if depth == 0:
                        break
                k += 1
            end = k + 1
            if body[t.start():t.end()].startswith("\\newcommand"):
                end = body.find("\n", t.start())        # a definition is one line here
                end = len(body) if end < 0 else end
            todo_spans.append((t.start(), end))

        # The author's rule names the em-dash AND the dash construction. `---` is caught below;
        # ` -- ` with spaces on both sides is the same construction spelled with an en-dash, and
        # it survived the em-dash pass. Compound terms (p53--cell-cycle), math ranges ($0$--$4$)
        # and page ranges carry no surrounding spaces, so they are untouched by this.
        for m in re.finditer(r"(?<!-) -- (?!-)", body):
            ls = body.rfind("\n", 0, m.start()) + 1
            if body[ls:m.start()].lstrip().startswith("%"):
                continue
            if any(a <= m.start() < b for a, b in todo_spans):  # unresolved placeholder text
                continue
            n_dash += 1
            fails.append("%s uses the dash construction: %r"
                         % (rel, " ".join(body[max(0, m.start()-44):m.end()+40].split())))

        for m in re.finditer(r"---", body):
            n_dash += 1
            fails.append("%s:%d uses an em-dash sentence: %r"
                         % (rel, body[:m.start()].count("\n") + 1,
                            " ".join(body[max(0, m.start()-40):m.start()+50].split())))
        for m in re.finditer(r"\\texttt\{([^}]*)\}", body):
            tok = m.group(1)
            if "@" in tok:            # the corresponding author's email is front matter
                continue
            if "/" in tok or "\\_" in tok or tok.endswith((".json", ".csv", ".py", ".md")):
                n_path += 1
                fails.append("%s:%d puts a path or code identifier in the prose: %r"
                             % (rel, body[:m.start()].count("\n") + 1, tok))
        # A semicolon is judged by the JOB IT IS DOING, not by its presence.
        for m in re.finditer(";", body):
            ls = body.rfind("\n", 0, m.start()) + 1
            le = body.find("\n", m.start())
            line = body[ls:le if le != -1 else len(body)]
            if line.lstrip().startswith("%"):                       # a source comment
                continue
            if "&" in line or line.rstrip().endswith("\\\\"):        # a table row
                continue
            if "\\keywords" in body[max(0, m.start() - 220):m.start()]:
                continue                                            # keyword separator
            # a list separator: the next thing in the source is another item, or the list ends.
            # `\item` usually begins the NEXT line, so looking at this one is not enough.
            if re.match(r"\s*(\\item|\\end\{(itemize|enumerate|description)\})",
                        body[m.end():m.end() + 60]):
                continue
            if body[:m.start()].count("$") % 2 == 1:                # inside math
                continue
            if re.search(r"(Table|Fig|Supplementary|Supplementary Table)[^;]{0,40}$",
                         body[max(0, m.start() - 60):m.start()]):   # several citations in one paren
                continue
            n_semi += 1
            fails.append("%s:%d joins clauses with a semicolon: %r"
                         % (rel, body[:m.start()].count("\n") + 1,
                            " ".join(body[max(0, m.start()-46):m.start()+44].split())))

        # An immediately repeated phrase. Four of these survived the path-removal pass, one of
        # them with no space at the join ("...patientthe up-revised patient..."), which is exactly
        # the shape a doubled-WORD check misses.
        flat = " ".join(body.split())
        # case-insensitive: a repeat whose first copy was sentence-initial ("At this cohort
        # sizeat this cohort size") is invisible to a case-sensitive backreference
        for m in re.finditer(r"\b((?:[\w,;:]+ ){1,5}[\w,;:]+)\s*\1\b", flat, re.I):
            # an enumeration legitimately repeats its own items ("stages II, III, II, III, II"),
            # so a repeat only counts when it contains a real word rather than numerals alone
            if not re.search(r"[A-Za-z]{4,}", m.group(1)):
                continue
            n_dup += 1
            fails.append("%s repeats a phrase immediately: %r" % (rel, m.group()[:80]))

        for m in re.finditer(r"(\S)(?:[ ]+|[ ]*\n[ ]*)([,.])", body):
            ls = body.rfind("\n", 0, m.start()) + 1
            if body[ls:m.start()].lstrip().startswith("%"):         # a source comment
                continue
            if body[m.end() - 1:m.end() + 2] == "...":              # an ellipsis, not a full stop
                continue
            n_space += 1
            fails.append("%s:%d leaves a space before %r, which is rewriter residue: %r"
                         % (rel, body[:m.start()].count("\n") + 1, m.group(2),
                            " ".join(body[max(0, m.start()-40):m.end()+30].split())))

    print("  ok  %-38s %d documents, 0 em-dashes, 0 clause semicolons, 0 paths, "
          "0 residue, 0 repeats" % ("manuscript prose", len(PROSE_DOCS)))

    # --- FIGURE GEOMETRY. This checked the insertion scale of figures included by
    # bib-submission/main.tex and supplementary.tex against the OUP class's two-column textwidth
    # (488.5pt) and the internal article build's textwidth (453.6pt) -- both BiB-specific and both
    # retired with the BiB build on 2026-09-27 (RETIRED_DOCS).
    print("  --  %-38s retired with the BiB build" % "figure geometry")

    # --- The clinical-baseline census. The manuscript's central claim is a claim about other
    # people's papers, and it was stated as "every one of the eleven" while the paper named ten at
    # most and half of those read report no clinical baseline at all. These pin what was read.
    cen = json.load(open(os.path.join(ROOT, "development", "clinical-baseline-census.json")))
    cs = cen["summary"]
    with_b = [p for p in cen["papers"] if p.get("reports_clinical_baseline") is True]
    check_grade = all("grade" in " ".join(p["covariates"]).lower() for p in with_b)
    if not check_grade:
        fails.append("a paper reporting a clinical baseline does not use grade; the claim that "
                     "every one of them does no longer holds")
    if any(p["uses_stage"] for p in with_b):
        fails.append("a paper reporting a clinical baseline DOES use stage; the paper's central "
                     "claim is falsified")
    if cs["verified_report_a_clinical_baseline"] != len(with_b):
        fails.append("census summary and rows disagree on how many report a baseline")
    print("  ok  %-38s %d report one (all grade, none stage), %d report none, %d unverified"
          % ("clinical-baseline census", cs["verified_report_a_clinical_baseline"],
             cs["verified_report_none"], cs["not_verified"]))

    # --- The parameter counts are quoted verbatim in the manuscript as having been "printed by
    # each method's own trainer". Until the results were pulled back we held only OUR TRANSCRIPTION
    # of that claim. Now the trainers' own output files are in the repository, so read them.
    import glob as _g
    TRAINER = {
        "survpath": _g.glob(os.path.join(ROOT, "experiments", "*survpath*", "**",
                                         "model_tcga_blca__*.txt"), recursive=True),
        "pibd": _g.glob(os.path.join(ROOT, "experiments", "*pibd*", "**",
                                     "model_parameters.txt"), recursive=True),
    }
    pc = load("parameter-counts.json")
    for who, paths in TRAINER.items():
        if not paths:
            fails.append("no trainer output file found for %s; the parameter count in the "
                         "manuscript is a transcription with nothing behind it" % who)
            continue
        seen = set()
        for path in paths:
            for line in open(path, errors="replace"):
                if "Total number of parameters" in line:
                    seen.add(int(line.split(":")[1].strip()))
                    break
        if len(seen) != 1:
            fails.append("%s trainer files disagree on the parameter count: %s" % (who, seen))
        elif seen.pop() != pc[who]["parameters"]:
            fails.append("%s: the trainer's own file and parameter-counts.json disagree" % who)
        else:
            print("  ok  %-38s %s, from %d trainer file(s)"
                  % ("%s parameter count, at source" % who,
                     "{:,}".format(pc[who]["parameters"]), len(paths)))

    # --- CROSS-DOCUMENT references. This bound bib-submission/main.tex's hand-written pointers to
    # bib-submission/supplementary.aux's own rendered labels -- both retired with the BiB build on
    # 2026-09-27 (RETIRED_DOCS plus its derived .aux). check_nc_si() binds the NC build's own
    # cross-document references against paper/nc-submission/supplementary.aux separately.
    print("  --  %-38s retired with the BiB build" % "supplementary cross-references")
    # The count is stated on TWO denominators since 2026-09-11: the entries on the released folds
    # (verified here or by definition) and all nine. A sentence that mixes the two was the defect.
    def _verified(e):
        return e["folds"].startswith("verified") or e["folds"] in ("definitional",
                                                                    "as reported in SurvPath")
    _ver = [e for e in pub["entries"] if _verified(e)]
    n_below_v = sum(1 for e in _ver if e["cindex"] < stage_c)
    n_ver = len(_ver)
    for path in SCAN:
        body = " ".join(open(path).read().split())
        # "every one of the nine" is a different claim (a value above all of them), not a count
        for m in re.finditer(r"(?<!every )\b(\w+) of the (\w+) published entries( on verified folds)?",
                             body, re.I):
            got_below, got_pub = m.group(1).lower(), m.group(2).lower()
            want = (WORD.get(n_below_v), WORD.get(n_ver)) if m.group(3) else \
                   (WORD.get(n_below), WORD.get(n_pub))
            if (got_below, got_pub) != want:
                fails.append("%s says '%s of the %s published entries%s'; recomputing from the "
                             "results files gives '%s of the %s'"
                             % (os.path.relpath(path, ROOT), got_below, got_pub,
                                m.group(3) or "", want[0], want[1]))
        for m in re.finditer(r"(\w+) of all (\w+)", body):
            if m.group(2).lower() in WORD.values() and m.group(1).lower() in WORD.values():
                if (m.group(1).lower(), m.group(2).lower()) != (WORD.get(n_below), WORD.get(n_pub)):
                    fails.append("%s says '%s of all %s'; the results files give '%s of all %s'"
                                 % (os.path.relpath(path, ROOT), m.group(1), m.group(2),
                                    WORD.get(n_below), WORD.get(n_pub)))
        print("  ok  %-38s %s of %s verified, %s of %s in all (%s)"
              % ("published-entry count", WORD.get(n_below_v), WORD.get(n_ver),
                 WORD.get(n_below), WORD.get(n_pub), os.path.relpath(path, ROOT)))

    for path in SCAN:
        body = open(path).read()
        # +0.0050 also appears as a CONFIDENCE-INTERVAL BOUND in the multiplicity table, which has
        # nothing to do with amendment A1. Only flag it where the surrounding text is actually
        # about the amendment -- the first version of this check fired on the CI bound.
        for m in re.finditer(r"\+0\.0050\$", body):
            window = body[max(0, m.start() - 700):m.end() + 700]
            about_a1 = ("defect" in window) or ("amendment" in window.lower())
            if about_a1 and "0.0048" not in window:
                line = body[:m.start()].count("\n") + 1
                fails.append("%s:%d gives amendment A1's +0.0050 without the net -0.0048 nearby; "
                             "that is the half of a two-part change that flatters us"
                             % (os.path.relpath(path, ROOT), line))

    for _doc in SCAN:
        _txt = open(_doc).read()
        for m in re.finditer(r"0\.002\s*(?:to|--|-)\s*0\.012", _txt):
            window = _txt[max(0, m.start() - 400):m.end() + 400]
            if "0.029" not in window:
                line = _txt[:m.start()].count("\n") + 1
                fails.append("%s:%d states the published gaps as 0.002-0.012 without naming the "
                             "0.029 one nearby -- that phrasing is what made the claim wrong"
                             % (os.path.relpath(_doc, ROOT), line))
    # Lines inside a Supplementary table that check_nc_si compares cell by cell with its result
    # file are exempt: a stale value there fails that comparison, and a current value can equal an
    # old token by coincidence (amendment A2 gave DeepMISL 0.5012 and stomach's D +0.0318).
    _bound = set()
    _sp = os.path.join(ROOT, "paper", "nc-submission", "supplementary.tex")
    _ss = open(_sp).read()
    for _lab in set(re.findall(r'expect\("(tab:[^"]+)"', open(os.path.abspath(__file__)).read())):
        _m = re.search(r"\\begin\{tabular\}\{[^}]*\}((?:(?!\\begin\{tabular\}).)*?)\\end\{tabular\}"
                       r"(?:(?!\\begin\{table\}).)*?\\label\{%s\}" % re.escape(_lab), _ss, re.S)
        if _m:
            _bound.update(range(_ss.count("\n", 0, _m.start(1)) + 1, _ss.count("\n", 0, _m.end(1)) + 2))
    for token, why in SUPERSEDED.items():
        for path in SCAN:
            body = open(path).read()
            for i, line in enumerate(body.split("\n"), 1):
                if path == _sp and i in _bound:
                    continue
                if token in line:
                    fails.append("SUPERSEDED VALUE %-9s at %s:%d -- %s"
                                 % (token, os.path.relpath(path, ROOT), i, why))

    # ---- THE CHIMERA EMBARGO. The second cohort's terms of use forbid publishing results obtained
    #      from it until the organisers' challenge journal paper and baseline journal paper are
    #      published (grand-challenge rules and data page, fetched 2026-09-11; only arXiv 2609.09510
    #      existed). The paragraph that carried its composition and two-arm results was therefore
    #      removed, and this gate replaces the ones that bound those numbers: every measured value
    #      from that cohort, and its name, must be ABSENT from every manuscript until permission is
    #      on record in research-control/chimera-external-20260911/.
    _perm = os.path.join(ROOT, "research-control", "chimera-external-20260911", "permission.json")
    _permitted = os.path.isfile(_perm) and json.load(open(_perm)).get("publication_permitted") is True
    if not _permitted:
        # The tokens themselves are embargoed values, so they live beside the protocol in
        # research-control/, which no public artifact contains, and not in this shipped file.
        _tokp = os.path.join(ROOT, "research-control", "chimera-external-20260911", "embargo-tokens.json")
        _EMB = json.load(open(_tokp))["tokens"] if os.path.isfile(_tokp) else ["CHIMERA"]
        # the BiB cover letter (paper/cover_letter.tex) was dropped from this scan on 2026-09-27
        # (RETIRED_DOCS); the NC cover letter is already in SCAN.
        for _doc in SCAN:
            if not os.path.isfile(_doc):
                continue
            _b = open(_doc).read()
            for _tok in _EMB:
                if _tok in _b:
                    fails.append("EMBARGO: %r from the CHIMERA cohort appears in %s"
                                 % (_tok, os.path.relpath(_doc, ROOT)))
        print("  ok  %-38s no CHIMERA name or value in %d documents"
              % ("CHIMERA embargo", len(SCAN)))

    # ---- the encoder-parity arm, read from its own result file. Its known-answer check is what
    #      licenses the number at all, so the gate refuses the prose if that check did not pass.
    _ep = os.path.join(ROOT, "experiments", "20260820-encoder-parity", "results",
                       "encoder-parity.json")
    if os.path.isfile(_ep):
        _e = json.load(open(_ep))
        if not _e["known_answer_check"]["reproduced"]:
            fails.append("the encoder-parity run's known-answer check did not reproduce, so its "
                         "number must not appear in the manuscript")
        # the membership check below read only bib-submission/main.tex; retired with the BiB
        # build on 2026-09-27 (RETIRED_DOCS). The known-answer check above is unaffected.
        print("  ok  encoder parity on CHIEF               full %.4f, slide %.4f, known-answer "
              "reproduced" % (_e["chief"]["full_method"], _e["chief"]["slide_alone"]))

    # ---- any prose that claims how many quantities the check script recomputes must agree with
    #      what the script reports. This drifted twice: the repository README said 27 against 42,
    #      and so did the cover letter, which is the copy an editor reads. The builder fixed the
    #      README by generating the number; prose that a person writes needs a gate instead.
    _vr = os.path.join(ROOT, "verify_release.py")
    if os.path.isfile(_vr):
        _r = subprocess.run([sys.executable, _vr], cwd=ROOT, capture_output=True, text=True)
        _m = re.search(r"(\d+) checks, (\d+) passed", _r.stdout or "")
        if _m:
            _n = _m.group(1)
            # the BiB cover letter (paper/cover_letter.tex) was dropped here on 2026-09-27
            # (RETIRED_DOCS); the NC cover letter is already in SCAN.
            for _p in list(SCAN):
                if not os.path.isfile(_p):
                    continue
                for _q in re.finditer(r"recomputes\s+(\d+)\s+(?:reported\s+)?quantit", open(_p).read()):
                    if _q.group(1) != _n:
                        fails.append("%s claims the check script recomputes %s quantities; it "
                                     "reports %s" % (os.path.relpath(_p, ROOT), _q.group(1), _n))
            print("  ok  recomputed-quantity count           prose agrees with the script at %s" % _n)

    # ---- the published identifiers. This bound release-identifiers.json's DOIs and repository URL,
    # and its unresolved-placeholder scan, against bib-submission/main.tex alone -- retired with the
    # BiB build on 2026-09-27 (RETIRED_DOCS).
    print("  --  %-38s retired with the BiB build" % "published identifiers")

    n_nc = check_nc(fails)
    print("  ok  %-38s %d quantities bound to their result files" % ("Nature Communications build", n_nc))

    print()
    if fails:
        print("FAIL (%d)" % len(fails))
        for f in fails:
            print("   " + f)
        return 1
    # len(CHECKS) is no longer meaningful here: those 29 literal bindings are retired with the
    # BiB build (see the "CHECKS table" line above) rather than verified this run. n_nc is what
    # this run actually verified end to end.
    print("PASS -- %d numbers verified, no superseded value present" % n_nc)
    return 0


NC_TEX = os.path.join(ROOT, "paper", "nc-submission", "main.tex")


def _written(value, literal):
    """Format a source value the way the manuscript wrote it: sign, decimals, thousands comma,
    percent (a literal ending in %% is the value times 100) and millions (a literal ending in M)."""
    lit = literal
    if lit.endswith("M"):
        return "%.1fM" % (value / 1e6)
    if lit.endswith("%"):
        return "%.1f%%" % (value * 100)
    if "," in lit and "." not in lit:
        return "{:,}".format(int(round(value)))
    body = lit.lstrip("+-")
    dp = len(body.split(".")[1]) if "." in body else 0
    s = "%.*f" % (dp, abs(value))
    if lit.startswith("+"):
        return ("+" if value >= 0 else "-") + s
    if lit.startswith("-") or value < 0:
        return "-" + s
    return s


def check_nc(fails):
    """Every number the Nature Communications build states, bound to the field it came from.

    Written 2026-09-11 when the NC manuscript was assembled. Each row is (label, source value,
    the number as written, a context the number must sit in). The context is copied from the
    manuscript, so a number cannot drift away from the thing it describes, and deleting the
    sentence is a failure rather than a silent pass. A missing source field is a failure too."""
    if not os.path.isfile(NC_TEX):
        fails.append("the NC manuscript is missing")
        return 0
    nc = " ".join(re.sub(r"^\s*%.*$", "", open(NC_TEX).read(), flags=re.M).split())
    PH = os.path.join(ROOT, "experiments", "20260911-blca-posthoc", "results")
    GX = os.path.join(ROOT, "experiments", "20260911-geo-external", "results")

    def j(*parts):
        return json.load(open(os.path.join(*parts)))

    a1, pp, pc = load("amendment-A1-clinical-provenance.json"), load("pibd-parity.json"), \
        load("parameter-counts.json")
    why, rg, cr = load("why-grade-fails.json"), load("regime-evidence.json"), \
        load("cold-panel-round1-measurements.json")
    bio = load("biology.json")["B5_within_molecular_subtype"]
    b3 = load("biology.json")["B3_what_the_image_arm_tracks"]["axes"]
    b4 = load("biology.json")["B4_what_the_omics_arm_tracks"]["axes"]
    c2 = load5("c2.json")["C5_gated_fusion"]
    proto = json.load(open(os.path.join(ROOT, "development", "benchmark-protocol.json")))
    unc = proto["metrics"]["uncertainty"]
    cen = json.load(open(os.path.join(ROOT, "development", "clinical-baseline-census.json")))
    cenp = {p["id"]: p for p in cen["papers"]}
    sn = j(ROOT, "experiments", "20260818-selection-null", "results", "selection-null.json")
    ep = j(ROOT, "experiments", "20260820-encoder-parity", "results", "encoder-parity.json")
    uf = j(PH, "unified-fusion-and-added-value.json")
    fu, av = uf["fusion_on_identical_inputs"], uf["added_value_inflation"]["constructions"]
    rs = j(PH, "resplit-and-site-cv.json")
    cal = j(PH, "calibration-and-decision-curve.json")
    ui = j(PH, "utility-intervals.json")
    am = j(PH, "modality-atlas-amended.json")
    dt = j(PH, "double-tied-probe.json")
    geo = j(GX, "geo-external.json")["cohorts"]
    holm = pp["holm_family_of_six"]["comparisons"]
    sgv = cr["stage_vs_grade_inference"]["paired_case_bootstrap"]
    tie = am["conditional_probe"]["clinically_tied_q05"]
    dtt = dt["transcriptome_and_clinically_tied_q40"]["amended_clinical_block"]
    g32, g31, g48 = geo["GSE32894"], geo["GSE31684"], geo["GSE48075"]
    bh = why["cohorts"]["blca"]
    abl = load("ablation-generalisation.json")["B_ablation_bladder"]["paired_vs_full"]
    # the five-study validation, promoted from replication by the dated charter amendment
    fc = j(ROOT, "experiments", "20260912-five-cohort-intervals", "results",
           "five-cohort-intervals.json")
    fcc, fcs = fc["cohorts"], fc["summary"]
    sg = j(ROOT, "experiments", "20260912-subgroup-fairness", "results", "subgroup-fairness.json")
    fl = j(ROOT, "experiments", "20260912-five-cohort-floors", "results",
           "five-cohort-floors.json")["cohorts"]
    enc = j(ROOT, "experiments", "20260913-encoder-sensitivity", "results",
            "five-cohort-gigassl.json")["cohorts"]
    fus = j(ROOT, "experiments", "20260913-five-cohort-fusion", "results",
            "five-cohort-fusion.json")["cohorts"]
    # the 2026-09-27 analyses (experiments/20260927-field-inflation): four further architectures,
    # three further fusions, and the inflation D against a score without information
    FI = os.path.join(ROOT, "experiments", "20260927-field-inflation", "analysis-results")
    fa = j(FI, "field-inflation-a.json")["models"]
    fbd = j(FI, "five-study-extensions-bd.json")
    nb = j(FI, "inflation-null-blca.json")
    ng = j(FI, "inflation-null-geo.json")["constructions"]
    n5 = j(FI, "inflation-null-five.json")["constructions"]
    nbc = nb["constructions"]
    ARCH = ("coattn", "abmil_wsi_pathways", "transmil_wsi_pathways", "deepmisl_wsi_pathways")
    if fbd["summary"]["D_comparisons"] != 15 or fbd["summary"]["D_fitted_fusion_ahead_interval_excluding_zero"] \
            or fbd["summary"]["D_rank_average_ahead_interval_excluding_zero"]:
        fails.append("NC says 15 further fusion comparisons and none separating; the run says %s"
                     % fbd["summary"])
    if not ("66 simplex points" in fbd["grid"] and "step 0.1" in fbd["grid"]):
        fails.append("NC describes a 66-point grid of step 0.1; the run says %r" % fbd["grid"])
    if (nb["no_signal"]["draws_full_data"], nb["no_signal"]["draws_per_resample"],
            nb["bootstrap"]["replicates"]) != (2000, 20, 6000):
        fails.append("NC describes 2,000 draws, 20 per resample and 6,000 resamples; the run says %s"
                     % nb["no_signal"])
    # which excess intervals exclude zero, stated in the Discussion and SI Note 8
    _ex = {k: v["excess_over_no_signal"]["ci95"][0] > 0
           for k, v in list(nbc.items()) + list(ng.items()) + [("five_" + c, v) for c, v in n5.items()]}
    _want_ex = {"ModRank": True, "SurvPath": False, "PIBD": True, "coattn": False,
                "abmil_wsi_pathways": False, "transmil_wsi_pathways": False,
                "deepmisl_wsi_pathways": False, "GSE32894": True, "GSE31684": False,
                "five_blca": True, "five_hnsc": False, "five_stad": False}
    if _ex != _want_ex:
        fails.append("NC names where the excess over a score without information is resolved; "
                     "the runs give %s" % _ex)
    # the architectures: D positive with an interval excluding zero, no added value over stage
    # resolved, and none beyond the benchmark (SI Note 8)
    # Since amendment A2 the four split: three with D resolved and no excess, and DeepMISL, whose
    # risks were constant within most validation folds, with D unresolved and an excess below zero.
    ARCH3 = ARCH[:3]
    for _m in ARCH:
        _v = fa[_m]
        _ok = _v["added_over_stage_boot"]["ci95"][0] <= 0 and _v["added_over_stage_boot"]["ci95"][1] >= 0
        if _m in ARCH3:
            _ok = _ok and _v["D"] > 0 and _v["D_boot"]["ci95"][0] > 0 \
                and nbc[_m]["excess_over_no_signal"]["ci95"][0] <= 0 \
                and _v["retraining_fold_pairs_with_constant_risk"] == 0
        else:
            _ok = _ok and _v["D_boot"]["ci95"][0] <= 0 <= _v["D_boot"]["ci95"][1] \
                and nbc[_m]["excess_over_no_signal"]["ci95"][1] < 0 \
                and (_v["retraining_fold_pairs"], _v["retraining_fold_pairs_with_constant_risk"],
                     _v["retraining_fold_pairs_with_fewer_than_half_distinct_risks"]) == (25, 13, 24)
        if not _ok:
            fails.append("NC's description of %s in SI Note 8 no longer matches the run" % _m)
    if "%.3f" % max(fa[_m]["D_p_holm"] for _m in ARCH3) != "0.022":
        fails.append("SI Note 8 says the largest Holm q over the three resolved architectures is 0.022")
    _dm = fa["deepmisl_wsi_pathways"]
    for _ph in ("DeepMISL returned one risk for the whole validation fold in 13 of its 25 retraining--fold pairs",
                "Its $D$ of $%+.4f$ ($[%+.4f, %+.4f]$)" % ((_dm["D"],) + tuple(_dm["D_boot"]["ci95"]))):
        if _ph not in " ".join(open(os.path.join(ROOT, "paper", "nc-submission", "supplementary.tex"))
                               .read().split()):
            fails.append("SI Note 8 should carry %r" % _ph)
    # analysis C: the reference gap in external cohorts (Results, Discussion, SI Note 9)
    _WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven",
              8: "eight", 9: "nine", 10: "ten"}
    gs = j(FI, "reference-gap-summary.json")
    gx = j(FI, "reference-gap-external.json")
    _ext_n = sum(v["n"] for v in gx["known_answers_geo"].values()) + sum(
        v["n"] for v in gx["cohorts"].values() if v.get("admitted"))
    _big = max(gs["rows"], key=lambda r_: r_["gap"])
    _g0 = j(FI, "reference-gap-gate0.json")["cohorts"][_big["cohort"]]["grade_levels"]
    _share = 100.0 * max(_g0.values()) / sum(_g0.values())
    for _ph in ("all six independent bladder cohorts that record both (%s patients)" % format(_ext_n, ","),
                "with intervals excluding zero in %s" % _WORDS[gs["external_only"]["interval_excludes_zero"]],
                "The largest gap, $%+.4f$, was in a muscle-invasive cohort" % _big["gap"],
                "%.1f\\%% of tumours share one grade" % _share,
                "(Spearman $%.2f$)" % gs["spearman_gap_vs_grade_entropy"]["rho"],
                "outperformed grade in all %s" % _WORDS[gs["gap_positive"]]):
        if _ph.replace("\\%", "\\%") not in " ".join(nc.split()):
            fails.append("NC's external reference-gap sentence should carry %r" % _ph)
    if not (gs["external_only"]["n"] == 6 and gs["external_only"]["gap_positive"] == 6 and gs["n"] == 9):
        fails.append("NC says six external cohorts, all positive, nine in all; the summary says %s" % gs)
    # cold panel round 3. Two statements survived amendment A2 while the results under them moved,
    # because nothing bound their DIRECTION. Each is now derived from its result file: the screened
    # variant's position against the stage block and ModRank, and the best value published on the
    # released folds outside bladder (PIBD ships the same split files) against ModRank.
    _gg = j(ROOT, "experiments", "20260911-geo-external", "results", "geo-gated.json")["cohorts"]
    _g48 = _gg["GSE48075"]
    _gx48 = geo["GSE48075"]["cv"]
    _dir = ("fell below both the stage block and ModRank"
            if _g48["cv"]["gated"]["mean"] < min(_gx48["clinical_stage"]["mean"], _gx48["modrank"]["mean"])
            else None)
    if _dir is None:
        fails.append("NC says the screened variant fell below the stage block and ModRank in GSE48075; "
                     "geo-gated.json no longer says so")
    _abl = load("ablation-generalisation.json")["A_generalisation"]["cohorts"]
    _s5p = j(ROOT, "development", "s5-results", "stage5.json")["cohorts"]
    _above = [c_ for c_ in ("blca", "brca", "coadread", "hnsc", "stad")
              if _s5p[c_]["OURS_wsi_omics_age_sex_stage"] > _abl[c_]["best_published_verified_folds"]]
    if _above != ["blca"]:
        fails.append("NC says ModRank is above the best value published on the released folds only in "
                     "bladder; the files say %s" % _above)
    _mp = j(ROOT, "experiments", "20260818-reporting-dump", "results", "metric-parity.json")
    for _ph in ("%s (%.4f against %.4f and %.4f)" % (_dir or "?", _g48["cv"]["gated"]["mean"],
                                                   _gx48["clinical_stage"]["mean"], _gx48["modrank"]["mean"]),
                "(%.3f in breast, %.3f in colorectal, %.3f in head and neck and %.3f in stomach)"
                % tuple(_abl[c_]["best_published_verified_folds"] for c_ in ("brca", "coadread", "hnsc", "stad")),
                "the two conventions give %.4f and %.4f" % (_mp["mean_of_per_fold"], _mp["pooled_out_of_fold"])):
        if _ph not in " ".join(nc.split()):
            fails.append("NC should carry %r" % _ph)
    # amendment A2 (Methods, SI Note 10): the chain frozen -> A1 -> A2 is read from the before-after
    # table, whose before column is the last commit with pre-A2 results, and from the A1 file
    _ba = j(ROOT, "experiments", "20260928-amendment-a2", "results", "before-after.json")
    _ba_p = [r_ for r_ in _ba["rows"] if r_["quantity"] == "ModRank, mean over five seeds"][0]
    if abs(_ba_p["after"] - a1["primary"]["value"]) > 5e-5:
        fails.append("the before-after table's A2 primary is not the A1 file's primary")
    for _ph in ("moved the primary from %.4f to %.4f" % (a1["frozen_run_for_comparison"]["primary"], _ba_p["before"]),
                "which moved the primary to %.4f" % _ba_p["after"],
                "tied scores share their average rank"):
        if _ph not in " ".join(nc.split()):
            fails.append("NC's amendment sentences should carry %r" % _ph)
    # ten comparisons, none separating. The manuscript rests a contrast on that, so it is asserted
    # rather than described: if any interval ever excludes zero the sentence has to change.
    _sepf = [(c, k) for c, v in fus.items()
             for k in ("ModRank_minus_concatenated", "ModRank_minus_stacked")
             if v[k]["ci95"][0] > 0 or v[k]["ci95"][1] < 0]
    if _sepf:
        fails.append("NC says none of the ten fusion comparisons separates from zero; these do: %s"
                     % _sepf)
    # the claim is that the SAME three separate and the SAME two do not, so both halves are asserted
    _sep = {k for k, v in enc.items()
            if v["differences"]["ours_minus_clinical_stage"]["ci95"][0] > 0}
    if _sep != {"blca", "brca", "hnsc"}:
        fails.append("NC says the same three studies separate under the second encoder; the run "
                     "gives %s" % sorted(_sep))
    # every study's paired difference is positive in all 24 re-partitions; the manuscript says so in
    # two places and the claim is worth an assertion rather than a bound literal
    if any(v["floor"]["partitions_with_a_positive_delta"] != 24 for v in fl.values()):
        fails.append("NC says the difference is positive in 24 of 24 partitions in every study; "
                     "the run says %s" % {k: v["floor"]["partitions_with_a_positive_delta"]
                                          for k, v in fl.items()})

    # the Discussion says both the sex gap and the age split have intervals covering zero
    for _k, _v in (("sex", sg["splits"]["sex"]["ours_men_minus_women"]),
                   ("age", sg["splits"]["age"]["ours_at_minus_above"])):
        if not _v["ci95"][0] < 0 < _v["ci95"][1]:
            fails.append("NC says the %s difference's interval covers zero; it is %s" % (_k, _v["ci95"]))

    # published entries below the stage-based clinical arm: stated in Results and, since
    # 2026-09-27, in the Discussion as the consequence of the correction. Neither was bound before.
    _pub = load("published-benchmark-table.json")["entries"]
    _ver = [e_ for e_ in _pub if e_["folds"].startswith(("verified", "definitional", "as reported in SurvPath"))]
    _stage = a1["arms_seed0"]["clinical_age_sex_stage_from_DIMAF_file"]
    _nv, _na = sum(e_["cindex"] < _stage for e_ in _ver), sum(e_["cindex"] < _stage for e_ in _pub)
    if (len(_ver), _nv, len(_pub), _na) != (4, 2, 9, 5):
        fails.append("NC says 2 of 4 verified and 5 of 9 published entries fall below the stage arm; "
                     "the table gives %d of %d and %d of %d" % (_nv, len(_ver), _na, len(_pub)))
    for _ph in ("exceeds two of the four published entries on verified folds and five of all nine",
                "Two of the four published entries on verified folds, and five of all nine,"):
        if _ph not in " ".join(nc.split()):
            fails.append("NC no longer carries %r" % _ph)

    # the one-column swap's range across the five studies, stated in Results as point estimates,
    # the same values SI Table 2 prints (stage5.json). It said +0.2212, the bootstrap mean from the
    # interval run, beside a table printing +0.2208, until 2026-09-27.
    _s5 = j(ROOT, "development", "s5-results", "stage5.json")["cohorts"]
    _sw = sorted(v["stage_minus_grade_same_pipeline"] for v in _s5.values())
    _want = "by $+%.4f$ to $+%.4f$" % (_sw[0], _sw[-1])
    if len(_sw) != 5 or _want not in " ".join(nc.split()):
        fails.append("NC's one-column-swap range should read %r (stage5.json)" % _want)

    # the second encoder covers all but nine of the five studies' patients (breast 8, stomach 1)
    _gg = j(ROOT, "experiments", "20260913-encoder-sensitivity", "results", "five-cohort-gigassl.json")["cohorts"]
    if fcs["total_cases"] - sum(v["n_cases"] for v in _gg.values()) != 9 or \
            "(all but nine of the 2,241)" not in " ".join(nc.split()):
        fails.append("NC says the second encoder covers all but nine of 2,241 patients; the run covers %d"
                     % sum(v["n_cases"] for v in _gg.values()))

    # Figure 1b's legend states what the transcriptome glyph lights (2026-09-27); the glyph itself
    # recomputes and asserts the same count in paper/figures/pathway_glyph.py
    _b3 = j(ROOT, "experiments", "20260817-blca-confirm", "results", "biology.json")
    _ph = "lights the %d of %d pathway groups" % (
        _b3["B3_what_the_image_arm_tracks"]["pathways_with_BH_q_below_0.05"],
        _b3["B2_survival_association"]["pathways_tested"])
    if _ph not in " ".join(nc.split()):
        fails.append("NC Figure 1 legend should say %r (biology.json)" % _ph)

    def fcd(cohort, key):
        return fcc[cohort]["differences"][key]

    ROWS = [
        # the clinical reference
        ("grade modal share", bh["grade"]["modal_share"], "94.4%",
         "two levels with 94.4\\% of cases in one of them"),
        ("grade entropy", bh["grade"]["normalised_entropy"], "0.311", "(normalised entropy 0.311)"),
        ("stage entropy", bh["stage"]["normalised_entropy"], "0.665", "five levels and 0.665 for stage"),
        ("DIMAF's reported clinical C", cenp["dimaf"]["blca_values"]["clinical"], "0.519",
         "bladder values are 0.519 and"),
        ("SurvPath/MMP reported clinical C", cenp["survpath"]["blca_values"]["age_sex_grade"], "0.570",
         "and 0.570 concordance"),
        # the abstract, rounded to three places
        ("abstract: grade reference", a1["arms_seed0"]["clinical_age_sex_GRADE_from_DIMAF_file"], '0.566',
         'raised its concordance from 0.566 to 0.664'),
        ("abstract: stage reference", a1["arms_seed0"]["clinical_age_sex_stage_from_DIMAF_file"], "0.664",
         'raised its concordance from 0.566 to 0.664'),
        ("abstract: primary", a1["primary"]["value"], "0.721", "It reached 0.721, the highest"),
        ("clinical with grade", a1["arms_seed0"]["clinical_age_sex_GRADE_from_DIMAF_file"], '0.5664',
         'The grade construction reaches 0.5664'),
        ("clinical with stage", a1["arms_seed0"]["clinical_age_sex_stage_from_DIMAF_file"], '0.6637',
         'the stage construction 0.6637'),
        ("stage minus grade", a1["arms_seed0"]["stage_minus_grade"], "+0.0972", "a gap of +0.0972"),
        ("stage minus grade, CI low", sgv["ci95"][0], "+0.047", '$[+0.047, +0.150]$'),
        ("stage minus grade, CI high", sgv["ci95"][1], '+0.150', '$[+0.047, +0.150]$'),
        # the regime
        ("noise fit in-sample", rg["capacity_sweep"]["arms"][0]["corrected_null"], "0.7474",
         "reaches 0.7474 concordance at sixteen parameters"),
        ("fusion oracle", rg["fusion"]["linear_fusion_oracle"]["value"], "0.7268", "(0.7268 against 0.7191)"),
        ("equal weight, deployed", rg["fusion"]["deployed_equal_weight"], "0.7191", "(0.7268 against 0.7191)"),
        ("what the oracle buys", rg["fusion"]["linear_fusion_oracle"]["value"]
         - rg["fusion"]["deployed_equal_weight"], "+0.008", "improves on equal weighting by only $+0.008$"),
        ("gated fusion", c2["gated_on_clinical_tertile"], "0.7493", "(0.7493 against 0.7502)"),
        ("its permuted control", c2["gated_on_PERMUTED_tertile_CONTROL"], "0.7502", "(0.7493 against 0.7502)"),
        ("split-reseed SD", unc["paired_split_reseed_sd"], "0.0073", "a standard deviation of 0.0073"),
        ("split-reseed bar", unc["reseed_bar_2sd"], "0.0145", "twice that, 0.0145, is the margin"),
        # The leave-one-out ablation, added 2026-09-12 after an external review found the subset
        # panel's numbers stated nowhere in the prose. The keys of paired_vs_full name the subset
        # that is KEPT, so "omics + clin" is the slide's removal; getting that backwards is the
        # error this comment exists to prevent.
        # the prose states each removal as the signed change in concordance (since 2026-09-27; it
        # said "costs 0.0218" beside a negative interval before), so the bound quantity is the difference
        ("drop the slide", abl["omics + clin"]["mean"], '-0.0221', 'by $-0.0221$ for the slide'),
        ("drop the slide, CI", abl["omics + clin"]["ci95"][0], '-0.0484', '$[-0.0484, +0.0039]$'),
        ("drop the slide, p", abl["omics + clin"]["p_two_sided"], '0.097', '$p=0.097$'),
        ("drop the transcriptome", abl["wsi + clin"]["mean"], '-0.0224',
         '$-0.0224$ for the transcriptome'),
        ("drop the transcriptome, CI", abl["wsi + clin"]["ci95"][0], '-0.0495', '$[-0.0495, +0.0055]$'),
        ("drop the transcriptome, p", abl["wsi + clin"]["p_two_sided"], '0.107', '$p=0.107$'),
        ("drop the clinical block", abl["wsi + omics"]["mean"], '-0.0390',
         '$-0.0390$ for the clinical block'),
        ("drop the clinical block, CI", abl["wsi + omics"]["ci95"][1], '-0.0086',
         '$[-0.0696, -0.0086]$'),
        ("drop the clinical block, p", abl["wsi + omics"]["p_two_sided"], '0.011', '$p=0.011$'),
        # --- the five studies, promoted to validation by the 2026-09-12 charter amendment. Three of
        # five separate from the corrected clinical arm and two do not, and the rows below bind both
        # halves: a later edit that keeps the three and drops the two fails here.
        ("five studies, cases", fcs["total_cases"], "2,241", "Across 2,241 patients and 411 events"),
        ("five studies, events", fcs["total_events"], "411", "2,241 patients and 411 events"),
        ("blca vs clinical", fcd("blca", "ours_minus_clinical_stage")["mean"], '+0.0784',
         'In bladder the margin is $+0.0784$'),
        ("blca vs clinical, CI low", fcd("blca", "ours_minus_clinical_stage")["ci95"][0], '+0.0298',
         '$[+0.0298, +0.1267]$'),
        ("blca vs clinical, CI high", fcd("blca", "ours_minus_clinical_stage")["ci95"][1], '+0.1267',
         '$[+0.0298, +0.1267]$'),
        ("blca vs clinical, p", fcd("blca", "ours_minus_clinical_stage")["p_two_sided"], "0.001",
         '($[+0.0298, +0.1267]$, $p=0.001$)'),
        ("brca vs clinical", fcd("brca", "ours_minus_clinical_stage")["mean"], '+0.1201',
         'In breast it is $+0.1201$ ($[+0.0519, +0.1897]$, $p<0.001$)'),
        ("brca vs clinical, CI low", fcd("brca", "ours_minus_clinical_stage")["ci95"][0], '+0.0519',
         '$[+0.0519, +0.1897]$'),
        ("brca vs clinical, CI high", fcd("brca", "ours_minus_clinical_stage")["ci95"][1], '+0.1897',
         '$[+0.0519, +0.1897]$'),
        ("hnsc vs clinical", fcd("hnsc", "ours_minus_clinical_stage")["mean"], '+0.0725',
         'in head and neck $+0.0725$ ($[+0.0206, +0.1231]$, $p=0.008$)'),
        ("hnsc vs clinical, CI low", fcd("hnsc", "ours_minus_clinical_stage")["ci95"][0], "+0.0206",
         "$[+0.0206, +0.1231]$"),
        ("hnsc vs clinical, CI high", fcd("hnsc", "ours_minus_clinical_stage")["ci95"][1], "+0.1231",
         "$[+0.0206, +0.1231]$"),
        ("hnsc vs clinical, p", fcd("hnsc", "ours_minus_clinical_stage")["p_two_sided"], '0.008',
         '($[+0.0206, +0.1231]$, $p=0.008$)'),
        # the two that do NOT separate, bound so they cannot quietly leave the paragraph
        ("stad vs clinical", fcd("stad", "ours_minus_clinical_stage")["mean"], '+0.0452',
         'Stomach reaches $+0.0452$'),
        ("stad vs clinical, CI low", fcd("stad", "ours_minus_clinical_stage")["ci95"][0], '-0.0155',
         '$[-0.0155, +0.1070]$'),
        ("coadread vs clinical", fcd("coadread", "ours_minus_clinical_stage")["mean"], '+0.0231',
         'colorectal $+0.0231$ ($[-0.0701, +0.1192]$)'),
        ("coadread vs clinical, CI low", fcd("coadread", "ours_minus_clinical_stage")["ci95"][0],
         '-0.0701', '$[-0.0701, +0.1192]$'),
        ("coadread events", fcc["coadread"]["events"], "37", "the study with 37 events"),
        ("coadread grade control", fcd("coadread", "ours_minus_grade_control")["mean"], '+0.0894',
         'zero, at $+0.0894$ with an interval of $[+0.0437, +0.1429]$'),
        ("coadread grade control, CI low", fcd("coadread", "ours_minus_grade_control")["ci95"][0],
         '+0.0437', '$[+0.0437, +0.1429]$'),
        ("brca three minus two modalities", fcd("brca", "ours_minus_wsi_omics")["mean"], '-0.0232',
         'negative in breast ($-0.0232$, $[-0.0749, +0.0277]$)'),
        ("brca three minus two, CI low", fcd("brca", "ours_minus_wsi_omics")["ci95"][0], '-0.0749',
         '$[-0.0749, +0.0277]$'),
        ("brca best published, any protocol", fcc["brca"]["published_best_any_protocol"][1], "0.794",
         "0.794 for breast and 0.832 for colorectal"),
        ("coadread best published, any protocol",
         fcc["coadread"]["published_best_any_protocol"][1], "0.832",
         "0.794 for breast and 0.832 for colorectal"),
        ("brca ModRank", fcc["brca"]["points"]["OURS_wsi_omics_age_sex_stage"], '0.6881',
         'against 0.6881 and 0.7140 here'),
        ("coadread ModRank", fcc["coadread"]["points"]["OURS_wsi_omics_age_sex_stage"], '0.7140',
         '0.6881 and 0.7140 here'),
        ("blca under the five-study protocol",
         fcc["blca"]["points"]["OURS_wsi_omics_age_sex_stage"], '0.7003',
         'Under it bladder reaches 0.7003, against 0.7214'),
        # --- the same five studies under a second slide encoder. The three that separate are bound
        # here; that the other two still cover zero is asserted below, because an edit that kept the
        # three and quietly dropped the two would otherwise pass.
        ("second encoder, blca", enc["blca"]["differences"]["ours_minus_clinical_stage"]["mean"],
         '+0.0574', '$+0.0574$ ($[+0.0044, +0.1101]$) in bladder'),
        ("second encoder, blca CI low",
         enc["blca"]["differences"]["ours_minus_clinical_stage"]["ci95"][0], '+0.0044',
         '$[+0.0044, +0.1101]$'),
        ("second encoder, brca", enc["brca"]["differences"]["ours_minus_clinical_stage"]["mean"],
         '+0.0849', '$+0.0849$ ($[+0.0163, +0.1531]$) in breast'),
        ("second encoder, brca CI low",
         enc["brca"]["differences"]["ours_minus_clinical_stage"]["ci95"][0], '+0.0163',
         '$[+0.0163, +0.1531]$'),
        ("second encoder, hnsc", enc["hnsc"]["differences"]["ours_minus_clinical_stage"]["mean"],
         '+0.0576', '$+0.0576$ ($[+0.0070, +0.1094]$) in head and neck'),
        ("second encoder, hnsc CI low",
         enc["hnsc"]["differences"]["ours_minus_clinical_stage"]["ci95"][0], '+0.0070',
         '$[+0.0070, +0.1094]$'),
        # the fusion alternatives. Breast is the one the main text names, because it is the study
        # where the rank average trails and the reason the title does not claim otherwise.
        ("breast concatenation", fus["brca"]["points"]["concatenated"], '0.7578',
         'the concatenated model reaches 0.7578 against 0.6881'),
        ("breast, ModRank minus concatenation", fus["brca"]["ModRank_minus_concatenated"]["mean"],
         '-0.0696', '($-0.0696$, $[-0.1399, +0.0010]$, $p=0.054$)'),
        ("breast, that CI low", fus["brca"]["ModRank_minus_concatenated"]["ci95"][0], '-0.1399',
         '$[-0.1399, +0.0010]$'),
        ("breast, that CI high", fus["brca"]["ModRank_minus_concatenated"]["ci95"][1], '+0.0010',
         '$[-0.1399, +0.0010]$'),
        ("breast, that p", fus["brca"]["ModRank_minus_concatenated"]["p_two_sided"], '0.054',
         '$p=0.054$)'),
        # The five "from X to Y" sentences are in the Supplementary Note, so they are checked
        # against that document in check_nc_si. Binding them here was the same error made earlier
        # today with the re-partition floors: every row in this list is matched against main.tex
        # alone, and text that lives in the supplementary fails as "no longer carried".
        # --- subgroup performance, the two TRIPOD rows that read "Not addressed" until 2026-09-12.
        # The men-minus-women difference is bound alongside its interval on purpose: the point
        # estimate alone would read as a finding, and the interval is what says it is not one.
        ("women, ModRank", sg["splits"]["sex"]["women"]["ours"], '0.6638',
         '(0.6638 against 0.7399)'),
        ("men, ModRank", sg["splits"]["sex"]["men"]["ours"], '0.7399', '(0.6638 against 0.7399)'),
        ("women, clinical arm", sg["splits"]["sex"]["women"]["clinical_stage"], '0.5664',
         '(0.5664 against 0.7014)'),
        ("men, clinical arm", sg["splits"]["sex"]["men"]["clinical_stage"], '0.7014',
         '(0.5664 against 0.7014)'),
        ("women in the cohort", sg["splits"]["sex"]["women"]["n"], "90", "with 90\nwomen and 31 events"),
        ("events among women", sg["splits"]["sex"]["women"]["events"], "31", "90\nwomen and 31 events"),
        # Since 2026-09-27 the Discussion states only the interval of the sex gap: the reported
        # difference is the bootstrap mean (+0.0811), which next to 0.7406 and 0.6616 reads as an
        # arithmetic slip (their difference is 0.0790). Full values stay in SI Table 11.
        ("ModRank sex gap, CI low", sg["splits"]["sex"]["ours_men_minus_women"]["ci95"][0], '-0.0244',
         '$[-0.0244, +0.1844]$'),
        ("ModRank sex gap, CI high", sg["splits"]["sex"]["ours_men_minus_women"]["ci95"][1], '+0.1844',
         '$[-0.0244, +0.1844]$'),
        # the age split is now stated as "small and its interval covers zero"; asserted below
        # race: obtained from the archive, and only one stratum clears the floors. The count is
        # bound because it is the reason no comparison is reported, not a descriptive aside.
        ("white patients", sg["splits"]["race"]["white"]["n"], "285",
         "Of the 359 patients, 285 are white"),
        ("white stratum, ModRank over clinical",
         sg["splits"]["race"]["white"]["ours_minus_clinical"]["mean"], '+0.0594',
         '$+0.0594$ ($[+0.0052, +0.1134]$)'),
        ("white stratum, CI low",
         sg["splits"]["race"]["white"]["ours_minus_clinical"]["ci95"][0], '+0.0052',
         '$[+0.0052, +0.1134]$'),
        ("white stratum, CI high",
         sg["splits"]["race"]["white"]["ours_minus_clinical"]["ci95"][1], '+0.1134',
         '$[+0.0052, +0.1134]$'),
        # --- per-study re-partition floors. The two that disagree with the bootstrap are bound
        # explicitly: colorectal fails this bar and stomach clears it while its interval covers
        # zero, and an edit that keeps one and drops the other fails here.
        ("coadread margin", fl["coadread"]["released_folds"]["margin"], '+0.0238',
         '$+0.0238$ against a bar of 0.0267'),
        ("coadread bar", fl["coadread"]["floor"]["bar_2sd"], '0.0267', 'against a bar of 0.0267'),
        ("stad margin", fl["stad"]["released_folds"]["margin"], '+0.0449',
         '$+0.0449$ against 0.0325'),
        ("stad bar", fl["stad"]["floor"]["bar_2sd"], '0.0325', '$+0.0449$ against 0.0325'),
        # The five floors themselves are written in the Supplementary Note, not here, so they are
        # checked by check_nc_si against that document. Binding them in this list was an error:
        # every row here is compared against main.tex alone, and eight rows failed as "no longer
        # carried" when the text had never been there.
        # the confirmatory result
        ("primary", a1["primary"]["value"], '0.7214', 'ModRank reaches 0.7214 concordance'),
        ("primary SD", a1["primary"]["sd_over_seeds"], '0.0043', 'standard deviation 0.0043'),
        ("over the best verified entry", a1["comparisons"]["vs_DIMAF_published"]["gap"], '+0.0424',
         '$+0.0424$ above the best entry whose folds we verified'),
        ("vs SurvPath", holm["ours vs SurvPath+clinical (seed-matched)"]["delta"], '+0.0338',
         'exceeds SurvPath by $+0.0338$ ($[-0.0053, +0.0727]$, $p=0.087$'),
        ("vs SurvPath, CI low", holm["ours vs SurvPath+clinical (seed-matched)"]["ci95"][0], '-0.0053',
         '($[-0.0053, +0.0727]$'),
        ("vs SurvPath, CI high", holm["ours vs SurvPath+clinical (seed-matched)"]["ci95"][1], '+0.0727',
         '($[-0.0053, +0.0727]$'),
        ("vs SurvPath, p", holm["ours vs SurvPath+clinical (seed-matched)"]["p_raw"], '0.087', '$p=0.087$'),
        ("vs PIBD", holm["ours vs PIBD+clinical (best-val checkpoint)"]["delta"], '+0.0252',
         'PIBD by $+0.0252$ ($[-0.0173, +0.0676]$, $p=0.238$)'),
        ("vs PIBD, CI low", holm["ours vs PIBD+clinical (best-val checkpoint)"]["ci95"][0], '-0.0173',
         '($[-0.0173, +0.0676]$'),
        ("vs PIBD, CI high", holm["ours vs PIBD+clinical (best-val checkpoint)"]["ci95"][1], '+0.0676',
         '($[-0.0173, +0.0676]$'),
        ("vs PIBD, p", holm["ours vs PIBD+clinical (best-val checkpoint)"]["p_raw"], '0.238', '$p=0.238$'),
        ("our coefficients", pc["ours"]["parameters"], "1,049", "with 1,049 parameters against 24.7"),
        ("SurvPath parameters", pc["survpath"]["parameters"], "24.7M", "against 24.7 and 26.9 million"),
        ("PIBD parameters", pc["pibd"]["parameters"], "26.9M", "against 24.7 and 26.9 million"),
        # fusion on identical inputs
        ("concatenated, mean", fu["concatenated_ridge_cox"]["mean"], '0.6861', '$0.6861 \\pm 0.0064$'),
        ("concatenated, SD", fu["concatenated_ridge_cox"]["sd"], '0.0064', '$0.6861 \\pm 0.0064$'),
        ("stacked, mean", fu["stacked_learned_weights"]["mean"], '0.7108', '$0.7108 \\pm 0.0032$'),
        ("stacked, SD", fu["stacked_learned_weights"]["sd"], '0.0032', '$0.7108 \\pm 0.0032$'),
        ("vs concatenated", fu["paired"]["ours_vs_concatenated"]["mean"], '+0.0342',
         'exceeds the concatenation by $+0.0342$ ($[+0.0015, +0.0681]$, $p=0.041$)'),
        ("vs concatenated, CI low", fu["paired"]["ours_vs_concatenated"]["ci95"][0], '+0.0015',
         '($[+0.0015, +0.0681]$'),
        ("vs concatenated, CI high", fu["paired"]["ours_vs_concatenated"]["ci95"][1], '+0.0681',
         '($[+0.0015, +0.0681]$'),
        ("vs concatenated, p", fu["paired"]["ours_vs_concatenated"]["p_two_sided"], '0.041', '$p=0.041$'),
        ("vs stacked", fu["paired"]["ours_vs_stacked"]["mean"], '+0.0097',
         'the stacking by $+0.0097$ ($[-0.0067, +0.0263]$)'),
        ("vs stacked, CI low", fu["paired"]["ours_vs_stacked"]["ci95"][0], '-0.0067', '($[-0.0067, +0.0263]$)'),
        ("vs stacked, CI high", fu["paired"]["ours_vs_stacked"]["ci95"][1], '+0.0263', '($[-0.0067, +0.0263]$)'),
        # amendment A2's tuned-ridge stacking (Results, Methods, SI Note 10)
        ("tuned stacking, mean", fu["stacked_tuned_ridge"]["mean"], "0.7147", "$0.7147 \\pm 0.0048$ with a tuned one"),
        ("tuned stacking, SD", fu["stacked_tuned_ridge"]["sd"], "0.0048", "$0.7147 \\pm 0.0048$ with a tuned one"),
        ("vs tuned stacking", fu["paired"]["ours_vs_stacked_tuned_ridge"]["mean"], "+0.0055",
         "the tuned stacking by $+0.0055$ ($[-0.0071, +0.0180]$)"),
        ("vs tuned stacking, CI low", fu["paired"]["ours_vs_stacked_tuned_ridge"]["ci95"][0], "-0.0071",
         "the tuned stacking by $+0.0055$ ($[-0.0071, +0.0180]$)"),
        ("vs tuned stacking, CI high", fu["paired"]["ours_vs_stacked_tuned_ridge"]["ci95"][1], "+0.0180",
         "the tuned stacking by $+0.0055$ ($[-0.0071, +0.0180]$)"),
        ("resplits vs concatenated", rs["resplits"]["ours_minus_concat"]["mean"], '+0.0471',
         '($+0.0471 \\pm 0.0151$)'),
        ("resplits vs concatenated, SD", rs["resplits"]["ours_minus_concat"]["sd"], '0.0151',
         '($+0.0471 \\pm 0.0151$)'),
        ("resplits vs stacked", rs["resplits"]["ours_minus_stacked"]["mean"], '+0.0107',
         '($+0.0107 \\pm 0.0062$)'),
        ("resplits vs stacked, SD", rs["resplits"]["ours_minus_stacked"]["sd"], '0.0062',
         '($+0.0107 \\pm 0.0062$)'),
        ("resplits vs clinical", rs["resplits"]["ours_minus_clinical"]["mean"], "+0.0549",
         '($+0.0549 \\pm 0.0082$) in 24 of 24'),
        ("resplits vs clinical, SD", rs["resplits"]["ours_minus_clinical"]["sd"], '0.0082',
         '($+0.0549 \\pm 0.0082$) in 24 of 24'),
        ("site-grouped sites", rs["site_grouped_cv"]["sites_total"], "33", "(33 sites, five folds)"),
        # the Discussion restates the two re-partition margins and the breast lean in prose (A4,
        # 2026-09-27), so each restatement is bound where it is written, not only where it first appears
        ("Discussion: resplits vs stacked", rs["resplits"]["ours_minus_stacked"]["mean"], '0.0107',
         'the stacking by a mean of 0.0107, which the benchmark cannot resolve'),
        ("Discussion: resplits vs concatenated", rs["resplits"]["ours_minus_concat"]["mean"], '0.0471',
         'the concatenation by 0.0471, which it can'),
        ("Discussion: breast concatenation lead", -fus["brca"]["ModRank_minus_concatenated"]["mean"],
         '0.0696', 'the concatenation leads it by 0.0696 without separating'),
        ("site-grouped ModRank", rs["site_grouped_cv"]["pooled"]["ours"], '0.7204',
         'ModRank holds 0.7204 while the stacking falls to 0.6855 and the concatenation to 0.6482'),
        ("site-grouped stacked", rs["site_grouped_cv"]["pooled"]["stacked"], "0.6855",
         "the stacking falls to 0.6855"),
        ("site-grouped concatenated", rs["site_grouped_cv"]["pooled"]["concat"], '0.6482',
         'the concatenation to 0.6482'),
        ("encoder parity", ep["chief"]["full_method"], '0.7004',
         'SurvPath itself consumed, ModRank reaches 0.7004'),
        # Both values are seed 0, which is the only seed the CHIEF parity run scored, so the
        # difference is stated on that seed rather than against the 5-seed primary. Written as
        # "near 0.022" until 2026-09-12, which rounded 0.0216 up and collided with the slide's
        # leave-one-out cost of 0.0218, a different quantity that happens to look the same.
        ("encoder contribution", a1["arms_seed0"]["OURS"] - ep["chief"]["full_method"], '0.0220',
         "the encoder's contribution at 0.0220 on the seed both values share"),
        ("corrected bar", sn["A_selection_null"]["bars"]["empirical_max_q95"]["bar"], '0.7708',
         'the bar of 0.7708 measured under an outcome-permuted null'),
        # added value and its inflation
        ("ModRank over grade", av["ours"]["added_over_weak"], '+0.1225',
         '$+0.1225$ ($[+0.0695, +0.1733]$) for ModRank'),
        ("ModRank over grade, CI low", av["ours"]["added_over_weak_boot"]["ci95"][0], '+0.0695',
         '($[+0.0695, +0.1733]$)'),
        ("ModRank over grade, CI high", av["ours"]["added_over_weak_boot"]["ci95"][1], '+0.1733',
         '($[+0.0695, +0.1733]$)'),
        ("SurvPath over grade", av["survpath"]["added_over_weak"], '+0.0656', '$+0.0656$ for SurvPath'),
        ("PIBD over grade", av["pibd_best_val"]["added_over_weak"], '+0.0842', '$+0.0842$ for PIBD'),
        ("ModRank over stage", av["ours"]["added_over_stage"], '+0.0588',
         'fall to $+0.0588$, $+0.0249$ and $+0.0336$'),
        ("SurvPath over stage", av["survpath"]["added_over_stage"], '+0.0249',
         'fall to $+0.0588$, $+0.0249$ and $+0.0336$'),
        ("PIBD over stage", av["pibd_best_val"]["added_over_stage"], '+0.0336',
         'fall to $+0.0588$, $+0.0249$ and $+0.0336$'),
        ("D ModRank", av["ours"]["D"], '+0.0637', 'so $D$ is $+0.0637$, $+0.0407$ and $+0.0506$'),
        ("D SurvPath", av["survpath"]["D"], '+0.0407', 'so $D$ is $+0.0637$, $+0.0407$ and $+0.0506$'),
        ("D PIBD", av["pibd_best_val"]["D"], '+0.0506', 'so $D$ is $+0.0637$, $+0.0407$ and $+0.0506$'),
        ("no-signal D, two arms", sum(nbc[k]["no_signal_D"]["mean"] for k in ("SurvPath", "PIBD") + ARCH) / 6,
         "0.029", "has $D = +0.029$ on average"),
        ("no-signal D, ModRank's three arms", nbc["ModRank"]["no_signal_D"]["mean"], '+0.0416',
         "$+0.0416$ in ModRank's three-arm form"),
        ("excess ModRank", nbc["ModRank"]["excess_over_no_signal"]["mean"], '+0.0214',
         '$+0.0214$ ($[+0.0050, +0.0388]$) for ModRank'),
        ("excess ModRank, CI low", nbc["ModRank"]["excess_over_no_signal"]["ci95"][0], '+0.0050', '($[+0.0050, +0.0388]$)'),
        ("excess ModRank, CI high", nbc["ModRank"]["excess_over_no_signal"]["ci95"][1], '+0.0388', '($[+0.0050, +0.0388]$)'),
        ("excess PIBD", nbc["PIBD"]["excess_over_no_signal"]["mean"], '+0.0235',
         '$+0.0235$ ($[+0.0072, +0.0404]$) for PIBD'),
        ("excess PIBD, CI low", nbc["PIBD"]["excess_over_no_signal"]["ci95"][0], '+0.0072', '($[+0.0072, +0.0404]$)'),
        ("excess PIBD, CI high", nbc["PIBD"]["excess_over_no_signal"]["ci95"][1], '+0.0404', '($[+0.0072, +0.0404]$)'),
        ("excess SurvPath", nbc["SurvPath"]["excess_over_no_signal"]["mean"], '+0.0104',
         '$+0.0104$ ($[-0.0047, +0.0261]$) for SurvPath'),
        ("excess SurvPath, CI low", nbc["SurvPath"]["excess_over_no_signal"]["ci95"][0], '-0.0047', '($[-0.0047, +0.0261]$)'),
        ("excess SurvPath, CI high", nbc["SurvPath"]["excess_over_no_signal"]["ci95"][1], '+0.0261', '($[-0.0047, +0.0261]$)'),
        ("architectures alone, lowest", min(fa[m]["C_alone"] for m in ARCH), '0.52', 'reach 0.52 to 0.60 on their own'),
        ("architectures alone, highest", max(fa[m]["C_alone"] for m in ARCH), "0.60", 'reach 0.52 to 0.60 on their own'),
        ("architectures D, lowest", min(fa[m]["D"] for m in ARCH[:3]), "+0.030", "Three have $D$ of $+0.030$ to $+0.033$"),
        ("architectures D, highest", max(fa[m]["D"] for m in ARCH[:3]), "+0.033", "Three have $D$ of $+0.030$ to $+0.033$"),
        ("DeepMISL D", fa["deepmisl_wsi_pathways"]["D"], "+0.005", "so its $D$ is $+0.005$"),
        ("DeepMISL constant folds", fa["deepmisl_wsi_pathways"]["retraining_fold_pairs_with_constant_risk"],
         "13", "fold in 13 of its 25 retraining--fold pairs"),
        ("ModRank vs clinical, Holm", holm["ours vs clinical alone"]["p_holm"], '0.067', '(Holm $q=0.067$)'),
        # calibration and utility
        ("mean predicted 2-year risk", cal["ours"]["by_horizon"]["24.0"]["mean_predicted_risk"], "0.342",
         "was 0.342 against an observed 0.341"),
        ("observed 2-year risk", cal["ours"]["by_horizon"]["24.0"]["observed_km_risk"], "0.341",
         "against an observed 0.341"),
        ("calibration slope", cal["ours"]["calibration_slope"], '0.771',
         'the calibration slope was 0.771 (clinical model 0.848)'),
        ("clinical calibration slope", cal["clinical"]["calibration_slope"], '0.848', '(clinical model 0.848)'),
        ("IPA 1 year", cal["ours"]["by_horizon"]["12.0"]["ipa"], "0.098", "was 0.098 at one year"),
        ("IPA 2 years", cal["ours"]["by_horizon"]["24.0"]["ipa"], '0.151', 'and 0.151 at two years'),
        ("clinical IPA 1 year", cal["clinical"]["by_horizon"]["12.0"]["ipa"], "0.047", "against 0.047 and 0.079"),
        ("clinical IPA 2 years", cal["clinical"]["by_horizon"]["24.0"]["ipa"], "0.079", "against 0.047 and 0.079"),
        ("IPA 3 years", cal["ours"]["by_horizon"]["36.0"]["ipa"], "0.105", "equal at three years (0.105)"),
        ("clinical IPA 3 years", cal["clinical"]["by_horizon"]["36.0"]["ipa"], "0.105",
         "equal at three years (0.105)"),
        ("net benefit at 40%", cal["ours"]["net_benefit_24m"]["0.4"], '0.138', '0.138 against 0.090'),
        ("clinical net benefit at 40%", cal["clinical"]["net_benefit_24m"]["0.4"], "0.090",
         '0.138 against 0.090'),
        ("treat-all net benefit at 40%", cal["treat_all_net_benefit_24m"]["0.4"], "-0.098",
         "and $-0.098$ at 40\\%"),
        ("IPA difference at 2 years", ui["ipa_modrank_minus_clinical"]["24.0"]["mean"], "+0.0724",
         'the two-year difference of $+0.0724$ has a 95\\% interval of $[-0.0056, +0.1504]$'),
        ("its CI low", ui["ipa_modrank_minus_clinical"]["24.0"]["ci95"][0], '-0.0056',
         '$[-0.0056, +0.1504]$'),
        ("its CI high", ui["ipa_modrank_minus_clinical"]["24.0"]["ci95"][1], '+0.1504',
         '$[-0.0056, +0.1504]$'),
        ("net benefit difference at 40%", ui["net_benefit_24m_modrank_minus_clinical"]["0.4"]["mean"],
         '+0.0474', 'At 40\\% the difference from the clinical model is $+0.0474$ ($[+0.0040, +0.0907]$)'),
        ("its CI low", ui["net_benefit_24m_modrank_minus_clinical"]["0.4"]["ci95"][0], '+0.0040',
         '($[+0.0040, +0.0907]$)'),
        ("its CI high", ui["net_benefit_24m_modrank_minus_clinical"]["0.4"]["ci95"][1], '+0.0907',
         '($[+0.0040, +0.0907]$)'),
        # the modalities
        ("slide-clinical Spearman", am["redundancy"]["spearman_titan_vs_clinical"], '0.2443',
         'weakly correlated (Spearman 0.2443)'),
        ("incumbent-clinical Spearman", am["redundancy"]["spearman_survpath_vs_clinical"], '0.1927',
         'and the clinical score (0.1927)'),
        ("clinically tied pairs", tie["pairs"], '1,211', '(1,211 of 24,219)'),
        ("comparable pairs", tie["pairs_total"], "24,219", '(1,211 of 24,219)'),
        ("clinical arm on the tied pairs", tie["arms"]["clinical"], '0.4975', 'at chance (0.4975)'),
        ("slide arm on the tied pairs", tie["arms"]["titan"], '0.6214', 'the slide arm reaches 0.6214'),
        ("doubly tied pairs", dtt["pairs"], '3,892', 'On the 3,892 pairs that neither'),
        ("slide on the doubly tied pairs", dtt["slide"], '0.6193', 'the slide arm reaches 0.6193'),
        ("transcriptome on them", dtt["omics"], '0.5540', 'sit at 0.5540 and 0.5633'),
        ("clinical on them", dtt["clinical"], '0.5633', 'sit at 0.5540 and 0.5633'),
        ("transcriptome arm vs basal axis", b4["basal"]["spearman"], '0.3601',
         'axes (Spearman 0.3601, $-0.372$ and 0.3673)'),
        ("transcriptome arm vs luminal axis", b4["luminal"]["spearman"], "-0.372",
         'axes (Spearman 0.3601, $-0.372$ and 0.3673)'),
        ("transcriptome arm vs EMT axis", b4["EMT"]["spearman"], '0.3673',
         'axes (Spearman 0.3601, $-0.372$ and 0.3673)'),
        ("slide arm, largest axis correlation", max(abs(v["spearman"]) for v in b3.values()), '0.2236',
         'is at most 0.2236 in magnitude'),
        ("ModRank, basal-leaning half", bio["basal_leaning"]["OURS"], "0.7109", '(0.7109 in the basal-leaning half and 0.6983'),
        ("ModRank, luminal-leaning half", bio["luminal_leaning"]["OURS"], '0.6983', '(0.7109 in the basal-leaning half and 0.6983'),
        # the external cohorts
        ("GSE32894 n", g32["n"], "224", "a Swedish cohort of 224 patients with 25 deaths"),
        ("GSE32894 events", g32["events"], "25", "224 patients with 25 deaths"),
        ("GSE32894 D", g32["D"]["point"], '+0.0586', '$D = +0.0586$ ($[+0.0179, +0.0984]$, $p=0.009$)'),
        ("GSE32894 D, CI low", g32["D"]["ci95"][0], '+0.0179', '($[+0.0179, +0.0984]$'),
        ("GSE32894 D, CI high", g32["D"]["ci95"][1], '+0.0984', '($[+0.0179, +0.0984]$'),
        ("GSE32894 D, p", g32["D"]["p_two_sided"], '0.009', '$p=0.009$'),
        ("GSE31684 n", g31["n"], "93", "cohort of 93 patients with 38 deaths"),
        ("GSE31684 events", g31["events"], "38", "93 patients with 38 deaths"),
        ("GSE32894 no-signal D", ng["GSE32894"]["no_signal_D"]["mean"], '+0.0197', '$+0.0197$ in this cohort'),
        ("GSE32894 excess", ng["GSE32894"]["excess_over_no_signal"]["mean"], '+0.0397',
         'excess was $+0.0397$ ($[+0.0011, +0.0780]$, $p=0.045$)'),
        ("GSE32894 excess, CI low", ng["GSE32894"]["excess_over_no_signal"]["ci95"][0], '+0.0011', '($[+0.0011,'),
        ("GSE32894 excess, CI high", ng["GSE32894"]["excess_over_no_signal"]["ci95"][1], '+0.0780', '+0.0780]$, $p=0.045$)'),
        ("GSE32894 excess, p", ng["GSE32894"]["excess_over_no_signal"]["p_two_sided"], '0.045', '$p=0.045$)'),
        ("GSE31684 D", g31["D"]["point"], '+0.0324', '$D$ was $+0.0324$ against $+0.0149$'),
        ("GSE31684 no-signal D", ng["GSE31684"]["no_signal_D"]["mean"], '+0.0149', '$D$ was $+0.0324$ against $+0.0149$'),
        ("GSE31684 excess", ng["GSE31684"]["excess_over_no_signal"]["mean"], '+0.0093',
         'not resolved ($+0.0093$, $[-0.0211, +0.0390]$)'),
        ("GSE31684 excess, CI low", ng["GSE31684"]["excess_over_no_signal"]["ci95"][0], '-0.0211', '$[-0.0211, +0.0390]$)'),
        ("GSE31684 excess, CI high", ng["GSE31684"]["excess_over_no_signal"]["ci95"][1], '+0.0390', '$[-0.0211, +0.0390]$)'),
        ("GSE31684 D without pre-cystectomy chemotherapy",
         g31["sensitivity_without_pre_cystectomy_chemotherapy"]["D"]["point"], '+0.1081',
         'gave $D = +0.1081$ ($[+0.0404, +0.1757]$)'),
        ("the same, CI low", g31["sensitivity_without_pre_cystectomy_chemotherapy"]["D"]["ci95"][0],
         '+0.0404', '($[+0.0404, +0.1757]$)'),
        ("the same, CI high", g31["sensitivity_without_pre_cystectomy_chemotherapy"]["D"]["ci95"][1],
         '+0.1757', '($[+0.0404, +0.1757]$)'),
        ("patients excluded by that sensitivity", g31["n"]
         - g31["sensitivity_without_pre_cystectomy_chemotherapy"]["n"], "3",
         "excludes the three patients given"),
        ("GSE32894 ModRank", g32["cv"]["modrank"]["mean"], '0.8792',
         '(0.8792 against 0.8672 for the stage block and 0.8748 for stacking)'),
        ("GSE32894 stage block", g32["cv"]["clinical_stage"]["mean"], "0.8672", "0.8672 for the stage block"),
        ("GSE32894 stacked", g32["cv"]["stacked"]["mean"], '0.8748', '0.8748 for stacking'),
        ("GSE32894 vs transcriptome, Holm", g32["paired_vs_modrank"]["transcriptome"]["p_holm"], '0.003',
         'correction ($p=0.003$ and $p=0.006$)'),
        ("GSE32894 vs concatenated, Holm", g32["paired_vs_modrank"]["concatenated"]["p_holm"], "0.006",
         'correction ($p=0.003$ and $p=0.006$)'),
        ("GSE48075 n", g48["n"], "73", "a third cohort of 73 patients"),
        ("GSE48075 ModRank", g48["cv"]["modrank"]["mean"], '0.6389', '(0.6389 against 0.6334)'),
        ("GSE48075 stage block", g48["cv"]["clinical_stage"]["mean"], '0.6334', '(0.6389 against 0.6334)'),
        ("GSE31684 transcriptome", g31["cv"]["transcriptome"]["mean"], '0.4747', 'were uninformative (0.4747)'),
        ("GSE31684 ModRank", g31["cv"]["modrank"]["mean"], '0.5896', '(0.5896 against 0.6591)'),
        ("GSE31684 stage block", g31["cv"]["clinical_stage"]["mean"], '0.6591', '(0.5896 against 0.6591)'),
        # Methods
        ("deployable variant", cr["transductive_vs_inductive"]["primary_inductive"], '0.7245',
         'reaches 0.7245 against the transductive 0.7214'),
        ("deployable variant, p",
         cr["transductive_vs_inductive"]["paired_case_bootstrap_inductive_minus_transductive"]["p_two_sided"],
         '0.59', '(paired $p=0.59$)'),
        ("frozen-run primary", a1["frozen_run_for_comparison"]["primary"], "0.7260",
         "moved the primary from 0.7260 to 0.7212"),
    ]
    n = 0
    for label, value, lit, ctx in ROWS:
        n += 1
        got = _written(float(value), lit)
        want = lit if not lit.endswith("M") else lit
        if got != want:
            fails.append("NC: %s is %s at source, written as %s" % (label, got, lit))
        if " ".join(ctx.split()) not in nc:
            fails.append("NC: the manuscript no longer carries %r (%s)" % (ctx, label))

    # derived claims, recomputed rather than read
    if "the gain over a grade-based reference is two to three times the gain over a stage-based one" in nc:
        for key in ("ours", "survpath", "pibd_best_val"):
            ratio = av[key]["added_over_weak"] / av[key]["added_over_stage"]
            if not 2.0 <= ratio <= 3.0:
                fails.append("NC Discussion says two to three times; for %s it is %.2f" % (key, ratio))
    else:
        fails.append("NC Discussion no longer states the grade-to-stage added-value ratio")
    sc = rs["site_grouped_cv"]["pooled"]
    if "Holding out whole hospitals widens both gaps" in nc and not (
            sc["ours"] - sc["stacked"] > rs["resplits"]["ours_minus_stacked"]["mean"]
            and sc["ours"] - sc["concat"] > rs["resplits"]["ours_minus_concat"]["mean"]):
        fails.append("NC says held-out hospitals widen both gaps; the site-grouped file disagrees")
    if "selected it among 206 candidates" in nc:
        _sr = j(ROOT, "experiments", "20260818-selection-null", "results", "selection-rank.json")
        if _sr["full_cohort_out_of_fold_candidates"] != 206:
            fails.append("NC says 206 candidates; the ledger says %d" % _sr["full_cohort_out_of_fold_candidates"])
    if "it reaches 0.7009, which is still above every published entry" in nc:
        _top = max(e_["cindex"] for e_ in load("published-benchmark-table.json")["entries"])
        if not ep["chief"]["full_method"] > _top:
            fails.append("NC says 0.7009 is above every published entry; the best is %.3f" % _top)
    for key in ("ours", "survpath", "pibd_best_val"):
        if not av[key]["added_over_stage"] < 0.5 * av[key]["added_over_weak"]:
            fails.append("NC says a stage-based reference MORE THAN HALVES the added value; for %s it "
                         "is %.4f against %.4f" % (key, av[key]["added_over_stage"],
                                                   av[key]["added_over_weak"]))
    for key in ("ours_minus_concat", "ours_minus_stacked", "ours_minus_clinical",
                "D_added_value_inflation"):
        if rs["resplits"][key]["fraction_positive"] != 1.0 or rs["resplits"]["partitions"] != 24:
            fails.append("NC says 24 of 24 re-partitions for %s; the file says %s of %d"
                         % (key, rs["resplits"][key]["fraction_positive"], rs["resplits"]["partitions"]))
    th = sorted(float(k) for k in cal["ours"]["net_benefit_24m"])
    ahead = [x for x in th if cal["ours"]["net_benefit_24m"][str(round(x, 2))]
             > max(cal["clinical"]["net_benefit_24m"][str(round(x, 2))],
                   cal["treat_all_net_benefit_24m"][str(round(x, 2))])]
    if "at every threshold examined from 20 to 60\\%" in nc:
        want = [x for x in th if 0.20 - 1e-9 <= x <= 0.60 + 1e-9]
        if not set(want) <= set(ahead):
            fails.append("NC says ModRank's net benefit leads at every threshold from 20 to 60%%; it "
                         "leads only at %s" % ahead)
        if max(th) < 0.60 - 1e-9 or [x for x in ahead if x < 0.20 - 1e-9]:
            fails.append("NC's 20-60%% range no longer matches the grid or the lead: %s" % ahead)
    else:
        fails.append("NC no longer states the net-benefit threshold range")
    _nbx = [th for th, v in ui["net_benefit_24m_modrank_minus_clinical"].items()
            if v["ci95"][0] > 0 or v["ci95"][1] < 0]
    if _nbx != ["0.4"] and "at the other four thresholds tested its interval includes zero" in nc:
        fails.append("NC says only the 40%% net-benefit interval excludes zero; it is %s" % _nbx)
    if any(v["ci95"][0] > 0 for v in ui["ipa_modrank_minus_clinical"].values()):
        fails.append("NC treats every IPA difference as unresolved; one interval now excludes zero")
    for key in ("GSE31684", "GSE32894", "GSE48075"):
        if geo[key]["decision"]["clears_every_comparator"] is not False:
            fails.append("NC says the two-arm rule clears no cohort's bar; %s clears it" % key)

    # the journal's hard limits, measured on the source rather than remembered
    raw = open(NC_TEX).read()
    _t = re.search(r"\\title\{(.*?)\}\n", raw, re.S).group(1)
    if len(_t.split()) > 15:
        fails.append("NC title is %d words; the journal allows 15" % len(_t.split()))
    _ab = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", raw, re.S).group(1)
    _ab = re.sub(r"\\noindent|\\TODO\{[^}]*\}", " ", _ab)
    if len(_ab.split()) > 150:
        fails.append("NC abstract is %d words; the journal allows 150" % len(_ab.split()))
    if re.search(r"\\cite|\b(TCGA|GEO|IPA|DSS|CI|HR)\b", _ab):
        fails.append("NC abstract carries a reference or an abbreviation")
    # prose rules on the NC build, which the lists above were written before
    if re.search(r"---", nc):
        fails.append("NC uses an em-dash")
    for m in re.finditer(r"\\texttt\{([^}]*)\}", nc):
        if "@" not in m.group(1):
            fails.append("NC puts a code identifier in the prose: %r" % m.group(1))
    # CROSS-DOCUMENT POINTERS. The main text and the Supplementary Information are separate
    # compiles, so each typed number is bound to the label it must render as, read from the other
    # document's own .aux. A pointer that no longer appears is a failure, not a skip.
    si_aux = os.path.join(ROOT, "paper", "nc-submission", "supplementary.aux")
    main_aux = os.path.join(ROOT, "paper", "nc-submission", "main.aux")
    if not (os.path.isfile(si_aux) and os.path.isfile(main_aux)):
        fails.append("NC: compile both nc-submission documents before the cross-document check")
    else:
        rs_si = dict(re.findall(r"\\newlabel\{([^}]*)\}\{\{([^}]*)\}", open(si_aux).read()))
        rs_main = dict(re.findall(r"\\newlabel\{([^}]*)\}\{\{([^}]*)\}", open(main_aux).read()))
        EXPECT = {"Supplementary Note~1": "note:census", "Supplementary Table~1": "tab:census",
                  "Supplementary Table~2": "tab:entropy", "Supplementary Note~2": "note:protocol",
                  "Supplementary Table~3": "tab:components", "Supplementary Fig.~1": "fig:campaign",
                  "Supplementary Table~4": "tab:perseed", "Supplementary Table~5": "tab:ties",
                  "Supplementary Note~3": "note:arms", "Supplementary Fig.~2": "fig:landscape",
                  "Supplementary Table~6": "tab:subtype", "Supplementary Fig.~3": "fig:biology",
                  "Supplementary Figs.~4": "fig:cohort", "and~5 draw": "fig:cases",
                  # Everything from here was inserted during 2026-09-12 and 13 and moved the numbers
                  # after it, twice for the tables. This half is written in LABEL order with the
                  # rendered number beside each, because the previous layout was in insertion order
                  # and that is how "Supplementary Table~7" came to appear twice in one dict, once
                  # for the five-study table and once for the fusion table. A dict keeps the last,
                  # so the first pointer stopped being checked and nothing said so. Sorted by label,
                  # a repeated number is visible on the page.
                  "Supplementary Note~4": "note:fivecohort",     # note 4
                  "Supplementary Note~5": "note:geo",            # note 5
                  "Supplementary Note~6": "note:calibration",    # note 6
                  "Supplementary Note~7": "note:subgroup",       # note 7
                  "Supplementary Note~8": "note:nulld",          # note 8, 2026-09-27
                  "Supplementary Note~9": "note:gapext",         # note 9, 2026-09-27
                  "Supplementary Table~7": "tab:fusion",         # table 7
                  "Supplementary Table~8": "tab:fivecohort",     # table 8
                  "Supplementary Table~9": "tab:geo",            # table 9
                  "Supplementary Table~10": "tab:calibration",   # table 10
                  "Supplementary Table~11": "tab:subgroup",      # table 11
                  "Supplementary Tables~12": "tab:architectures",  # table 12, 2026-09-27
                  "and~13)": "tab:nulld",                        # table 13, 2026-09-27
                  "Supplementary Table~14": "tab:gapext",        # table 14, 2026-09-27
                  "Supplementary Note~10": "note:amendA2",       # note 10, amendment A2
                  "Supplementary Table~15": "tab:amendA2",       # table 15, amendment A2
                  "Supplementary Table~16": "tab:tripod"}        # table 16
        # A duplicated KEY cannot be caught by inspecting this dict, because Python collapses it
        # before anything runs: writing "Supplementary Table~7" twice silently drops the first
        # pointer, which is what happened on 2026-09-13. What can be caught is the consequence.
        # Every note and table the Supplementary defines must be pointed at by exactly one phrase,
        # so a collapsed entry shows up as a label nobody checks.
        _si_items = {k for k in rs_si if k.startswith(("note:", "tab:"))}
        _unchecked = _si_items - set(EXPECT.values())
        if _unchecked:
            fails.append("NC: the Supplementary defines %s, which no main-text pointer is checked "
                         "against" % sorted(_unchecked))
        if len(set(EXPECT.values())) != len(EXPECT):
            fails.append("NC: two phrases in the pointer table map to the same label")
        for phrase, label in EXPECT.items():
            want = re.search(r"(\d+)", phrase).group(1)
            if rs_si.get(label) != want:
                fails.append("NC: %s renders as %s in the Supplementary Information, but the main "
                             "text points at %s" % (label, rs_si.get(label), want))
            if phrase not in nc:
                fails.append("NC: the main text no longer carries %r (%s)" % (phrase, label))
        valid = {k.split(":")[0]: set() for k in rs_si}
        for k, v in rs_si.items():
            valid[k.split(":")[0]].add(v)
        for m in re.finditer(r"Supplementary (Table|Note|Figs?\.)~(\d+)", nc):
            kind = {"Table": "tab", "Note": "note"}.get(m.group(1), "fig")
            if m.group(2) not in valid.get(kind, set()):
                fails.append("NC: the main text cites %s, which the Supplementary Information does "
                             "not contain" % m.group(0))
        si = " ".join(open(os.path.join(ROOT, "paper", "nc-submission", "supplementary.tex")).read()
                      .split())
        for phrase, label in (("Figure~%sd", "fig:regime"), ("Figure~%sc", "fig:robust")):
            p_ = phrase % rs_main.get(label, "?")
            if p_ not in si:
                fails.append("NC: the Supplementary Information should point at %r (%s)" % (p_, label))
        for m in re.finditer(r"Figure~(\d+)", si):
            if m.group(1) not in {v for k, v in rs_main.items() if k.startswith("fig:")}:
                fails.append("NC: the Supplementary Information cites main Figure %s, which does "
                             "not exist" % m.group(1))

    n += check_nc_si(fails)

    # FIGURE WIDTH, G45: every figure the NC build includes is 180 mm wide (510.2 pt, tolerance
    # 0.2 mm). Since 2026-09-27 the main text names each figure by its path relative to the .tex
    # (the kit's figure gates do not read \\graphicspath); a graphicspath, if one returns, is still
    # honoured, and the .tex's own directory is always tried first.
    _gpm = re.search(r"\\graphicspath\{((?:\{[^}]*\})+)\}", open(NC_TEX).read())
    _dirs = [os.path.dirname(NC_TEX)] + (
        [os.path.normpath(os.path.join(os.path.dirname(NC_TEX), d_))
         for d_ in re.findall(r"\{([^}]*)\}", _gpm.group(1))] if _gpm else [])
    for _doc in (NC_TEX, os.path.join(ROOT, "paper", "nc-submission", "supplementary.tex")):
        for _m in re.finditer(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", open(_doc).read()):
            if _doc != NC_TEX:
                continue            # the SI carries earlier figures at their own width; checked by eye
            _hit = next((os.path.join(d_, _m.group(1)) for d_ in _dirs
                         if os.path.isfile(os.path.join(d_, _m.group(1)))), None)
            if _hit is None:
                fails.append("NC includes %s, which does not resolve to a file" % _m.group(1))
                continue
            _o = subprocess.run(["pdfinfo", _hit], capture_output=True, text=True).stdout
            _w = float([l_ for l_ in _o.split("\n") if "Page size" in l_][0].split()[2])
            if abs(_w - 180.0 / 25.4 * 72) > 0.2 / 25.4 * 72:
                fails.append("NC figure %s is %.1f pt wide; 180 mm is 510.2 pt" % (os.path.relpath(_hit, ROOT), _w))

    # display items numbered in citation order
    src = open(NC_TEX).read()
    for pre in ("fig:", "tab:"):
        declared = re.findall(r"\\label\{(%s[^}]+)\}" % pre, src)
        first = {}
        for m in re.finditer(r"\\ref\{(%s[^}]+)\}" % pre, src):
            first.setdefault(m.group(1), m.start())
        if [l for l in declared if l not in first]:
            fails.append("NC never cites %s" % [l for l in declared if l not in first])
        if sorted(declared, key=lambda l: first.get(l, 10 ** 9)) != declared:
            fails.append("NC display items are not numbered in citation order")
    return n


def _si_rows(src, label):
    """The rows of the tabular that carries \\label{label}, as lists of cleaned cells."""
    # the body may not contain another tabular: without that, the lazy match began at the FIRST
    # tabular in the document and ran on to this label, so a row was looked up in every earlier
    # table as well (found 2026-09-27, amendment A2)
    m = re.search(r"\\begin\{tabular\}\{[^}]*\}((?:(?!\\begin\{tabular\}).)*?)\\end\{tabular\}"
                  r"(?:(?!\\begin\{table\}).)*?\\label\{%s\}" % re.escape(label), src, re.S)
    if not m:
        return None
    body = re.sub(r"\\(toprule|midrule|bottomrule|addlinespace)|\\cmidrule\([^)]*\)\{[^}]*\}", "", m.group(1))
    rows = []
    for r in body.split("\\\\"):
        r = r.strip()
        if not r or "\\multicolumn{3}{c}{grade}" in r:
            continue
        cells = [re.sub(r"\$|\\textbf\{([^}]*)\}|\\emph\{([^}]*)\}", lambda k: k.group(1) or k.group(2) or "",
                        c).replace("\\%", "%").strip() for c in r.split("&")]
        rows.append(cells)
    return rows


def check_nc_si(fails):
    """The Supplementary Information's tables, each row rebuilt from its result file and compared
    cell by cell. A table that cannot be found is a failure."""
    sp = os.path.join(ROOT, "paper", "nc-submission", "supplementary.tex")
    src = open(sp).read()
    PH = os.path.join(ROOT, "experiments", "20260911-blca-posthoc", "results")
    n = 0

    def expect(label, want):
        nonlocal n
        got = _si_rows(src, label)
        if got is None:
            fails.append("NC SI: no table labelled %s" % label)
            return
        for w in want:
            n += 1
            if w not in got:
                fails.append("NC SI %s: expected the row %s, not found among %s"
                             % (label, w, [g[:3] for g in got]))

    why = load("why-grade-fails.json")["cohorts"]
    NAME = {"blca": "bladder", "brca": "breast", "coadread": "colorectal", "hnsc": "head and neck",
            "stad": "stomach"}
    rows = []
    for k, nm in NAME.items():
        c = why[k]
        g, s = c["grade"], c["stage"]
        stage = [str(s["n_levels"]), "%.1f%%" % (100 * s["modal_share"]), "%.3f" % s["normalised_entropy"]]
        if g is None:
            continue                     # the no-grade rows carry a multicolumn note, checked below
        rows.append([nm, str(c["n"]), str(g["n_levels"]), "%.1f%%" % (100 * g["modal_share"]),
                     "%.3f" % g["normalised_entropy"]] + stage + ["%+.4f" % c["stage_minus_grade"]])
    expect("tab:entropy", rows)
    for k in ("brca", "coadread"):
        c = why[k]
        s = c["stage"]
        frag = "%s & %d & \\multicolumn{3}{c}{no usable grade in the released file} & %d & %.1f\\%% & %.3f & $%+.4f$" % (
            NAME[k], c["n"], s["n_levels"], 100 * s["modal_share"], s["normalised_entropy"],
            c["stage_minus_grade"])
        n += 1
        if frag not in src:
            fails.append("NC SI tab:entropy: the %s row does not read %r" % (k, frag))

    comp = load("component-slate.json")["components"]
    got = _si_rows(src, "tab:components")
    if got is None:
        fails.append("NC SI: no table labelled tab:components")
        got = []
    for c in comp:
        n += 1
        ctl = "n/a" if c["control"] is None else "%+.4f" % c["control"]
        row = [r for r in got if r and r[0] == c["id"].replace("$", "")]
        want = ["%+.4f" % c["predicted"], "%+.4f" % c["observed"], ctl, c["status"]]
        if not row or row[0][2:] != want:
            fails.append("NC SI tab:components: %s reads %s, component-slate.json gives %s"
                         % (c["id"], row[0][2:] if row else None, want))

    a1 = load("amendment-A1-clinical-provenance.json")
    expect("tab:perseed", [[str(r["seed"]), "%.4f" % r["OURS"], "%.4f" % r["wsi_titan"],
                            "%.4f" % r["omics_combine"], "%.4f" % r["clinical"]] for r in a1["per_seed"]]
           + [["mean", "%.4f" % a1["primary"]["value"], "", "", ""],
              ["standard deviation", "%.4f" % a1["primary"]["sd_over_seeds"], "", "", ""]])

    cp = json.load(open(os.path.join(PH, "modality-atlas-amended.json")))["conditional_probe"]
    TIE = [("all comparable pairs", "all_pairs"), ("closest 40% on the clinical model", "clinically_tied_q40"),
           ("closest 20%", "clinically_tied_q20"), ("closest 10%", "clinically_tied_q10"),
           ("closest 5%", "clinically_tied_q05")]
    expect("tab:ties", [[lab, "{:,}".format(cp[k]["pairs"]), "%.4f" % cp[k]["arms"]["clinical"],
                         "%.4f" % cp[k]["arms"]["titan"], "%.4f" % cp[k]["arms"]["survpath"]]
                        for lab, k in TIE])

    b5 = load("biology.json")["B5_within_molecular_subtype"]
    expect("tab:subtype", [["%s-leaning" % h, "%d / %d" % (b5["%s_leaning" % h]["n"], b5["%s_leaning" % h]["events"]),
                            "%.4f" % b5["%s_leaning" % h]["clinical"], "%.4f" % b5["%s_leaning" % h]["image"],
                            "%.4f" % b5["%s_leaning" % h]["omics"], "%.4f" % b5["%s_leaning" % h]["OURS"]]
                           for h in ("basal", "luminal")])

    geo = json.load(open(os.path.join(ROOT, "experiments", "20260911-geo-external", "results",
                                      "geo-external.json")))["cohorts"]
    C = ("GSE31684", "GSE32894", "GSE48075")
    geo_rows = [["patients / deaths from the disease"] + ["%d / %d" % (geo[c]["n"], geo[c]["events"]) for c in C],
                ["pathway groups kept"] + [str(geo[c]["pathways_kept"]) for c in C]]
    for lab, key in (("two-arm ModRank", "modrank"), ("clinical stage block", "clinical_stage"),
                     ("transcriptome", "transcriptome"), ("concatenated", "concatenated"),
                     ("stacked", "stacked")):
        geo_rows.append([lab] + ["%.4f" % geo[c]["cv"][key]["mean"] for c in C])
    geo_rows.append(["bar (twice the re-partition SD)"] + ["%.4f" % geo[c]["floor"]["bar"] for c in C])
    geo_rows.append(["ModRank clears every comparator"] +
                    ["yes" if geo[c]["decision"]["clears_every_comparator"] else "no" for c in C])
    geo_rows.append(["D from a grade-based reference"] +
                    ["%+.4f" % geo[c]["D"]["point"] if geo[c].get("D") else "not defined" for c in C])
    expect("tab:geo", geo_rows)

    # The five-study table was hand-written on 2026-09-12 and every cell was checked against its
    # source by hand, once. That is exactly the check this file exists to replace: a hand check
    # passes at the moment it is made and says nothing about the next edit. Rebuilt here from the
    # run's own output so a changed cell fails.
    fcj = json.load(open(os.path.join(ROOT, "experiments", "20260912-five-cohort-intervals",
                                      "results", "five-cohort-intervals.json")))["cohorts"]
    fc_rows = []
    for key, lab in (("blca", "bladder"), ("brca", "breast"), ("coadread", "colorectal"),
                     ("hnsc", "head and neck"), ("stad", "stomach")):
        v, p, d = fcj[key], fcj[key]["points"], fcj[key]["differences"]
        fc_rows.append([lab, str(v["n_cases"]), str(v["events"]),
                        "%.4f" % p["age_sex_stage"], "%.4f" % p["OURS_wsi_omics_age_sex_stage"],
                        "%+.4f" % d["ours_minus_clinical_stage"]["mean"],
                        "%+.4f" % d["ours_minus_grade_control"]["mean"],
                        "%+.4f" % d["ours_minus_wsi_omics"]["mean"]])
    expect("tab:fivecohort", fc_rows)

    sgj = json.load(open(os.path.join(ROOT, "experiments", "20260912-subgroup-fairness",
                                      "results", "subgroup-fairness.json")))["splits"]
    # Every row of the table, the race row included. It was added by hand on 2026-09-12 with two
    # cells that had never been printed anywhere, and this list did not cover it, so nothing would
    # have caught them. A rebuilt row is only a check for the rows it is given.
    sg_rows = []
    for axis, key, lab in (("sex", "men", "men"), ("sex", "women", "women"),
                           ("age", "at or below", "age at or below 68"),
                           ("age", "above the median", "age above 68"),
                           ("race", "white", "recorded white")):
        k = key if key in sgj[axis] else next(x for x in sgj[axis] if x.startswith(key))
        v = sgj[axis][k]
        sg_rows.append([lab, str(v["n"]), str(v["events"]), "%.4f" % v["clinical_stage"],
                        "%.4f" % v["slide"], "%.4f" % v["omics"], "%.4f" % v["ours"]])
    expect("tab:subgroup", sg_rows)

    fusj = json.load(open(os.path.join(ROOT, "experiments", "20260913-five-cohort-fusion",
                                       "results", "five-cohort-fusion.json")))["cohorts"]
    _FI = os.path.join(ROOT, "experiments", "20260927-field-inflation", "analysis-results")
    _jl = lambda n_: json.load(open(os.path.join(_FI, n_)))
    fbd, fa = _jl("five-study-extensions-bd.json"), _jl("field-inflation-a.json")["models"]
    nbc, ng = _jl("inflation-null-blca.json")["constructions"], _jl("inflation-null-geo.json")["constructions"]
    n5 = _jl("inflation-null-five.json")["constructions"]
    fu_rows = []
    _fd = fbd["fusions_D"]
    for key, lab in (("blca", "bladder"), ("brca", "breast"), ("coadread", "colorectal"),
                     ("hnsc", "head and neck"), ("stad", "stomach")):
        v, w = fusj[key], _fd[key]
        fu_rows.append([lab, "rank average (ModRank)", "%.4f" % v["points"]["ModRank"], "", ""])
        for nm, pt, d in (("concatenated", v["points"]["concatenated"], v["ModRank_minus_concatenated"]),
                          ("stacked", v["points"]["stacked"], v["ModRank_minus_stacked"]),
                          ("simplex weights", w["points"]["simplex"], w["ModRank_minus_simplex"]),
                          ("stacked with interactions", w["points"]["interaction"], w["ModRank_minus_interaction"]),
                          ("weights by clinical tertile", w["points"]["gated"], w["ModRank_minus_gated"]),
                          ("stacked, tuned ridge", v["stacked_tuned_ridge"]["point"],
                           v["stacked_tuned_ridge"]["ModRank_minus_stacked_tuned_ridge"])):
            _m = "$0.0000$" if abs(d["mean"]) < 5e-5 else "$%+.4f$" % d["mean"]
            fu_rows.append(["", nm, "%.4f" % pt, _m, "$[%+.4f, %+.4f]$" % tuple(d["ci95"])])
    expect("tab:fusion", [[c.replace("$", "") for c in r_] for r_ in fu_rows])

    ar_rows = []
    for key, lab, pub in (("coattn", "MCAT", "0.598"), ("abmil_wsi_pathways", "ABMIL with pathways", "0.562"),
                          ("transmil_wsi_pathways", "TransMIL with pathways", "0.630"),
                          ("deepmisl_wsi_pathways", "DeepMISL with pathways", "0.518")):
        v = fa[key]
        ps = v["per_seed_C_alone"]
        mu = sum(ps) / len(ps)
        sdv = (sum((x - mu) ** 2 for x in ps) / (len(ps) - 1)) ** 0.5
        ar_rows.append([lab, "%.4f" % v["C_alone"], "$%.4f \\pm %.4f$" % (mu, sdv), pub,
                        "$%+.4f$" % v["added_over_weak"], "$%+.4f$" % v["added_over_stage"]])
    expect("tab:architectures", [[c.replace("$", "") for c in r_] for r_ in ar_rows])

    nd_rows = []
    for src_, key, lab in ((nbc, "ModRank", "ModRank"), (nbc, "SurvPath", "SurvPath + clinical"),
                           (nbc, "PIBD", "PIBD + clinical"), (nbc, "coattn", "MCAT + clinical"),
                           (nbc, "abmil_wsi_pathways", "ABMIL + clinical"),
                           (nbc, "transmil_wsi_pathways", "TransMIL + clinical"),
                           (nbc, "deepmisl_wsi_pathways", "DeepMISL + clinical"),
                           (n5, "blca", "bladder"), (n5, "hnsc", "head and neck"), (n5, "stad", "stomach"),
                           (ng, "GSE32894", "GSE32894 (Lund)"), (ng, "GSE31684", "GSE31684 (cystectomy)")):
        v = src_[key]
        x, q = v["excess_over_no_signal"], v["no_signal_D"]["q025_q975"]
        nd_rows.append([lab, "$%+.4f$" % v["D"], "$%+.4f$" % v["no_signal_D"]["mean"],
                        "$[%+.4f, %+.4f]$" % tuple(q), "$%+.4f$" % x["point"], "$%+.4f$" % x["mean"],
                        "$[%+.4f, %+.4f]$" % tuple(x["ci95"])])
    expect("tab:nulld", [[c.replace("$", "") for c in r_] for r_ in nd_rows])

    gs = _jl("reference-gap-summary.json")
    gx = _jl("reference-gap-external.json")
    s5c = json.load(open(os.path.join(ROOT, "experiments", "20260911-five-study-repro", "results",
                                      "stage5.json")))["cohorts"]
    SET = {"TCGA-BLCA": ("benchmark", "DSS"), "TCGA-HNSC": ("benchmark", "DSS"), "TCGA-STAD": ("benchmark", "DSS"),
           "GSE31684": ("cystectomy", "DSS"), "GSE32894": ("mixed", "DSS"), "GSE19915": ("mixed", "DSS"),
           "E-MTAB-1803": ("MIBC", "OS"), "GSE13507": ("mixed", "CSS"), "E-MTAB-4321": ("NMIBC", "PFS")}
    gx_rows = []
    for r_ in gs["rows"]:
        c_ = r_["cohort"]
        if c_.startswith("TCGA-"):
            v_ = s5c[c_[5:].lower()]
            n_, ev_, cs_, cg_ = v_["n_cases"], v_["events"], v_["arms"]["age_sex_stage"], v_["arms"]["age_sex_grade"]
        else:
            v_ = gx["known_answers_geo"].get(c_) or gx["cohorts"][c_]
            n_, ev_, cs_, cg_ = v_["n"], v_["events"], v_["C_stage_model"], v_["C_grade_model"]
        gx_rows.append([c_, SET[c_][0], SET[c_][1], "n / events".replace("n", "%d" % n_, 1).replace("events", "%d" % ev_),
                        "%.3f" % r_["grade_entropy"], "%.4f" % cs_, "%.4f" % cg_, "%+.4f" % r_["gap"],
                        "[%+.4f, %+.4f]" % tuple(r_["ci95"])])
    expect("tab:gapext", gx_rows)

    cal = json.load(open(os.path.join(PH, "calibration-and-decision-curve.json")))
    cal_rows = []
    for h, lab in (("12.0", "12 months"), ("24.0", "24 months"), ("36.0", "36 months")):
        o, c = cal["ours"]["by_horizon"][h], cal["clinical"]["by_horizon"][h]
        cal_rows.append([lab, "%.4f / %.4f" % (o["mean_predicted_risk"], o["observed_km_risk"]),
                         "%.4f" % o["brier"], "%.3f" % o["ipa"],
                         "%.4f / %.4f" % (c["mean_predicted_risk"], c["observed_km_risk"]),
                         "%.4f" % c["brier"], "%.3f" % c["ipa"]])
    expect("tab:calibration", cal_rows)
    at = cal["at_risk"]
    for frag in ("At risk: %d, %d and %d patients at 12, 24 and 36 months" % (at["12.0"], at["24.0"], at["36.0"]),
                 "where %d patients remain at risk against %d at one year" % (at["36.0"], at["12.0"])):
        n += 1
        if frag not in " ".join(src.split()):
            fails.append("NC SI: %r does not match the calibration file" % frag)
    # the five per-study re-partition floors, which live in the Supplementary Note rather than the
    # main text, with the sentence that reports each pair rebuilt from the run's own output
    flo = json.load(open(os.path.join(ROOT, "experiments", "20260912-five-cohort-floors",
                                      "results", "five-cohort-floors.json")))["cohorts"]
    NAMES = (("blca", "bladder"), ("brca", "breast"), ("coadread", "colorectal"),
             ("hnsc", "head and neck"), ("stad", "stomach"))
    frag = "are %s and %s in %s, %s and %s in %s, %s and %s in %s, %s and %s in %s, and %s and %s in %s." % tuple(
        x for k, lab in NAMES
        for x in ("%.4f" % flo[k]["floor"]["sd_of_paired_delta"],
                  "%.4f" % flo[k]["floor"]["bar_2sd"], lab))
    n += 1
    if frag not in " ".join(src.split()):
        fails.append("NC SI: the five re-partition floors do not read %r" % frag)
    n += 1
    if not all(v["floor"]["partitions_with_a_positive_delta"] == 24 for v in flo.values()):
        fails.append("NC SI claims 24 of 24 partitions positive in every study; the run disagrees")

    # The second encoder: three studies fall and two rise. Both halves of that sentence are rebuilt
    # from the two runs, because the first draft claimed all five fell and nothing was watching.
    t5 = json.load(open(os.path.join(ROOT, "experiments", "20260912-five-cohort-intervals",
                                     "results", "five-cohort-intervals.json")))["cohorts"]
    g5 = json.load(open(os.path.join(ROOT, "experiments", "20260913-encoder-sensitivity",
                                     "results", "five-cohort-gigassl.json")))["cohorts"]

    def pair(c):
        return ("%.4f" % t5[c]["points"]["OURS_wsi_omics_age_sex_stage"],
                "%.4f" % g5[c]["points"]["OURS_wsi_omics_age_sex_stage"])
    frag = ("bladder from %s to %s, breast from %s to %s and head and neck from %s to %s, and it "
            "rises slightly in the two that do not, stomach from %s to %s and colorectal from "
            "%s to %s." % (pair("blca") + pair("brca") + pair("hnsc") + pair("stad")
                           + pair("coadread")))
    n += 1
    if frag not in " ".join(src.split()):
        fails.append("NC SI: the second-encoder values do not read %r" % frag)
    n += 1
    _fell = [c for c in ("blca", "brca", "hnsc") if g5[c]["points"]["OURS_wsi_omics_age_sex_stage"]
             < t5[c]["points"]["OURS_wsi_omics_age_sex_stage"]]
    _rose = [c for c in ("stad", "coadread") if g5[c]["points"]["OURS_wsi_omics_age_sex_stage"]
             > t5[c]["points"]["OURS_wsi_omics_age_sex_stage"]]
    if len(_fell) != 3 or len(_rose) != 2:
        fails.append("NC SI says three studies fall and two rise under the second encoder; the runs "
                     "give fell=%s rose=%s" % (_fell, _rose))

    ui = json.load(open(os.path.join(PH, "utility-intervals.json")))
    flat = " ".join(src.split())

    def iv(d):
        return "$%+.4f$ ($[%+.4f, %+.4f]$)" % (d["mean"], d["ci95"][0], d["ci95"][1])
    for h, when in (("12.0", "at one year"), ("24.0", "at two years"), ("36.0", "at three years")):
        n += 1
        frag = "%s %s" % (iv(ui["ipa_modrank_minus_clinical"][h]), when)
        if frag not in flat:
            fails.append("NC SI: the IPA difference %s does not read %r" % (when, frag))
    for th in ("0.2", "0.3", "0.4", "0.5", "0.6"):
        n += 1
        frag = iv(ui["net_benefit_24m_modrank_minus_clinical"][th])
        if frag not in flat:
            fails.append("NC SI: the net-benefit difference at %s does not read %r" % (th, frag))
    sl = {k: cal["recalibration_slopes_per_fold"][k] for k in ("ours", "clinical")}
    frag = "ranged from %.2f to %.2f across the five folds for ModRank and from %.2f to %.2f" % (
        min(sl["ours"]), max(sl["ours"]), min(sl["clinical"]), max(sl["clinical"]))
    n += 1
    if frag not in " ".join(src.split()):
        fails.append("NC SI: the recalibration coefficients do not read %r" % frag)

    # amendment A2 (Note 10 and Table 15). The table is bound cell by cell to the before-after
    # file; the note's prose to the result files each number comes from; the tied-time counts are
    # recomputed from the per-case table rather than quoted.
    ba = json.load(open(os.path.join(ROOT, "experiments", "20260928-amendment-a2", "results",
                                     "before-after.json")))["rows"]

    def _ba_fmt(r_, v):
        if v is None:
            return "new"
        return ("%+.4f" % v) if r_["signed"] else "%.4f" % v
    expect("tab:amendA2", [[r_["quantity"], _ba_fmt(r_, r_["before"]), _ba_fmt(r_, r_["after"])]
                           for r_ in ba])
    _pc = json.load(open(os.path.join(PH, "percase-canonical-vectors.json")))["cases"]
    _tc = collections.Counter(round(c_["months"], 6) for c_ in _pc)
    _tied = [k_ for k_, v_ in _tc.items() if v_ > 1]
    _tev = [k_ for k_ in _tied if any(c_["event"] == 1 and round(c_["months"], 6) == k_ for c_ in _pc)]
    fuj = json.load(open(os.path.join(PH, "unified-fusion-and-added-value.json")))["fusion_on_identical_inputs"]
    _t = fuj["stacked_tuned_ridge"]
    _tp = fuj["paired"]["ours_vs_stacked_tuned_ridge"]
    _dmf = fa["deepmisl_wsi_pathways"]
    _baq = {r_["quantity"]: r_ for r_ in ba}
    _gg = json.load(open(os.path.join(ROOT, "experiments", "20260911-geo-external", "results",
                                      "geo-gated.json")))["cohorts"]
    _g48 = _gg["GSE48075"]
    _geo48 = json.load(open(os.path.join(ROOT, "experiments", "20260911-geo-external", "results",
                                         "geo-external.json")))["cohorts"]["GSE48075"]["cv"]
    _sx = json.load(open(os.path.join(ROOT, "experiments", "20260912-subgroup-fairness", "results",
                                      "subgroup-fairness.json")))["splits"]["sex"]["clinical_gap_minus_ours_gap"]
    _sn2 = json.load(open(os.path.join(ROOT, "experiments", "20260818-selection-null", "results",
                                       "selection-null.json")))["A_selection_null"]["bars"]["empirical_max_q95"]["bar"]
    _bar_before = _baq["selection-null bar"]["before"]
    for _ph in ("The development cohort has %d tied times, %d of them with an event, covering %d patients"
                % (len(_tied), len(_tev), sum(_tc[k_] for k_ in _tied)),
                "it reaches $%.4f \\pm %.4f$, and ModRank exceeds it by $%+.4f$ ($[%+.4f, %+.4f]$)"
                % (_t["mean"], _t["sd"], _tp["mean"], _tp["ci95"][0], _tp["ci95"][1]),
                "in 13 of its 25 retraining--fold pairs",
                "its $D$ fell from $%+.3f$ to $%+.3f$" % (_baq["D, DeepMISL with pathways"]["before"],
                                                         _baq["D, DeepMISL with pathways"]["after"]),
                "their gap rose from $%+.3f$ to $%+.3f$" % (_baq["GSE19915, stage minus grade reference"]["before"],
                                                            _baq["GSE19915, stage minus grade reference"]["after"]),
                "The selection bar moved from %.4f to %.4f" % (_bar_before, _sn2),
                "dropped the transcriptome arm in %d of 25 folds and the stage block in %d"
                % (_g48["gate_pass_counts"]["transcriptome_fail"], _g48["gate_pass_counts"]["clinical_stage_fail"]),
                "(%.4f against %.4f for the block and %.4f for ModRank)"
                % (_g48["cv"]["gated"]["mean"], _geo48["clinical_stage"]["mean"], _geo48["modrank"]["mean"]),
                "equalled ModRank (%.4f)" % _gg["GSE32894"]["cv"]["gated"]["mean"],
                "and ModRank's is $%+.4f$ $[%+.4f, %+.4f]$" % (_sx["mean"], _sx["ci95"][0], _sx["ci95"][1])):
        n += 1
        if _ph not in flat:
            fails.append("NC SI should carry %r" % _ph)
    if _bar_before >= 0.7716 + 5e-5 or _bar_before <= 0.7716 - 5e-5:
        fails.append("the before-after table's pre-A2 selection bar is not 0.7716")
    return n


if __name__ == "__main__":
    sys.exit(main())
