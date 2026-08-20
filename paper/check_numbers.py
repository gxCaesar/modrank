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

import json
import math
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "experiments", "20260817-blca-confirm", "results")
TEX = os.path.join(ROOT, "paper", "main.tex")


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
    tex = open(TEX).read()
    fails = []

    for label, fname, path, expected, tol in CHECKS:
        try:
            got = dig(load(fname), path)
        except (KeyError, IndexError, FileNotFoundError) as exc:
            fails.append("%-38s UNRESOLVABLE  %s" % (label, exc))
            continue
        if abs(float(got) - float(expected)) > tol:
            fails.append("%-38s file=%s  paper=%s  (%s -> %s)"
                         % (label, got, expected, fname, path))
        else:
            print("  ok  %-38s %s" % (label, got))

    # --- the contradiction check: a superseded value must not survive anywhere in the manuscript
    SUPERSEDED = {
        "0.0318": "margin vs SurvPath under the OLD mixed construction",
        "0.0334": "margin vs SurvPath under the OLD seed-0 numerator",
        "0.116":  "p vs SurvPath, superseded",
        "0.263":  "p vs PIBD, superseded",
        "0.6982": "PIBD+clinical under the OLD seed-0 clinical block",
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
    # methods.tex is \input by the internal manuscript, so a claim can live there and be invisible
    # to a scan of main.tex alone. It was: the first run of this check missed it entirely.
    SCAN = (TEX,
            os.path.join(ROOT, "paper", "methods.tex"),
            os.path.join(ROOT, "paper", "bib-submission", "main.tex"),
            os.path.join(ROOT, "paper", "bib-submission", "supplementary.tex"))
    # --- A DERIVED value the manuscript states and no file stores: the margin over the incumbent
    # minus the selection-inflation term. It was written as +0.0183 when the arithmetic gives
    # +0.0178, and nothing caught it because neither number is a stored field.
    a1 = load("amendment-A1-clinical-provenance.json")
    gap = a1["comparisons"]["vs_DIMAF_published"]["gap"]
    clears = round(gap - 0.0244, 4)
    for doc in (os.path.join(ROOT, "paper", "bib-submission", "main.tex"),
                os.path.join(ROOT, "paper", "main.tex")):
        body = open(doc).read()
        if "clears" in body and "0.0244" in body:
            want = "$+%.4f$" % clears
            if want not in body:
                fails.append("%s states the margin clears selection inflation by something other "
                             "than %s (%.4f minus 0.0244)"
                             % (os.path.relpath(doc, ROOT), want, gap))
    print("  ok  %-38s +%.4f = %.4f - 0.0244" % ("margin over selection inflation", clears, gap))

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
    for _doc in (os.path.join(ROOT, "paper", "bib-submission", "main.tex"),
                 os.path.join(ROOT, "paper", "bib-submission", "supplementary.tex"),
                 os.path.join(ROOT, "paper", "main.tex"),
                 os.path.join(ROOT, "paper", "methods.tex")):
        _b = " ".join(open(_doc).read().split())
        _rel = os.path.relpath(_doc, ROOT)
        if re.search(r"0\.679\s*\+\s*0\.0244", _b):
            fails.append("%s states the frozen decision rule with the ledger-final term "
                         "(0.679 + 0.0244); the protocol froze 0.679 + %.4f" % (_rel, _infl))
        # a document that states the term at all must state the frozen N beside it
        if re.search(r"sqrt\{2\\ln N\}|sqrt\{2 \\ln N\}|selection.inflation", _b, re.I):
            if str(_n_frozen) not in _b:
                fails.append("%s states the selection-inflation term without the frozen N=%d"
                             % (_rel, _n_frozen))
    print("  ok  %-38s frozen N=%d, inflation %.4f, target %.4f"
          % ("pre-registered selection term", _n_frozen, _infl, _target))

    # --- A CROSS-DOCUMENT DISPLAY-ITEM POINTER. The body and the supplement are separate compiled
    # documents, so a \ref between them renders as ?? and the only workable pointer is a typed
    # number. A typed number is checked by nobody: the supplement sent readers to "Figure 5e" for
    # the twelve-case series, which is Figure 7e, and Figure 5 has a real panel e about something
    # else, so the reader landed somewhere plausible and wrong. Found 2026-08-20. The number is
    # therefore bound to the main document's own .aux, and the figure is identified by what its
    # caption says rather than by a label spelled twice.
    _mainaux = os.path.join(ROOT, "paper", "bib-submission", "main.aux")
    _mainsrc = os.path.join(ROOT, "paper", "bib-submission", "main.tex")
    if os.path.exists(_mainaux):
        _num = {m.group(1): m.group(2) for m in
                re.finditer(r"\\newlabel\{([^}]*)\}\{\{([^}]*)\}\{\d+\}", open(_mainaux).read())}
        _mb = open(_mainsrc).read()
        _owner = None
        for _m in re.finditer(r"\\begin\{figure\*?\}(.*?)\\end\{figure\*?\}", _mb, re.S):
            if "twelve cases" in _m.group(1):
                _lab = re.search(r"\\label\{([^}]*)\}", _m.group(1))
                _owner = _lab.group(1) if _lab else None
        if _owner is None or _owner not in _num:
            fails.append("no main-text figure caption names the twelve-case series, so the "
                         "supplement's pointer cannot be checked")
        else:
            _want = "Figure~%se of the main paper" % _num[_owner]
            _sb = " ".join(open(os.path.join(ROOT, "paper", "bib-submission",
                                             "supplementary.tex")).read().split())
            if _want not in _sb:
                fails.append("the supplement does not point the twelve-case series at %r "
                             "(%s is Figure %s in the main document)"
                             % (_want, _owner, _num[_owner]))
            else:
                print("  ok  %-38s twelve-case series -> Figure %se"
                      % ("cross-document figure pointer", _num[_owner]))

    # --- THE COVER LETTER gets the two prose rules the loop above applies to the manuscript and
    # skipped for it. Its register is different (it is addressed to an editor, so the second person
    # is correct there and is not checked), but a path or a clause semicolon is wrong in a letter
    # for the same reason it is wrong in the paper.
    _cl = os.path.join(ROOT, "paper", "cover_letter.tex")
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
        print("  ok  %-38s no path, no clause semicolon" % "cover letter prose")

    # --- THE MODALITY FIGURE'S NUMBERS. Its panels come from the S5 evidence directory rather
    # than the confirmatory results, so the table above never reaches them. Every quantity the
    # new Results subsection states in prose is recomputed here from the frozen source.
    PROSE_MODALITY_DOCS = [os.path.join(ROOT, "paper", "bib-submission", "main.tex"),
                           os.path.join(ROOT, "paper", "main.tex")]
    _atlas, _c6 = load5("atlas.json"), load5("c6.json")
    _tie = _atlas["conditional_probe"]["clinically_tied_q05"]
    _img = [_c6["single_arms"][k] for k in ("wsi_titan", "wsi_chief_mean", "wsi_chief_dispersion")]
    _om = [_c6["single_arms"][k] for k in ("omics_xena", "omics_combine", "omics_hallmarks")]
    MOD_CHECKS = [
        ("clinical alone", _c6["single_arms"]["clinical"], "0.6856"),
        ("slide alone", _c6["single_arms"]["wsi_titan"], "0.6596"),
        ("pathway means alone", _c6["single_arms"]["omics_combine"], "0.6510"),
        ("spread within the image family", max(_img) - min(_img), "0.104"),
        ("spread within the omics family", max(_om) - min(_om), "0.069"),
        ("rho slide vs clinical", _atlas["redundancy"]["spearman_titan_vs_clinical"], "0.213"),
        ("rho incumbent vs clinical", _atlas["redundancy"]["spearman_survpath_vs_clinical"],
         "0.166"),
        ("clinical on the tightest ties", _tie["arms"]["clinical"], "0.5012"),
        ("slide on the tightest ties", _tie["arms"]["titan"], "0.6120"),
    ]
    for label, value, written in MOD_CHECKS:
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
    print("  ok  %-38s %d quantities, both documents" % ("modality figure", len(MOD_CHECKS)))

    # --- THE REPORTING DUMP. Two supplementary figures and one Results paragraph are built on
    # it, and it is a SECOND execution of the frozen arm construction, so the first thing checked
    # is that it reproduced the frozen primary. Everything else here is recomputed from its own
    # per-case table rather than trusted from a summary field.
    dump = json.load(open(os.path.join(ROOT, "experiments", "20260818-reporting-dump", "results",
                                       "reporting-dump.json")))
    if not dump.get("frozen_primary_matched") or abs(dump["primary"] - 0.7212) > 5e-5:
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
    for label, val, written in (("per-fold minimum", min(_pf), "0.678"),
                                ("per-fold maximum", max(_pf), "0.782"),
                                ("per-fold range", max(_pf) - min(_pf), "0.103")):
        if "%.3f" % val != written:
            fails.append("supplement: %s is %.4f recomputed, written as %s"
                         % (label, val, written))
        for doc in PROSE_MODALITY_DOCS + [os.path.join(ROOT, "paper", "bib-submission",
                                                       "supplementary.tex")]:
            if written not in open(doc).read():
                fails.append("%s never states %s (%s)" % (os.path.relpath(doc, ROOT), written,
                                                          label))
    print("  ok  %-38s primary reproduced, %d folds recomputed"
          % ("reporting dump", len(_folds)))

    # --- METRIC PARITY, answering round 2's G10. The paper now states both conventions on the
    # same predictions, so both are recomputed from the dump's own per-case table rather than
    # read back from the summary that produced them.
    _mp = os.path.join(ROOT, "experiments", "20260818-reporting-dump", "results",
                       "metric-parity.json")
    if os.path.exists(_mp):
        mp = json.load(open(_mp))
        for label, key, written in (("pooled out-of-fold", "pooled_out_of_fold", "0.7239"),
                                    ("mean of per-fold", "mean_of_per_fold", "0.7272"),
                                    ("their difference", "difference_mean_minus_pooled",
                                     "+0.0034")):
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
        SEL = [("null single-candidate sd", A["null_single_candidate"]["sd"], "0.0397"),
               ("null single-candidate mean", A["null_single_candidate"]["mean"], "0.5011"),
               ("null max q95", A["null_max_over_the_family"]["q95"], "0.5926"),
               ("gaussian bar with the null sigma", A["bars"]["gaussian_with_the_null_sigma"]["bar"],
                "0.8087"),
               ("empirical max bar", A["bars"]["empirical_max_q95"]["bar"], "0.7716"),
               ("slide-arm span across encoders", B["slide_arm_span"], "0.1372"),
               ("full-method span across encoders", B["full_method_span"], "0.0521"),
               ("full method, weakest encoder", B["full_method_min"], "0.6704"),
               ("full method, strongest encoder", B["full_method_max"], "0.7225")]
        for label, value, written in SEL:
            if "%.4f" % value != written:
                fails.append("selection null: %s is %.4f at source, written as %s"
                             % (label, value, written))
            for doc in PROSE_MODALITY_DOCS:
                if written not in open(doc).read():
                    fails.append("%s never states %s (%s)"
                                 % (os.path.relpath(doc, ROOT), written, label))
        # the two claims the paragraphs rest on, asserted rather than assumed
        if 0.7212 >= A["bars"]["empirical_max_q95"]["bar"]:
            fails.append("the paper says 0.7212 does NOT clear the corrected bar; it now does")
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
            # The supplement enumerates the strata in prose. The gate binds each count to its
            # stratum by proximity rather than by an exact phrase, so the sentence may be
            # rewritten freely and a count may not drift away from the thing it counts.
            _sup = " ".join(open(os.path.join(ROOT, "paper", "bib-submission",
                                              "supplementary.tex")).read().split())
            for key, word in (("full_cohort_out_of_fold", "full-cohort out-of-fold"),
                              ("subgroup_slice", "subgroup slice"),
                              ("in_sample_oracle", "in-sample oracle"),
                              ("permuted_or_null_control", "permuted or null control"),
                              ("cross_cohort_transfer_probe", "cross-cohort transfer probe")):
                count = sr["strata"][key]
                if not re.search(r"\b%d\b[^.]{0,60}?%s" % (count, re.escape(word)), _sup):
                    fails.append("the supplement does not bind %d to \"%s\"; a stratification "
                                 "whose parts are not all written down does not close for a "
                                 "reader" % (count, word))
        print("  ok  %-38s corrected bar %.4f, %d of 7 encoders clear 0.679"
              % ("selection null and encoder swap",
                 A["bars"]["empirical_max_q95"]["bar"],
                 B["encoders_whose_full_method_still_beats_the_incumbent"]))

    # --- DISPLAY-ITEM ORDER in the submission document. A copy-editor checks this and it broke
    # twice today, once when a figure was inserted and once when a table's float sat 30,000
    # characters after its first citation. Gated on the BiB main only: the internal draft is
    # Nature-shaped, its abstract cites Figure 2 first, and that ordering is deliberate there.
    _sub = open(os.path.join(ROOT, "paper", "bib-submission", "main.tex")).read()
    for pre, what in (("fig:", "figures"), ("tab:", "tables")):
        declared = re.findall(r"\\label\{(%s[^}]+)\}" % pre, _sub)
        firstcite = {}
        for m in re.finditer(r"\\ref\{(%s[^}]+)\}" % pre, _sub):
            firstcite.setdefault(m.group(1), m.start())
        never = [l for l in declared if l not in firstcite]
        if never:
            fails.append("submission document never cites %s" % never)
        byc = sorted(declared, key=lambda l: firstcite.get(l, 10 ** 9))
        if byc != declared:
            fails.append("submission %s are not cited in numbering order: declared %s, cited %s"
                         % (what, declared, byc))
    print("  ok  %-38s all cited, in numbering order" % "display-item order")

    # --- THE COVER LETTER. It is the one document an editor reads beside the paper, and a letter
    # that disagrees with its manuscript is the inconsistency they are guaranteed to see. Every
    # number in it must appear in the submission document, and the letter carries the same prose
    # rules as the paper.
    _cl = os.path.join(ROOT, "paper", "cover_letter.tex")
    if os.path.exists(_cl):
        letter = open(_cl).read()
        sub = open(os.path.join(ROOT, "paper", "bib-submission", "main.tex")).read()
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
    PROSE_DOCS = [os.path.join(ROOT, "paper", "bib-submission", "main.tex"),
                  os.path.join(ROOT, "paper", "bib-submission", "supplementary.tex"),
                  os.path.join(ROOT, "paper", "main.tex"),
                  os.path.join(ROOT, "paper", "methods.tex")]
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
    PROSE_DOCS = [os.path.join(ROOT, "paper", "bib-submission", "main.tex"),
                  os.path.join(ROOT, "paper", "bib-submission", "supplementary.tex"),
                  os.path.join(ROOT, "paper", "main.tex"),
                  os.path.join(ROOT, "paper", "methods.tex")]
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

    # --- FIGURE GEOMETRY. Authoring at final size only means something if the insertion scale is
    # exactly 1. The main manuscript is TWO-COLUMN under the OUP class -- textwidth 488.5pt, not
    # the 453.6pt of the internal article build -- so figures drawn at 6.30in were being enlarged
    # 1.077x and a 9pt label was landing at 9.69pt. The supplementary is a plain article whose
    # textwidth really is 453.6pt, so its figures are a different target.
    import subprocess as _sp
    import glob as _g
    FIGDIR = os.path.join(ROOT, "paper", "figures")
    TARGET = {"main": 488.5, "supp": 453.6}
    # Only figures a manuscript actually INCLUDES. Globbing the directory also measured a
    # superseded draft that nothing references, which is noise; this way an included-but-missing
    # figure fails too, which the glob could never catch.
    included = []
    for doc, key in ((os.path.join(ROOT, "paper", "bib-submission", "main.tex"), "main"),
                     (os.path.join(ROOT, "paper", "bib-submission", "supplementary.tex"), "supp")):
        for m in re.finditer(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", open(doc).read()):
            included.append((m.group(1), key))
    for base, key in sorted(set(included)):
        fn = os.path.join(FIGDIR, base)
        want = TARGET[key]
        if not os.path.exists(fn):
            fails.append("%s is included by the manuscript but does not exist" % base); continue
        try:
            out = _sp.run(["pdfinfo", fn], capture_output=True, text=True).stdout
            w = float([l for l in out.split("\n") if "Page size" in l][0].split()[2])
        except Exception as exc:
            fails.append("cannot read %s to check its insertion scale: %s" % (base, exc)); continue
        if abs(w - want) > 0.6:
            fails.append("%s is %.1fpt wide but is inserted into a %.1fpt text block; the scale is "
                         "%.4f, so source pt is not rendered pt" % (base, w, want, want / w))
    print("  ok  %-38s %d included figures, insertion scale 1.0000"
          % ("figure geometry", len(set(included))))

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

    # --- CROSS-DOCUMENT references. The supplementary is a separate compile, so its numbers cannot
    # be \ref'd from the main text and are written by hand. Three of six were wrong the first time,
    # because the labels were positional -- `tab:s4` rendered as Table S2. This reads the
    # supplementary's own .aux and checks every hand-written number against what it actually is.
    aux = os.path.join(ROOT, "paper", "bib-submission", "supplementary.aux")
    bibmain = os.path.join(ROOT, "paper", "bib-submission", "main.tex")
    if not os.path.exists(aux):
        fails.append("supplementary.aux is missing; cross-document numbers cannot be checked "
                     "(compile paper/bib-submission/supplementary.tex first)")
    else:
        rendered = dict(re.findall(r"\\newlabel\{([^}]*)\}\{\{(S\d+)\}", open(aux).read()))
        EXPECT = {"Supplementary Table~S1": "tab:sup-census",
                  "Supplementary Table~S3": "tab:sup-components",
                  "Supplementary Table~S4": "tab:sup-perseed",
                  # the three thin campaign figures merged into one on 2026-08-18, so the
                  # main text now points at panels of S1 rather than at S1, S2 and S3
                  "Supplementary Figure~S1a,~b": "fig:sup-campaign",
                  "Supplementary Figure~S1c,~d": "fig:sup-campaign",
                  "Supplementary Figure~S1e,~f": "fig:sup-campaign",
                  "Supplementary Figure~S3": "fig:sup-perfold"}
        body = open(bibmain).read()
        # A MISSING phrase is a failure, not a skip. The first version of this check skipped what
        # it could not find, so mutating "Table~S3" to "Table~S2" made the check pass by making the
        # thing it was checking disappear -- which is the exact defect it exists to catch.
        for phrase, label in EXPECT.items():
            want = re.search(r"S\d+", phrase).group()   # the phrase may carry panel letters
            got = rendered.get(label)
            if got != want:
                fails.append("%s renders as %s in the supplementary, but the main text is written "
                             "against %s" % (label, got, want))
            elif phrase not in body:
                fails.append("main.tex no longer cites %r (%s); either it was renumbered by hand "
                             "or the pointer was dropped" % (phrase, label))
        # and nothing may point at a supplementary item that does not exist
        valid = set(rendered.values())
        for m in re.finditer(r"Supplementary (?:Table|Figure)~(S\d+)", body):
            if m.group(1) not in valid:
                fails.append("main.tex cites Supplementary %s, which the supplementary does not "
                             "contain" % m.group(1))
        print("  ok  %-38s %d mappings, %d valid targets"
              % ("supplementary cross-references", len(EXPECT), len(valid)))
    for path in SCAN:
        body = open(path).read()
        for m in re.finditer(r"(\w+) of the (\w+) published entries", body):
            got_below, got_pub = m.group(1), m.group(2)
            if got_below != WORD.get(n_below) or got_pub != WORD.get(n_pub):
                line = body[:m.start()].count("\n") + 1
                fails.append("%s:%d says '%s of the %s published entries'; recomputing from the "
                             "results files gives '%s of the %s'"
                             % (os.path.relpath(path, ROOT), line, got_below, got_pub,
                                WORD.get(n_below), WORD.get(n_pub)))
        print("  ok  %-38s %s of the %s (%s)"
              % ("published-entry count", WORD.get(n_below), WORD.get(n_pub),
                 os.path.relpath(path, ROOT)))

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

    for m in re.finditer(r"0\.002\s*(?:to|--|-)\s*0\.012", tex):
        window = tex[max(0, m.start() - 400):m.end() + 400]
        if "0.029" not in window:
            line = tex[:m.start()].count("\n") + 1
            fails.append("main.tex:%d states the published gaps as 0.002-0.012 without naming the "
                         "0.029 one nearby -- that phrasing is what made the claim wrong" % line)
    for token, why in SUPERSEDED.items():
        for path in SCAN:
            body = open(path).read()
            for i, line in enumerate(body.split("\n"), 1):
                if token in line:
                    fails.append("SUPERSEDED VALUE %-9s at %s:%d -- %s"
                                 % (token, os.path.relpath(path, ROOT), i, why))

    print()
    if fails:
        print("FAIL (%d)" % len(fails))
        for f in fails:
            print("   " + f)
        return 1
    print("PASS -- %d numbers verified, no superseded value present" % len(CHECKS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
