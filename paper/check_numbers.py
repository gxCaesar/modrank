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
import subprocess
import sys

import numpy as np

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
    # methods.tex is \input by the internal manuscript, so a claim can live there and be invisible
    # to a scan of main.tex alone. It was: the first run of this check missed it entirely.
    SCAN = (TEX,
            os.path.join(ROOT, "paper", "methods.tex"),
            os.path.join(ROOT, "paper", "bib-submission", "main.tex"),
            os.path.join(ROOT, "paper", "bib-submission", "supplementary.tex"),
            os.path.join(ROOT, "paper", "nc-submission", "main.tex"),
            os.path.join(ROOT, "paper", "nc-submission", "supplementary.tex"))
    # --- STALE CLAIMS, banned as phrases. Found by an external review on 2026-09-11: after the
    # selection correction was shown to be mis-specified (0.0397, not 0.0073), four sentences still
    # said the result "clears the selection-inflation term", and the Limitations quoted a p that
    # matched no source. The derived +0.0178 margin that used to be checked here was one of those
    # sentences, so the check that guarded its arithmetic now guards its absence. Each phrase below
    # is a claim the evidence no longer supports; none of them may come back by a later edit.
    STALE = {
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
    _gr = a1["arms_seed0"]["OURS_with_grade_instead_of_stage"]
    if "%.4f" % _gr != "0.6856":
        fails.append("ModRank with the grade arm is %.4f at source, the prose says 0.6856" % _gr)
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
    _bio = load("biology.json")["B5_within_molecular_subtype"]
    _p_axis = "$p=%.4f$" % float("%.2g" % _bio["subtype_axis_itself"]["logrank_p"])
    for _doc in (os.path.join(ROOT, "paper", "bib-submission", "main.tex"), TEX):
        _src = open(_doc).read()
        _tb = re.search(r"\\label\{tab:subtype\}.*?\\end\{tabular\}", _src, re.S)
        if not _tb:
            fails.append("%s has no subtype table" % os.path.relpath(_doc, ROOT))
            continue
        _rows = {}
        for _row in _tb.group(0).split("\\\\"):
            _cells = [re.sub(r"\\textbf\{([^}]*)\}", r"\1", c).strip() for c in _row.split("&")]
            if len(_cells) == 6 and _cells[0].endswith("-leaning"):
                _rows[_cells[0].split()[-1]] = _cells[1:]
        for _k, _lab in (("basal_leaning", "basal-leaning"), ("luminal_leaning", "luminal-leaning")):
            _g = _bio[_k]
            _exp = ["%d / %d" % (_g["n"], _g["events"])] + ["%.4f" % _g[a] for a in
                                                             ("clinical", "image", "omics", "OURS")]
            if _rows.get(_lab) != _exp:
                fails.append("%s: the %s row reads %s; biology.json gives %s"
                             % (os.path.relpath(_doc, ROOT), _lab, _rows.get(_lab), _exp))
        _para = _src[_src.find("And the two modalities are informative"):_tb.start()]
        for _v in (["%.4f" % _bio[k][a] for k in ("basal_leaning", "luminal_leaning")
                    for a in ("image", "clinical", "OURS")] + [_p_axis]):
            if _v not in _para:
                fails.append("%s: the subtype paragraph does not state %s from biology.json"
                             % (os.path.relpath(_doc, ROOT), _v))
    print("  ok  %-38s rows and prose bound to biology.json B5" % "within-subtype table")

    # --- THE LIMITATIONS' p FOR THE CLINICAL-ALONE COMPARISON, bound to the Holm family.
    _hc = load("pibd-parity.json")["holm_family_of_six"]["comparisons"]["ours vs clinical alone"]
    for _doc in (os.path.join(ROOT, "paper", "bib-submission", "main.tex"), TEX):
        _b = " ".join(open(_doc).read().split())
        _need = "clinical arm alone at $p=%.4f$, which is $%.4f$ after Holm" % (_hc["p_raw"],
                                                                                _hc["p_holm"])
        if _need not in _b:
            fails.append("%s: the Limitations do not state the clinical-alone comparison as %r"
                         % (os.path.relpath(_doc, ROOT), _need))
    print("  ok  %-38s p %.4f, Holm %.4f" % ("clinical-alone p in Limitations", _hc["p_raw"],
                                             _hc["p_holm"]))

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
    for _cl in (os.path.join(ROOT, "paper", "cover_letter.tex"),
                os.path.join(ROOT, "paper", "nc-submission", "cover_letter.tex")):
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
    PROSE_MODALITY_DOCS = [os.path.join(ROOT, "paper", "bib-submission", "main.tex"),
                           os.path.join(ROOT, "paper", "main.tex")]
    _atlas, _c6 = load5("atlas.json"), load5("c6.json")
    _am = json.load(open(os.path.join(ROOT, "experiments", "20260911-blca-posthoc", "results",
                                      "modality-atlas-amended.json")))
    if not _am["known_answer"]["reproduced"]:
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
         ["clinical_age_sex_stage_from_DIMAF_file"], "0.6638"),
        ("slide alone", _c6["single_arms"]["wsi_titan"], "0.6596"),
        ("pathway means alone", _c6["single_arms"]["omics_combine"], "0.6510"),
        ("spread within the image family", max(_img) - min(_img), "0.104"),
        ("spread within the omics family", max(_om) - min(_om), "0.069"),
        # redundancy and tied pairs from the AMENDED atlas (s19b), which first reproduces the
        # development atlas on its own pre-amendment block and then swaps only the clinical block
        ("rho slide vs clinical, amended", _am["redundancy"]["spearman_titan_vs_clinical"], "0.2445"),
        ("rho incumbent vs clinical, amended", _am["redundancy"]["spearman_survpath_vs_clinical"],
         "0.1924"),
        ("clinical on the tightest ties, amended", _tie["arms"]["clinical"], "0.4959"),
        ("slide on the tightest ties, amended", _tie["arms"]["titan"], "0.6220"),
        ("tightest-tie pair count, amended", _tie["pairs"], "1,213"),
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
    # the two correlations the supplement's pathway paragraph states, recomputed from the dump.
    # Until 2026-09-11 it said the slide and transcriptome arms share Spearman 0.48; that is the
    # slide arm against the INCUMBENT's joint score. The two arms share 0.28.
    _sl = np.array([c_["seed_mean"]["wsi_titan"] for c_ in dump["cases"]])
    _om = np.array([c_["seed_mean"]["omics_combine"] for c_ in dump["cases"]])
    _rk = lambda v: np.argsort(np.argsort(v)).astype(float)
    _rho = float(np.corrcoef(_rk(_sl), _rk(_om))[0, 1])
    _rbt = float(np.corrcoef([p_["rho_slide"] for p_ in dump["pathways"]],
                             [p_["rho_omics"] for p_ in dump["pathways"]])[0, 1])
    _supb = " ".join(open(os.path.join(ROOT, "paper", "bib-submission", "supplementary.tex")).read().split())
    for _frag in ("sharing Spearman %.2f between their out-of-fold scores" % _rho,
                  "profiles correlate at $r=%.2f$" % _rbt):
        if _frag not in _supb:
            fails.append("BiB supplement: %r does not match the reporting dump" % _frag)


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
    for _cl, _subdoc in ((os.path.join(ROOT, "paper", "cover_letter.tex"),
                          os.path.join(ROOT, "paper", "bib-submission", "main.tex")),
                         (os.path.join(ROOT, "paper", "nc-submission", "cover_letter.tex"),
                          os.path.join(ROOT, "paper", "nc-submission", "main.tex"))):
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
    PROSE_DOCS = [os.path.join(ROOT, "paper", "bib-submission", "main.tex"),
                  os.path.join(ROOT, "paper", "bib-submission", "supplementary.tex"),
                  os.path.join(ROOT, "paper", "main.tex"),
                  os.path.join(ROOT, "paper", "methods.tex"),
                  os.path.join(ROOT, "paper", "nc-submission", "main.tex"),
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
    PROSE_DOCS = [os.path.join(ROOT, "paper", "bib-submission", "main.tex"),
                  os.path.join(ROOT, "paper", "bib-submission", "supplementary.tex"),
                  os.path.join(ROOT, "paper", "main.tex"),
                  os.path.join(ROOT, "paper", "methods.tex"),
                  os.path.join(ROOT, "paper", "nc-submission", "main.tex"),
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
        for _doc in SCAN + (os.path.join(ROOT, "paper", "cover_letter.tex"),):
            if not os.path.isfile(_doc):
                continue
            _b = open(_doc).read()
            for _tok in _EMB:
                if _tok in _b:
                    fails.append("EMBARGO: %r from the CHIMERA cohort appears in %s"
                                 % (_tok, os.path.relpath(_doc, ROOT)))
        print("  ok  %-38s no CHIMERA name or value in %d documents"
              % ("CHIMERA embargo", len(SCAN) + 1))

    # ---- the encoder-parity arm, read from its own result file. Its known-answer check is what
    #      licenses the number at all, so the gate refuses the prose if that check did not pass.
    _ep = os.path.join(ROOT, "experiments", "20260820-encoder-parity", "results",
                       "encoder-parity.json")
    if os.path.isfile(_ep):
        _e = json.load(open(_ep))
        _sub = open(os.path.join(ROOT, "paper", "bib-submission", "main.tex")).read()
        if not _e["known_answer_check"]["reproduced"]:
            fails.append("the encoder-parity run's known-answer check did not reproduce, so its "
                         "number must not appear in the manuscript")
        for _k, _v in (("chief full method", _e["chief"]["full_method"]),
                       ("chief slide arm", _e["chief"]["slide_alone"])):
            if ("%.4f" % _v) not in _sub:
                fails.append("%s is %.4f in the result file and does not appear in the manuscript"
                             % (_k, _v))
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
            for _p in [os.path.join(ROOT, "paper", "cover_letter.tex")] + list(SCAN):
                if not os.path.isfile(_p):
                    continue
                for _q in re.finditer(r"recomputes\s+(\d+)\s+(?:reported\s+)?quantit", open(_p).read()):
                    if _q.group(1) != _n:
                        fails.append("%s claims the check script recomputes %s quantities; it "
                                     "reports %s" % (os.path.relpath(_p, ROOT), _q.group(1), _n))
            print("  ok  recomputed-quantity count           prose agrees with the script at %s" % _n)

    # ---- the published identifiers, bound to the one file the builder reads them from.
    #      CITATION.cff, .zenodo.json and this manuscript all take these strings from
    #      release-identifiers.json, and a manuscript quoting a different DOI than the archive
    #      actually carries is the kind of error nobody re-reads for.
    _idp = os.path.join(ROOT, "tools_release", "release-identifiers.json")
    if os.path.isfile(_idp):
        _ident = json.load(open(_idp))
        _sub = open(os.path.join(ROOT, "paper", "bib-submission", "main.tex")).read()
        _named = 0
        for _key, _label in (("concept_doi", "concept DOI"), ("version_doi", "version DOI"),
                             ("repository_url", "repository URL")):
            _val = _ident.get(_key)
            if not _val:
                continue
            _named += 1
            if _val not in _sub:
                fails.append("the %s in the identifiers file (%s) does not appear in the "
                             "manuscript" % (_label, _val))
        # A placeholder that survived is the failure this gate was added to close.
        _open = [i + 1 for i, ln in enumerate(_sub.split("\n"), 0)
                 if re.search(r"\\TODO(?:port|author)?\{", ln) and "newcommand" not in ln]
        if _ident.get("concept_doi") and _open:
            fails.append("the identifiers exist and the manuscript still carries an unresolved "
                         "placeholder at line %s" % ", ".join(str(n) for n in _open))
        if _named:
            print("  ok  published identifiers               %d bound to the identifiers file, "
                  "%d placeholder%s left"
                  % (_named, len(_open), "" if len(_open) == 1 else "s"))

    n_nc = check_nc(fails)
    print("  ok  %-38s %d quantities bound to their result files" % ("Nature Communications build", n_nc))

    print()
    if fails:
        print("FAIL (%d)" % len(fails))
        for f in fails:
            print("   " + f)
        return 1
    print("PASS -- %d numbers verified, no superseded value present" % len(CHECKS))
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
        ("abstract: grade reference", a1["arms_seed0"]["clinical_age_sex_GRADE_from_DIMAF_file"], "0.567",
         "raised its concordance from 0.567 to 0.664"),
        ("abstract: stage reference", a1["arms_seed0"]["clinical_age_sex_stage_from_DIMAF_file"], "0.664",
         "raised its concordance from 0.567 to 0.664"),
        ("abstract: primary", a1["primary"]["value"], "0.721", "It reached 0.721, the highest value"),
        ("abstract: re-partitions", rs["resplits"]["partitions"], "24", "in all 24 re-partitions"),
        ("clinical with grade", a1["arms_seed0"]["clinical_age_sex_GRADE_from_DIMAF_file"], "0.5666",
         "The grade construction reaches 0.5666"),
        ("clinical with stage", a1["arms_seed0"]["clinical_age_sex_stage_from_DIMAF_file"], "0.6638",
         "the stage construction 0.6638"),
        ("stage minus grade", a1["arms_seed0"]["stage_minus_grade"], "+0.0972", "a gap of +0.0972"),
        ("stage minus grade, CI low", sgv["ci95"][0], "+0.047", "$[+0.047, +0.149]$"),
        ("stage minus grade, CI high", sgv["ci95"][1], "+0.149", "$[+0.047, +0.149]$"),
        # the regime
        ("noise fit in-sample", rg["capacity_sweep"]["arms"][0]["corrected_null"], "0.7474",
         "reaches 0.7474 concordance at sixteen parameters"),
        ("fusion oracle", rg["fusion"]["linear_fusion_oracle"]["value"], "0.7268", "(0.7268 against 0.7191)"),
        ("equal weight, deployed", rg["fusion"]["deployed_equal_weight"], "0.7191", "(0.7268 against 0.7191)"),
        ("what the oracle buys", rg["fusion"]["linear_fusion_oracle"]["value"]
         - rg["fusion"]["deployed_equal_weight"], "+0.008", "is worth only $+0.008$ over equal weighting"),
        ("gated fusion", c2["gated_on_clinical_tertile"], "0.7493", "(0.7493 against 0.7502)"),
        ("its permuted control", c2["gated_on_PERMUTED_tertile_CONTROL"], "0.7502", "(0.7493 against 0.7502)"),
        ("split-reseed SD", unc["paired_split_reseed_sd"], "0.0073", "a standard deviation of 0.0073"),
        ("split-reseed bar", unc["reseed_bar_2sd"], "0.0145", "a margin must exceed 0.0145"),
        # the confirmatory result
        ("primary", a1["primary"]["value"], "0.7212", "ModRank reaches 0.7212 concordance"),
        ("primary SD", a1["primary"]["sd_over_seeds"], "0.0048", "standard deviation 0.0048"),
        ("over the best verified entry", a1["comparisons"]["vs_DIMAF_published"]["gap"], "+0.0422",
         "$+0.0422$ above the best entry whose folds we verified"),
        ("vs SurvPath", holm["ours vs SurvPath+clinical (seed-matched)"]["delta"], "+0.0340",
         "exceeds SurvPath by $+0.0340$ ($[-0.0056, +0.0731]$, $p=0.084$"),
        ("vs SurvPath, CI low", holm["ours vs SurvPath+clinical (seed-matched)"]["ci95"][0], "-0.0056",
         "($[-0.0056, +0.0731]$"),
        ("vs SurvPath, CI high", holm["ours vs SurvPath+clinical (seed-matched)"]["ci95"][1], "+0.0731",
         "($[-0.0056, +0.0731]$"),
        ("vs SurvPath, p", holm["ours vs SurvPath+clinical (seed-matched)"]["p_raw"], "0.084", "$p=0.084$"),
        ("vs PIBD", holm["ours vs PIBD+clinical (best-val checkpoint)"]["delta"], "+0.0249",
         "PIBD by $+0.0249$ ($[-0.0181, +0.0675]$, $p=0.241$)"),
        ("vs PIBD, CI low", holm["ours vs PIBD+clinical (best-val checkpoint)"]["ci95"][0], "-0.0181",
         "($[-0.0181, +0.0675]$"),
        ("vs PIBD, CI high", holm["ours vs PIBD+clinical (best-val checkpoint)"]["ci95"][1], "+0.0675",
         "($[-0.0181, +0.0675]$"),
        ("vs PIBD, p", holm["ours vs PIBD+clinical (best-val checkpoint)"]["p_raw"], "0.241", "$p=0.241$"),
        ("our coefficients", pc["ours"]["parameters"], "1,049", "with 1,049 parameters against 24.7"),
        ("SurvPath parameters", pc["survpath"]["parameters"], "24.7M", "against 24.7 and 26.9 million"),
        ("PIBD parameters", pc["pibd"]["parameters"], "26.9M", "against 24.7 and 26.9 million"),
        # fusion on identical inputs
        ("concatenated, mean", fu["concatenated_ridge_cox"]["mean"], "0.6864", "$0.6864 \\pm 0.0063$"),
        ("concatenated, SD", fu["concatenated_ridge_cox"]["sd"], "0.0063", "$0.6864 \\pm 0.0063$"),
        ("stacked, mean", fu["stacked_learned_weights"]["mean"], "0.7106", "$0.7106 \\pm 0.0044$"),
        ("stacked, SD", fu["stacked_learned_weights"]["sd"], "0.0044", "$0.7106 \\pm 0.0044$"),
        ("vs concatenated", fu["paired"]["ours_vs_concatenated"]["mean"], "+0.0340",
         "exceeds the concatenation by $+0.0340$ ($[+0.0014, +0.0674]$, $p=0.042$)"),
        ("vs concatenated, CI low", fu["paired"]["ours_vs_concatenated"]["ci95"][0], "+0.0014",
         "($[+0.0014, +0.0674]$"),
        ("vs concatenated, CI high", fu["paired"]["ours_vs_concatenated"]["ci95"][1], "+0.0674",
         "($[+0.0014, +0.0674]$"),
        ("vs concatenated, p", fu["paired"]["ours_vs_concatenated"]["p_two_sided"], "0.042", "$p=0.042$"),
        ("vs stacked", fu["paired"]["ours_vs_stacked"]["mean"], "+0.0108",
         "the stacking by $+0.0108$ ($[-0.0052, +0.0269]$)"),
        ("vs stacked, CI low", fu["paired"]["ours_vs_stacked"]["ci95"][0], "-0.0052", "($[-0.0052, +0.0269]$)"),
        ("vs stacked, CI high", fu["paired"]["ours_vs_stacked"]["ci95"][1], "+0.0269", "($[-0.0052, +0.0269]$)"),
        ("resplits vs concatenated", rs["resplits"]["ours_minus_concat"]["mean"], "+0.0464",
         "($+0.0464 \\pm 0.0158$)"),
        ("resplits vs concatenated, SD", rs["resplits"]["ours_minus_concat"]["sd"], "0.0158",
         "($+0.0464 \\pm 0.0158$)"),
        ("resplits vs stacked", rs["resplits"]["ours_minus_stacked"]["mean"], "+0.0111",
         "($+0.0111 \\pm 0.0065$)"),
        ("resplits vs stacked, SD", rs["resplits"]["ours_minus_stacked"]["sd"], "0.0065",
         "($+0.0111 \\pm 0.0065$)"),
        ("resplits vs clinical", rs["resplits"]["ours_minus_clinical"]["mean"], "+0.0549",
         "($+0.0549 \\pm 0.0083$) in 24 of 24"),
        ("resplits vs clinical, SD", rs["resplits"]["ours_minus_clinical"]["sd"], "0.0083",
         "($+0.0549 \\pm 0.0083$) in 24 of 24"),
        ("site-grouped sites", rs["site_grouped_cv"]["sites_total"], "33", "(33 sites, five folds)"),
        ("site-grouped ModRank", rs["site_grouped_cv"]["pooled"]["ours"], "0.7210",
         "ModRank holds 0.7210 while the stacking falls to 0.6855 and the concatenation to 0.6487"),
        ("site-grouped stacked", rs["site_grouped_cv"]["pooled"]["stacked"], "0.6855",
         "the stacking falls to 0.6855"),
        ("site-grouped concatenated", rs["site_grouped_cv"]["pooled"]["concat"], "0.6487",
         "the concatenation to 0.6487"),
        ("encoder parity", ep["chief"]["full_method"], "0.7009",
         "SurvPath itself consumed, ModRank reaches 0.7009"),
        ("encoder contribution", a1["arms_seed0"]["OURS"] - ep["chief"]["full_method"], "0.022",
         "the encoder's contribution near 0.022"),
        ("corrected bar", sn["A_selection_null"]["bars"]["empirical_max_q95"]["bar"], "0.7716",
         "the correctly measured bar of 0.7716"),
        # added value and its inflation
        ("ModRank over grade", av["ours"]["added_over_weak"], "+0.1204",
         "$+0.1204$ ($[+0.0684, +0.1708]$) for ModRank"),
        ("ModRank over grade, CI low", av["ours"]["added_over_weak_boot"]["ci95"][0], "+0.0684",
         "($[+0.0684, +0.1708]$)"),
        ("ModRank over grade, CI high", av["ours"]["added_over_weak_boot"]["ci95"][1], "+0.1708",
         "($[+0.0684, +0.1708]$)"),
        ("SurvPath over grade", av["survpath"]["added_over_weak"], "+0.0670", "$+0.0670$ for SurvPath"),
        ("PIBD over grade", av["pibd_best_val"]["added_over_weak"], "+0.0841", "$+0.0841$ for PIBD"),
        ("ModRank over stage", av["ours"]["added_over_stage"], "+0.0594",
         "fall to $+0.0594$, $+0.0253$ and $+0.0345$"),
        ("SurvPath over stage", av["survpath"]["added_over_stage"], "+0.0253",
         "fall to $+0.0594$, $+0.0253$ and $+0.0345$"),
        ("PIBD over stage", av["pibd_best_val"]["added_over_stage"], "+0.0345",
         "fall to $+0.0594$, $+0.0253$ and $+0.0345$"),
        ("D ModRank", av["ours"]["D"], "+0.0610", "$D = +0.0610$ ($[+0.0284, +0.0955]$)"),
        ("D ModRank, CI low", av["ours"]["D_boot"]["ci95"][0], "+0.0284", "($[+0.0284, +0.0955]$)"),
        ("D ModRank, CI high", av["ours"]["D_boot"]["ci95"][1], "+0.0955", "($[+0.0284, +0.0955]$)"),
        ("D SurvPath", av["survpath"]["D"], "+0.0417", "$+0.0417$ ($[+0.0196, +0.0668]$)"),
        ("D SurvPath, CI low", av["survpath"]["D_boot"]["ci95"][0], "+0.0196", "($[+0.0196, +0.0668]$)"),
        ("D SurvPath, CI high", av["survpath"]["D_boot"]["ci95"][1], "+0.0668", "($[+0.0196, +0.0668]$)"),
        ("D PIBD", av["pibd_best_val"]["D"], "+0.0496", "$+0.0496$ ($[+0.0238, +0.0769]$)"),
        ("D PIBD, CI low", av["pibd_best_val"]["D_boot"]["ci95"][0], "+0.0238", "($[+0.0238, +0.0769]$)"),
        ("D PIBD, CI high", av["pibd_best_val"]["D_boot"]["ci95"][1], "+0.0769", "($[+0.0238, +0.0769]$)"),
        ("D over resplits", rs["resplits"]["D_added_value_inflation"]["mean"], "+0.0628",
         "positive in 24 of 24 re-partitions ($+0.0628 \\pm 0.0086$)"),
        ("D over resplits, SD", rs["resplits"]["D_added_value_inflation"]["sd"], "0.0086",
         "($+0.0628 \\pm 0.0086$)"),
        ("ModRank vs clinical, Holm", holm["ours vs clinical alone"]["p_holm"], "0.065", "(Holm $q=0.065$)"),
        # calibration and utility
        ("mean predicted 2-year risk", cal["ours"]["by_horizon"]["24.0"]["mean_predicted_risk"], "0.342",
         "was 0.342 against an observed 0.341"),
        ("observed 2-year risk", cal["ours"]["by_horizon"]["24.0"]["observed_km_risk"], "0.341",
         "against an observed 0.341"),
        ("calibration slope", cal["ours"]["calibration_slope"], "0.770",
         "the calibration slope was 0.770 (clinical model 0.849)"),
        ("clinical calibration slope", cal["clinical"]["calibration_slope"], "0.849", "(clinical model 0.849)"),
        ("IPA 1 year", cal["ours"]["by_horizon"]["12.0"]["ipa"], "0.098", "was 0.098 at one year"),
        ("IPA 2 years", cal["ours"]["by_horizon"]["24.0"]["ipa"], "0.152", "and 0.152 at two years"),
        ("clinical IPA 1 year", cal["clinical"]["by_horizon"]["12.0"]["ipa"], "0.047", "against 0.047 and 0.079"),
        ("clinical IPA 2 years", cal["clinical"]["by_horizon"]["24.0"]["ipa"], "0.079", "against 0.047 and 0.079"),
        ("IPA 3 years", cal["ours"]["by_horizon"]["36.0"]["ipa"], "0.105", "equal at three years (0.105)"),
        ("clinical IPA 3 years", cal["clinical"]["by_horizon"]["36.0"]["ipa"], "0.105",
         "equal at three years (0.105)"),
        ("net benefit at 40%", cal["ours"]["net_benefit_24m"]["0.4"], "0.136", "0.136 against 0.090"),
        ("clinical net benefit at 40%", cal["clinical"]["net_benefit_24m"]["0.4"], "0.090",
         "0.136 against 0.090"),
        ("treat-all net benefit at 40%", cal["treat_all_net_benefit_24m"]["0.4"], "-0.098",
         "and $-0.098$ at 40\\%"),
        ("IPA difference at 2 years", ui["ipa_modrank_minus_clinical"]["24.0"]["mean"], "+0.0724",
         "the two-year difference of $+0.0724$ has a 95\\% interval of $[-0.0054, +0.1502]$"),
        ("its CI low", ui["ipa_modrank_minus_clinical"]["24.0"]["ci95"][0], "-0.0054",
         "$[-0.0054, +0.1502]$"),
        ("its CI high", ui["ipa_modrank_minus_clinical"]["24.0"]["ci95"][1], "+0.1502",
         "$[-0.0054, +0.1502]$"),
        ("net benefit difference at 40%", ui["net_benefit_24m_modrank_minus_clinical"]["0.4"]["mean"],
         "+0.0454", "At 40\\% the difference from the clinical model is $+0.0454$ ($[+0.0006, +0.0915]$)"),
        ("its CI low", ui["net_benefit_24m_modrank_minus_clinical"]["0.4"]["ci95"][0], "+0.0006",
         "($[+0.0006, +0.0915]$)"),
        ("its CI high", ui["net_benefit_24m_modrank_minus_clinical"]["0.4"]["ci95"][1], "+0.0915",
         "($[+0.0006, +0.0915]$)"),
        # the modalities
        ("slide-clinical Spearman", am["redundancy"]["spearman_titan_vs_clinical"], "0.2445",
         "slide and clinical scores is 0.2445"),
        ("incumbent-clinical Spearman", am["redundancy"]["spearman_survpath_vs_clinical"], "0.1924",
         "and the clinical score 0.1924"),
        ("clinically tied pairs", tie["pairs"], "1,213", "(1,213 of 24,219)"),
        ("comparable pairs", tie["pairs_total"], "24,219", "(1,213 of 24,219)"),
        ("clinical arm on the tied pairs", tie["arms"]["clinical"], "0.4959", "at chance (0.4959)"),
        ("slide arm on the tied pairs", tie["arms"]["titan"], "0.6220", "the slide arm reaches 0.6220"),
        ("doubly tied pairs", dtt["pairs"], "3,893", "On the 3,893 pairs that neither"),
        ("slide on the doubly tied pairs", dtt["slide"], "0.6212", "the slide arm reaches 0.6212"),
        ("transcriptome on them", dtt["omics"], "0.5543", "sit at 0.5543 and 0.5634"),
        ("clinical on them", dtt["clinical"], "0.5634", "sit at 0.5543 and 0.5634"),
        ("transcriptome arm vs basal axis", b4["basal"]["spearman"], "0.3595",
         "axes (Spearman 0.3595, $-0.372$ and 0.3672)"),
        ("transcriptome arm vs luminal axis", b4["luminal"]["spearman"], "-0.372",
         "axes (Spearman 0.3595, $-0.372$ and 0.3672)"),
        ("transcriptome arm vs EMT axis", b4["EMT"]["spearman"], "0.3672",
         "axes (Spearman 0.3595, $-0.372$ and 0.3672)"),
        ("slide arm, largest axis correlation", max(abs(v["spearman"]) for v in b3.values()), "0.2245",
         "is at most 0.2245 in magnitude"),
        ("ModRank, basal-leaning half", bio["basal_leaning"]["OURS"], "0.7109", "(0.7109 and 0.6997"),
        ("ModRank, luminal-leaning half", bio["luminal_leaning"]["OURS"], "0.6997", "(0.7109 and 0.6997"),
        # the external cohorts
        ("GSE32894 n", g32["n"], "224", "a Swedish cohort of 224 patients with 25 deaths"),
        ("GSE32894 events", g32["events"], "25", "224 patients with 25 deaths"),
        ("GSE32894 D", g32["D"]["point"], "+0.0581", "$D = +0.0581$ ($[+0.0175, +0.0971]$, $p=0.007$)"),
        ("GSE32894 D, CI low", g32["D"]["ci95"][0], "+0.0175", "($[+0.0175, +0.0971]$"),
        ("GSE32894 D, CI high", g32["D"]["ci95"][1], "+0.0971", "($[+0.0175, +0.0971]$"),
        ("GSE32894 D, p", g32["D"]["p_two_sided"], "0.007", "$p=0.007$"),
        ("GSE31684 n", g31["n"], "93", "cohort of 93 patients with 38 deaths"),
        ("GSE31684 events", g31["events"], "38", "93 patients with 38 deaths"),
        ("GSE31684 D", g31["D"]["point"], "+0.0253", "without reaching significance ($+0.0253$)"),
        ("GSE31684 D without pre-cystectomy chemotherapy",
         g31["sensitivity_without_pre_cystectomy_chemotherapy"]["D"]["point"], "+0.1094",
         "gave $+0.1094$ ($[+0.0424, +0.1774]$)"),
        ("the same, CI low", g31["sensitivity_without_pre_cystectomy_chemotherapy"]["D"]["ci95"][0],
         "+0.0424", "($[+0.0424, +0.1774]$)"),
        ("the same, CI high", g31["sensitivity_without_pre_cystectomy_chemotherapy"]["D"]["ci95"][1],
         "+0.1774", "($[+0.0424, +0.1774]$)"),
        ("patients excluded by that sensitivity", g31["n"]
         - g31["sensitivity_without_pre_cystectomy_chemotherapy"]["n"], "3",
         "excludes the three patients given"),
        ("GSE32894 ModRank", g32["cv"]["modrank"]["mean"], "0.8779",
         "(0.8779 against 0.8672 for the stage block and 0.8747 for stacking)"),
        ("GSE32894 stage block", g32["cv"]["clinical_stage"]["mean"], "0.8672", "0.8672 for the stage block"),
        ("GSE32894 stacked", g32["cv"]["stacked"]["mean"], "0.8747", "0.8747 for stacking"),
        ("GSE32894 vs transcriptome, Holm", g32["paired_vs_modrank"]["transcriptome"]["p_holm"], "0.004",
         "correction ($p=0.004$ for both)"),
        ("GSE32894 vs concatenated, Holm", g32["paired_vs_modrank"]["concatenated"]["p_holm"], "0.004",
         "correction ($p=0.004$ for both)"),
        ("GSE48075 n", g48["n"], "73", "a third cohort of 73 patients"),
        ("GSE48075 ModRank", g48["cv"]["modrank"]["mean"], "0.6384", "(0.6384 against 0.6322)"),
        ("GSE48075 stage block", g48["cv"]["clinical_stage"]["mean"], "0.6322", "(0.6384 against 0.6322)"),
        ("GSE31684 transcriptome", g31["cv"]["transcriptome"]["mean"], "0.4742", "were uninformative (0.4742)"),
        ("GSE31684 ModRank", g31["cv"]["modrank"]["mean"], "0.5898", "(0.5898 against 0.6590)"),
        ("GSE31684 stage block", g31["cv"]["clinical_stage"]["mean"], "0.6590", "(0.5898 against 0.6590)"),
        # Methods
        ("deployable variant", cr["transductive_vs_inductive"]["primary_inductive"], "0.7247",
         "reaches 0.7247 against the transductive 0.7212"),
        ("deployable variant, p",
         cr["transductive_vs_inductive"]["paired_case_bootstrap_inductive_minus_transductive"]["p_two_sided"],
         "0.69", "(paired $p=0.69$)"),
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
                  "Supplementary Table~4": "tab:perseed", "Supplementary Note~5": "note:calibration",
                  "Supplementary Table~8": "tab:calibration", "Supplementary Table~5": "tab:ties",
                  "Supplementary Note~3": "note:arms", "Supplementary Fig.~2": "fig:landscape",
                  "Supplementary Table~6": "tab:subtype", "Supplementary Fig.~3": "fig:biology",
                  "Supplementary Figs.~4": "fig:cohort", "and~5 draw": "fig:cases",
                  "Supplementary Note~4": "note:geo", "Supplementary Table~7": "tab:geo",
                  "Supplementary Table~9": "tab:tripod"}
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
    # 0.2 mm), resolved through the manuscript's own graphicspath order.
    _gp = re.search(r"\\graphicspath\{((?:\{[^}]*\})+)\}", open(NC_TEX).read()).group(1)
    _dirs = [os.path.normpath(os.path.join(os.path.dirname(NC_TEX), d_)) for d_ in re.findall(r"\{([^}]*)\}", _gp)]
    for _doc in (NC_TEX, os.path.join(ROOT, "paper", "nc-submission", "supplementary.tex")):
        for _m in re.finditer(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", open(_doc).read()):
            if _doc != NC_TEX:
                continue            # the SI carries earlier figures at their own width; checked by eye
            _hit = next((os.path.join(d_, _m.group(1)) for d_ in _dirs
                         if os.path.isfile(os.path.join(d_, _m.group(1)))), None)
            if _hit is None:
                fails.append("NC includes %s, which exists in no graphicspath directory" % _m.group(1))
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
    m = re.search(r"\\begin\{tabular\}\{[^}]*\}(.*?)\\end\{tabular\}(?:(?!\\begin\{table\}).)*?"
                  r"\\label\{%s\}" % re.escape(label), src, re.S)
    if not m:
        return None
    body = re.sub(r"\\(toprule|midrule|bottomrule)|\\cmidrule\([^)]*\)\{[^}]*\}", "", m.group(1))
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
    return n


if __name__ == "__main__":
    sys.exit(main())
