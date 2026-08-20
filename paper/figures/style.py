"""Figure style for this manuscript, with the two assertions wired into the build.

The generating rule: lightness is an ORDERED channel, so it is spent only on ordered quantities.
Two palettes follow and mixing them is the error to avoid -- a pale colour must never mean "a
category" in one panel and "high on the scale" in the next.

GEOMETRY IS NOW LOCKED. The venue was decided at J8 on 2026-08-17: Briefings in Bioinformatics,
Problem Solving Protocol. BiB states no figure width, no maximum figure count and no dpi rule of
its own, so what binds is (a) OUP's general artwork rules -- vector PDF with all fonts embedded,
sans-serif lettering of uniform size, which this module already satisfies -- and (b) the page the
figure has to fit on. Figures are authored at 6.30 in, the manuscript's text width, and inserted
at \textwidth so insertion scale is 1.0000 and source pt equals rendered pt.

ONE CONFLICT IS RECORDED AND DELIBERATELY NOT RESOLVED HERE. OUP's general guidance says "Colour
figures should be supplied in CMYK not RGB colours"; BiB has been online-only open access since
January 2024 and states no colour-space rule of its own. Converting a luminance-ordered ladder to
CMYK would move the luminance gaps the ladder exists for, so this is a decision to take with the
editorial office, not one to let a colour-management default take silently. Everything below is
sRGB. See briefings_in_bioinformatics_lock.json -> j8_completion_20260817.figures_technical.
"""

from __future__ import annotations

import matplotlib as mpl
import numpy as np

# -- nominal: method family, arm, cohort. Chroma fixed at C* 35, lightness free, hues spread wide,
#    plus exactly one pure neutral for the baseline category.
NOMINAL = ["#6D98D3", "#BC6A79", "#C98E66", "#89A467", "#3CAC9C"]
NEUTRAL = "#8C8C8C"

# -- ordered ladder: anything with a direction. CIELCh, chroma near-constant, lightness alone makes
#    the ladder. Dark -> pale; the palest is fills only.
LADDER = ["#667DB8", "#CB7F9F", "#A4A8B6", "#CAAEDF", "#9CD0F3", "#F3CFED", "#DEEFFF"]

INK = "#28303F"

# Risk tertiles are an ORDERED quantity, so they take three entries from the ladder -- chosen by
# maximising the minimum adjacent luminance gap rather than by eye. Computed in
# figure-source-data.json -> F4; the gap is 0.272, five times the 0.05 greyscale floor.
RISK = {"high": "#667DB8", "middle": "#CAAEDF", "low": "#DEEFFF"}

# Named roles, so a colour is never chosen at a call site.
ROLE = {"ours": "#6D98D3", "clinical": "#89A467", "slide": "#C98E66",
        "omics": "#3CAC9C", "competitor": "#BC6A79", "published": "#8C8C8C"}


def _luminance(hexes):
    """Relative luminance, sRGB, for the greyscale-survival assertion."""
    out = []
    for h in hexes:
        r, g, b = (int(h[i:i + 2], 16) / 255.0 for i in (1, 3, 5))
        f = lambda c: c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
        out.append(0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b))
    return np.array(out)


def assert_palettes():
    """Two mirror-image checks. Both have caught real errors and both run at import time.

    1. The ORDERED ladder must survive greyscale: min adjacent luminance gap >= 0.05.
    2. The NOMINAL set must NOT smuggle an ordering into an unordered variable: luminance spread
       <= 0.22, and it must share no colour with the ladder.
    """
    lad = _luminance(LADDER)
    gaps = np.abs(np.diff(lad))
    assert gaps.min() >= 0.05, (
        "ordered ladder fails greyscale: min adjacent luminance gap %.4f < 0.05" % gaps.min())

    nom = _luminance(NOMINAL + [NEUTRAL])
    spread = nom.max() - nom.min()
    assert spread <= 0.22, (
        "nominal palette luminance spread %.4f > 0.22 -- lightness is ordering an unordered "
        "variable" % spread)

    shared = set(c.lower() for c in NOMINAL + [NEUTRAL]) & set(c.lower() for c in LADDER)
    assert not shared, "nominal and ordered palettes share %s" % shared

    # 3. the risk triple must itself survive greyscale, and must come FROM the ladder -- a risk
    #    tertile is ordered, so encoding it in a nominal hue would be the same error in reverse.
    trip = [RISK["high"], RISK["middle"], RISK["low"]]
    assert all(c in LADDER for c in trip), "risk colours must come from the ordered ladder"
    tl = _luminance(trip)
    assert np.abs(np.diff(tl)).min() >= 0.05, (
        "risk triple fails greyscale: min gap %.4f" % np.abs(np.diff(tl)).min())
    assert (np.diff(tl) > 0).all(), "risk must run dark (high) to pale (low), monotonically"


def apply():
    assert_palettes()
    mpl.rcParams.update({
        # Arial everywhere INCLUDING math. mathtext has its own font stack and ignores
        # font.sans-serif entirely, so a figure that is Arial everywhere a person looks still
        # renders every $...$ in DejaVu unless these four are set.
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold", "mathtext.sf": "Arial",
        # every size pinned, or the 10pt default leaks in through anything not set
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
        "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
        "xtick.color": INK, "ytick.color": INK,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "legend.frameon": False,
        "figure.dpi": 200,
        # live text, editable by a typesetter. NOTE: this emits Identity-H and would VIOLATE an
        # AAAI submission; it is required by the journal families this manuscript targets. Same
        # rcParam, opposite sign -- a module shared between a conference and a journal submission
        # must not share it.
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.transparent": False,
    })


def save(fig, path):
    """Save at the authored size. bbox_inches=None does NOT disable tight bbox.

    `savefig` reads None as "argument not given" and falls back to rcParams["savefig.bbox"], so the
    trap survives being passed explicitly. The rc_context is what actually turns it off, and with
    it native size equals figsize and source pt equals rendered pt.
    """
    import matplotlib.pyplot as plt  # noqa: F401
    with mpl.rc_context({"savefig.bbox": None, "savefig.pad_inches": 0}):
        fig.savefig(path)
        # SVG beside every PDF. OUP accepts both and a typesetter can edit either, but they are
        # NOT interchangeable here: pdf.fonttype 42 embeds a subset in the PDF, while SVG keeps
        # text as text with no font travelling with it, so the SVG renders correctly only where
        # Arial exists. Both are shipped and the PDF is the one \includegraphics reads.
        fig.savefig(path[:-4] + ".svg" if path.endswith(".pdf") else path + ".svg")
    return path
