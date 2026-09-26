"""Which of the 275 pathway groups the slide arm tracks, for the transcriptome glyph in Figures 1b
and 8a.

The glyph lights the pathways whose Spearman correlation with the slide arm's out-of-fold score has
a Benjamini-Hochberg q below 0.05. The per-pathway correlations are in the reporting dump; the
p values follow from the t distribution with n - 2 degrees of freedom, and the count this gives is
asserted equal to the count the biology results record (59), so the glyph cannot drift from the
number printed beside it.
"""

from __future__ import annotations

import json
import os

import numpy as np
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))


def slide_tracked():
    """Boolean mask over the 275 pathways in the reporting dump's order."""
    dump = json.load(open(os.path.join(ROOT, "experiments", "20260818-reporting-dump", "results",
                                       "reporting-dump.json")))
    bio = json.load(open(os.path.join(ROOT, "experiments", "20260817-blca-confirm", "results",
                                      "biology.json")))
    n = dump["n"]
    rho = np.array([p["rho_slide"] for p in dump["pathways"]])
    t = rho * np.sqrt((n - 2) / (1 - rho ** 2))
    pv = 2 * stats.t.sf(np.abs(t), n - 2)
    order = np.argsort(pv)
    m = len(pv)
    q = np.empty(m)
    q[order] = np.minimum.accumulate((pv[order] * m / np.arange(1, m + 1))[::-1])[::-1]
    mask = q < 0.05
    want = bio["B3_what_the_image_arm_tracks"]["pathways_with_BH_q_below_0.05"]
    assert m == bio["B2_survival_association"]["pathways_tested"], "pathway count drifted"
    assert int(mask.sum()) == want, "recomputed %d slide-tracked pathways, results say %d" % (
        int(mask.sum()), want)
    return mask
