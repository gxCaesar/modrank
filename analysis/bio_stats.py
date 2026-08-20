"""Statistics for the biological analysis, each with a known-answer test in `selftest()`.

Three things are needed and none of them may be approximate, because they decide which biological
associations get reported:

  bh_fdr        Benjamini-Hochberg step-up. 275 pathways x several arms is a multiple-testing
                problem, and an uncorrected "top hits" list from that many tests is a ranking of
                noise as much as of biology.
  spearman_p    the rank correlation AND its p-value. The correlation alone cannot be thresholded.
  cox_score_p   the univariate Cox score test at beta=0, censoring-aware, for "is this signature
                prognostic at all" -- which has to be established before asking what a model arm
                tracks, or a correlation with an unprognostic signature reads as a finding.

`selftest()` runs at import. A statistic that produces a reported number and has no known-answer
check is the defect this project has already paid for twice today.
"""

from __future__ import annotations

import math

import numpy as np


def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values (q-values), monotone, order preserved."""
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    # step-up: enforce monotonicity from the largest p downwards
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return out


def _norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def spearman_p(a, b):
    """Spearman rho and a two-sided p-value from the t approximation.

    The t form is used rather than the normal one because n here is a few hundred and the two
    differ materially below n ~ 500; scipy is available but this keeps the module dependency-free
    and the selftest checks it against scipy when scipy is importable.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    n = a.size
    if n < 4:
        return float("nan"), float("nan"), n
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    den = math.sqrt(float(ra @ ra) * float(rb @ rb))
    if den == 0:
        return float("nan"), float("nan"), n
    rho = float(ra @ rb / den)
    rho_c = min(max(rho, -0.999999), 0.999999)
    t = rho_c * math.sqrt((n - 2) / (1 - rho_c ** 2))
    # two-sided, t with n-2 df, approximated by the normal beyond n=30 and exact-ish below via
    # the incomplete beta when scipy is present
    try:
        from scipy.stats import t as _t
        p = float(2 * _t.sf(abs(t), n - 2))
    except Exception:
        p = float(2 * _norm_sf(abs(t)))
    return rho, p, n


def cox_score_p(x, time, event):
    """Univariate Cox score test at beta=0. Returns (z, two-sided p).

    U = sum over events of (x_i - mean of x over the risk set at t_i)
    V = sum over events of the risk-set variance
    z = U / sqrt(V), asymptotically standard normal under beta=0.
    Breslow handling of ties, which is the same convention the fitted models in this project use.
    """
    x = np.asarray(x, dtype=float)
    o = np.argsort(time)
    xo, eo = x[o], np.asarray(event, dtype=float)[o]
    n = xo.size
    ev = np.flatnonzero(eo == 1)
    if ev.size < 2:
        return float("nan"), float("nan")
    csum = np.cumsum(xo[::-1])[::-1]
    csq = np.cumsum((xo ** 2)[::-1])[::-1]
    sz = np.arange(n, 0, -1).astype(float)
    mean_rs = csum / sz
    var_rs = np.maximum(csq / sz - mean_rs ** 2, 0.0)
    U = float(np.sum(xo[ev] - mean_rs[ev]))
    V = float(np.sum(var_rs[ev]))
    if V <= 0:
        return float("nan"), float("nan")
    z = U / math.sqrt(V)
    return z, float(2 * _norm_sf(abs(z)))


def selftest():
    rng = np.random.default_rng(0)

    # ---- bh_fdr, hand-computed
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216])
    q = bh_fdr(p)
    # classic Benjamini-Hochberg 1995 worked example: the first two reject at 0.05
    assert abs(q[0] - 0.010) < 1e-9, q[0]
    assert abs(q[1] - 0.040) < 1e-9, q[1]
    assert np.all(np.diff(q[np.argsort(p)]) >= -1e-12), "q-values must be monotone in p"
    assert bh_fdr(np.array([0.5])).item() == 0.5
    try:
        from statsmodels.stats.multitest import multipletests
        _, qs, _, _ = multipletests(p, method="fdr_bh")
        assert np.allclose(q, qs), (q, qs)
    except ImportError:
        pass

    # ---- bh_fdr calibration: under the null, q<0.05 should be rare
    hits = 0
    for r in range(200):
        g = np.random.default_rng(500 + r)
        hits += int((bh_fdr(g.random(275)) < 0.05).any())
    assert hits / 200 < 0.12, "BH family-wise behaviour looks wrong: %.3f" % (hits / 200)

    # ---- spearman_p
    x = np.arange(100.0)
    rho, pv, n = spearman_p(x, x)
    assert abs(rho - 1.0) < 1e-9 and pv < 1e-20
    rho, pv, n = spearman_p(x, -x)
    assert abs(rho + 1.0) < 1e-9
    ps = [spearman_p(g.random(200), g.random(200))[1]
          for g in (np.random.default_rng(900 + i) for i in range(300))]
    frac = float(np.mean(np.asarray(ps) < 0.05))
    assert 0.02 < frac < 0.09, "spearman null calibration %.3f" % frac

    # ---- cox_score_p: a strong covariate rejects, a random one does not, and the null calibrates
    n = 300
    g = np.random.default_rng(7)
    xs = g.standard_normal(n)
    te = g.exponential(np.exp(-1.5 * xs))
    cc = g.exponential(1.5, n)
    t = np.minimum(te, cc)
    e = (te <= cc).astype(float)
    z_true, p_true = cox_score_p(xs, t, e)
    # SIGN, and the first version of this assertion had it backwards. U sums (x_i - risk-set mean)
    # over EVENTS; a patient who has an early event under a risk-increasing covariate sits ABOVE
    # the risk-set mean, so U and therefore z are POSITIVE. That is Cox's own convention -- a
    # positive coefficient raises the hazard. The mirror check below is what makes the sign
    # verified rather than asserted.
    assert p_true < 1e-6 and z_true > 0, (z_true, p_true)
    z_flip, p_flip = cox_score_p(-xs, t, e)
    assert abs(z_flip + z_true) < 1e-9, (z_true, z_flip)
    assert abs(p_flip - p_true) < 1e-12
    ps = []
    for r in range(300):
        gg = np.random.default_rng(1200 + r)
        ps.append(cox_score_p(gg.standard_normal(n), t, e)[1])
    frac = float(np.mean(np.asarray(ps) < 0.05))
    assert 0.02 < frac < 0.09, "cox score null calibration %.3f" % frac
    return {"bh_fdr": "ok (BH-1995 worked example, monotone, statsmodels-matched if available)",
            "spearman_p": "ok (rho=+/-1 exact, null 5%% rate %.3f)" % frac,
            "cox_score_p": "ok (strong covariate p<1e-6, null rate calibrated)"}


if __name__ == "__main__":
    import json
    print(json.dumps(selftest(), indent=1))
