"""Amendment A2: the Breslow objective and average-rank percentiles in analysis/blca_common.py.

Known answers for the two functions every reported concordance passes through: an O(n^2) Breslow
brute force, invariance to row order, equality with the pre-A2 objective when no times tie, the
analytic gradient against finite differences, and ranks that do not depend on sort order.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "analysis"))
import blca_common as bc  # noqa: E402


def _data(seed=0, n=120, ties=True):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, 3))
    t = np.round(rng.exponential(10.0, n), 0 if ties else 6) + 0.5
    e = (rng.random(n) < 0.6).astype(float)
    return X, t, e


def _brute(b, X, t, e, alpha):
    eta = X @ b
    ll = sum(eta[i] - np.log(np.exp(eta[t >= t[i]]).sum()) for i in np.flatnonzero(e == 1))
    return -ll + 0.5 * alpha * float(b @ b)


def _obj(b, X, t, e, alpha, a2):
    o = np.argsort(t, kind="stable")
    return bc.cox_objective(b, X[o], t[o], e[o], alpha, a2=a2)


def test_breslow_equals_brute_force_with_ties():
    X, t, e = _data()
    assert len(np.unique(t)) < len(t)
    b = np.array([0.3, -0.2, 0.1])
    assert abs(_obj(b, X, t, e, 1.0, True)[0] - _brute(b, X, t, e, 1.0)) < 1e-9


def test_breslow_is_invariant_to_row_order_and_the_old_objective_is_not():
    X, t, e = _data()
    b = np.array([0.3, -0.2, 0.1])
    p = np.random.default_rng(5).permutation(len(t))
    new = [_obj(b, X[q], t[q], e[q], 1.0, True)[0] for q in (np.arange(len(t)), p)]
    assert abs(new[0] - new[1]) < 1e-9
    old = [bc.cox_objective(b, X[q][np.argsort(t[q])], np.sort(t[q]), e[q][np.argsort(t[q])], 1.0, a2=False)[0]
           for q in (np.arange(len(t)), p)]
    assert abs(old[0] - old[1]) > 1e-6          # the defect A2 removes


def test_no_ties_old_and_new_agree():
    X, t, e = _data(ties=False)
    assert len(np.unique(t)) == len(t)
    b = np.array([0.3, -0.2, 0.1])
    assert abs(_obj(b, X, t, e, 1.0, True)[0] - _obj(b, X, t, e, 1.0, False)[0]) < 1e-9


def test_gradient_matches_finite_differences():
    X, t, e = _data()
    b = np.array([0.3, -0.2, 0.1])
    g = _obj(b, X, t, e, 1.0, True)[1]
    h = 1e-6
    fd = [(_obj(b + h * d, X, t, e, 1.0, True)[0] - _obj(b - h * d, X, t, e, 1.0, True)[0]) / (2 * h)
          for d in np.eye(3)]
    assert np.max(np.abs(g - np.array(fd))) < 1e-5


def test_average_ranks():
    v = np.array([0.5, 0.1, 0.5, 0.9, 0.1, 0.5])
    assert bc.avg_rank(v).tolist() == [3.0, 0.5, 3.0, 5.0, 0.5, 3.0]
    w = np.random.default_rng(2).random(50)
    assert np.array_equal(bc.avg_rank(w), np.argsort(np.argsort(w)).astype(float))


def test_pct_is_order_independent_under_a2():
    if not bc.A2:
        return
    v = np.array([3, 1, 3, 2, 3, 1, 2], float)
    p = np.random.default_rng(3).permutation(len(v))
    assert np.allclose(bc.pct(v)[p], bc.pct(v[p]))


def test_score_test_breslow_risk_sets():
    import s5_omics_arm
    X, t, e = _data()
    ev = np.flatnonzero(e == 1)
    # brute force: U and V at beta = 0 with full risk sets {j : t_j >= t_i}
    U = np.zeros(X.shape[1]); V = np.zeros(X.shape[1])
    for i in ev:
        R = t >= t[i]
        U += X[i] - X[R].mean(0)
        V += np.maximum((X[R] ** 2).mean(0) - X[R].mean(0) ** 2, 0.0)
    want = np.abs(U) / np.sqrt(np.maximum(V, 1e-12))
    if bc.A2:
        assert np.allclose(s5_omics_arm.score_test(X, t, e), want, atol=1e-10)
        p = np.random.default_rng(9).permutation(len(t))
        assert np.allclose(s5_omics_arm.score_test(X[p], t[p], e[p]), want, atol=1e-10)


def test_ranks_switch():
    v = np.array([2.0, 1.0, 2.0, 3.0])
    if bc.A2:
        assert bc.ranks(v).tolist() == [1.5, 0.0, 1.5, 3.0]


def test_tuned_stacking_ridge_selection():
    # the added comparator: ridge from the grid by pooled held-out concordance, first strict maximum
    rng = np.random.default_rng(4)
    n = 150
    Z = rng.random((n, 3))
    t = np.round(rng.exponential(10.0, n) * np.exp(-1.5 * Z[:, 0]), 0) + 0.5
    e = (rng.random(n) < 0.6).astype(float)
    inner = np.array_split(rng.permutation(n), 3)
    got_ridge, got_w = bc.stack_weights_tuned(Z, t, e, inner)
    scores = []
    for al in bc.STACK_RIDGES:
        ip = np.zeros(n)
        for ite in inner:
            itr = np.setdiff1d(np.arange(n), ite)
            mu = Z[itr].mean(0)
            ip[ite] = (Z[ite] - mu) @ bc.cox_fit(Z[itr] - mu, t[itr], e[itr], al)
        scores.append(bc.cindex(ip, t, e))
    assert bc.STACK_RIDGES == (0.1, 1.0, 10.0, 100.0)
    assert got_ridge == bc.STACK_RIDGES[int(np.argmax(scores))]
    assert np.allclose(got_w, bc.cox_fit(Z - Z.mean(0), t, e, got_ridge))
    assert got_w[0] > 0


if __name__ == "__main__":
    # runnable without pytest: every test in this file, in order, with a count
    import inspect
    tests = [(n, f) for n, f in sorted(globals().items()) if n.startswith("test_") and inspect.isfunction(f)]
    for n, f in tests:
        f()
        print("ok  %s" % n)
    print("%d tests passed (BLCA_A2=%s)" % (len(tests), "1" if bc.A2 else "0"))
