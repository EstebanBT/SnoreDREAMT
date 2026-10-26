"""Faithful Python ports of CircStat 2012a functions (P. Berens, J Stat Softw 2009; ~/toolbox/CircStat2012a).

Only the functions needed for Step 5 are ported, following the MATLAB code line by line (including its tie handling
and assumption checks). Angles in radians.
"""
import warnings

import numpy as np
from scipy import stats


def circ_r(alpha, w=None, d=0.0):
    alpha = np.asarray(alpha, float)
    w = np.ones_like(alpha) if w is None else np.asarray(w, float)
    r = np.abs(np.sum(w * np.exp(1j * alpha))) / np.sum(w)
    if d != 0:
        r = d / 2 / np.sin(d / 2) * r
    return r


def circ_mean(alpha, w=None):
    alpha = np.asarray(alpha, float)
    w = np.ones_like(alpha) if w is None else np.asarray(w, float)
    ok = ~(np.isnan(alpha) | np.isnan(w))
    return np.angle(np.sum(w[ok] * np.exp(1j * alpha[ok])))


def circ_dist(x, y):
    """Pairwise difference x_i - y_i around the circle (angle(exp(1i*x)./exp(1i*y)))."""
    return np.angle(np.exp(1j * np.asarray(x, float)) / np.exp(1j * np.asarray(y, float)))


def circ_dist2(x, y=None):
    """All pairwise differences x_i - y_j (matrix)."""
    x = np.asarray(x, float).ravel()
    y = x if y is None else np.asarray(y, float).ravel()
    return np.angle(np.exp(1j * x)[:, None] / np.exp(1j * y)[None, :])


def circ_median(alpha):
    beta = np.mod(np.asarray(alpha, float).ravel(), 2 * np.pi)
    n = len(beta)
    dd = circ_dist2(beta, beta)
    m1 = np.sum(dd >= 0, axis=0)
    m2 = np.sum(dd <= 0, axis=0)
    dm = np.abs(m1 - m2)
    if n % 2 == 1:
        idx = [int(np.argmin(dm))]
        m = dm[idx[0]]
    else:
        m = dm.min()
        idx = list(np.flatnonzero(dm == m)[:2])
    if m > 1:
        warnings.warn('Ties detected.')
    md = circ_mean(beta[idx])
    if abs(circ_dist(circ_mean(beta), md)) > abs(circ_dist(circ_mean(beta), md + np.pi)):
        md = np.mod(md + np.pi, 2 * np.pi)
    return md


def circ_kappa(alpha=None, w=None, R=None):
    """Kappa of a von Mises approximation; pass R directly to mirror circ_kappa(rw) with a scalar (N = 1)."""
    if R is None:
        alpha = np.asarray(alpha, float).ravel()
        N = len(alpha)
        R = circ_r(alpha, w) if N > 1 else alpha[0]
    else:
        N = 1
    if R < 0.53:
        kappa = 2 * R + R ** 3 + 5 * R ** 5 / 6
    elif R < 0.85:
        kappa = -.4 + 1.39 * R + 0.43 / (1 - R)
    else:
        kappa = 1 / (R ** 3 - 4 * R ** 2 + 3 * R)
    if 1 < N < 15:
        kappa = max(kappa - 2 * (N * kappa) ** -1, 0) if kappa < 2 else (N - 1) ** 3 * kappa / (N ** 3 + N)
    return kappa


def _ww_assumption_ok(rw, n_mean):
    """CircStat circ_wwtest checkAssumption: True if the test is applicable."""
    if n_mean >= 11:
        return rw >= .45
    if 7 <= n_mean < 11:
        return rw >= .5
    if 5 <= n_mean < 7:
        return rw >= .55
    return False


def circ_wwtest(alpha1, alpha2):
    """Two-sample Watson-Williams test for equal mean directions. Returns (pval, F, assumption_ok)."""
    alpha = np.concatenate([np.asarray(alpha1, float).ravel(), np.asarray(alpha2, float).ravel()])
    idx = np.concatenate([np.ones(len(alpha1)), 2 * np.ones(len(alpha2))])
    w = np.ones_like(alpha)
    s, n = 2, w.sum()
    pn = np.array([np.sum(idx == u) for u in (1, 2)], float)
    pr = np.array([circ_r(alpha[idx == u]) for u in (1, 2)])
    r = circ_r(alpha, w)
    rw = np.sum(pn * pr) / n
    ok = _ww_assumption_ok(rw, pn.mean())
    kk = circ_kappa(R=rw)
    beta = 1 + 3 / (8 * kk)
    A = np.sum(pr * pn) - r * n
    B = n - np.sum(pr * pn)
    F = beta * (n - s) * A / (s - 1) / B
    return 1 - stats.f.cdf(F, s - 1, n - s), F, ok


def circ_cmtest(alpha1, alpha2):
    """Two-sample common median test (non-parametric). Returns (pval, grand median, P statistic)."""
    alpha = np.concatenate([np.asarray(alpha1, float).ravel(), np.asarray(alpha2, float).ravel()])
    idx = np.concatenate([np.ones(len(alpha1)), 2 * np.ones(len(alpha2))])
    N = len(idx)
    med = circ_median(alpha)
    n = np.array([np.sum(idx == u) for u in (1, 2)], float)
    m = np.array([np.sum(circ_dist(alpha[idx == u], med) < 0) for u in (1, 2)], float)
    if np.any(n < 10):
        warnings.warn('Test not applicable. Sample size in at least one group to small.')
    M = m.sum()
    P = (N ** 2 / (M * (N - M))) * np.sum(m ** 2 / n) - N * M / (N - M)
    return 1 - stats.chi2.cdf(P, 1), med, P


def circ_confmean(alpha, xi=0.05, w=None, d=0.0):
    alpha = np.asarray(alpha, float)
    w = np.ones_like(alpha) if w is None else np.asarray(w, float)
    r = circ_r(alpha, w, d)
    n = w.sum()
    R = n * r
    c2 = stats.chi2.ppf(1 - xi, 1)
    if np.sqrt(c2 / 2 / n) < r < .9:
        t = np.sqrt((2 * n * (2 * R ** 2 - n * c2)) / (4 * n - c2))
    elif r >= .9:
        t = np.sqrt(n ** 2 - (n ** 2 - R ** 2) * np.exp(c2 / n))
    else:
        return np.nan                     # MATLAB: warning 'Requirements for confidence levels not met.'
    return np.arccos(t / R)


def circ_mtest(alpha, dir_, xi=0.05):
    """One-sample test for mean direction: h = True if dir_ lies outside the (1 - xi) CI of the mean."""
    mu = circ_mean(alpha)
    t = circ_confmean(alpha, xi)
    h = bool(abs(circ_dist(dir_, mu)) > t) if np.isfinite(t) else False
    return h, mu, mu + t, mu - t


def circ_medtest(alpha, md):
    """Exact sign test for the median direction being md.

    Deviation from CircStat 2012a: when n1 == n2 the two binomial tails [0:min] and [max:n] overlap in MATLAB and the
    shared term is counted twice, giving p > 1; the result is capped at 1 here."""
    alpha = np.asarray(alpha, float).ravel()
    n = len(alpha)
    d = circ_dist(alpha, md)
    n1, n2 = int(np.sum(d < 0)), int(np.sum(d > 0))
    ks = np.concatenate([np.arange(0, min(n1, n2) + 1), np.arange(max(n1, n2), n + 1)])
    return float(min(1.0, np.sum(stats.binom.pmf(ks, n, 0.5))))
