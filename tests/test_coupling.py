import numpy as np

from snoredreamt import conditions, coupling


def test_rayleigh_formula_and_null_calibration():
    a = np.array([0.1, 0.2, 0.3, 2, 3])
    p, z = coupling.rayleigh(a)
    n, R = len(a), abs(np.exp(1j * a).sum())
    assert np.isclose(z, R ** 2 / n)                      # Rayleigh z = R^2 / n
    # Zar (1999) approximation used by CircStat
    assert np.isclose(p, np.exp(np.sqrt(1 + 4 * n + 4 * (n ** 2 - R ** 2)) - (1 + 2 * n)))
    # under H0 the p-values are ~uniform: about 5 % below 0.05
    rng = np.random.default_rng(1)
    ps = np.array([coupling.rayleigh(rng.uniform(-np.pi, np.pi, 30))[0] for _ in range(4000)])
    assert 0.035 < np.mean(ps < 0.05) < 0.065


def test_uniform_and_concentrated():
    rng = np.random.default_rng(0)
    assert coupling.rayleigh(rng.uniform(-np.pi, np.pi, 2000))[0] > 0.01
    conc = rng.vonmises(0.5, 4, 200)
    assert coupling.rayleigh(conc)[0] < 1e-10
    assert abs(coupling.circ_mean(conc) - 0.5) < 0.1
    assert coupling.vtest(conc, 0.5)[0] < 1e-10
    assert coupling.vtest(conc, 0.5 + np.pi)[0] > 0.99
    assert abs(coupling.ppc(rng.uniform(-np.pi, np.pi, 5000))) < 0.01


def test_paint_and_runs():
    m = np.zeros(20, bool)
    m[5:8] = True
    s, e = conditions.runs(m)
    assert s.tolist() == [5] and e.tolist() == [7]
    p = conditions.paint(20, s, e, 2, 3)
    assert np.flatnonzero(p).tolist() == list(range(3, 11))


def test_peth_percentages():
    h, edges = coupling.peth([100, 105, 300], [100, 300], fs=100, half_s=0.5, bin_s=0.05)
    assert np.isclose(h.sum(), 100)
