import numpy as np

from snoredreamt import circstat as cs


def test_median_basic():
    a = np.array([0.1, 0.2, 0.3, 0.4, 2.0])
    assert np.isclose(cs.circ_median(a), 0.3)
    # wraps around 0: median near 2*pi - 0.05 or -0.05
    b = np.array([-0.2, -0.1, 0.05, 0.1, 6.0])
    assert abs(cs.circ_dist(cs.circ_median(b), 0.0)) < 0.2


def test_ww_and_cm_detect_shift_and_not_null():
    rng = np.random.default_rng(0)
    a1, a2 = rng.vonmises(0, 4, 200), rng.vonmises(0.8, 4, 200)
    p_ww, F, ok = cs.circ_wwtest(a1, a2)
    assert ok and p_ww < 1e-6
    assert cs.circ_cmtest(a1, a2)[0] < 1e-4
    ps = [cs.circ_cmtest(rng.vonmises(0, 1, 60), rng.vonmises(0, 1, 60))[0] for _ in range(400)]
    assert 0.02 < np.mean(np.array(ps) < 0.05) < 0.09          # roughly nominal under H0


def test_ww_assumption_flag_for_weak_concentration():
    rng = np.random.default_rng(1)
    _, _, ok = cs.circ_wwtest(rng.uniform(-np.pi, np.pi, 100), rng.uniform(-np.pi, np.pi, 100))
    assert not ok


def test_mtest_and_medtest():
    rng = np.random.default_rng(2)
    a = rng.vonmises(0.6, 3, 80)
    h, mu, ul, ll = cs.circ_mtest(a, 0.0)
    assert h and ll < mu < ul
    assert not cs.circ_mtest(a, 0.6)[0]
    assert cs.circ_medtest(a, 0.0) < 0.01
    assert cs.circ_medtest(a, 0.6) > 0.05


def test_corrcl_and_k_groups():
    rng = np.random.default_rng(3)
    x = rng.normal(size=200)
    a = np.angle(np.exp(1j * (0.8 * x + rng.vonmises(0, 4, 200))))     # phase depends on x
    rho, p = cs.circ_corrcl(a, x)
    assert rho > 0.5 and p < 1e-6
    rho0, p0 = cs.circ_corrcl(rng.uniform(-np.pi, np.pi, 200), x)
    assert p0 > 0.001
    g = [rng.vonmises(m, 3, 60) for m in (0, 0, 1.0)]
    assert cs.circ_cmtest_k(g)[0] < 1e-3 and cs.circ_wwtest_k(g)[0] < 1e-3
    assert np.isclose(cs.circ_cmtest(g[0], g[1])[0], cs.circ_cmtest_k(g[:2])[0])
