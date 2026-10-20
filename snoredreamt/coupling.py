"""Circular statistics (CircStat semantics) and coupling measures."""
import numpy as np
from scipy import stats


def circ_mean(a):
    return np.angle(np.mean(np.exp(1j * np.asarray(a))))


def circ_r(a):
    return np.abs(np.mean(np.exp(1j * np.asarray(a))))


def rayleigh(a):
    """CircStat circ_rtest: returns (p, z)."""
    a = np.asarray(a)
    n = len(a)
    R = n * circ_r(a)
    z = R ** 2 / n
    p = np.exp(np.sqrt(1 + 4 * n + 4 * (n ** 2 - R ** 2)) - (1 + 2 * n))
    return min(p, 1.0), z


def vtest(a, direction=0.0):
    """CircStat circ_vtest: H0 uniform vs H1 non-uniform with mean `direction`. Returns (p, V)."""
    a = np.asarray(a)
    n = len(a)
    R = n * circ_r(a)
    V = R * np.cos(circ_mean(a) - direction)
    u = V * np.sqrt(2 / n)
    return 1 - stats.norm.cdf(u), V


def ppc(a):
    """Pairwise phase consistency (Vinck 2010), unbiased by n."""
    a = np.asarray(a)
    n = len(a)
    if n < 2:
        return np.nan
    s = np.abs(np.sum(np.exp(1j * a))) ** 2
    return (s - n) / (n * (n - 1))


def equalizer(phase_samples, max_n=200_000, rng=None):
    """Phase equalisation: map phases through their empirical CDF over the analysis mask so that time is uniformly
    distributed over phase (otherwise e.g. cycle-based phase over-represents the longer expiratory half and
    Rayleigh tests on event phases are biased). The anchor is kept: phase 0 (detected peak) maps to 0.
    Returns a function phase -> equalised phase."""
    x = np.asarray(phase_samples, float)
    x = x[np.isfinite(x)]
    if len(x) > max_n:
        x = (rng or np.random.default_rng(0)).choice(x, max_n, replace=False)
    srt = np.sort(x)
    F = lambda p: np.searchsorted(srt, p, side='right') / len(srt)  # noqa: E731
    f0 = F(0.0)
    return lambda p: np.angle(np.exp(1j * 2 * np.pi * (F(np.asarray(p, float)) - f0)))


def null_mvl(phase, mask_idx, n, n_draws, rng):
    """MVL of n random time points drawn from mask_idx (samples where phase is defined), n_draws times."""
    draws = rng.choice(mask_idx, size=(n_draws, n), replace=True)
    return np.abs(np.mean(np.exp(1j * phase[draws].astype(float)), axis=1))


def phase_rate(event_phase, phase, mask_idx, fs, n_bins=18):
    """Occupancy-normalised event rate (events/min) per phase bin."""
    edges = np.linspace(-np.pi, np.pi, n_bins + 1)
    occ = np.histogram(phase[mask_idx].astype(float), edges)[0] / fs / 60
    cnt = np.histogram(event_phase, edges)[0]
    with np.errstate(invalid='ignore', divide='ignore'):
        return cnt / occ, occ, cnt


def peth(event_times, ref_times, fs, half_s=1.5, bin_s=0.05):
    """Schreiner 2023 PETH: events relative to reference times (e.g. inhalation peaks), % of events in window."""
    edges = np.arange(-half_s, half_s + bin_s / 2, bin_s)
    ev = np.sort(np.asarray(event_times))
    lags = []
    for r in ref_times:
        lo, hi = np.searchsorted(ev, [r - half_s * fs, r + half_s * fs])
        lags.append((ev[lo:hi] - r) / fs)
    lags = np.concatenate(lags) if lags else np.array([])
    h = np.histogram(lags, edges)[0].astype(float)
    return 100 * h / h.sum() if h.sum() else h, edges
