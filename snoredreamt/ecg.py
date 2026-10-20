"""R-peak detection (Pan-Tompkins-style) for 100 Hz ECG."""
import numpy as np
from scipy import signal, stats


def detect_rpeaks(ecg, fs, band=(5.0, 25.0), min_rr_s=0.33, win_s=10.0):
    """Return R-peak sample indices.

    Band-pass -> derivative^2 -> 150 ms moving integration -> peaks above 0.3 x the local (10 s) 95th percentile,
    >= min_rr apart; each peak refined to the extremum of the band-passed ECG (polarity from skewness) within 80 ms.
    """
    x = np.nan_to_num(np.asarray(ecg, float))
    sos = signal.butter(3, band, btype='band', fs=fs, output='sos')
    bp = signal.sosfiltfilt(sos, x)
    sign = 1.0 if stats.skew(bp[: min(len(bp), 3600 * int(fs))]) >= 0 else -1.0
    bp = sign * bp
    integ = np.convolve(np.gradient(bp) ** 2, np.ones(int(0.15 * fs)) / (0.15 * fs), 'same')
    w = int(win_s * fs)
    nw = int(np.ceil(len(integ) / w))
    ref = np.array([np.percentile(integ[i * w:(i + 1) * w], 95) for i in range(nw)])
    ref = np.repeat(ref, w)[: len(integ)]
    cand, _ = signal.find_peaks(integ, height=0.3 * ref, distance=int(min_rr_s * fs))
    r = int(0.08 * fs)
    peaks = np.array([max(c - r, 0) + np.argmax(bp[max(c - r, 0):c + r + 1]) for c in cand], dtype=np.int64)
    return np.unique(peaks), sign


def instantaneous_hr(rpeaks, fs, n, max_rr_s=2.0):
    """Heart rate (bpm) interpolated to every sample; NaN outside plausible RR intervals."""
    rr = np.diff(rpeaks) / fs
    t = rpeaks[1:]
    ok = (rr > 0.3) & (rr < max_rr_s)
    hr = np.full(n, np.nan)
    if ok.sum() > 2:
        hr = np.interp(np.arange(n), t[ok], 60 / rr[ok], left=np.nan, right=np.nan)
    return hr
