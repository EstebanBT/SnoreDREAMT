"""Python port of the slythm SO and spindle detectors (canonical copy: /home/esteban/BreathREM/+slythm, 2026-04-07).

Faithful to detectSOs.m / detectSOs_splitStages.m / detectSpindles.m with filtType 'but',
criterionCenter 'median', criterionVar 'scaleCenter', no padding, no merging, no double threshold.

Indices are 0-based samples (MATLAB indices minus 1). Deliberate deviations from the MATLAB code:
- Per-stage SO thresholds (splitStages) are used only for stages with >= `min_stage_minutes` of valid data;
  sparser stages (DREAMT often has <10 min N3) borrow the N2 threshold instead of a noisy own threshold.
- Events are returned sorted in time (splitStages concatenates stage blocks unsorted).
- Spindle `duration` is (end - start + 1) / fs, consistent with the length criterion; MATLAB stores (end - start) / fs.
- Detection can be run on other stages (REM/wake) with thresholds taken from NREM, for specificity checks.
- Optional spectral-peak criterion for spindles (config `spectral_peak`), added because slythm's 1.5x median
  threshold was not specific in DREAMT (wake/N1 rates >= N2 from alpha/beta activity).
Not ported: thresType 'average' (buggy in MATLAB), findEvtFree, searchIndivEvtPeak, raw_trough.
"""
import numpy as np
import pandas as pd
from scipy import signal


# ----------------------------------------------------------------------------- filtering
def _butter_stable(order, wn, btype):
    """FieldTrip 'reduce' instability fix: lower the order until all poles are inside the unit circle."""
    while True:
        b, a = signal.butter(order, wn, btype=btype)
        if np.all(np.abs(np.roots(a)) < 1) or order == 1:
            return b, a
        order -= 1


def _filtfilt_matlab(b, a, x):
    """MATLAB filtfilt: odd reflection with 3*(nfilt-1) samples, steady-state initial conditions."""
    return signal.filtfilt(b, a, x, padtype='odd', padlen=3 * (max(len(a), len(b)) - 1))


def ft_highpass(x, fs, f, order):
    """ft_preproc_highpassfilter(x, fs, f, order, 'but', 'twopass', 'reduce'): demeaned, mean not restored."""
    b, a = _butter_stable(order, f / (fs / 2), 'high')
    return _filtfilt_matlab(b, a, x - np.nanmean(x))


def ft_lowpass(x, fs, f, order):
    """ft_preproc_lowpassfilter(x, fs, f, order, 'but', 'twopass', 'reduce'): demeaned, mean restored."""
    b, a = _butter_stable(order, f / (fs / 2), 'low')
    m = np.nanmean(x)
    return _filtfilt_matlab(b, a, x - m) + m


def bandpass_but(x, fs, band, orders):
    return ft_lowpass(ft_highpass(np.asarray(x, float), fs, band[0], orders[0]), fs, band[1], orders[1])


def envelope_rms(x, n):
    """MATLAB envelope(x, n, 'rms') upper envelope: moving RMS, window shrinks at edges.
    Even n: n/2-1 samples back, n/2 forward; odd n: (n-1)/2 each side."""
    m = np.mean(x)
    xc = x - m
    back, fwd = (n // 2 - 1, n // 2) if n % 2 == 0 else (n // 2, n // 2)
    c = np.concatenate([[0.0], np.cumsum(xc * xc)])
    i = np.arange(len(x))
    lo = np.maximum(i - back, 0)
    hi = np.minimum(i + fwd, len(x) - 1)
    return np.sqrt((c[hi + 1] - c[lo]) / (hi - lo + 1)) + m


# ----------------------------------------------------------------------------- helpers
def _window_all(mask_cum, start, end):
    """True where mask is True for every sample in [start, end] (inclusive). mask_cum = cumsum of ~mask with leading 0."""
    return (mask_cum[end + 1] - mask_cum[start]) == 0


def _bad_cumsum(ok):
    return np.concatenate([[0], np.cumsum(~ok)])


def _runs(mask):
    d = np.diff(np.concatenate([[0], mask.astype(np.int8), [0]]))
    return np.flatnonzero(d > 0), np.flatnonzero(d < 0) - 1


# ----------------------------------------------------------------------------- slow oscillations
def so_candidates(xf):
    """Zero-crossing candidates on the band-passed signal (detectSOs.m l.219-262)."""
    d = np.diff(np.sign(xf))
    p2n = np.flatnonzero(d < 0)
    n2p = np.flatnonzero(d > 0)
    if len(p2n) and len(n2p) and n2p[0] < p2n[0]:
        n2p = n2p[1:]
    k = min(len(p2n) - 1, len(n2p))
    p2n, nxt, n2p = p2n[:k], p2n[1:k + 1], n2p[:k]
    keep = (n2p >= p2n) & (nxt >= n2p)          # guard against non-alternating crossings (exact zeros)
    p2n, nxt, n2p = p2n[keep], nxt[keep], n2p[keep]
    tr_pos = np.array([s + np.argmin(xf[s:e + 1]) for s, e in zip(p2n, n2p)], dtype=np.int64)
    pk_pos = np.array([s + np.argmax(xf[s:e + 1]) for s, e in zip(n2p, nxt)], dtype=np.int64)
    c = pd.DataFrame({'start': p2n, 'pivot': n2p, 'end': nxt, 'trough': tr_pos, 'peak': pk_pos})
    c['trough_amp'] = xf[tr_pos]
    c['peak_amp'] = xf[pk_pos]
    c['pk2pk'] = np.abs(c.trough_amp) + c.peak_amp
    return c


def detect_sos(x, stage, valid, fs, p, detect_stages=None, thresholds=None):
    """Detect slow oscillations.

    x: raw EEG (µV); stage: int codes per sample; valid: bool artifact-free mask; p: config['detect']['so'].
    detect_stages: stages to detect in (default p['stages']); thresholds: {stage: (thr_trough, thr_pk2pk)}
    to reuse (e.g. NREM thresholds for a REM specificity run). Returns (events DataFrame, thresholds dict).
    """
    xf = bandpass_but(x, fs, p['bpfreq'], p['filt_order'])
    c = so_candidates(xf)
    c['duration'] = (c.end - c.start) / fs
    len_ok = (c.duration >= p['criterion_len'][0]) & (c.duration <= p['criterion_len'][1])
    scale = p['criterion_scale']

    if thresholds is None:
        thresholds, minutes = {}, {}
        for s in p['stages']:
            ok = (stage == s) & valid
            minutes[s] = ok.sum() / fs / 60
            inwin = _window_all(_bad_cumsum(ok), c.start.to_numpy(), c.end.to_numpy())
            pool = c[inwin & len_ok]
            if len(pool):
                thresholds[s] = (scale * np.median(np.abs(pool.trough_amp)), scale * np.median(pool.pk2pk))
        ref = p['stages'][0]   # N2
        for s in list(thresholds):
            if s != ref and minutes[s] < p.get('min_stage_minutes', 0) and ref in thresholds:
                thresholds[s] = thresholds[ref]

    out = []
    for s in (detect_stages or p['stages']):
        thr = thresholds.get(s, thresholds.get(p['stages'][0]))
        if thr is None:
            continue
        ok = (stage == s) & valid
        inwin = _window_all(_bad_cumsum(ok), c.start.to_numpy(), c.end.to_numpy())
        sel = inwin & len_ok & (np.abs(c.trough_amp) >= thr[0]) & (c.pk2pk >= thr[1])
        out.append(c[sel].assign(stage=s))
    ev = pd.concat(out).sort_values('start').reset_index(drop=True) if out else c.iloc[:0].assign(stage=[])
    return ev, thresholds


# ----------------------------------------------------------------------------- spindles
def detect_spindles(x, stage, valid, fs, p, detect_stages=None, threshold=None):
    """Detect spindles (RMS envelope above median*scale of in-stage, valid envelope samples).

    Returns (events DataFrame, threshold). With `threshold` given (e.g. from NREM) it is reused, for REM checks.
    """
    xf = bandpass_but(x, fs, p['bpfreq'], p['filt_order'])
    env = envelope_rms(xf, int(round(p['envelope_win_s'] * fs)))
    if threshold is None:
        pool_ok = np.isin(stage, p['stages']) & valid
        threshold = p['criterion_scale'] * np.median(env[pool_ok])
    ok = np.isin(stage, detect_stages or p['stages']) & valid
    sta, end = _runs((env >= threshold) & ok)
    n = end - sta + 1
    keep = (n >= round(p['criterion_len'][0] * fs)) & (n <= round(p['criterion_len'][1] * fs))
    sta, end = sta[keep], end[keep]
    rows = []
    for s, e in zip(sta, end):
        seg = xf[s:e + 1]
        pk = signal.find_peaks(seg)[0]
        rows.append((s, e, s + np.argmax(seg), s + np.argmin(seg), seg.max(), seg.min(),
                     s + np.argmax(env[s:e + 1]), env[s:e + 1].max(), env[s:e + 1].mean(),
                     fs / np.mean(np.diff(pk)) if len(pk) > 1 else np.nan))
    ev = pd.DataFrame(rows, columns=['start', 'end', 'peak', 'trough', 'peak_amp', 'trough_amp',
                                     'env_peak', 'env_max', 'env_mean', 'freq'])
    ev['mid'] = np.round((ev.start + ev.end) / 2).astype(np.int64)
    ev['duration'] = (ev.end - ev.start + 1) / fs
    ev['stage'] = stage[ev.start.to_numpy()] if len(ev) else []
    if p.get('spectral_peak'):
        ev['peak_ok'] = spectral_peak(np.asarray(x, float), ev, fs, **p['spectral_peak'])
        ev = ev[ev.peak_ok].reset_index(drop=True)
    return ev, threshold


def spectral_peak(x, ev, fs, sigma=(12, 16), below=(8, 11), above=(18, 25), min_win_s=1.0, nfft=256):
    """Not in slythm (added for DREAMT specificity): the event's sigma-band power density must exceed both the
    adjacent alpha and beta bands. Periodogram of the raw EEG (Hann), window = event, at least `min_win_s`, centred."""
    half_min = int(min_win_s * fs / 2)
    ok = np.zeros(len(ev), bool)
    for k, (s, e) in enumerate(zip(ev.start.to_numpy(), ev.end.to_numpy())):
        c, half = (s + e) // 2, max((e - s + 1) // 2, half_min)
        seg = x[max(c - half, 0): c + half]
        seg = (seg - seg.mean()) * np.hanning(len(seg))
        P = np.abs(np.fft.rfft(seg, max(nfft, len(seg)))) ** 2
        ff = np.fft.rfftfreq(max(nfft, len(seg)), 1 / fs)
        pw = [P[(ff >= lo) & (ff <= hi)].mean() for lo, hi in (sigma, below, above)]
        ok[k] = pw[0] > pw[1] and pw[0] > pw[2]
    return ok


# ----------------------------------------------------------------------------- complexes
def find_complexes(sos, spindles, fs, max_lag_s, spindle_time='trough'):
    """SO-spindle complexes (replicateSchreinerCoupling.m l.244-257): an SO whose trough is followed by a
    spindle time point within (0, max_lag) s (strict, as in the MATLAB code). Returns the SO rows with the
    first matching spindle's index."""
    st = np.sort(spindles[spindle_time].to_numpy())
    order = np.argsort(spindles[spindle_time].to_numpy())
    tr = sos.trough.to_numpy()
    i = np.searchsorted(st, tr, side='right')            # first spindle strictly after the trough
    has = (i < len(st))
    lag = np.full(len(tr), np.inf)
    lag[has] = st[i[has]] - tr[has]
    is_cx = lag < max_lag_s * fs
    cx = sos[is_cx].copy()
    cx['spindle_idx'] = spindles.index.to_numpy()[order[i[is_cx]]]
    cx['spindle_lag_s'] = lag[is_cx] / fs
    return cx.reset_index(drop=True)
