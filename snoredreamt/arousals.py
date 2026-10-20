"""Rule-based EEG arousal detection (single EEG + chin EMG + optional ECG), ported from MATLAB.

Copyright (C) 2017-2019 Isaac Fernández-Varela
Copyright (C) 2017-2019 Diego Alvarez-Estevez
Python port and 100 Hz / DREAMT adaptations (C) 2026 the SnoreDREAMT authors.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details: <http://www.gnu.org/licenses/>.

This module is a derivative work of https://github.com/bigsasi/arousals-detection (GPL-3; local copy in
third_party/arousals-detection), which implements:

[1] I. Fernández-Varela, D. Alvarez-Estevez, E. Hernández-Pereira, V. Moret-Bonillo, "A simple and robust method
    for the automatic scoring of EEG arousals in polysomnographic recordings", Computers in Biology and Medicine
    87:77-86, 2017.
[2] D. Alvarez-Estevez, I. Fernández-Varela, "Large-scale validation of an automatic EEG arousal detection
    algorithm using different heterogeneous databases", Sleep Medicine 57:6-14, 2019.

Ported files: arousalEvents/{arousalDetection, findArousalEvents, findBetaArousals, findAlphaArousals,
removeSpindles, removeAwakeEvents, updateWithPowerBeta, updateWithPowerAlpha, updateWithOverlappingEmg,
updateWithNonContinuousOverlappingEmg, updateWithNonOverlappingEmg, updateWithAmplitude}.m,
events/{findPowerEvents, findAmplitudeEvents, getSignalEvents, getEventIndexes, mergeEventsCloserThan,
mergeEvents, updateAlphaBetaOverlapping}.m, signalProcessing/{signalPower, signalAmplitude, signalCycles,
baseline}.m, artifacts/ecg/{ecgReliabilityCheck, ecgConditioning}.m, libs/Polyman/*.m, loadFile.m (filter
order) and the parameter values of config.ini / conf_filters.xml. The code (not the papers) is followed where
they differ, e.g. alpha trigger 1.5 x baseline (paper 2.5), alpha power-extension threshold 4000 (i.e. disabled;
config.ini "Default was 4"), max 20 s per epoch (paper 15 s), min alpha duration 4 s (paper 4.5 s), beta events
merged across gaps < 1 s while alpha events are not merged (config names AlphaUnion/BetaUnion are applied
crosswise in findPowerEvents.m), and the EMG overlap reference median includes the event itself (paper: excludes).
Off-by-one index quirks of the MATLAB code (1-based indices built from 0-based times) are reproduced on purpose.

Deviations from the MATLAB code
-------------------------------
1. Beta band: MATLAB uses 16 Hz .. fs/2 (64-128 Hz in the original data). At 100 Hz we use 16-45 Hz
   (`beta_band`), i.e. below Nyquist with margin for the resampling roll-off. The power-event code and the
   bin-selection rules of signalPower.m are otherwise unchanged (3 s Hamming window, 0.2 s step).
2. Mains: the original notches 50/60 Hz (and 100 Hz) on EEG and EMG. At fs = 100 Hz a 60 Hz notch is impossible
   (60 > Nyquist) and 50 Hz is the Nyquist bin, so instead we notch the *aliases* of the 60 Hz harmonics seen in
   DREAMT: a narrow line at 20 Hz (120/180 Hz aliased) is present in C4-M1 and/or CHIN of most subjects and a
   40 Hz line (60/240 Hz aliased) in a few (sleep spectra of every third subject checked). The same Polyman 2nd-order notch
   (bandwidth 1 Hz) is applied at `notch_freqs` = (20, 40) Hz to EEG and EMG. ecgReliabilityCheck.m's hard-coded
   50 Hz notch on the ECG is skipped (it would sit on the Nyquist bin, ahead of an 8-16 Hz band-pass).
3. EMG high-pass 15 Hz (Polyman 1st order) as in conf_filters.xml; with fs = 100 Hz only 15-50 Hz of EMG
   (plus aliased power) remains.
4. updateWithNonOverlappingEmg.m uses 0.125 s amplitude windows (0.12 s / 0.48 s blocks at 125 and 250 Hz so that
   windows are whole samples). 0.125 s is 12.5 samples at 100 Hz, so the 125/250 Hz branch (0.12 s / 0.48 s) is
   used whenever 0.125 * fs is not an integer.
5. R-peaks: libs/qrs/nqrsdetect.m is not ported; `snoredreamt.ecg.detect_rpeaks` (Pan-Tompkins style, zero-phase)
   is used, or `rpeaks` passed by the caller. Beats whose 0.5 s window is cut by the recording edge are skipped
   (MATLAB would error).
6. Hypnogram: the per-sample DREAMT stage vector replaces the 30 s hypnogram. Wake/REM look-ups take the stage at
   the event's start/end sample, which respects DREAMT's epoch grid (offset from sample 0 by up to 30 s); the
   ">20 s within one epoch" duration rule uses the same scored epoch grid (MATLAB: grid anchored at t = 0).
   Events whose start or end lies in unscored samples (stage -1, lights-on/preparation) are removed like wake
   (`remove_unscored`; MATLAB only knows W = 0).
7. MATLAB's updateWithPower{Beta,Alpha}.m set `endSample = startSample + round(added_duration * fs)` (the added
   duration, not the total) - a bug that later EMG updates partly overwrite. It is reproduced by default
   (`matlab_endsample_quirk=True`) because the published thresholds were tuned with it.
8. signalCycles.m turns zero first differences into NaN (0/0) and its zero-filling loop is a no-op; reproduced.
9. Edge guards (index clipping at the recording end, empty alpha/beta event sets) where MATLAB would error.
10. `valid` (artifact mask) is not part of the original algorithm: it never changes detection; it only adds the
    column `valid_frac` (fraction of event samples that are valid) so callers can filter.
11. MATLAB round() (half away from zero) is reproduced with `_mround`; numpy's banker's rounding is never used
    for index arithmetic.

Not in the original either (known limitations): no 10 s of prior *sleep* check after wake (only the 10 s gap to
the previous arousal), no theta band, single EEG derivation. Onsets are the centre of the first 3 s power window
above threshold: on synthetic bursts they are ~1 s early for strong arousals and up to ~0.3 s (rarely ~2 s) late for
weak ones. On DREAMT the detector counts ~0.38x the clinical Arousal Index (median ratio; Spearman 0.65, n = 98).
"""
import numpy as np
import pandas as pd
from scipy import signal as sps

DEFAULT_PARAMS = {
    # [PowerEvents]  (config.ini)
    'power_alpha_threshold': 1.5,
    'power_alpha_union': 1.0,          # sic: applied to the BETA event vector in findPowerEvents.m
    'power_beta_threshold': 2.0,
    'power_beta_union': 0.0,           # sic: applied to the ALPHA event vector
    'power_window_step': 0.2,
    'power_window_size': 3.0,
    'power_background_time': 10.0,
    # [AmplitudeEvents] (chin EMG)
    'amp_threshold': 2.5,
    'amp_window_step': 0.1,
    'amp_window_size': 0.2,
    'amp_background_time': 30.0,
    'amp_merge_near': 0.5,
    # [UpdateWith]
    'update_with_alpha_events': False,
    'update_with_alpha_events_union_time': 2.0,
    'update_with_beta_events': True,
    'update_with_beta_events_union_time': 1.0,
    # [UpdateWithPower]
    'beta_power_threshold': 3.0,
    'beta_power_limit_time': 3,
    'beta_power_increase_time': 1,
    'beta_power_previous_seconds': 10,
    'alpha_power_threshold': 4000.0,   # config.ini: "Default was 4" -> extension effectively disabled
    'alpha_power_limit_time': 1,
    'alpha_power_increase_time': 1,
    'alpha_power_previous_seconds': 3,
    # [UpdateWithAmplitude]
    'amplitude_update_threshold': 4.0,
    # [UpdateWithFutureEmg]
    'future_emg_limit': 30,
    'future_emg_threshold': 2.5,
    # [UpdateWithEmg]
    'emg_update_window_size': 0.1,
    'emg_update_window_step': 0.1,
    'emg_update_threshold': 1.4,
    # constants hard-coded in the MATLAB functions
    'alpha_band': (8.0, 12.0),
    'beta_band': (16.0, 45.0),         # MATLAB: (16, fs/2); see deviation 1
    'spindle_bandpass': (1.0, 32.0),
    'min_beta_duration': 3.0,
    'min_alpha_duration': 4.0,
    'max_duration_per_epoch': 20.0,
    'alpha_min_emg_s': 1.0,
    'min_stable_sleep': 10.0,
    'edge_exclusion_s': 60.0,
    'max_update_s': 15.0,
    'emg_amp_floor_uv': 15.0,          # updateWithNonOverlappingEmg.m: max(previousAmp, 15)
    # preprocessing (conf_filters.xml + 100 Hz adaptations)
    'emg_highpass': 15.0,
    'notch_freqs': (20.0, 40.0),       # aliases of 60 Hz harmonics at fs = 100 Hz; see deviation 2
    'notch_bandwidth': 1.0,
    'ecg_artifact_removal': True,
    # port options
    'matlab_endsample_quirk': True,
    'remove_unscored': True,
}

SLEEP_STAGES = (1, 2, 3, 5)


# ----------------------------------------------------------------------------------------------- helpers
def _mround(x):
    """MATLAB round (half away from zero) -> int."""
    return int(np.sign(x) * np.floor(abs(x) + 0.5))


def _mcolon(a, step, b):
    """Values of MATLAB a:step:b."""
    if b < a:
        return np.zeros(0)
    n = int(np.floor((b - a) / step + 1e-10)) + 1
    return a + step * np.arange(n)


def _sl(a, i, j):
    """MATLAB a(i:j) with 1-based inclusive indices (clipped to the array)."""
    return a[max(i - 1, 0):max(j, 0)]


# --------------------------------------------------------------------------------------- Polyman filters
def polyman_bandpass(x, fs, gain, hp, lp):
    r, s = 2 * np.pi * lp / fs, 2 * np.pi * hp / fs
    t = 4 + 2 * (r + s) + r * s
    b = [gain * 2 * r / t, 0, -gain * 2 * r / t]
    a = [1, -(8 - 2 * r * s) / t, -(-4 + 2 * (r + s) - r * s) / t]
    return sps.lfilter(b, a, x)


def polyman_highpass(x, fs, gain, hp):
    r = 2 * np.pi * hp * ((1 / fs) / 2)
    t = gain / (r + 1)
    return sps.lfilter([t, -t], [1, (r - 1) / (r + 1)], x)


def polyman_notch(x, fs, center, bandwidth, gain=1.0):
    if not (0 < center <= fs / 2) or not (0 < bandwidth <= fs / 4):
        raise ValueError('notch outside the Nyquist interval')
    r = 2 * np.pi * center / fs
    s = 2 * np.pi * (center - bandwidth / 2) / fs
    t = 1 - np.sqrt((np.cos(r) - np.cos(s)) ** 2 + (np.sin(r) - np.sin(s)) ** 2) * np.sqrt(np.e ** 2 - 1)
    cx, cy = gain * np.cos(r), gain * np.sin(r)
    b = [gain, -2 * cx, cx ** 2 + cy ** 2]
    a = [1, -2 * t * cx, (t * cx) ** 2 + (t * cy) ** 2]
    return sps.lfilter(b, a, x)


# ------------------------------------------------------------------------------------ signal processing
def signal_amplitude(x, rate, window_size, window_step):
    """signalAmplitude.m: peak-to-peak per window; trailing windows run to the END of the signal."""
    x = np.asarray(x, float)
    step, size = rate * window_step, rate * window_size
    n = len(x)
    cw = max(int(np.floor((n - size) / step)), 0)
    nsz = int(np.floor(size))
    if cw:
        starts = np.floor(np.arange(cw) * step).astype(np.int64)
        if nsz == int(round(step)) and abs(step - nsz) < 1e-9:
            w = x[:cw * nsz].reshape(cw, nsz)
        else:
            w = np.lib.stride_tricks.sliding_window_view(x, nsz)[starts]
        amps = list(w.max(1) - w.min(1))
    else:
        amps = []
    for cs in _mcolon(cw * step + 1, step, n):
        seg = x[_mround(cs) - 1:]
        amps.append(seg.max() - seg.min())
    return np.asarray(amps, float)


def signal_power(x, window_size, rate, window_step, lcut, hcut, block=10000):
    """signalPower.m: Hamming-windowed periodogram band power per window (same bin rules and padding)."""
    x = np.asarray(x, float)
    step_samples = rate * window_step
    N = int(round(window_size * rate))
    n = len(x)
    ham = sps.windows.hamming(N, sym=True)
    n_unique = int(np.ceil((N + 1) / 2))
    limit_pts = n_unique - 1 if N % 2 == 0 else n_unique
    fstep = rate / N
    f1 = [int(np.ceil(lc / fstep)) + 1 for lc in lcut]
    f2 = [min(int(np.ceil(hc / fstep)) + 1, n_unique) for hc in hcut]
    num_windows = len(_mcolon(1, step_samples, n))
    windows_limit = len(_mcolon(1, step_samples, n - N))
    starts = np.floor(1 + np.arange(windows_limit) * step_samples).astype(np.int64) - 1
    view = np.lib.stride_tricks.sliding_window_view(x, N)
    out = [np.empty(num_windows) for _ in lcut]
    for b0 in range(0, windows_limit, block):
        w = view[starts[b0:b0 + block]] * ham
        mx = np.abs(np.fft.rfft(w, N, axis=1)) ** 2 / N        # (windows, n_unique), 0-based bins
        for k in range(len(lcut)):
            a, b = f1[k], f2[k]                                # 1-based bins as in MATLAB
            if a > 1 and b <= limit_pts:
                t = mx[:, a - 1:b].sum(1) * 2
            elif a > 1 and b > limit_pts:
                t = mx[:, a - 1:limit_pts].sum(1) * 2 + mx[:, limit_pts:b].sum(1)
            elif b < limit_pts:
                t = mx[:, 0] + mx[:, 1:b].sum(1) * 2
            else:
                t = mx[:, 0] + mx[:, 1:limit_pts].sum(1) * 2 + mx[:, limit_pts:b].sum(1)
            out[k][b0:b0 + len(t)] = t / N
    for k in range(len(lcut)):
        out[k][windows_limit:] = out[k][windows_limit - 1]
    return out


def baseline(background_time, window_step, p):
    """baseline.m (circular buffer) in closed form: mean of the previous `size` values, initialised with the
    mean of the first `length(1:step:T)` values."""
    p = np.asarray(p, float)
    size = _mround(background_time / window_step)
    nb = len(_mcolon(1, window_step, background_time))
    m = p[:nb].mean()
    k = np.arange(1, len(p) + 1)                                # 1-based window index
    S = np.concatenate([[0.0], np.cumsum(p)])
    lo = np.maximum(k - size - 1, 0)
    return (S[k - 1] - S[lo] + np.maximum(0, size + 1 - k) * m) / size


def signal_cycles(x):
    """signalCycles.m: diff(sign(diff(x))) with 0/0 -> NaN (MATLAB behaviour)."""
    d = np.diff(x)
    with np.errstate(invalid='ignore', divide='ignore'):
        s = np.abs(d) / d
    return np.diff(s)


def merge_events_closer_than(ev, union_time, window_step):
    ev = ev.copy()
    if union_time > 0:
        on = np.flatnonzero(ev)
        for a, b in zip(on[:-1], on[1:]):
            if b - a < union_time / window_step:
                ev[a:b + 1] = 1
    return ev


def get_signal_events(p, threshold, union_time, window_step, base):
    ev = (p > threshold * base).astype(np.int8)
    ev[0] = ev[-1] = 0
    return merge_events_closer_than(ev, union_time, window_step)


def get_event_indexes(ev):
    """getEventIndexes.m -> 1-based (first 1, last 1) of each run."""
    d = np.diff(ev.astype(np.int8))
    init = np.flatnonzero(d == 1) + 2
    ends = np.flatnonzero(d == -1) + 1
    if len(init) and len(ends):
        if ends[0] < init[0]:
            ends = ends[1:]
        if len(ends) and ends[-1] < init[-1]:
            init = init[:-1]
    return init, ends


def update_alpha_beta_overlapping(ea, eb, union):
    """updateAlphaBetaOverlapping.m (1-based arithmetic kept); returns (events, n_updates)."""
    ea = ea.copy()
    union = int(round(union))
    ini, end = get_event_indexes(ea)
    if not len(ini):
        return ea, 0
    L = len(eb)
    B = lambda i: eb[i - 1]                                    # noqa: E731  MATLAB eb(i)
    updates = 0
    for k in range(len(ini)):
        i0, e0 = int(ini[k]), int(end[k])
        old = e0 - i0
        si = _sl(eb, max(i0 - union, 1), max(e0 - 1, 1))
        if si.any() and np.flatnonzero(si)[0] == 0 and i0 - union - 1 >= 1:
            k1 = 0
            while i0 - union - k1 - 1 >= 1 and B(i0 - union - k1 - 1):
                k1 += 1
                if i0 - union - k1 == 1:
                    break
            i0 = i0 - union - k1
        si = _sl(eb, min(e0 + 1, L), min(e0 + union, L))
        if si.any() and np.flatnonzero(si)[-1] == len(si) - 1 and e0 + union + 1 <= L:
            k1 = 0
            while e0 + union + k1 + 1 <= L and B(e0 + union + k1 + 1):
                k1 += 1
                if e0 + union + k1 == L:
                    break
            e0 = e0 + union + k1
        if e0 - i0 != old:
            ea[i0 - 1:e0] = 1
            updates += 1
    return ea, updates


# ---------------------------------------------------------------------------------- ECG artifact removal
def _mean_filter_matlab(x, n):
    y = np.convolve(x, np.ones(n) / n)
    y[:n] = y[n - 1]
    y[-n:] = y[-n]
    h = n // 2
    return y[h - 1:len(y) - h]


def ecg_reliability(ecg, sr):
    """ecgReliabilityCheck.m (without the 50 Hz notch): per-sample 1 = reliable ECG."""
    ecg = np.asarray(ecg, float)
    L = len(ecg)
    rel = np.ones(L)
    f = polyman_bandpass(ecg, sr, 1, 8, 16)
    f = np.abs(np.diff(f))
    f = _mean_filter_matlab(f, _mround(0.08 * sr))
    fn = _mean_filter_matlab(f, _mround(0.5 * sr))
    if len(fn) < L:
        fn = np.concatenate([fn, np.full(L - len(fn), fn[-1])])
    if len(f) < L:
        f = np.concatenate([f, np.full(L - len(f), f[-1])])
    f, fn = f[:L], fn[:L]
    with np.errstate(invalid='ignore', divide='ignore'):
        f = f / fn
    # peak regions: above the mean of non-overlapping 10 s windows
    ws = int(10 * sr)
    nw = (len(f) - ws) // ws
    peak = np.zeros(L, bool)
    if nw > 0:
        W = f[:nw * ws].reshape(nw, ws)
        peak[:nw * ws] = (W > W.mean(1, keepdims=True)).ravel()
    d = np.diff(peak.astype(np.int8))
    ini = list(np.flatnonzero(d == 1) + 2)
    ends = list(np.flatnonzero(d == -1) + 1)
    if peak[0]:
        ini = [1] + ini
    if peak[-1]:
        ends = ends + [L]
    if len(ini) > len(ends):
        ini = ini[:-1]
    if len(ends) > len(ini):
        ends = ends[1:]
    ini, ends = np.asarray(ini, int), np.asarray(ends, int)
    peaks = np.array([f[a - 1:b].max() for a, b in zip(ini, ends)]) if len(ini) else np.zeros(0)
    idx_peaks = np.array([a + np.argmax(f[a - 1:b]) + 1 for a, b in zip(ini, ends)]) if len(ini) else np.zeros(0)
    # 30 s windows: CV of the 20 largest peaks (40 bpm) > 0.05 -> unreliable
    ws = int(30 * sr)
    nw = (L - ws) // ws
    n_pick = _mround(40 * 30 / 60)
    for k in range(1, nw + 1):
        off = ((k - 1) * 30 + 1) * sr                          # sic: window shifted by 1 s (MATLAB)
        m = peaks[(idx_peaks >= off) & (idx_peaks <= off + ws - 1)]
        top = np.sort(m)[::-1][:min(len(m), n_pick)]
        ratio = np.std(top, ddof=1) / np.mean(top) if len(top) > 1 else (0.0 if len(top) else np.nan)
        rel[(k - 1) * ws:(k - 1) * ws + ws] = 0.0 if ratio > 0.05 else 1.0
    return rel


def ecg_conditioning(x, fs, rel, rpeaks, length_s=0.5, up=0.01):
    """ecgConditioning.m: subtract a running (alpha = 0.01) Hann-tapered template around each R peak."""
    y = np.asarray(x, float).copy()
    half = _mround(int(np.floor(fs * length_s)) / 2)
    tmpl = np.zeros(2 * half)
    taper = sps.windows.hann(2 * half, sym=True)
    n = len(y)
    for r in np.asarray(rpeaks, np.int64):
        a, b = r - half, r + half
        if a < 0 or b > n:
            continue                                           # deviation 5 (MATLAB errors at the edges)
        win = y[a:b]
        if rel[a:b].sum() < (b - a):
            tmpl = (1 - up) * tmpl
        else:
            tmpl = (1 - up) * tmpl + up * win
            y[a:b] = win - taper * tmpl
    return y


# --------------------------------------------------------------------------------------- core algorithm
class _Ctx:
    pass


def _new_event(start_idx, end_idx, er, rate, alpha):
    start = (start_idx - 1) / er
    dur = (end_idx - start_idx) / er
    ss = _mround(start * rate)
    es = ss + _mround(dur * rate) if alpha else _mround((start + dur) * rate)
    return {'startIdx': int(start_idx), 'endIdx': int(end_idx), 'start': start, 'duration': dur,
            'startSample': ss, 'endSample': es, 'reason': ''}


def _update_with_power(ev, c, power, kind):
    p = c.p
    pr = c.er
    thr = p[f'{kind}_power_threshold']
    limit_t = p[f'{kind}_power_limit_time'] * pr
    inc = int(p[f'{kind}_power_increase_time'] * pr)
    prev_n = int(p[f'{kind}_power_previous_seconds'] * pr)
    end_idx = ev['startIdx'] + _mround(ev['duration'] * pr)
    prev_power = _sl(power, ev['startIdx'] - prev_n, ev['startIdx'] - 1).mean()
    x = c.eeg
    prev_amp = np.abs(_sl(x, ev['startSample'] - int(10 * c.rate), ev['startSample'] - 1)).max()
    ev_amp = np.abs(_sl(x, ev['startSample'], ev['endSample'])).max()
    if ev_amp > 8 * prev_amp:
        return ev
    last = new_end = end_idx
    while new_end < end_idx + 30 * pr:
        seg = _sl(power, last + 1, last + inc)
        if not len(seg):
            break
        if seg.mean() > thr * prev_power:
            new_end = last + inc
        last += inc
        if last >= new_end + limit_t:
            break
    if new_end > end_idx:
        add = (new_end - end_idx) / pr
        if add > p['max_update_s']:
            return ev
        ev['duration'] += add
        ev['endSample'] = ev['startSample'] + _mround((add if p['matlab_endsample_quirk'] else ev['duration'])
                                                      * c.rate)
    return ev


def _update_with_overlapping_emg(ev, c):
    p = c.p
    amp_rate = 1 / p['emg_update_window_step']
    d = int(15 * c.rate)
    aux = np.sort(signal_amplitude(_sl(c.emg, ev['startSample'] - d, ev['endSample'] + d), c.rate,
                                   p['emg_update_window_size'], p['emg_update_window_step']))
    thr = p['emg_update_threshold'] * aux[max(_mround(0.5 * len(aux)), 1) - 1]
    sA = _mround(ev['start'] * amp_rate)
    eA = _mround((ev['start'] + ev['duration']) * amp_rate)
    au = c.emg_amp_update
    if (_sl(au, eA - 10, eA) > thr).sum() >= 8:
        below = np.flatnonzero(au[eA - 1:] <= thr)
        k = eA + below[0] if len(below) else len(au)
    else:
        k = eA
    dur = (k - sA) / amp_rate
    if dur > ev['duration'] + 0.1:
        ev['duration'] = dur
        ev['endSample'] = ev['startSample'] + _mround(dur * c.rate)
    return ev


def _update_with_noncontinuous_emg(ev, c):
    rate = c.rate
    cut = 0.52 if rate in (125, 250) else 0.5
    cs = _mround(cut * rate)
    if ev['duration'] < cut:
        return ev
    limit = int(np.floor(60 / cut))
    ev_limit = int(np.floor(ev['duration'] / cut))
    if ev_limit >= limit:
        return ev                                              # MATLAB: empty loop -> no update
    s0 = ev['startSample'] - int(3 * rate)                     # 1-based first sample of previousWindowAxe
    wl = int(3 * rate)
    x = c.emg
    prev = x[s0 - 1:s0 - 1 + wl]
    prev_diff = np.mean(np.abs(x[s0 - 1 - cs:s0 - 1 - cs + wl] - prev))
    idx = np.arange(1, limit + 1)
    st = s0 - 1 + idx * cs
    ok = st + wl <= len(x)
    view = np.lib.stride_tricks.sliding_window_view(x, wl)
    diffs = np.full(limit, np.nan)
    diffs[ok] = np.abs(view[st[ok]] - prev).mean(1)
    ev_diff = diffs[:ev_limit]
    if np.all(ev_diff < 2 * prev_diff):
        return ev
    m = ev_diff.mean()
    sd = ev_diff.std(ddof=1) if ev_limit > 1 else 0.0
    post = np.full(limit, np.nan)
    i = ev_limit
    for i in range(ev_limit + 1, limit + 1):
        post[i - 1] = diffs[i - 1]
        if post[i - 1] < m - sd:
            break
        if np.sum(post < 1.1 * (m - sd)) > 2:
            break
    dur = ev['duration'] + (i - ev_limit - 6) * cut
    if dur > ev['duration']:
        if dur > c.p['max_update_s']:
            return ev
        ev['endSample'] = ev['startSample'] + _mround(dur * rate)
        ev['duration'] = dur
    return ev


def _update_with_nonoverlapping_emg(ev, c):
    p, rate, x = c.p, c.rate, c.emg
    lim = p['future_emg_limit']
    if rate in (125, 250) or abs(0.125 * rate - round(0.125 * rate)) > 1e-9:
        wsize, aws = 0.48, 0.12                                # deviation 4
    else:
        wsize, aws = 0.5, 0.125
    an = _mround(aws * rate)
    aww = int(np.floor(wsize / aws + 1e-9))
    prev = np.abs(_sl(x, ev['startSample'] - int(10 * rate), ev['startSample'] - int(1 * rate) - 1))
    prev = prev[:len(prev) // an * an]
    prev_amp = max(np.median(prev.reshape(-1, an).max(1)), p['emg_amp_floor_uv'])
    if ev['endSample'] + lim * rate - 1 > len(x):
        return ev
    fut = np.abs(_sl(x, ev['endSample'], ev['endSample'] + int(lim * rate) - 1))
    fut = fut[:len(fut) // an * an]
    fut_amp = fut.reshape(-1, an).max(1)
    nw = len(fut_amp) // aww
    big = new_end = 0
    i = 0
    for i in range(1, nw + 1):
        amp = np.median(fut_amp[(i - 1) * aww:i * aww])
        if amp > prev_amp * p['future_emg_threshold']:
            big += 1
        else:
            if big >= 2:
                new_end = i - 1
            if i - 1 - new_end >= 4 / wsize:
                break
            big = 0
    if i == nw and big > 0:
        new_end = nw
    if new_end > 0:
        dur = ev['duration'] + new_end * wsize
        ev['endSample'] = ev['startSample'] + _mround(dur * rate)
        ev['duration'] = dur
    return ev


def _update_with_amplitude(ev, c):
    thr = c.p['amplitude_update_threshold']
    x, rate = c.eeg, int(c.rate)
    ss = ev['startSample']
    es = min(len(x), ev['endSample'])
    prev = signal_amplitude(_sl(x, ss - int(c.er * rate), ss), rate, 1, 1)   # sic: eventsRate * rate = 5 s
    prev_amp = prev[prev > 0].mean() if np.any(prev > 0) else np.inf
    seg = _sl(x, ss, es)
    if seg.max() - seg.min() > thr * prev_amp:
        for sample in range(es + rate, len(x) - rate + 1, rate):
            a = signal_amplitude(_sl(x, sample, sample + rate), rate, 1, 1)
            if (a[a > 0].mean() if np.any(a > 0) else 0) > thr * prev_amp:
                es = sample
            else:
                break
    dur = (es - ss) / rate
    if dur > ev['duration'] + 1:
        ev['duration'] = dur                                   # sic: endSample not updated
    return ev


def _emg_activity(ev, c, min_s=None):
    seg = _sl(c.emg_events, _mround(ev['start'] * c.emg_er), _mround((ev['start'] + ev['duration']) * c.emg_er))
    if min_s is None:
        return bool(seg.any())
    return seg.sum() > min_s * c.emg_er


def _too_long(ev, c):
    lim = c.p['max_duration_per_epoch']
    if ev['duration'] <= lim:
        return False
    off = c.epoch_offset
    end = ev['start'] + ev['duration']
    e1 = np.ceil((ev['start'] - off) / 30)
    e2 = np.ceil((end - off) / 30)
    if e1 != e2:
        first = e1 * 30 + off - ev['start']
        second = end - (e1 * 30 + off)
        return first > lim or second > lim
    return True


def _find_spindles(ev, c):
    rate = c.rate
    first = _mround((ev['start'] - 0.5) * rate)
    last = _mround((ev['duration'] + ev['start']) * rate)
    nxt = _mround(rate * 0.125)
    delta = _mround(rate * 0.125)
    ws = _mround(rate * 0.5)
    C = c.cyc_cum                                              # C[k] = sum(cycles(1:k) > 0)
    n = len(c.eeg)
    out = []
    i = first
    while i < last - ws:
        j = i + ws
        num = C[j - 2] - C[i - 1]
        freq = num / 0.5
        while 12 <= freq <= 14 and j + delta < n:
            j += delta
            num = C[j - 2] - C[i - 1]
            freq = num / ((j - i) / rate)
        if (j - delta - i) > 0.5 * rate:
            out.append((i / rate, (j - delta - i) / rate))
            i = j + nxt
        else:
            i += delta
    return out


def _spindle_free(ev, c):
    sp = _find_spindles(ev, c)
    if not sp:
        return True
    starts = [ev['start']] + [s + d for s, d in sp]
    durs = [sp[0][0] - ev['start']] + [sp[k + 1][0] - starts[k + 1] for k in range(len(sp) - 1)]
    durs.append(ev['duration'] + ev['start'] - starts[-1])
    return any(d >= 3 for d in durs)


def _candidates(ev_vec, c, alpha):
    ini, end = get_event_indexes(ev_vec)
    evs = [_new_event(a, b, c.er, c.rate, alpha) for a, b in zip(ini, end)]
    edge = c.p['edge_exclusion_s']
    return [e for e in evs if e['start'] >= edge and e['startSample'] <= len(c.eeg) - edge * c.rate
            and e['startIdx'] != e['endIdx']]


def _update_chain(ev, c, power, kind):
    ev = _update_with_power(ev, c, power, kind)
    ev = _update_with_overlapping_emg(ev, c)
    ev = _update_with_noncontinuous_emg(ev, c)
    ev = _update_with_nonoverlapping_emg(ev, c)
    return _update_with_amplitude(ev, c)


def _find_beta_arousals(c):
    out = []
    for ev in _candidates(c.events_beta, c, alpha=False):
        ev = _update_chain(ev, c, c.power_beta, 'beta')
        emg_act = _emg_activity(ev, c)
        ev['endIdx'] = ev['startIdx'] + int(np.floor(ev['duration'] * c.er))
        valid = True
        if ev['duration'] < c.p['min_beta_duration']:
            valid, ev['reason'] = False, 'too short'
        if _too_long(ev, c):
            valid, ev['reason'] = False, 'too long'
        if valid and not emg_act and not _spindle_free(ev, c):
            valid, ev['reason'] = False, 'spindles'
        ev['kind'], ev['ok'] = 'beta', valid
        out.append(ev)
    return out


def _find_alpha_arousals(c):
    out = []
    for ev in _candidates(c.events_alpha, c, alpha=True):
        emg_act = False
        if ev['duration'] > 1:
            ev = _update_chain(ev, c, c.power_alpha, 'alpha')
            emg_act = _emg_activity(ev, c, c.p['alpha_min_emg_s'])
        valid = True
        if ev['duration'] < c.p['min_alpha_duration']:
            valid, ev['reason'] = False, 'too short'
        if _too_long(ev, c):
            valid, ev['reason'] = False, 'too long'
        if valid and not emg_act:
            valid, ev['reason'] = False, 'no EMG'
        ev['kind'], ev['ok'] = 'alpha', valid
        out.append(ev)
    return out


def _stage_at(c, t):
    return int(c.stage[min(max(int(np.floor(t * c.rate)), 0), len(c.stage) - 1)])


def _epoch_offset(stage, fs):
    ch = np.flatnonzero(np.diff(stage) != 0) + 1
    if not len(ch):
        return 0.0
    v, n = np.unique(ch % int(30 * fs), return_counts=True)
    return float(v[np.argmax(n)] / fs)


def preprocess(eeg, emg, ecg, fs, rpeaks=None, params=None):
    """Notch (alias lines), EMG high-pass and ECG-artifact removal (loadFile.m filter chain)."""
    p = {**DEFAULT_PARAMS, **(params or {})}
    eeg = np.nan_to_num(np.asarray(eeg, float))
    emg = np.nan_to_num(np.asarray(emg, float))
    for f in p['notch_freqs'] or ():
        if f < fs / 2:
            eeg = polyman_notch(eeg, fs, f, p['notch_bandwidth'])
            emg = polyman_notch(emg, fs, f, p['notch_bandwidth'])
    emg = polyman_highpass(emg, fs, 1, p['emg_highpass'])
    if ecg is not None and p['ecg_artifact_removal']:
        ecg = np.nan_to_num(np.asarray(ecg, float))
        if rpeaks is None:
            from snoredreamt.ecg import detect_rpeaks
            rpeaks, _ = detect_rpeaks(ecg, fs)
        rel = ecg_reliability(ecg, fs)
        eeg = ecg_conditioning(eeg, fs, rel, rpeaks)
        emg = ecg_conditioning(emg, fs, rel, rpeaks)
    return eeg, emg


def detect_arousals(eeg, emg, ecg, stage, fs, valid=None, rpeaks=None, params=None, return_all=False):
    """Detect EEG arousals.

    eeg, emg, ecg : 1-D arrays in µV at `fs` (ecg may be None -> no ECG artifact removal).
    stage : per-sample hypnogram (0 W, 1-3 N1-N3, 5 REM, -1 unscored).
    valid : optional per-sample bool artifact mask (True = clean); only used for the `valid_frac` column.
    rpeaks : optional R-peak sample indices (0-based); detected with snoredreamt.ecg.detect_rpeaks if None.
    return_all : also return rejected candidates (column `reason` says why; `ok` False).

    Returns a DataFrame with start/end (0-based samples), onset_s, duration_s, kind ('alpha'/'beta'),
    stage (at onset), stage_end, valid_frac.
    """
    p = {**DEFAULT_PARAMS, **(params or {})}
    if p['beta_band'][1] >= fs / 2 or p['alpha_band'][1] >= fs / 2:
        raise ValueError('band edges must be below Nyquist')
    stage = np.asarray(stage)
    x, y = preprocess(eeg, emg, ecg, fs, rpeaks, p)
    c = _Ctx()
    c.p, c.rate, c.eeg, c.emg, c.stage = p, fs, x, y, stage
    c.epoch_offset = _epoch_offset(stage, fs)
    # findPowerEvents.m
    step = p['power_window_step']
    c.er = 1 / step
    pa, pb = signal_power(x, p['power_window_size'], fs, step, [p['alpha_band'][0], p['beta_band'][0]],
                          [p['alpha_band'][1], p['beta_band'][1]])
    nsh = _mround(0.5 * p['power_window_size'] / step)
    for arr in (pa, pb):
        arr[nsh - 1:] = arr[:len(arr) - nsh + 1].copy()
    ba = baseline(p['power_background_time'], step, pa)
    bb = baseline(p['power_background_time'], step, pb)
    ea = get_signal_events(pa, p['power_alpha_threshold'], 0, step, ba)
    eb = get_signal_events(pb, p['power_beta_threshold'], 0, step, bb)
    if p['update_with_alpha_events']:
        u = 1
        while u:
            eb, u = update_alpha_beta_overlapping(eb, ea, p['update_with_alpha_events_union_time'] / step)
    if p['update_with_beta_events']:
        u = 1
        while u:
            ea, u = update_alpha_beta_overlapping(ea, eb, p['update_with_beta_events_union_time'] / step)
    ea = merge_events_closer_than(ea, p['power_beta_union'], step)    # sic (crossed names)
    eb = merge_events_closer_than(eb, p['power_alpha_union'], step)
    c.power_alpha, c.power_beta, c.events_alpha, c.events_beta = pa, pb, ea, eb
    # findAmplitudeEvents.m (EMG)
    astep = p['amp_window_step']
    amp = signal_amplitude(y, fs, p['amp_window_size'], astep)
    bamp = baseline(p['amp_background_time'], astep, amp)
    nsh = _mround(0.5 * p['amp_window_size'] / astep)
    if nsh > 1:
        amp[nsh - 1:] = amp[:len(amp) - nsh + 1].copy()
        bamp[nsh - 1:] = bamp[:len(bamp) - nsh + 1].copy()
    eemg = get_signal_events(amp, p['amp_threshold'], p['amp_merge_near'], astep, bamp)
    iso = np.zeros(len(eemg), bool)
    iso[1:-2] = (eemg[:-3] == 0) & (eemg[1:-2] == 1) & (eemg[2:-1] == 0)
    eemg[iso] = 0
    c.emg_events, c.emg_er = eemg, 1 / astep
    c.emg_amp_update = signal_amplitude(y, fs, p['emg_update_window_size'], p['emg_update_window_step'])
    sm = polyman_bandpass(x, fs, 1, *p['spindle_bandpass'])
    c.cyc_cum = np.concatenate([[0], np.cumsum(signal_cycles(sm) > 0)])
    # findArousalEvents.m
    beta = _find_beta_arousals(c) if eb.any() else []
    for ev in beta:                                            # removeBetaArousals
        if ev['ok']:
            a = _mround(ev['start'] * c.er) - int(10 * c.er)
            b = _mround((ev['start'] + ev['duration']) * c.er)
            c.events_alpha[max(a - 1, 0):b] = 0
    alpha = _find_alpha_arousals(c) if c.events_alpha.any() else []
    # mergeEvents(alpha, beta): stable by start, ties -> beta first
    cand = sorted([e for e in alpha + beta if e['ok']], key=lambda e: (e['start'], e['kind'] == 'alpha'))
    rejected = [e for e in alpha + beta if not e['ok']]
    kept = []
    for ev in cand:                                            # removeAwakeEvents + REM rule
        s0, s1 = _stage_at(c, ev['start']), _stage_at(c, ev['start'] + ev['duration'])
        if s0 == 0 or s1 == 0:
            ev['reason'] = 'awake'
        elif p['remove_unscored'] and (s0 < 0 or s1 < 0):
            ev['reason'] = 'unscored'
        elif s0 == 5 and s1 == 5 and not _emg_activity(ev, c):
            ev['reason'] = 'rem'
        if ev['reason']:
            ev['ok'] = False
            rejected.append(ev)
        else:
            kept.append(ev)
    keep = np.ones(len(kept), bool)                            # removeIfClose
    for i in range(1, len(kept)):
        prev_end = kept[i - 1]['start'] + kept[i - 1]['duration']
        if keep[i - 1] and kept[i]['start'] <= prev_end + p['min_stable_sleep']:
            keep[i] = False
    final = [e for e, k in zip(kept, keep) if k]
    for e, k in zip(kept, keep):
        if not k:
            e['ok'], e['reason'] = False, 'within 10 s of previous'
            rejected.append(e)
    rem = np.zeros(len(final), bool)                           # removeIncluded (on the removeIfClose output)
    for i in range(len(final) - 1):
        if final[i]['start'] + final[i]['duration'] > final[i + 1]['start']:
            rem[i + 1] = True
    for e, r in zip(final, rem):
        if r:
            e['ok'], e['reason'] = False, 'included'
            rejected.append(e)
    final = [e for e, r in zip(final, rem) if not r]
    return _to_frame(final + (rejected if return_all else []), c, valid)


def _to_frame(evs, c, valid):
    cols = ['start', 'end', 'onset_s', 'duration_s', 'kind', 'stage', 'stage_end', 'valid_frac', 'ok', 'reason']
    if not evs:
        return pd.DataFrame(columns=cols)
    n = len(c.eeg)
    rows = []
    for e in evs:
        s = min(_mround(e['start'] * c.rate), n - 1)
        en = min(_mround((e['start'] + e['duration']) * c.rate), n)
        vf = float(np.mean(valid[s:max(en, s + 1)])) if valid is not None else np.nan
        rows.append({'start': s, 'end': en, 'onset_s': e['start'], 'duration_s': e['duration'], 'kind': e['kind'],
                     'stage': _stage_at(c, e['start']), 'stage_end': _stage_at(c, e['start'] + e['duration']),
                     'valid_frac': vf, 'ok': e['ok'], 'reason': e['reason']})
    return pd.DataFrame(rows, columns=cols).sort_values('start', kind='stable').reset_index(drop=True)
