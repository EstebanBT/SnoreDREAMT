"""Respiratory cycle detection and phase.

Conventions (after polarity calibration, so that inspiration is positive-going):
- Belts (ABDOMEN/THORAX): peak = end of inspiration (maximal expansion), trough = end of expiration.
- Airflow (FLOW thermistor, PTAF): peak = inhalation peak (maximal inspiratory flow), as BreathMetrics 'inhale peak'.
Phase 0 = peak, +-pi = trough. Hilbert phase is offset-corrected so that its circular mean at detected peaks is 0.

Two phase estimators:
- `hilbert_phase`: Butterworth band-pass at the subject's breathing peak frequency +- half-width, two-pass, angle of
  the analytic signal (user's BreathTMR pipeline / Schreiner 2023).
- `cycle_phase`: piecewise-linear phase between detected peaks and troughs (port of slythm.doubleInterpolation, with
  the segment-length bug fixed): 0 -> pi from peak to trough, -pi -> 0 from trough to next peak. NaN where breath
  detection fails or the breath amplitude is below `min_amp_frac` of the subject's median clean amplitude.
"""
import numpy as np
from scipy import signal


def lowpass(x, fs, f=1.0, order=4):
    sos = signal.butter(order, f, fs=fs, output='sos')
    return signal.sosfiltfilt(sos, np.nan_to_num(np.asarray(x, float)))


def highpass(x, fs, f=0.05, order=2):
    sos = signal.butter(order, f, btype='high', fs=fs, output='sos')
    return signal.sosfiltfilt(sos, x)


def peak_frequency(x, fs, mask, fmin=0.1, fmax=0.6, seg_s=60):
    """Breathing peak frequency from the mean spectrum of 60 s segments fully inside `mask`."""
    w = int(seg_s * fs)
    starts = [i for i in range(0, len(x) - w, w) if mask[i:i + w].all()]
    if len(starts) < 3:
        return np.nan
    segs = np.stack([x[i:i + w] for i in starts])
    f, p = signal.welch(segs, fs, nperseg=w, axis=1)
    p = p.mean(0)
    sel = (f >= fmin) & (f <= fmax)
    return f[sel][np.argmax(p[sel])]


def detect_breaths(x, fs, ref_mask, min_breath_s=1.5, prom_frac=0.3):
    """Peaks and troughs of a (polarity-corrected, low-passed) respiratory signal.

    Prominence threshold = prom_frac x median peak-to-trough amplitude estimated in `ref_mask` (clean sleep).
    Returns arrays peaks, troughs (alternating, trough i lies between peak i and peak i+1) and per-breath amplitude.
    """
    # robust breath amplitude: median of 10 s ranges in the reference mask
    w = int(10 * fs)
    idx = [i for i in range(0, len(x) - w, w) if ref_mask[i:i + w].all()]
    amp_ref = np.median([np.ptp(x[i:i + w]) for i in idx]) if idx else np.ptp(x) / 10
    dist = int(min_breath_s * fs)
    pk, _ = signal.find_peaks(x, prominence=prom_frac * amp_ref, distance=dist)
    tr = np.array([a + np.argmin(x[a:b]) for a, b in zip(pk[:-1], pk[1:])], dtype=np.int64)
    pk = pk[:-1] if len(pk) else pk
    amp = x[pk] - x[tr] if len(pk) else np.array([])
    return pk.astype(np.int64), tr, amp, amp_ref


def cycle_phase(n, peaks, troughs, amp=None, amp_ref=None, min_amp_frac=0.25, max_cycle_s=15, fs=100):
    """Piecewise-linear phase (radians) at every sample; NaN where undefined.

    Segment peak_i -> trough_i: 0 -> pi; trough_i -> peak_{i+1}: -pi -> 0. Breaths with amplitude below
    min_amp_frac x amp_ref, or cycles longer than max_cycle_s, are NaN (apnea / detection failure).
    """
    ph = np.full(n, np.nan)
    for i in range(len(peaks) - 1):
        p0, t0, p1 = peaks[i], troughs[i], peaks[i + 1]
        if p1 - p0 > max_cycle_s * fs:
            continue
        if amp is not None and amp_ref is not None and amp[i] < min_amp_frac * amp_ref:
            continue
        ph[p0:t0] = np.linspace(0, np.pi, t0 - p0, endpoint=False)
        ph[t0:p1] = np.linspace(-np.pi, 0, p1 - t0, endpoint=False)
    return ph


def hilbert_phase(x, fs, f0, half_width=0.1, order=4):
    """Narrow-band Hilbert phase around f0 (two-pass Butterworth, FieldTrip default order 4)."""
    sos = signal.butter(order, [max(f0 - half_width, 0.02), f0 + half_width], btype='band', fs=fs, output='sos')
    return np.angle(signal.hilbert(signal.sosfiltfilt(sos, x - np.mean(x))))


def align_phase(ph, peaks):
    """Rotate phase so that its circular mean at `peaks` is 0. Returns (aligned phase, offset removed)."""
    v = ph[peaks]
    v = v[np.isfinite(v)]
    off = np.angle(np.mean(np.exp(1j * v))) if len(v) else 0.0
    return np.angle(np.exp(1j * (ph - off))), off


# ------------------------------------------------------------------ polarity calibration
def _breath_vectors(y, ph, starts, ends):
    """Per-breath phase vector of the anomaly of y: sum((y - breath mean) * exp(i*phase)) / n."""
    v = []
    for a, b in zip(starts, ends):
        yy, pp = y[a:b], ph[a:b]
        ok = np.isfinite(yy) & np.isfinite(pp)
        if ok.sum() > 10:
            v.append(np.mean((yy[ok] - yy[ok].mean()) * np.exp(1j * pp[ok])))
    return np.array(v)


def _vote(v, project):
    """Project per-breath vectors on the expected direction; return (vote, z, preferred angle in degrees)."""
    if len(v) < 50:
        return 0, np.nan, np.nan
    s = project(v)
    z = np.mean(s) / (np.std(s) / np.sqrt(len(s)))
    return (int(np.sign(z)) if abs(z) > 3 else 0), float(z), float(np.degrees(np.angle(np.mean(v))))


def polarity_anchors(x, fs, ref_mask, hr=None, snore_env=None):
    """Evidence that `x` (low-passed belt signal) has inspiration positive-going: +1 (as is), -1 (flip), 0 (unclear).

    Uses the cycle phase of x as given (0 = peak, -pi..0 = rising limb):
    - snore: snore energy is centred in inspiration -> expected near -90 deg (rising limb).
    - rsa: heart rate peaks around end-inspiration (respiratory sinus arrhythmia) -> expected near 0 deg.
    - timing: inspiration (rising limb) shorter than expiration; weakest anchor, fails with flow limitation (snoring).
    """
    pk, tr, amp, amp_ref = detect_breaths(x, fs, ref_mask)
    ok = ref_mask[pk[:-1]] & ref_mask[tr[:-1]] & ref_mask[pk[1:]]
    out = {'n_breaths': int(ok.sum())}
    if ok.sum() < 50:
        return dict(out, timing=0, rsa=0, snore=0)
    p0, t0, p1 = pk[:-1][ok], tr[:-1][ok], pk[1:][ok]
    ratio = np.median((p1 - t0) / (t0 - p0))
    out.update(timing=int(np.sign(1 - ratio)) if abs(np.log(ratio)) > 0.1 else 0, ie_ratio=float(ratio))
    ph = cycle_phase(len(x), pk, tr, amp, amp_ref, fs=fs)
    starts, ends = p0, p1
    if hr is not None:
        out['rsa'], out['rsa_z'], out['rsa_deg'] = _vote(_breath_vectors(hr, ph, starts, ends), np.real)
    else:
        out['rsa'] = 0
    if snore_env is not None:
        out['snore'], out['snore_z'], out['snore_deg'] = _vote(_breath_vectors(snore_env, ph, starts, ends),
                                                              lambda v: -np.imag(v))
    else:
        out['snore'] = 0
    return out


def decide_polarity(a, strong_z=5.0):
    """Combine anchors: a strong snore anchor wins, then a strong RSA anchor, then majority. Flags conflicts."""
    votes = {k: a.get(k, 0) for k in ['snore', 'rsa', 'timing']}
    zs = {'snore': abs(a.get('snore_z', 0) or 0), 'rsa': abs(a.get('rsa_z', 0) or 0)}
    if votes['snore'] and zs['snore'] >= strong_z:
        sign = votes['snore']
    elif votes['rsa'] and zs['rsa'] >= strong_z:
        sign = votes['rsa']
    else:
        sign = int(np.sign(sum(votes.values())))
    nonzero = [v for v in votes.values() if v]
    conflict = bool(nonzero) and (min(nonzero) != max(nonzero))
    return sign, conflict


def joint_belt_polarity(a_abd, a_tho, r_abd_tho, strong_z=5.0, min_r=0.3):
    """Pool anchor evidence across belts. In normal breathing thorax and abdomen move in phase, so if their
    clean-breathing correlation r_abd_tho is clearly positive they share recording polarity (negative: opposite).
    Returns (abdomen sign, thorax sign, conflict, pooled z dict)."""
    rel = int(np.sign(r_abd_tho)) if abs(r_abd_tho) >= min_r else 0

    def pooled(key):
        za = a_abd.get(key + '_z', np.nan)
        zt = a_tho.get(key + '_z', np.nan)
        za = 0.0 if not np.isfinite(za) else za
        zt = 0.0 if not np.isfinite(zt) or rel == 0 else rel * zt
        return (za + zt) / (np.sqrt(2) if rel else 1)

    z = {'snore': pooled('snore'), 'rsa': pooled('rsa')}
    timing = a_abd.get('timing', 0) + (rel * a_tho.get('timing', 0) if rel else 0)
    if abs(z['snore']) >= strong_z:
        s = int(np.sign(z['snore']))
    elif abs(z['rsa']) >= strong_z:
        s = int(np.sign(z['rsa']))
    else:
        s = int(np.sign(np.sign(z['snore']) * (abs(z['snore']) > 3) + np.sign(z['rsa']) * (abs(z['rsa']) > 3)
                        + np.sign(timing)))
    votes = [int(np.sign(v)) for v in (z['snore'], z['rsa']) if abs(v) > 3]
    conflict = len(set(votes)) > 1
    tho = s * rel if rel else resp_sign_single(a_tho, strong_z)
    return s, tho, conflict, z


def resp_sign_single(a, strong_z=5.0):
    return decide_polarity(a, strong_z)[0]


def airflow_polarity(flow, belt, fs, mask, max_lag_s=2.0):
    """Airflow ~ d(volume)/dt: sign and lag of the strongest correlation between airflow and the derivative of a
    polarity-calibrated belt signal, within lags [-0.5, max_lag] s (thermistors lag the belts).
    Returns (sign, lag_s, r)."""
    dv = np.gradient(belt)
    a, b = flow[mask], dv[mask]
    a = (a - a.mean()) / a.std()
    b = (b - b.mean()) / b.std()
    lags = np.arange(int(-0.5 * fs), int(max_lag_s * fs) + 1, 5)
    r = np.array([np.mean(a[l:] * b[:len(b) - l]) if l >= 0 else np.mean(a[:l] * b[-l:]) for l in lags])
    i = np.argmax(np.abs(r))
    return int(np.sign(r[i])), lags[i] / fs, float(r[i])
