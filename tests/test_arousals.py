import numpy as np
import pytest
from scipy import signal

from snoredreamt import arousals as A

FS = 100


# ---------------------------------------------------------------- literal (loop) transcriptions of MATLAB code
def matlab_baseline(T, step, p):
    size = A._mround(T / step)
    nb = len(A._mcolon(1, step, T))
    bg = np.full(size, p[:nb].mean())
    nxt = 1
    out = np.empty(len(p))
    for k in range(1, len(p) + 1):
        out[k - 1] = bg.sum() / size
        nxt += 1
        bg[nxt - 1] = p[k - 1]
        if nxt == size:
            nxt = 0
    return out


def matlab_signal_amplitude(x, rate, size, step):
    ss, sz = rate * step, rate * size
    cw = int(np.floor((len(x) - sz) / ss))
    amps = []
    for w in range(cw):
        seg = x[int(np.floor(w * ss)):int(np.floor(w * ss)) + int(np.floor(sz))]
        amps.append(seg.max() - seg.min())
    cs = cw * ss + 1
    while cs <= len(x):
        seg = x[A._mround(cs) - 1:]
        amps.append(seg.max() - seg.min())
        cs += ss
    return np.array(amps)


def naive_band_power(x, rate, size, step, lo, hi):
    N = int(size * rate)
    ham = signal.windows.hamming(N, sym=True)
    fstep = rate / N
    f1, f2 = int(np.ceil(lo / fstep)) + 1, min(int(np.ceil(hi / fstep)) + 1, N // 2 + 1)
    out = []
    for s in range(0, len(x) - N, int(rate * step)):
        mx = np.abs(np.fft.fft(x[s:s + N] * ham)) ** 2 / N
        out.append(mx[f1 - 1:f2].sum() * 2 / N)          # f1 > 1 and f2 <= N/2 branch
    return np.array(out)


def test_baseline_matches_matlab_loop():
    p = np.random.default_rng(0).gamma(2, size=800)
    for T, step in [(10, 0.2), (30, 0.1)]:
        assert np.allclose(A.baseline(T, step, p), matlab_baseline(T, step, p))


@pytest.mark.parametrize('size,step', [(0.2, 0.1), (0.1, 0.1), (1, 1)])
def test_signal_amplitude_matches_matlab_loop(size, step):
    x = np.random.default_rng(1).normal(size=1234)
    assert np.allclose(A.signal_amplitude(x, FS, size, step), matlab_signal_amplitude(x, FS, size, step))


def test_signal_power_matches_naive():
    x = np.random.default_rng(2).normal(size=60 * FS)
    pa, pb = A.signal_power(x, 3, FS, 0.2, [8, 16], [12, 45])
    na = naive_band_power(x, FS, 3, 0.2, 8, 12)
    nb = naive_band_power(x, FS, 3, 0.2, 16, 45)
    assert np.allclose(pa[:len(na)], na) and np.allclose(pb[:len(nb)], nb)
    assert len(pa) == len(A._mcolon(1, 20, len(x)))           # padded to MATLAB's numWindows


def test_event_indexes_one_based():
    ev = np.array([0, 1, 1, 0, 0, 1, 0], np.int8)
    ini, end = A.get_event_indexes(ev)
    assert ini.tolist() == [2, 6] and end.tolist() == [3, 6]


def test_merge_closer_than():
    ev = np.array([0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0], np.int8)
    out = A.merge_events_closer_than(ev, 1.0, 0.2)            # gap < 5 windows merged
    assert out[1:6].all() and not out[6:11].any()


def test_notch_removes_alias_line():
    t = np.arange(120 * FS) / FS
    x = A.polyman_notch(np.sin(2 * np.pi * 20 * t), FS, 20, 1)
    assert np.abs(x[30 * FS:]).max() < 0.02


def test_find_spindles_on_13hz_burst():
    c = A._Ctx()
    rng = np.random.default_rng(3)
    x = rng.normal(0, 1, 40 * FS)
    t = np.arange(len(x)) / FS
    x[10 * FS:16 * FS] += 30 * np.sin(2 * np.pi * 13 * t[10 * FS:16 * FS])
    sm = A.polyman_bandpass(x, FS, 1, 1, 32)
    c.rate, c.eeg = FS, x
    c.cyc_cum = np.concatenate([[0], np.cumsum(A.signal_cycles(sm) > 0)])
    ev = {'start': 10.0, 'duration': 6.0}
    sp = A._find_spindles(ev, c)
    assert sp and sum(d for _, d in sp) > 3
    assert not A._spindle_free(ev, c)


# --------------------------------------------------------------------------------------- synthetic nights
def background(minutes=10, seed=0):
    rng = np.random.default_rng(seed)
    n = minutes * 60 * FS
    b, a = signal.butter(2, 4, fs=FS)
    eeg = signal.lfilter(b, a, rng.normal(0, 60, n)) + rng.normal(0, 2, n)   # delta-dominant sleep EEG
    emg = rng.normal(0, 3, n)
    return eeg, emg, rng


def add_arousal(eeg, emg, rng, t0, dur, alpha=20.0, beta=8.0, emg_gain=5.0):
    s, e = int(t0 * FS), int((t0 + dur) * FS)
    t = np.arange(e - s) / FS
    eeg[s:e] += alpha * np.sin(2 * np.pi * 10 * t) + beta * rng.normal(size=e - s)
    if emg_gain:
        emg[s:e] += rng.normal(0, 3 * emg_gain, e - s)


def detect(eeg, emg, stage, **kw):
    return A.detect_arousals(eeg, emg, None, stage, FS, **kw)


def test_detects_arousal_with_onset_precision():
    eeg, emg, rng = background()
    add_arousal(eeg, emg, rng, 300.0, 6.0)
    df = detect(eeg, emg, np.full(len(eeg), 2, np.int8))
    assert len(df) == 1
    ev = df.iloc[0]
    assert abs(ev.onset_s - 300.0) <= 1.0
    assert 3.0 <= ev.duration_s <= 10.0
    assert ev.stage == 2 and ev.start == pytest.approx(ev.onset_s * FS, abs=1)


def test_spindle_only_burst_rejected():
    eeg, emg, rng = background()
    t = np.arange(2 * FS) / FS
    for t0 in (200, 320, 440):
        s = t0 * FS
        eeg[s:s + 2 * FS] += 40 * np.sin(2 * np.pi * 13 * t) * np.hanning(2 * FS)
    df = detect(eeg, emg, np.full(len(eeg), 2, np.int8))
    assert len(df) == 0


def test_event_in_wake_removed():
    eeg, emg, rng = background()
    add_arousal(eeg, emg, rng, 300.0, 6.0)
    stage = np.full(len(eeg), 2, np.int8)
    stage[270 * FS:330 * FS] = 0
    df = detect(eeg, emg, stage, return_all=True)
    assert not df.ok.any()
    assert (df.reason == 'awake').any()


def test_short_event_removed():
    eeg, emg, rng = background()
    add_arousal(eeg, emg, rng, 300.0, 1.0, alpha=10.0, beta=4.0, emg_gain=0)
    df = detect(eeg, emg, np.full(len(eeg), 2, np.int8), return_all=True)
    assert not df.ok.any()
    assert (df.reason == 'too short').any()


def test_rem_event_without_emg_removed_and_with_emg_kept():
    eeg, emg, rng = background()
    # beta-dominant bursts (an alpha-rich burst without EMG is already rejected by the spindle rule)
    add_arousal(eeg, emg, rng, 200.0, 6.0, alpha=0, beta=15.0, emg_gain=0)
    add_arousal(eeg, emg, rng, 400.0, 6.0, alpha=0, beta=15.0, emg_gain=5.0)
    df = detect(eeg, emg, np.full(len(eeg), 5, np.int8), return_all=True)
    ok = df[df.ok]
    assert len(ok) == 1 and abs(ok.onset_s.iloc[0] - 400) <= 1.0
    assert ((df.reason == 'rem') & (df.onset_s.sub(200).abs() < 2)).any()


def test_close_events_second_removed():
    eeg, emg, rng = background()
    add_arousal(eeg, emg, rng, 300.0, 5.0)
    add_arousal(eeg, emg, rng, 312.0, 5.0)
    df = detect(eeg, emg, np.full(len(eeg), 2, np.int8), return_all=True)
    assert df.ok.sum() == 1 and df[df.ok].onset_s.iloc[0] < 305
    assert (df.reason == 'within 10 s of previous').any()


def test_valid_frac_column():
    eeg, emg, rng = background()
    add_arousal(eeg, emg, rng, 300.0, 6.0)
    valid = np.ones(len(eeg), bool)
    valid[302 * FS:304 * FS] = False
    df = detect(eeg, emg, np.full(len(eeg), 2, np.int8), valid=valid)
    assert 0 < df.valid_frac.iloc[0] < 1


def test_band_above_nyquist_rejected():
    eeg, emg, _ = background(minutes=3)
    with pytest.raises(ValueError):
        detect(eeg, emg, np.full(len(eeg), 2, np.int8), params={'beta_band': (16, 50)})
