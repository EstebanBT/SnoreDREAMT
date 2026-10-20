import numpy as np
import pandas as pd
import pytest

from snoredreamt import detectors as d
from snoredreamt.io import load_config

FS = 100
CFG = load_config()['detect']


def matlab_envelope_rms(x, n):
    """Literal transcription of MATLAB R2024a envelope.m movrmsSame (1-based loops)."""
    m = x.mean()
    xc = x - m
    nd2 = n >> 1
    n1, n2 = (nd2 - 1, nd2) if n % 2 == 0 else (nd2, nd2)
    y = np.empty_like(xc)
    for i in range(1, len(x) + 1):
        k1, k2 = max(i - n1, 1), min(i + n2, len(x))
        y[i - 1] = np.sqrt(np.sum(xc[k1 - 1:k2] ** 2) / (k2 - k1 + 1))
    return y + m


@pytest.mark.parametrize('n', [20, 21, 4, 5])
def test_envelope_rms_matches_matlab(n):
    x = np.random.default_rng(0).normal(size=500) + 0.3
    assert np.allclose(d.envelope_rms(x, n), matlab_envelope_rms(x, n))


def synthetic_night(seed=1, minutes=20):
    rng = np.random.default_rng(seed)
    n = minutes * 60 * FS
    x = rng.normal(0, 5, n)                       # background
    so_times = np.arange(30, minutes * 60 - 30, 20) * FS
    t = np.arange(-FS, FS) / FS                   # 2 s
    so = -80 * np.sin(np.pi * (t + 1)) * (np.abs(t) <= 1)   # one full 1 Hz-ish cycle, trough first
    for s in so_times:
        x[s - FS:s + FS] += np.where(t < 0, -80 * np.sin(np.pi * (t + 1)), 60 * np.sin(np.pi * t))
    sp_times = so_times[::2] + 5 * FS
    ts = np.arange(FS) / FS
    for s in sp_times:
        x[s:s + FS] += 25 * np.sin(2 * np.pi * 13 * ts) * np.hanning(FS)
    stage = np.full(n, 2, np.int8)
    return x, stage, np.ones(n, bool), so_times, sp_times


def test_detect_sos_finds_inserted_waves():
    x, stage, valid, so_times, _ = synthetic_night()
    ev, thr = d.detect_sos(x, stage, valid, FS, CFG['so'])
    troughs = ev.trough.to_numpy()
    hit = [np.min(np.abs(troughs - (s - FS // 2))) <= 15 for s in so_times]  # trough at -0.5 s
    assert np.mean(hit) > 0.9
    assert len(ev) < 2 * len(so_times)
    assert 2 in thr


def test_detect_spindles_finds_inserted_bursts():
    x, stage, valid, _, sp_times = synthetic_night()
    ev, thr = d.detect_spindles(x, stage, valid, FS, CFG['spindle'])
    hit = [((ev.start <= s + 50) & (ev.end >= s + 50)).any() for s in sp_times]
    assert np.mean(hit) > 0.9
    assert ev.freq.between(11, 15).mean() > 0.8


def test_stage_and_validity_masks_respected():
    x, stage, valid, so_times, _ = synthetic_night()
    stage[: len(stage) // 2] = 5                  # first half REM
    valid[-60 * FS:] = False                      # last minute artifact
    ev, _ = d.detect_sos(x, stage, valid, FS, CFG['so'])
    assert (ev.start >= len(stage) // 2).all()
    assert (ev.end < len(stage) - 60 * FS).all()


def test_complexes():
    sos = pd.DataFrame({'trough': [100, 500, 900]})
    sp = pd.DataFrame({'trough': [150, 700, 905]}, index=[10, 11, 12])
    cx = d.find_complexes(sos, sp, FS, 1.5)
    assert cx.trough.tolist() == [100, 900]
    assert cx.spindle_idx.tolist() == [10, 12]
