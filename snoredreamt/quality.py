"""Signal-quality masks (DREAMT has no artifact annotations).

EEG: per `eeg_win_s` window, robust z (median/MAD over scored-sleep windows) of log SD, < 0.5 Hz power (sweat/drift)
and 20-40 Hz power (muscle); a window is bad if any |z| exceeds `robust_z` or it is flat; bad windows padded.
Respiratory channels: per 30 s epoch, flat, saturated (>= 1 % of samples at the channel's absolute maximum),
or amplitude below `low_amp_frac` of the subject's median sleep-epoch SD.
"""
import numpy as np
from scipy import signal

SLEEP = (1, 2, 3, 5)


def _robust_z(v, ref):
    med = np.nanmedian(v[ref])
    mad = 1.4826 * np.nanmedian(np.abs(v[ref] - med))
    return (v - med) / (mad if mad > 0 else np.nan)


def _band_rms(x, fs, band, btype='band'):
    sos = signal.butter(4, band, btype=btype, fs=fs, output='sos')
    return signal.sosfiltfilt(sos, x)


def _to_samples(win_ok, w, n):
    out = np.zeros(n, bool)
    out[:len(win_ok) * w] = np.repeat(win_ok, w)
    return out


def eeg_mask(x, stage, fs, q):
    """Return (valid samples, per-window table dict)."""
    x = np.asarray(x, float)
    w = int(q['eeg_win_s'] * fs)
    nw = len(x) // w
    W = x[:nw * w].reshape(nw, w)
    sd = W.std(1)
    low = np.sqrt((_band_rms(x, fs, 0.5, 'low')[:nw * w].reshape(nw, w) ** 2).mean(1))
    hf = np.sqrt((_band_rms(x, fs, [20, 40])[:nw * w].reshape(nw, w) ** 2).mean(1))
    ref = np.isin(stage[:nw * w:w], SLEEP) & (sd > q['flat_std'])
    z = {'z_sd': _robust_z(np.log(sd + 1e-9), ref), 'z_low': _robust_z(np.log(low + 1e-9), ref),
         'z_hf': _robust_z(np.log(hf + 1e-9), ref)}
    bad = (sd <= q['flat_std']) | (np.abs(z['z_sd']) > q['robust_z']) | (z['z_low'] > q['robust_z']) | \
          (z['z_hf'] > q['robust_z'])
    pad = int(q['eeg_pad_s'] * fs)
    bad_s = ~_to_samples(~bad, w, len(x))
    if pad:
        bad_s = np.convolve(bad_s, np.ones(2 * pad + 1), 'same') > 0
    return ~bad_s, {'sd': sd, 'bad': bad, **z}


EEG_CHANNELS = ['F4-M1', 'C4-M1', 'O2-M1', 'Fp1-O2', 'T3 - CZ', 'CZ - T4']
RESP_CHANNELS = ['FLOW', 'PTAF', 'THORAX', 'ABDOMEN']


def compute_masks(d, cfg):
    """All masks for one subject: {channel: bool samples}, plus per-window/epoch tables."""
    fs, q, st = cfg['fs'], cfg['quality'], d['stage']
    masks, tables = {}, {}
    for ch in EEG_CHANNELS:
        masks[ch], tables[ch] = eeg_mask(d[ch], st, fs, q)
    for ch in RESP_CHANNELS:
        masks[ch], tables[ch] = resp_mask(d[ch], st, fs, q)
    return masks, tables


def save_masks(path, masks, tables, n):
    out = {'n': n}
    for ch, m in masks.items():
        out['mask__' + ch] = np.packbits(m)
    for ch, t in tables.items():
        for k, v in t.items():
            out[f'tab__{ch}__{k}'] = v
    np.savez_compressed(path, **out)


def load_masks(path):
    """Return {channel: bool sample mask}."""
    with np.load(path) as z:
        n = int(z['n'])
        return {k[6:]: np.unpackbits(z[k])[:n].astype(bool) for k in z.files if k.startswith('mask__')}


def resp_mask(x, stage, fs, q):
    """Return (valid samples, per-epoch bad reasons dict) for a respiratory channel."""
    x = np.asarray(x, float)
    w = int(q['epoch_s'] * fs)
    ne = len(x) // w
    E = x[:ne * w].reshape(ne, w)
    sd = E.std(1)
    amax = np.nanmax(np.abs(x))
    flat = sd <= q['flat_std']
    sat = (np.abs(E) >= 0.999 * amax).mean(1) >= q['saturation_frac']
    sleep = np.isin(stage[:ne * w:w], SLEEP) & ~flat
    low = sd < q['low_amp_frac'] * (np.median(sd[sleep]) if sleep.any() else np.inf)
    bad = flat | sat | low
    return _to_samples(~bad, w, len(x)), {'flat': flat, 'saturated': sat, 'low_amp': low, 'sd': sd}
