"""Breath-by-breath snoring detection from the 100 Hz SNORE channel (no snoring annotations exist in DREAMT).

Per breath (consecutive ABDOMEN peaks; polarity-free):
- loudness: RMS of SNORE band-passed 20-45 Hz over the breath, in dB re the subject's quiet reference, where the
  quiet reference = min(10th percentile of event-free NREM breaths, median of wake breaths). The wake reference
  handles subjects who snore through most of NREM (their NREM 10th percentile is itself snoring).
- phasicity: ratio of the louder to the quieter half-cycle RMS (snoring is generated in inspiration, so it is
  concentrated in one limb; continuous noise is not).
A breath snores if loudness >= `min_db` and phasicity >= `min_phasic`; episodes need >= `min_run` consecutive
snoring breaths. Absolute SNORE gain differs ~1000x between subjects, so everything is within-subject.
"""
import numpy as np
import pandas as pd
from scipy import signal


def band_rms_cumsum(x, fs, band):
    sos = signal.butter(4, band, btype='band', fs=fs, output='sos')
    y = signal.sosfiltfilt(sos, np.nan_to_num(np.asarray(x, float)))
    return np.concatenate([[0.0], np.cumsum(y * y)])


def _seg_rms(c, a, b):
    return np.sqrt((c[b] - c[a]) / np.maximum(b - a, 1))


def breath_features(snore, peaks, troughs, fs, band=(20, 45), ptaf=None):
    """Per-breath snore features between consecutive peaks (troughs[i] lies between peaks[i], peaks[i+1])."""
    c = band_rms_cumsum(snore, fs, band)
    p0, t0, p1 = peaks[:-1], troughs[:len(peaks) - 1], peaks[1:]
    whole = _seg_rms(c, p0, p1)
    limb1, limb2 = _seg_rms(c, p0, t0), _seg_rms(c, t0, p1)
    df = pd.DataFrame({'start': p0, 'mid': t0, 'end': p1, 'snore_rms': whole,
                       'phasic': np.maximum(limb1, limb2) / np.maximum(np.minimum(limb1, limb2), 1e-12),
                       'duration_s': (p1 - p0) / fs})
    if ptaf is not None:
        cp = band_rms_cumsum(ptaf, fs, band)
        df['ptaf_flutter'] = _seg_rms(cp, p0, p1)
    return df


def label_snoring(br, stage, quiet_mask, wake_mask, min_db=10.0, min_phasic=2.0, min_run=3):
    """Add columns snore_db, snore_breath, snore (episode-level, >= min_run consecutive breaths)."""
    q = br.snore_rms[quiet_mask[br.start.to_numpy()]].to_numpy()
    w = br.snore_rms[wake_mask[br.start.to_numpy()]].to_numpy()
    refs = [np.percentile(q, 10)] if len(q) >= 50 else []
    if len(w) >= 30:
        refs.append(np.median(w))
    ref = min(refs) if refs else np.percentile(br.snore_rms, 10)
    br = br.copy()
    br['snore_db'] = 20 * np.log10(np.maximum(br.snore_rms, 1e-12) / ref)
    br['snore_breath'] = (br.snore_db >= min_db) & (br.phasic >= min_phasic)
    # runs of consecutive snoring breaths (breaths are consecutive when end == next start)
    s = br.snore_breath.to_numpy()
    contiguous = np.concatenate([[False], br.start.to_numpy()[1:] == br.end.to_numpy()[:-1]])
    run_id = np.cumsum(~(s & np.concatenate([[False], s[:-1]]) & contiguous))
    run_len = pd.Series(s).groupby(run_id).transform('sum').to_numpy()
    br['snore'] = s & (run_len >= min_run)
    br.attrs['quiet_ref'] = float(ref)
    return br
