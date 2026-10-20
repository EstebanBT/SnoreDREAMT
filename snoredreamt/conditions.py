"""Sample-level analysis masks.

Step 2 'clean' NREM: scored N2/N3, outside apneas/hypopneas (any type) with `pre_event_s` before onset and
`post_event_s` after the end excluded, and outside detected arousals padded by `exclude_pad_s`.
Snoring is not excluded here (it is a Step 3 condition).
"""
import numpy as np


def runs(mask):
    d = np.diff(np.concatenate([[0], mask.astype(np.int8), [0]]))
    return np.flatnonzero(d > 0), np.flatnonzero(d < 0) - 1


def respiratory_events(d):
    return d['ev_oa'] | d['ev_ca'] | d['ev_hyp'] | d['ev_mixed']


def paint(n, starts, ends, pre, post):
    """Boolean mask True on [start - pre, end + post] for every interval (inclusive ends, sample units)."""
    diff = np.zeros(n + 1, np.int32)
    lo = np.clip(np.asarray(starts) - pre, 0, n)
    hi = np.clip(np.asarray(ends) + post + 1, 0, n)
    np.add.at(diff, lo, 1)
    np.add.at(diff, hi, -1)
    return np.cumsum(diff[:-1]) > 0


def event_exclusion(d, fs, pre_s, post_s):
    s, e = runs(respiratory_events(d))
    return paint(len(d['stage']), s, e, int(pre_s * fs), int(post_s * fs))


def arousal_exclusion(n, arousals, fs, pad_s):
    if arousals is None or len(arousals) == 0:
        return np.zeros(n, bool)
    return paint(n, arousals.start.to_numpy(), arousals.end.to_numpy(), int(pad_s[0] * fs), int(pad_s[1] * fs))


def clean_mask(d, cfg, arousals=None, stages=None, exclude_arousals=True):
    fs, c = cfg['fs'], cfg['clean']
    ok = np.isin(d['stage'], stages if stages is not None else c['stages'])
    ok &= ~event_exclusion(d, fs, c['pre_event_s'], c['post_event_s'])
    if exclude_arousals:
        ok &= ~arousal_exclusion(len(ok), arousals, fs, cfg['arousals']['exclude_pad_s'])
    return ok
