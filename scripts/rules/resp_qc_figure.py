"""Snakemake rule `resp_qc_figure`: visual polarity check for one subject.

Three 30-s clean N2 snippets: calibrated FLOW and ABDOMEN (inspiration should be upward), THORAX, SNORE (bursts
should sit on inspiration) and heart rate (should rise during inspiration); plus the anchor statistics.
"""
import json

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt import ecg, resp  # noqa: E402
from snoredreamt.io import load_subject  # noqa: E402
from snoredreamt.quality import load_masks  # noqa: E402

cfg = snakemake.config  # noqa: F821
fs = cfg['fs']
sid = snakemake.wildcards.sid  # noqa: F821
d = load_subject(sid, cfg)
m = load_masks(snakemake.input.quality)  # noqa: F821
info = json.load(open(snakemake.input.info))  # noqa: F821
st = d['stage']
ev = d['ev_oa'] | d['ev_ca'] | d['ev_hyp'] | d['ev_mixed']
clean = np.isin(st, [2, 3]) & ~ev & m['FLOW'] & m['ABDOMEN']
rp, _ = ecg.detect_rpeaks(d['ECG'], fs)
hr = ecg.instantaneous_hr(rp, fs, len(st))

w = 30 * fs
starts = [i for i in range(0, len(st) - w, w) if clean[i:i + w].all()]
picks = [starts[int(q * (len(starts) - 1))] for q in (0.2, 0.5, 0.8)] if len(starts) >= 3 else starts
fig, axs = plt.subplots(5, max(len(picks), 1), figsize=(18, 10), sharex=True, squeeze=False)
for j, s0 in enumerate(picks):
    sl = slice(s0, s0 + w)
    t = np.arange(w) / fs
    for i, ch in enumerate(['FLOW', 'ABDOMEN', 'THORAX']):
        x = resp.highpass(resp.lowpass(d[ch], fs), fs)[sl] * (info[ch]['sign'] or 1)
        axs[i, j].plot(t, x, 'k', lw=.8)
        axs[i, j].set_ylabel(f'{ch} x{info[ch]["sign"]:+d}')
    axs[3, j].plot(t, d['SNORE'][sl], 'k', lw=.3)
    axs[3, j].set_ylabel('SNORE')
    axs[4, j].plot(t, hr[sl], 'r')
    axs[4, j].set_ylabel('HR')
a, f = info['ABDOMEN'], info['FLOW']
fig.suptitle(f'{sid}: ABD sign {a["sign"]:+d} (pooled snore z {a.get("pooled_snore_z", np.nan):.1f}, rsa z '
             f'{a.get("pooled_rsa_z", np.nan):.1f}, r abd-tho {a.get("r_abd_tho", np.nan):.2f}, conflict {a["conflict"]}) | '
             f'FLOW sign {f["sign"]:+d} (snore z {f.get("snore_z", np.nan):.1f}, r with belts {f["r_with_dabdomen"]:.2f})'
             '\nExpected after calibration: inspiration upward in all traces, SNORE bursts and HR rise on the upstroke')
fig.tight_layout()
fig.savefig(snakemake.output[0], dpi=55)  # noqa: F821
