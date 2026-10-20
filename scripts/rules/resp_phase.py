"""Snakemake rule `resp_phase`: polarity calibration, breath detection and respiratory phase for one subject.

Belts (ABDOMEN, THORAX): polarity from snore / RSA / timing anchors (resp.polarity_anchors, resp.decide_polarity).
Airflow (FLOW, PTAF): polarity and lag from correlation with d/dt of the calibrated ABDOMEN signal.
Phases saved (float16 radians, NaN = undefined): {FLOW,ABDOMEN} x {hilbert,cycle}; 0 = peak (see resp.py).
"""
import json

import numpy as np
import pandas as pd
from scipy import signal

from snoredreamt import ecg, resp
from snoredreamt.io import load_subject
from snoredreamt.quality import load_masks

cfg = snakemake.config  # noqa: F821
p = cfg['resp']
fs = cfg['fs']
d = load_subject(snakemake.wildcards.sid, cfg)  # noqa: F821
masks = load_masks(snakemake.input.quality)  # noqa: F821
st = d['stage']
any_event = d['ev_oa'] | d['ev_ca'] | d['ev_hyp'] | d['ev_mixed']
nrem = np.isin(st, [2, 3])

rp, _ = ecg.detect_rpeaks(d['ECG'], fs, band=cfg['ecg']['band'], min_rr_s=cfg['ecg']['min_rr_s'])
hr = ecg.instantaneous_hr(rp, fs, len(st))
sos = signal.butter(4, p['snore_band'], btype='band', fs=fs, output='sos')
snore_env = resp.lowpass(np.abs(signal.sosfiltfilt(sos, d['SNORE'].astype(float))), fs, 2.0)

sig, info, raw, anchors = {}, {'sid': snakemake.wildcards.sid}, {}, {}  # noqa: F821
for ch in ['ABDOMEN', 'THORAX']:
    ref = nrem & ~any_event & masks[ch]
    raw[ch] = resp.highpass(resp.lowpass(d[ch], fs, p['lowpass_hz']), fs, p['highpass_hz'])
    anchors[ch] = resp.polarity_anchors(raw[ch], fs, ref, hr=hr, snore_env=snore_env)
both = nrem & ~any_event & masks['ABDOMEN'] & masks['THORAX']
r_at = float(np.corrcoef(raw['ABDOMEN'][both], raw['THORAX'][both])[0, 1]) if both.sum() > 600 * fs else 0.0
s_abd, s_tho, conflict, z = resp.joint_belt_polarity(anchors['ABDOMEN'], anchors['THORAX'], r_at, p['strong_z'])
for ch, s in [('ABDOMEN', s_abd), ('THORAX', s_tho)]:
    sig[ch] = raw[ch] * (s or 1)
    info[ch] = dict(anchors[ch], sign=s, conflict=conflict, pooled_snore_z=z['snore'], pooled_rsa_z=z['rsa'],
                    r_abd_tho=r_at)
# volume proxy: sum of z-scored calibrated belts (abdomen only if the thorax relation is unclear)
vol = sig['ABDOMEN'] / sig['ABDOMEN'][both].std()
if abs(r_at) >= 0.3:
    vol = vol + sig['THORAX'] / sig['THORAX'][both].std()

for ch in ['FLOW', 'PTAF']:
    ref = nrem & ~any_event & masks[ch] & masks['ABDOMEN']
    x = resp.highpass(resp.lowpass(d[ch], fs, p['lowpass_hz']), fs, p['highpass_hz'])
    if ref.sum() > 600 * fs:
        sign, lag, r = resp.airflow_polarity(x, vol, fs, ref)
    else:
        sign, lag, r = 0, np.nan, np.nan
    # direct anchor: snore energy should be centred on the inspiratory-flow peak (phase 0) of correctly signed flow
    pk, tr, amp, amp_ref = resp.detect_breaths(x, fs, ref)
    ph = resp.cycle_phase(len(x), pk, tr, amp, amp_ref, fs=fs)
    okb = ref[pk[:-1]] & ref[pk[1:]]
    v = resp._breath_vectors(snore_env, ph, pk[:-1][okb], pk[1:][okb])
    snore_vote, snore_z, _ = resp._vote(v, np.real)
    # The thermistor-belt relation is not consistent across DREAMT subjects (lag -1.5..+1.4 s in subjects with
    # unambiguous anchors), so a strong snore anchor (inspiration in any sensor) takes priority.
    corr_sign = sign
    if abs(snore_z) >= p['strong_z']:
        sign = snore_vote
    conflict = bool(corr_sign and snore_vote and abs(snore_z) >= p['strong_z'] and corr_sign != snore_vote)
    sig[ch] = x * (sign or 1)
    info[ch] = dict(sign=sign, lag_s=lag, r_with_dabdomen=r, snore_z=snore_z, conflict=conflict)

phases, breaths = {}, []
for ch in p['phase_channels']:
    ref = nrem & ~any_event & masks[ch]
    x = sig[ch]
    pk, tr, amp, amp_ref = resp.detect_breaths(x, fs, ref, p['min_breath_s'], p['prominence_frac'])
    valid_pk = masks[ch][pk]
    period = np.diff(pk)[ref[pk[:-1]]] / fs                 # breath periods in clean NREM
    max_cycle = min(p['max_cycle_s'], p['max_cycle_rel'] * np.median(period)) if len(period) else p['max_cycle_s']
    cyc = resp.cycle_phase(len(x), pk, tr, amp, amp_ref, p['min_amp_frac'], max_cycle, fs)
    cyc[~masks[ch]] = np.nan
    f0 = resp.peak_frequency(x, fs, ref)
    if not np.isfinite(f0):            # severe OSA: few event-free minutes -> use all valid NREM
        f0 = resp.peak_frequency(x, fs, nrem & masks[ch])
    if np.isfinite(f0):
        hil = resp.hilbert_phase(x, fs, f0, p['hilbert_half_width'])
        hil, off = resp.align_phase(hil, pk[valid_pk & ref[pk]])
        hil[~masks[ch]] = np.nan
    else:
        hil, off = np.full(len(x), np.nan), np.nan
    phases[f'{ch}_hilbert'] = hil.astype(np.float16)
    phases[f'{ch}_cycle'] = cyc.astype(np.float16)
    # narrow-band phase only where an actual breath is detected in the raw signal (no filter ringing in apneas)
    phases[f'{ch}_hilbertgated'] = np.where(np.isfinite(cyc), hil, np.nan).astype(np.float16)
    info[ch].update(f0=float(f0), hilbert_offset_rad=float(off), n_breaths=int(len(pk)), max_cycle_s=float(max_cycle),
                    amp_ref=float(amp_ref), cycle_defined_frac_nrem=float(np.isfinite(cyc[nrem]).mean()))
    breaths.append(pd.DataFrame({'channel': ch, 'peak': pk, 'trough': tr, 'amp': amp, 'valid': valid_pk,
                                 'stage': st[pk]}))

np.savez_compressed(snakemake.output.phase, **phases)  # noqa: F821
pd.concat(breaths).to_parquet(snakemake.output.breaths)  # noqa: F821
with open(snakemake.output.info, 'w') as f:  # noqa: F821
    json.dump(info, f, indent=1, default=float)
