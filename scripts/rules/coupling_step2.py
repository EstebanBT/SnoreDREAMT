"""Snakemake rule `coupling_step2`: respiration coupling of SOs, spindles, SO-spindle complexes and arousals in
clean NREM for one subject, for every phase variant ({FLOW, ABDOMEN} x {hilbert, cycle}).

Event time points: SO trough (F4), spindle onset (C4, as Schreiner 2023), complex SO trough (F4), arousal onset.
Clean mask: N2/N3 outside apneas/hypopneas (+ margins) and outside detected arousals (conditions.clean_mask).
Arousals: their own onsets in N1-N3 clean breathing (arousals are of course not excluded for themselves).
Negative control: REM 'SOs' detected with N2 thresholds, in clean REM.
"""
import json

import numpy as np
import pandas as pd

from snoredreamt import conditions, coupling, resp
from snoredreamt.io import load_subject

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
sid = snakemake.wildcards.sid  # noqa: F821
fs = cfg['fs']
P = cfg['coupling']
rng = np.random.default_rng(int(sid[1:]))

d = load_subject(sid, cfg)
n = len(d['stage'])
ph = np.load(inp.phase)
info = json.load(open(inp.info))
breaths = pd.read_parquet(inp.breaths)
ar = pd.read_parquet(inp.arousals)

masks = {
    'nrem_clean': conditions.clean_mask(d, cfg, ar),
    'sleep_clean_with_arousals': conditions.clean_mask(d, cfg, ar, stages=[1, 2, 3], exclude_arousals=False),
    'rem_clean': conditions.clean_mask(d, cfg, ar, stages=[5]),
}
sos = pd.read_parquet(inp.sos)
spi = pd.read_parquet(inp.spindles)
events = {
    'SO': (sos.trough[sos.nrem_event].to_numpy(), 'nrem_clean'),
    'spindle': (spi.start[spi.nrem_event].to_numpy(), 'nrem_clean'),
    'complex': (pd.read_parquet(inp.complexes).trough.to_numpy(), 'nrem_clean'),
    'arousal': (ar.start.to_numpy(), 'sleep_clean_with_arousals'),
    # negative control: REM 'SOs' detected with the N2 thresholds (should not lock to breathing unless artifactual)
    'SO_REM_control': (sos.trough[~sos.nrem_event & (sos.stage == 5)].to_numpy(), 'rem_clean'),
}

summ, rates, peths, phases, locked = [], [], [], [], []
# calibrated, z-scored breathing signals for event-locked averages
calib = {}
for ch in ['FLOW', 'ABDOMEN']:
    x = resp.highpass(resp.lowpass(d[ch], fs, cfg['resp']['lowpass_hz']), fs, cfg['resp']['highpass_hz'])
    x = x * (info[ch]['sign'] or 1)
    m = masks['nrem_clean']
    calib[ch] = (x - x[m].mean()) / x[m].std() if m.sum() > 100 else x * np.nan

half = int(P['locked_half_s'] * fs)
for ev_name, (t, mname) in events.items():
    mask = masks[mname]
    t = t[(t >= 0) & (t < n)]
    t = t[mask[t]]
    for ch in ['FLOW', 'ABDOMEN']:
        # event-locked breathing (independent of the phase estimator)
        tt = t[(t > half) & (t < n - half)]
        if len(tt):
            seg = calib[ch][tt[:, None] + np.arange(-half, half + 1)]
            locked.append(pd.DataFrame({'sid': sid, 'event': ev_name, 'channel': ch, 'n': len(tt),
                                        'lag_s': np.arange(-half, half + 1) / fs, 'value': np.nanmean(seg, 0)}))
        # PETH around inhalation peaks of this channel inside the mask
        pk = breaths.peak[(breaths.channel == ch) & breaths.valid].to_numpy()
        pk = pk[mask[pk]]
        h, edges = coupling.peth(t, pk, fs, P['peth_half_s'], P['peth_bin_s'])
        peths.append(pd.DataFrame({'sid': sid, 'event': ev_name, 'channel': ch, 'bin_center': edges[:-1] + P['peth_bin_s'] / 2,
                                   'pct': h, 'n_events_in_window': np.nan}))
        for method in ['hilbert', 'hilbertgated', 'cycle']:
            var = f'{ch}_{method}'
            defined = np.isfinite(ph[var])
            idx = np.flatnonzero(mask & defined)
            if len(idx) < 60 * fs:
                continue
            # occupancy non-uniformity of the raw phase over the mask (0 = time uniform across phase)
            occ_sub = rng.choice(idx, min(len(idx), 200_000), replace=False)
            occ_r = coupling.circ_r(ph[var][occ_sub].astype(float))
            phi = np.full(n, np.nan)
            if P['equalize'][method]:
                # equalised phase (uniform occupancy within the mask, 0 = detected peak); raw phase kept
                eq = coupling.equalizer(ph[var][idx], rng=rng)
                phi[defined] = eq(ph[var][defined])
            else:
                phi[defined] = ph[var][defined]
            te = t[defined[t]]
            a = phi[te]
            phases.append(pd.DataFrame({'sid': sid, 'event': ev_name, 'variant': var, 'sample': te, 'phase': a,
                                        'phase_raw': ph[var][te].astype(float), 'stage': d['stage'][te]}))
            row = dict(sid=sid, event=ev_name, variant=var, n=len(a), mask_min=len(idx) / fs / 60, occupancy_R=occ_r,
                       density=len(a) / (len(idx) / fs / 60) if len(idx) else np.nan)
            if len(a) >= 2:
                p_r, z_r = coupling.rayleigh(a)
                row.update(mean=coupling.circ_mean(a), mvl=coupling.circ_r(a), rayleigh_p=p_r, rayleigh_z=z_r,
                           ppc=coupling.ppc(a))
                if len(idx) > len(a):
                    null = coupling.null_mvl(phi, idx, len(a), P['n_null'], rng)
                    row.update(null_mvl_mean=null.mean(), mvl_z=(row['mvl'] - null.mean()) / null.std(),
                               mvl_p_null=(np.sum(null >= row['mvl']) + 1) / (len(null) + 1))
                r, occ, cnt = coupling.phase_rate(a, phi, idx, fs, P['n_bins'])
                rates.append(pd.DataFrame({'sid': sid, 'event': ev_name, 'variant': var, 'bin': np.arange(P['n_bins']),
                                           'rate': r, 'occupancy_min': occ, 'count': cnt}))
            summ.append(row)

pd.DataFrame(summ).to_parquet(out.summary)
pd.concat(rates).to_parquet(out.rates)
pd.concat(peths).to_parquet(out.peth)
pd.concat(phases).to_parquet(out.phases)
pd.concat(locked).to_parquet(out.locked)
