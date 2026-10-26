"""Snakemake rule `coupling_step3b`: event- and bin-level data for the power-improved Step 3 analyses, one subject.

Event sets (idea 2): SO (F4), SO_large (>= subject median |trough|), SO_cons (F4 SO with a C4/O2 SO trough within
+-0.3 s), SO_large_cons, spindle (C4 onset), spindle_cons (C4 spindle overlapping a spindle on another channel).
Phase variants (idea 2c): BEST_hilbertgated = FLOW gated Hilbert where defined, else ABDOMEN gated Hilbert rotated by
the subject's circular FLOW-ABDOMEN offset (estimated where both are defined in clean NREM); plus the two sources.
Conditions (idea 3a): clean, resp_event (pooled OA+CA+mixed+hypopnea), and the single conditions.
Breath predictors (idea 3b): snore loudness (dB) of the ABDOMEN breath containing the event, airflow (FLOW) and
effort (ABDOMEN) peak-to-peak amplitude of that breath relative to the subject's median clean N2 breath (log2).
Outputs: per-event table, per phase-bin counts + occupancy (idea 4), info json.
"""
import json

import numpy as np
import pandas as pd

from snoredreamt import coupling, resp
from snoredreamt.io import load_subject

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
sid = snakemake.wildcards.sid  # noqa: F821
fs = cfg['fs']
P = cfg['step3b']
codes = cfg['condition_codes']

z = np.load(inp.labels)
lab, st = z['labels'], z['stage']
n = len(lab)
ph = np.load(inp.phase)
d = load_subject(sid, cfg)


def load_events(kind, ch):
    e = pd.read_parquet(f"{cfg['paths']['derivatives']}/events/{sid}/{kind}_{ch}.parquet")
    return e[e.nrem_event].reset_index(drop=True)


# ------------------------------------------------------------------ event sets
so = load_events('sos', 'F4')
amp = np.abs(so.trough_amp.to_numpy())
large = amp >= np.median(amp) if len(amp) else np.zeros(0, bool)
c = P['consensus']['SO']
other_tr = np.sort(np.concatenate([load_events('sos', ch).trough.to_numpy() for ch in c['others']]))
if len(other_tr):
    pos = np.clip(np.searchsorted(other_tr, so.trough.to_numpy()), 1, len(other_tr) - 1)
    near = np.minimum(np.abs(other_tr[pos] - so.trough.to_numpy()), np.abs(other_tr[pos - 1] - so.trough.to_numpy()))
    so_cons = near <= c['max_lag_s'] * fs
else:
    so_cons = np.zeros(len(so), bool)
sp = load_events('spindles', P['consensus']['spindle']['primary'])
cover = np.zeros(n + 1, np.int32)
for ch in P['consensus']['spindle']['others']:
    o = load_events('spindles', ch)
    np.add.at(cover, o.start.to_numpy().astype(np.int64), 1)
    np.add.at(cover, o.end.to_numpy().astype(np.int64) + 1, -1)
cover = np.cumsum(cover[:-1]) > 0
cc = np.concatenate([[0], np.cumsum(cover)])
sp_cons = (cc[sp.end.to_numpy() + 1] - cc[sp.start.to_numpy()]) > 0
# complexes (re-built with the new pieces): F4 SO trough followed by a consensus C4 spindle envelope peak within
# (0, complex_max_lag_s]; three SO sources (all, large, consensus); phase time point = SO trough
spk = np.sort(sp.env_peak.to_numpy()[sp_cons].astype(np.int64))
def has_spindle_after(t):
    i = np.searchsorted(spk, t, side='right')
    ok = i < len(spk)
    lag = np.full(len(t), np.inf)
    lag[ok] = spk[i[ok]] - t[ok]
    return lag <= P['complex_max_lag_s'] * fs
cx = has_spindle_after(so.trough.to_numpy())
sets = {
    'complex': (so.trough.to_numpy()[cx], amp[cx]),
    'complex_large': (so.trough.to_numpy()[cx & large], amp[cx & large]),
    'complex_cons': (so.trough.to_numpy()[cx & so_cons], amp[cx & so_cons]),
    'SO': (so.trough.to_numpy(), amp),
    'SO_large': (so.trough.to_numpy()[large], amp[large]),
    'SO_cons': (so.trough.to_numpy()[so_cons], amp[so_cons]),
    'SO_large_cons': (so.trough.to_numpy()[large & so_cons], amp[large & so_cons]),
    'spindle': (sp.start.to_numpy(), sp.env_max.to_numpy()),
    'spindle_cons': (sp.start.to_numpy()[sp_cons], sp.env_max.to_numpy()[sp_cons]),
}

# ------------------------------------------------------------------ phase variants
clean_nrem = (lab == codes['clean']) & np.isin(st, [2, 3])
fl = ph['FLOW_hilbertgated'].astype(float)
ab = ph['ABDOMEN_hilbertgated'].astype(float)
both = clean_nrem & np.isfinite(fl) & np.isfinite(ab)
offset = coupling.circ_mean(fl[both] - ab[both]) if both.sum() > 60 * fs else np.nan
offset_r = coupling.circ_r(fl[both] - ab[both]) if both.sum() > 60 * fs else np.nan
best = fl.copy()
if np.isfinite(offset):
    fill = ~np.isfinite(fl) & np.isfinite(ab)
    best[fill] = np.angle(np.exp(1j * (ab[fill] + offset)))
PH = {'BEST_hilbertgated': best, 'ABDOMEN_hilbertgated': ab, 'FLOW_hilbertgated': fl}

# ------------------------------------------------------------------ breath predictors
br = pd.read_parquet(inp.breaths)        # ABDOMEN breaths with snore features (conditions rule)
fsig = resp.highpass(resp.lowpass(d['FLOW'], fs, cfg['resp']['lowpass_hz']), fs, cfg['resp']['highpass_hz'])
asig = resp.highpass(resp.lowpass(d['ABDOMEN'], fs, cfg['resp']['lowpass_hz']), fs, cfg['resp']['highpass_hz'])
bs, be = br.start.to_numpy(), br.end.to_numpy()
fl_ptp = np.array([np.ptp(fsig[a:b]) if b > a else np.nan for a, b in zip(bs, be)])
ab_ptp = np.array([np.ptp(asig[a:b]) if b > a else np.nan for a, b in zip(bs, be)])
ref = clean_nrem[bs] & (st[bs] == 2)
br['flow_log2'] = np.log2(fl_ptp / np.nanmedian(fl_ptp[ref])) if ref.any() else np.nan
br['effort_log2'] = np.log2(ab_ptp / np.nanmedian(ab_ptp[ref])) if ref.any() else np.nan


def breath_of(t):
    i = np.searchsorted(bs, t, side='right') - 1
    ok = (i >= 0) & (t < be[np.clip(i, 0, len(be) - 1)])
    return np.where(ok, i, -1)


# ------------------------------------------------------------------ snore-burst references (phase-matched exclusion)
from snoredreamt.snoring import burst_references  # noqa: E402
refs = {'burst': [], 'pseudo': []}
for stg in [2, 3]:
    r_ = burst_references(sid, cfg, stage=stg, default_relpos=P['default_burst_relpos'])
    if r_ is not None:
        refs['burst'].append(r_['burst'])
        refs['pseudo'].append(r_['pseudo'])
ref_all = np.sort(np.concatenate(refs['burst'] + refs['pseudo'])) if any(len(x) for x in refs['burst'] + refs['pseudo']) \
    else np.array([], dtype=np.int64)


def family(name):
    return 'complex' if name.startswith('complex') else ('SO' if name.startswith('SO') else 'spindle')


def evoked_flag(t, fam):
    lo, hi = P['evoked_windows_s'][fam]
    if len(ref_all) == 0:
        return np.zeros(len(t), bool)
    i = np.searchsorted(ref_all, t, side='right') - 1
    lag = np.where(i >= 0, (t - ref_all[np.clip(i, 0, None)]) / fs, np.inf)
    # the window can reach back past the most recent reference only if it is shorter than a breath: check 2 refs
    lag2 = np.where(i >= 1, (t - ref_all[np.clip(i - 1, 0, None)]) / fs, np.inf)
    return ((lag >= lo) & (lag < hi)) | ((lag2 >= lo) & (lag2 < hi))


evoked_mask = {}
for fam, (lo, hi) in P['evoked_windows_s'].items():
    m = np.zeros(n + 1, np.int32)
    if len(ref_all):
        np.add.at(m, np.clip(ref_all + int(lo * fs), 0, n), 1)
        np.add.at(m, np.clip(ref_all + int(hi * fs), 0, n), -1)
    evoked_mask[fam] = np.cumsum(m[:-1]) > 0

# ------------------------------------------------------------------ outputs
ev_rows = []
for name, (t, a) in sets.items():
    t = t[(t >= 0) & (t < n)]
    a = a[:len(t)]
    keep = (lab[t] >= 0) & np.isin(st[t], [2, 3])
    t, a = t[keep], a[keep]
    bi = breath_of(t)
    df = pd.DataFrame({'sid': sid, 'event_set': name, 'sample': t, 'stage': st[t], 'label': lab[t], 'amp': a,
                       'breath': bi})
    for var, phi in PH.items():
        df[f'phase_{var}'] = phi[t]
    for col in ['snore_db', 'flow_log2', 'effort_log2']:
        df[col] = np.where(bi >= 0, br[col].to_numpy()[np.clip(bi, 0, len(br) - 1)], np.nan)
    df['evoked_window'] = evoked_flag(t, family(name))
    ev_rows.append(df)
E = pd.concat(ev_rows)
E.to_parquet(out.events)

cond_masks = {'clean': lab == codes['clean'], 'resp_event': np.isin(lab, P['pooled_event_codes'])}
for name, code in codes.items():
    if name != 'clean':
        cond_masks[name] = lab == code
edges = np.linspace(-np.pi, np.pi, P['n_bins'] + 1)
bin_rows = []
for sname, stages in P['stage_sets'].items():
    in_stage = np.isin(st, stages)
    for cname, cm in cond_masks.items():
        m = cm & in_stage
        if m.sum() < 30 * fs:
            continue
        for var, phi in PH.items():
            for excl in ['none', 'evoked']:
                for name in sets:
                    ok = m & np.isfinite(phi)
                    if excl == 'evoked':
                        ok = ok & ~evoked_mask[family(name)]
                    occ = np.histogram(phi[ok], edges)[0] / fs / 60
                    g = E[(E.event_set == name) & np.isin(E.stage, stages)]
                    tt = g['sample'].to_numpy()
                    a = phi[tt[ok[tt]]]
                    cnt = np.histogram(a, edges)[0]
                    bin_rows.append(pd.DataFrame({'sid': sid, 'stage_set': sname, 'condition': cname, 'variant': var,
                                                  'event_set': name, 'exclusion': excl, 'bin': np.arange(P['n_bins']),
                                                  'bin_center': (edges[:-1] + edges[1:]) / 2, 'count': cnt,
                                                  'occupancy_min': occ}))
pd.concat(bin_rows).to_parquet(out.bins)
with open(out.info, 'w') as f:
    json.dump({'sid': sid, 'flow_abd_offset_rad': float(offset), 'flow_abd_offset_consistency_R': float(offset_r),
               'n_events': {k: int(len(v[0])) for k, v in sets.items()},
               'n_burst_refs': int(sum(len(x) for x in refs['burst'])),
               'n_pseudo_refs': int(sum(len(x) for x in refs['pseudo'])),
               'best_phase_from_flow_frac_nrem': float(np.isfinite(fl[np.isin(st, [2, 3])]).mean())}, f, indent=1)
