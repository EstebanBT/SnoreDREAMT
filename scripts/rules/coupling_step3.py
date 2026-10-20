"""Snakemake rule `coupling_step3`: respiration coupling per breathing condition vs matched clean controls, one subject.

For each condition (snoring, hypopnea, obstructive, central, mixed, recovery, and the negative control
'clean_split' = random clean runs with snoring-like durations):
  - condition runs: contiguous samples with that label in the selected stages;
  - matched controls: clean samples of the same stage within +-match_window_min of each run's centre, nearest
    first, equal total duration, without reuse within a condition; runs with < 50 % matched duration are dropped;
  - reference phase mu_ref: preferred phase of clean events OUTSIDE that condition's matched-control mask
    (independent data); projected coupling PC = mean cos(phase - mu_ref) (signed, unbiased by n, polarity-free);
  - also PPC, mean angle, MVL, n, density, phase-defined fraction, event amplitude.
A whole-night comparison (condition vs all other clean samples) is stored as sensitivity ('control' = 'all_clean').
"""
import numpy as np
import pandas as pd

from snoredreamt import conditions, coupling
from snoredreamt.io import load_subject

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
sid = snakemake.wildcards.sid  # noqa: F821
fs = cfg['fs']
P3 = cfg['step3']
rng = np.random.default_rng(int(sid[1:]) + 7)

lab_z = np.load(inp.labels)
lab, st = lab_z['labels'], lab_z['stage']
n = len(lab)
ph = np.load(inp.phase)
sos = pd.read_parquet(inp.sos)
spi = pd.read_parquet(inp.spindles)
cx = pd.read_parquet(inp.complexes)
sos = sos[sos.nrem_event]
spi = spi[spi.nrem_event]
big = np.abs(sos.trough_amp) >= np.median(np.abs(sos.trough_amp)) if len(sos) else np.zeros(0, bool)
events = {'SO': (sos.trough.to_numpy(), np.abs(sos.trough_amp.to_numpy())),
          # exploratory: coupling in clean NREM is carried by the larger SOs (scripts/explore_coupling_modulators.py)
          'SO_large': (sos.trough[big].to_numpy(), np.abs(sos.trough_amp[big].to_numpy())),
          'spindle': (spi.start.to_numpy(), spi.env_max.to_numpy()),
          'complex': (cx.trough.to_numpy(), np.abs(cx.trough_amp.to_numpy()))}
codes = cfg['condition_codes']
win = int(P3['match_window_min'] * 60 * fs)


def matched_controls(cond_mask, clean_mask):
    """Return (kept condition mask, matched control mask)."""
    s, e = conditions.runs(cond_mask)
    order = rng.permutation(len(s))
    used = np.zeros(n, bool)
    keep = np.zeros(n, bool)
    clean_idx = np.flatnonzero(clean_mask)
    for k in order:
        a, b = s[k], e[k]
        dur = b - a + 1
        c = (a + b) // 2
        lo, hi = np.searchsorted(clean_idx, [c - win, c + win])
        cand = clean_idx[lo:hi]
        cand = cand[~used[cand]]
        if len(cand) < P3['min_match_frac'] * dur:
            continue
        pick = cand[np.argsort(np.abs(cand - c), kind='stable')[:dur]]
        used[pick] = True
        keep[a:b + 1] = True
    return keep, used


def clean_split_pseudo(clean_mask, template_runs):
    """Negative control: random clean runs with the durations of the snoring runs (or 30 s if none)."""
    durs = template_runs if len(template_runs) else np.full(20, 30 * fs)
    pseudo = np.zeros(n, bool)
    s, e = conditions.runs(clean_mask)
    long_ok = [(a, b) for a, b in zip(s, e) if b - a > 3 * 60 * fs]
    if not long_ok:
        return pseudo
    for dlen in durs:
        a, b = long_ok[rng.integers(len(long_ok))]
        if b - a <= dlen:
            continue
        st0 = rng.integers(a, b - dlen)
        pseudo[st0:st0 + dlen] = True
    return pseudo


# phase variants; methods flagged in config coupling.equalize are equalised over the subject's clean NREM samples
clean_all = (lab == codes['clean']) & np.isin(st, [2, 3])
PH = {}
for var in P3['variants']:
    phi = ph[var].astype(float)
    if cfg['coupling']['equalize'][var.split('_')[1]]:
        idx = np.flatnonzero(clean_all & np.isfinite(phi))
        eq = coupling.equalizer(phi[idx], rng=rng)
        ok = np.isfinite(phi)
        phi[ok] = eq(phi[ok])
    PH[var] = phi

rows, ev_rows = [], []
for stage_set_name, stage_set in P3['stage_sets'].items():
    in_stage = np.isin(st, stage_set)
    clean = (lab == codes['clean']) & in_stage
    sn_s, sn_e = conditions.runs((lab == codes['snoring']) & in_stage)
    cond_masks = {name: (lab == code) & in_stage for name, code in codes.items() if name != 'clean'}
    pseudo = clean_split_pseudo(clean, (sn_e - sn_s + 1))
    cond_masks['clean_split'] = pseudo
    for cname, cmask in cond_masks.items():
        base_clean = clean & ~pseudo if cname == 'clean_split' else clean
        if cmask.sum() < P3['min_condition_s'] * fs:
            continue
        # matched within-stage: match each stage separately so controls share the run's stage
        keep = np.zeros(n, bool)
        ctrl = np.zeros(n, bool)
        for sv in stage_set:
            k, u = matched_controls(cmask & (st == sv), base_clean & (st == sv))
            keep |= k
            ctrl |= u
        for control_name, cond_m, ctrl_m in [('matched', keep, ctrl), ('all_clean', cmask, base_clean)]:
            ref_m = base_clean & ~ctrl_m if control_name == 'matched' else base_clean
            for var in P3['variants']:
                phi = PH[var]
                defined = np.isfinite(phi)
                for ev, (t, amp) in events.items():
                    t_ref = t[ref_m[t] & defined[t]]
                    holdout = np.ones(len(t), bool)       # events allowed on the control side
                    if control_name == 'all_clean':
                        # whole-night control: split clean events (even rank -> mu_ref, odd rank -> control side) so
                        # the control projection is not computed on the events that defined mu_ref (circularity)
                        is_ref = np.zeros(len(t), bool)
                        cand = np.flatnonzero(ref_m[t] & defined[t])
                        is_ref[cand[0::2]] = True
                        t_ref = t[is_ref]
                        holdout = ~is_ref
                    if len(t_ref) < P3['min_ref_events']:
                        continue
                    mu = coupling.circ_mean(phi[t_ref])
                    for side, m in [('condition', cond_m), ('control', ctrl_m)]:
                        sel = m[t] & defined[t] & (holdout if side == 'control' else True)
                        a = phi[t[sel]]
                        minutes = m.sum() / fs / 60
                        row = dict(sid=sid, stage_set=stage_set_name, condition=cname, control=control_name,
                                   variant=var, event=ev, side=side, n=len(a), minutes=minutes,
                                   minutes_phase_defined=(m & defined).sum() / fs / 60,
                                   density=m[t].sum() / minutes if minutes > 0 else np.nan,
                                   amp_median=np.median(amp[m[t]]) if m[t].any() else np.nan,
                                   mu_ref=mu, n_ref=len(t_ref))
                        if control_name == 'matched' and stage_set_name == 'N2':
                            ev_rows.append(pd.DataFrame({'sid': sid, 'condition': cname, 'variant': var, 'event': ev,
                                                         'side': side, 'sample': t[sel], 'phase': a}))
                        if len(a) >= 2:
                            row.update(pc=np.mean(np.cos(a - mu)), ppc=coupling.ppc(a), mean=coupling.circ_mean(a),
                                       mvl=coupling.circ_r(a), vec_x=np.mean(np.cos(a)), vec_y=np.mean(np.sin(a)))
                        rows.append(row)
pd.DataFrame(rows).to_parquet(out.table)
(pd.concat(ev_rows) if ev_rows else pd.DataFrame(columns=['sid', 'condition', 'variant', 'event', 'side', 'sample',
                                                         'phase'])).to_parquet(out.events)
