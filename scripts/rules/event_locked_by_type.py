"""Snakemake rule `event_locked_by_type`: SO / spindle / arousal rates around respiratory-event onsets and ends, split by
event type (hypopnea, obstructive, central, mixed apnea) and pooled.

Same computation as `event_locked_rates` (N2/N3 anchors, rates relative to the subject's clean-NREM rate, subjects with
>= 5 events of the type and >= 10 min clean NREM). Purpose: check that the 15 s post-event (recovery) window, chosen on
the pooled curve (70 % hypopneas), also holds for each event type. Descriptive only.
Summary: mean rate in 4-s windows after event end relative to clean NREM (95 % CI), and the first window whose CI
reaches the clean level. (A median-based settling criterion is unusable here: per-subject rates in 2-s bins are mostly
zero for arousals and for SOs after central apneas.) For obstructive apneas also without the subject contributing most
events.
"""
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt import conditions  # noqa: E402
from snoredreamt.io import load_subject  # noqa: E402

cfg = snakemake.config  # noqa: F821
out = snakemake.output  # noqa: F821
fs = cfg['fs']
DER = cfg['paths']['derivatives']
ARV = cfg['arousals']['exclusion']
half_pre, half_post, bw = 30, 60, 2.0
edges = np.arange(-half_pre, half_post + bw / 2, bw)
TYPES = {'ev_hyp': 'hypopnea', 'ev_oa': 'obstructive', 'ev_ca': 'central', 'ev_mixed': 'mixed'}
rows = []
for sid in snakemake.params.subjects:  # noqa: F821
    d = load_subject(sid, cfg)
    st = d['stage']
    ar = pd.read_parquet(f'{DER}/arousals/{sid}_{ARV}.parquet')
    sos = pd.read_parquet(f'{DER}/events/{sid}/sos_F4.parquet')
    spi = pd.read_parquet(f'{DER}/events/{sid}/spindles_C4.parquet')
    times = {'SO': sos.trough[sos.nrem_event].to_numpy(), 'spindle': spi.start[spi.nrem_event].to_numpy(),
             'arousal': ar.start.to_numpy()}
    clean = conditions.clean_mask(d, cfg, ar)
    clean_ar = conditions.clean_mask(d, cfg, ar, exclude_arousals=False)
    if clean.sum() / fs / 60 < 10:
        continue
    base = {ev: (np.sum((clean_ar if ev == 'arousal' else clean)[t]) /
                 ((clean_ar if ev == 'arousal' else clean).sum() / fs / 60)) if len(t) else np.nan
            for ev, t in times.items()}
    s, e = conditions.runs(conditions.respiratory_events(d))
    kind = np.array([next(TYPES[k] for k in TYPES if d[k][a:b + 1].any()) for a, b in zip(s, e)])
    for typ in list(TYPES.values()) + ['pooled']:
        sel = np.ones(len(s), bool) if typ == 'pooled' else kind == typ
        for anchor_name, anchor in [('onset', s[sel]), ('end', e[sel])]:
            anchor = anchor[np.isin(st[anchor], [2, 3])]
            if len(anchor) < 5:
                continue
            for ev, t in times.items():
                lags = np.concatenate([(t[(t >= a - half_pre * fs) & (t < a + half_post * fs)] - a) / fs
                                       for a in anchor])
                h = np.histogram(lags, edges)[0] / len(anchor) / (bw / 60)
                rows.append(pd.DataFrame({'sid': sid, 'type': typ, 'anchor': anchor_name, 'event': ev,
                                          'lag_s': edges[:-1] + bw / 2, 'rate': h, 'baseline': base[ev],
                                          'n_resp_events': len(anchor)}))
T = pd.concat(rows)
T['rel'] = T.rate / T.baseline
T.to_parquet(out.table)


WINDOWS = [(0, 4), (4, 8), (8, 12), (12, 16), (16, 30)]


def windows(sub, absolute):
    """Mean over subjects (95 % CI) of the rate in 4-s windows after event end, relative to clean NREM. SOs and
    spindles: mean of per-subject relative rates. Arousals: ratio of mean rate to mean clean rate (many subjects have no
    arousal in clean NREM, so per-subject ratios are undefined); CI by bootstrap over subjects."""
    out, rng = {}, np.random.default_rng(0)
    for a, b in WINDOWS:
        w = sub[(sub.lag_s > a) & (sub.lag_s < b)]
        if absolute:
            r = w.groupby('sid').rate.mean()
            bl = sub.groupby('sid').baseline.first().reindex(r.index)
            est = r.mean() / bl.mean()
            boot = []
            for _ in range(1000):
                ix = rng.integers(0, len(r), len(r))
                boot.append(r.iloc[ix].mean() / bl.iloc[ix].mean())
            lo, hi = np.percentile(boot, [2.5, 97.5])
        else:
            r = w[np.isfinite(w.rel)].groupby('sid').rel.mean()
            est = r.mean()
            se = r.std(ddof=1) / np.sqrt(len(r))
            lo, hi = est - 1.96 * se, est + 1.96 * se
        out[f'{a}-{b}s'] = f'{est:.2f} [{lo:.2f}, {hi:.2f}]'
        out[f'_{a}-{b}s_lo'] = lo
    back = next((f'{a}-{b}s' for a, b in WINDOWS if out[f'_{a}-{b}s_lo'] <= 1), 'not within 30 s')
    return {k: v for k, v in out.items() if not k.startswith('_')}, back


summ = []
E = T[T.anchor == 'end']
for typ in list(TYPES.values()) + ['pooled']:
    for ev in ['SO', 'spindle', 'arousal']:
        sub = E[(E.type == typ) & (E.event == ev)]
        if sub.sid.nunique() < 3:
            continue
        variants = [('all subjects', sub)]
        if typ == 'obstructive':
            top = E[E.type == typ].groupby('sid').n_resp_events.first().idxmax()
            variants.append((f'without {top}', sub[sub.sid != top]))
        for vname, v in variants:
            w, back = windows(v, absolute=(ev == 'arousal'))
            summ.append(dict(type=typ, event=ev, variant=vname, n_subjects=v.sid.nunique(),
                             n_events=int(v.groupby('sid').n_resp_events.first().sum()), **w,
                             first_window_reaching_clean_level=back))
S = pd.DataFrame(summ)
S.to_csv(out.summary, index=False)
with pd.option_context('display.width', 250, 'display.max_columns', 20):
    print(S.to_string(index=False))

# figure: rates after event END by type (decision-relevant for the recovery window)
plt.rcParams.update({'font.size': 12, 'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 200, 'savefig.bbox': 'tight'})
COL = {'SO': '#0072B2', 'spindle': '#D55E00', 'arousal': '#CC79A7'}
cols = ['hypopnea', 'obstructive', 'central', 'pooled']
fig, axs = plt.subplots(2, 4, figsize=(16, 6.5), sharex=True, sharey='row', gridspec_kw={'height_ratios': [1.3, 1]})
for j, typ in enumerate(cols):
    for i, (evs, col) in enumerate([(['SO', 'spindle'], 'rel'), (['arousal'], 'rate')]):
        ax = axs[i, j]
        for ev in evs:
            p = E[(E.type == typ) & (E.event == ev)].pivot_table(index='sid', columns='lag_s', values=col)
            mu = p.mean()
            se = p.std(ddof=1) / np.sqrt(p.notna().sum())
            ax.plot(p.columns, mu, color=COL[ev], lw=2, label=ev)
            ax.fill_between(p.columns, mu - se, mu + se, color=COL[ev], alpha=.2)
        ax.axvspan(0, 15, color='#E69F00', alpha=.15, lw=0)
        ax.axvline(0, color='k', lw=.8)
        if col == 'rel':
            ax.axhline(1, color='grey', lw=.8, ls='--')
        else:
            ax.axhline(E[(E.type == typ) & (E.event == 'arousal')].groupby('sid').baseline.first().mean(),
                       color='grey', lw=.8, ls='--')
        if i == 0:
            n = E[(E.type == typ)].sid.nunique()
            ax.set_title(f'{typ} (N = {n} subjects)')
        else:
            ax.set_xlabel('time from event end (s)')
axs[0, 0].set_ylabel('rate relative\nto clean NREM')
axs[1, 0].set_ylabel('arousals per min\n(dashed = clean NREM)')
axs[0, 0].legend(loc='upper right')
axs[1, 0].legend(loc='upper right')
fig.suptitle('Event rates after the end of each respiratory-event type (mean ± SEM; shaded = 15 s recovery window)',
             y=1.02)
fig.savefig(out.fig)
