"""Snakemake rule `event_locked_rates`: SO / spindle / arousal rates around respiratory-event onsets and ends (N2/N3),
relative to each subject's clean-NREM rate. Used to choose the post-event exclusion window from the data."""
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt import conditions  # noqa: E402
from snoredreamt.io import load_subject  # noqa: E402

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
fs = cfg['fs']
DER = cfg['paths']['derivatives']
half_pre, half_post, bw = 30, 60, 2.0
edges = np.arange(-half_pre, half_post + bw / 2, bw)
rows = []
for sid in snakemake.params.subjects:  # noqa: F821
    d = load_subject(sid, cfg)
    st = d['stage']
    ar = pd.read_parquet(f'{DER}/arousals/{sid}_primary.parquet')
    sos = pd.read_parquet(f'{DER}/events/{sid}/sos_F4.parquet')
    spi = pd.read_parquet(f'{DER}/events/{sid}/spindles_C4.parquet')
    times = {'SO': sos.trough[sos.nrem_event].to_numpy(), 'spindle': spi.start[spi.nrem_event].to_numpy(),
             'arousal': ar.start.to_numpy()}
    clean = conditions.clean_mask(d, cfg, ar)
    clean_ar = conditions.clean_mask(d, cfg, ar, exclude_arousals=False)   # baseline for arousals themselves
    cmin = clean.sum() / fs / 60
    s, e = conditions.runs(conditions.respiratory_events(d))
    for anchor_name, anchor in [('onset', s), ('end', e)]:
        anchor = anchor[np.isin(st[anchor], [2, 3])]
        if len(anchor) < 5 or cmin < 10:
            continue
        for ev, t in times.items():
            m = clean_ar if ev == 'arousal' else clean
            base = np.sum(m[t]) / (m.sum() / fs / 60) if len(t) else np.nan   # events/min in clean NREM
            lags = np.concatenate([(t[(t >= a - half_pre * fs) & (t < a + half_post * fs)] - a) / fs for a in anchor])
            h = np.histogram(lags, edges)[0] / len(anchor) / (bw / 60)  # events/min
            rows.append(pd.DataFrame({'sid': sid, 'anchor': anchor_name, 'event': ev, 'lag_s': edges[:-1] + bw / 2,
                                      'rate': h, 'baseline': base, 'n_resp_events': len(anchor)}))
T = pd.concat(rows)
T['rel'] = T.rate / T.baseline
T.to_parquet(out.table)

fig, axs = plt.subplots(2, 3, figsize=(15, 7), sharex=True)
for i, anchor in enumerate(['onset', 'end']):
    for j, ev in enumerate(['SO', 'spindle', 'arousal']):
        p = T[(T.anchor == anchor) & (T.event == ev) & np.isfinite(T.rel)].pivot(index='sid', columns='lag_s', values='rel')
        ax = axs[i, j]
        ax.plot(p.columns, p.median(), 'k', label='median')
        ax.fill_between(p.columns, p.quantile(.25), p.quantile(.75), color='grey', alpha=.3, label='IQR')
        ax.axhline(1, color='r', lw=.7)
        ax.axvline(0, color='b', lw=.7)
        for x in (10, 15, 20, 30):
            ax.axvline(x, color='orange', lw=.4, ls=':')
        ax.set_title(f'{ev} rate around respiratory-event {anchor} (N={len(p)} subjects)', fontsize=9)
        ax.set_ylabel('rate / clean-NREM rate')
        ax.set_ylim(0, None)
axs[1, 0].set_xlabel('time from event (s)')
axs[0, 0].legend(fontsize=7)
fig.suptitle('Event-locked SO / spindle / arousal rates around apnea/hypopnea onsets and ends (N2/N3); '
             'dotted = candidate post-event windows', fontsize=10)
fig.tight_layout()
fig.savefig(out.fig, dpi=80)

# where does the median relative rate after event end settle within +-20 % of the late (40-60 s) level?
for ev in ['SO', 'spindle', 'arousal']:
    p = T[(T.anchor == 'end') & (T.event == ev) & np.isfinite(T.rel)].pivot(index='sid', columns='lag_s', values='rel')
    med = p.median()
    late = med[(med.index > 40)].mean()
    post = med[med.index > 0]
    dev = (post - late).abs() > 0.2 * late
    settle = post.index[dev.to_numpy()].max() if dev.any() else 0
    print(f'{ev}: late level {late:.2f}; last post-end bin deviating >20% from it: {settle:.0f} s; '
          f'peak {post.max():.2f} at {post.idxmax():.0f} s')
