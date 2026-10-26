"""Are SOs / spindles evoked by snore bursts? Event rates around snore-burst peaks vs the equivalent breath moment of
quiet clean breaths (pseudo-burst at the subject's median relative burst position), N2. Generalises
snore_evoked_so.py to large SOs (trough), consensus spindles (onset and envelope peak). Used to choose exclusion
windows for the phase-matched exclusion in coupling_step3b.
"""
import json

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt.io import load_config  # noqa: E402
from snoredreamt.snoring import burst_references  # noqa: E402

cfg = load_config()
fs, DER, FIG = cfg['fs'], cfg['paths']['derivatives'], cfg['paths']['figures']
info = pd.read_csv(f'{DER}/step3/conditions_summary.csv')
sids = info.sid[(info.snoring_N2 >= 10) & (info.clean_N2 >= 10)].tolist()
edges = np.arange(-2, 4.0001, 0.1)
rows = []
for sid in sids:
    refs = burst_references(sid, cfg)
    if refs is None:
        continue
    burst, pseudo = refs['burst'], refs['pseudo']
    ev = pd.read_parquet(f'{DER}/step3b/{sid}/events.parquet')
    sp = pd.read_parquet(f'{DER}/events/{sid}/spindles_C4.parquet')
    sp = sp[sp.nrem_event].set_index('start')
    times = {'SO_large trough': ev[ev.event_set == 'SO_large']['sample'].to_numpy(),
             'spindle_cons onset': ev[ev.event_set == 'spindle_cons']['sample'].to_numpy()}
    times['spindle_cons env peak'] = sp.env_peak.reindex(times['spindle_cons onset']).dropna().to_numpy().astype(int)
    for name, t in times.items():
        t = np.sort(t)
        for rname, ref in [('snore burst peak', burst), ('quiet pseudo-burst', pseudo)]:
            lags = np.concatenate([(t[(t >= r - 2 * fs) & (t < r + 4 * fs)] - r) / fs for r in ref])
            h = np.histogram(lags, edges)[0] / len(ref) / 0.1
            rows.append(pd.DataFrame({'sid': sid, 'event': name, 'ref': rname, 'lag': edges[:-1] + 0.05, 'rate': h}))
R = pd.concat(rows)
R.to_parquet(f'{DER}/step3b/snore_evoked_events.parquet')
summary = {}
fig, axs = plt.subplots(1, 3, figsize=(17, 4.5))
for ax, (ev, g) in zip(axs, R.groupby('event')):
    piv = {k: gg.pivot(index='sid', columns='lag', values='rate') for k, gg in g.groupby('ref')}
    diff = piv['snore burst peak'] - piv['quiet pseudo-burst']
    for k, c in [('snore burst peak', 'tab:orange'), ('quiet pseudo-burst', 'grey')]:
        p = piv[k]
        ax.plot(p.columns, p.mean(), color=c, label=f'{k} (N={len(p)})')
        ax.fill_between(p.columns, p.mean() - p.sem(), p.mean() + p.sem(), color=c, alpha=.25)
    # 0.5 s sliding windows: where is the burst excess significant?
    sig = []
    for lo in np.arange(-1.0, 3.5, 0.25):
        w = (diff.columns >= lo) & (diff.columns < lo + 0.5)
        dd = diff.loc[:, w].mean(1)
        p = stats.wilcoxon(dd).pvalue
        sig.append(dict(window=f'{lo:.2f}-{lo + .5:.2f}', excess=dd.mean(), p=p))
        if p < 0.05 and dd.mean() > 0:
            ax.axvspan(lo, lo + 0.5, color='yellow', alpha=.15)
    summary[ev] = sig
    ax.axvline(0, color='k', lw=.6)
    ax.set_title(f'{ev} around burst peaks (yellow: excess p<.05, 0.5 s windows)', fontsize=8)
    ax.set_xlabel('time from (pseudo-)burst peak (s)')
    ax.set_ylabel('events per reference per s')
    ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(f'{FIG}/step3b_snore_evoked_events.png', dpi=80)
for ev, sig in summary.items():
    print(ev)
    print(pd.DataFrame(sig).round(4).to_string(index=False))
json.dump(summary, open(f'{DER}/step3b/snore_evoked_events.json', 'w'), indent=1, default=float)
