"""Confound check for the snore-loudness dose-response: are SOs *evoked* by snore bursts (auditory K-complexes)?

Per subject with snoring in N2:
- snoring breaths: reference = time of the snore-envelope maximum within the breath (burst peak);
- quiet clean breaths: reference = the same relative position within the breath (median burst position of that
  subject's snoring breaths), i.e. the equivalent respiratory moment without a burst.
PETH of large-SO troughs (F4, >= subject median amplitude) around the references, -2..+3 s, 100 ms bins, as rate
(events per reference per second). Evoked SOs predict an excess at ~0.3-1 s after burst peaks only.
Also: large-SO probability per breath as a function of within-subject snore loudness (dB).
"""
import json

import matplotlib
import numpy as np
import pandas as pd
from scipy import signal, stats

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt import resp  # noqa: E402
from snoredreamt.io import load_config, load_subject  # noqa: E402

cfg = load_config()
fs, DER, FIG = cfg['fs'], cfg['paths']['derivatives'], cfg['paths']['figures']
info = pd.read_csv(f'{DER}/step3/conditions_summary.csv')
sids = info.sid[(info.snoring_N2 >= 10) & (info.clean_N2 >= 10)].tolist()
edges = np.arange(-2, 3.0001, 0.1)
rows, perbreath = [], []
for sid in sids:
    d = load_subject(sid, cfg)
    st = d['stage']
    lab = np.load(f'{DER}/step3/{sid}/conditions_primary.npz')['labels']
    br = pd.read_parquet(f'{DER}/step3/{sid}/snore_breaths_primary.parquet')
    so = pd.read_parquet(f'{DER}/events/{sid}/sos_F4.parquet')
    so = so[so.nrem_event]
    so = so[np.abs(so.trough_amp) >= np.median(np.abs(so.trough_amp))]
    tr = np.sort(so.trough.to_numpy())
    sos_ = signal.butter(4, cfg['snoring']['band'], btype='band', fs=fs, output='sos')
    env = resp.lowpass(np.abs(signal.sosfiltfilt(sos_, d['SNORE'].astype(float))), fs, 5.0)
    br = br[(st[br.start.to_numpy()] == 2)]
    sn = br[br.snore & (lab[br.start.to_numpy()] == 1)]
    qu = br[~br.snore & (lab[br.start.to_numpy()] == 0)]
    if len(sn) < 50 or len(qu) < 50:
        continue
    peaks = np.array([a + np.argmax(env[a:b]) for a, b in zip(sn.start, sn.end)])
    relpos = np.median((peaks - sn.start.to_numpy()) / (sn.end - sn.start).to_numpy())
    qref = (qu.start + relpos * (qu.end - qu.start)).round().astype(int).to_numpy()
    for name, ref in [('snore burst peak', peaks), ('quiet breath, same breath position', qref)]:
        lags = np.concatenate([(tr[(tr >= r - 2 * fs) & (tr < r + 3 * fs)] - r) / fs for r in ref])
        h = np.histogram(lags, edges)[0] / len(ref) / 0.1
        rows.append(pd.DataFrame({'sid': sid, 'ref': name, 'lag': edges[:-1] + 0.05, 'rate': h, 'n_ref': len(ref)}))
    # per-breath large-SO occurrence vs snore loudness (snoring + quiet clean breaths)
    allb = pd.concat([sn.assign(kind='snore'), qu.assign(kind='quiet')])
    has = np.array([((tr >= a) & (tr < b)).any() for a, b in zip(allb.start, allb.end)])
    perbreath.append(pd.DataFrame({'sid': sid, 'snore_db': allb.snore_db.to_numpy(), 'has_so': has,
                                   'kind': allb.kind.to_numpy()}))
R = pd.concat(rows)
B = pd.concat(perbreath)
B['snore_db_c'] = B.snore_db - B.groupby('sid').snore_db.transform('mean')

piv = {k: g.pivot(index='sid', columns='lag', values='rate') for k, g in R.groupby('ref')}
diff = piv['snore burst peak'] - piv['quiet breath, same breath position']
win = (diff.columns > 0.3) & (diff.columns < 1.0)
early = diff.loc[:, win].mean(1)
res = {'n_subjects': len(diff), 'excess_rate_0.3-1.0s_mean': float(early.mean()),
       'excess_rate_0.3-1.0s_wilcoxon_p': float(stats.wilcoxon(early).pvalue),
       'snoring_rate_0.3-1.0s': float(piv['snore burst peak'].loc[:, win].mean(1).mean()),
       'quiet_rate_0.3-1.0s': float(piv['quiet breath, same breath position'].loc[:, win].mean(1).mean())}
import statsmodels.formula.api as smf  # noqa: E402
m = smf.logit('has_so ~ snore_db_c', B.assign(has_so=B.has_so.astype(int))).fit(
    disp=0, cov_type='cluster', cov_kwds={'groups': pd.factorize(B.sid)[0]})
res.update({'so_per_breath_logit_beta_per_10dB': float(10 * m.params['snore_db_c']),
            'so_per_breath_p': float(m.pvalues['snore_db_c']),
            'p_so_snore_breath': float(B.has_so[B.kind == 'snore'].mean()),
            'p_so_quiet_breath': float(B.has_so[B.kind == 'quiet'].mean())})
print(json.dumps(res, indent=1))
json.dump(res, open(f'{DER}/step3b/snore_evoked_so.json', 'w'), indent=1)

fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
for k, c in [('snore burst peak', 'tab:orange'), ('quiet breath, same breath position', 'grey')]:
    p = piv[k]
    axs[0].plot(p.columns, p.mean(), color=c, label=f'{k} (N={len(p)})')
    axs[0].fill_between(p.columns, p.mean() - p.sem(), p.mean() + p.sem(), color=c, alpha=.25)
axs[0].axvline(0, color='k', lw=.6)
axs[0].axvspan(0.3, 1.0, color='yellow', alpha=.15)
axs[0].set_xlabel('time from reference (s)')
axs[0].set_ylabel('large-SO troughs per reference per s')
axs[0].legend(fontsize=7)
axs[0].set_title('Large SOs around snore-burst peaks vs equivalent quiet-breath moment', fontsize=9)
q = pd.qcut(B.snore_db_c, 10, duplicates='drop')
mm = B.groupby(q, observed=True).has_so.mean()
axs[1].plot(B.groupby(q, observed=True).snore_db_c.mean(), mm, 'o-')
axs[1].set_xlabel('breath snore loudness (dB, within-subject centred)')
axs[1].set_ylabel('P(large SO trough within breath)')
axs[1].set_title('Large-SO occurrence per breath vs snore loudness', fontsize=9)
fig.tight_layout()
fig.savefig(f'{FIG}/step3b_snore_evoked_so.png', dpi=80)
