"""Exploratory (Step 3 support): where is respiration coupling concentrated in clean NREM?

Per event, y = cos(phase - mu_group_LOSO) (projection on the leave-one-subject-out group phase; polarity-unambiguous
subjects; ABDOMEN gated Hilbert and FLOW Hilbert). Modulators: SO amplitude quartile (within subject), stage,
distance to the nearest apnea/hypopnea (s), subject AHI tertile. Output: table + figure.
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt import conditions  # noqa: E402
from snoredreamt.coupling import circ_mean  # noqa: E402
from snoredreamt.io import load_config, load_subject, participants  # noqa: E402

cfg = load_config()
DER, FIG = Path(cfg['paths']['derivatives']), Path(cfg['paths']['figures'])
fs = cfg['fs']
Q = pd.read_csv(DER / 'step2/resp_summary.csv').set_index('sid')
amb = {'FLOW': set(Q.index[Q.flow_ambiguous]), 'ABDOMEN': set(Q.index[Q.abdomen_ambiguous])}
VARS = ['ABDOMEN_hilbertgated', 'FLOW_hilbert']
S2 = pd.concat(pd.read_parquet(DER / f'step2/{s}/summary.parquet') for s in Q.index)
S2 = S2[(S2.n >= 20) & S2.variant.isin(VARS)]
S2 = S2[[sid not in amb[v.split('_')[0]] for sid, v in zip(S2.sid, S2.variant)]]
ahi = participants(cfg)['AHI']

rows = []
for sid in Q.index:
    e = pd.read_parquet(DER / f'step2/{sid}/event_phases.parquet')
    e = e[e.variant.isin(VARS) & e.event.isin(['SO', 'spindle'])]
    if e.empty:
        continue
    d = load_subject(sid, cfg)
    ev = conditions.respiratory_events(d)
    idx = np.flatnonzero(ev)
    so = pd.read_parquet(DER / f'events/{sid}/sos_F4.parquet').set_index('trough')
    for (event, var), g in e.groupby(['event', 'variant']):
        if sid in amb[var.split('_')[0]]:
            continue
        others = S2[(S2.event == event) & (S2.variant == var) & (S2.sid != sid)]
        if len(others) < 10:
            continue
        mu = circ_mean(others['mean'].to_numpy())
        g = g.copy()
        g['y'] = np.cos(g.phase - mu)
        if len(idx):
            pos = np.searchsorted(idx, g['sample'].to_numpy())
            dl = np.abs(g['sample'].to_numpy() - idx[np.clip(pos - 1, 0, len(idx) - 1)])
            dr = np.abs(idx[np.clip(pos, 0, len(idx) - 1)] - g['sample'].to_numpy())
            g['dist_event_s'] = np.minimum(dl, dr) / fs
        else:
            g['dist_event_s'] = np.inf
        if event == 'SO':
            amp = np.abs(so.reindex(g['sample'].to_numpy()).trough_amp.to_numpy())
            q = np.clip(np.ceil(pd.Series(amp).rank(pct=True).to_numpy() * 4), 1, 4)
            g['amp_q'] = pd.Series(q).map({1: 'Q1 small', 2: 'Q2', 3: 'Q3', 4: 'Q4 large'}).to_numpy()
        g['ahi'] = ahi.get(sid, np.nan)
        rows.append(g)
E = pd.concat(rows)
E['dist_bin'] = pd.cut(E.dist_event_s, [0, 60, 300, np.inf], labels=['<1 min', '1-5 min', '>5 min']).astype(str)
E['ahi_group'] = pd.cut(E.ahi, [-1, 5, 15, 30, 1000], labels=['AHI<5', '5-15', '15-30', '>=30']).astype(str)
E['stage_name'] = E.stage.map({2: 'N2', 3: 'N3'})


def summarise(by):
    def subj_then_group(g):
        s = g.groupby('sid').y.mean()           # subject means, then group mean (subjects weighted equally)
        return pd.Series({'mean_y': s.mean(), 'sem_y': s.sem(), 'n_subjects': len(s), 'n_events': len(g),
                          'event_level_mean': g.y.mean(), 'event_level_sem': g.y.sem()})
    return E.groupby(['event', 'variant', by]).apply(subj_then_group).reset_index().rename(columns={by: 'level'}) \
        .assign(modulator=by)


T = pd.concat([summarise(b) for b in ['amp_q', 'stage_name', 'dist_bin', 'ahi_group']])
T = T[T.level.notna() & (T.level != 'nan')]
T.to_csv(DER / 'step3/explore_coupling_modulators.csv', index=False)
with pd.option_context('display.width', 200):
    print(T.round(4).to_string())

mods = ['amp_q', 'stage_name', 'dist_bin', 'ahi_group']
fig, axs = plt.subplots(2, len(mods), figsize=(17, 7), sharey='row')
for i, ev in enumerate(['SO', 'spindle']):
    for j, m in enumerate(mods):
        ax = axs[i, j]
        for k, var in enumerate(VARS):
            t = T[(T.event == ev) & (T.variant == var) & (T.modulator == m)]
            x = np.arange(len(t)) + (k - .5) * .2
            ax.errorbar(x, t.event_level_mean, t.event_level_sem, fmt='o', capsize=3, label=var)
            ax.set_xticks(np.arange(len(t)))
            ax.set_xticklabels(t.level, fontsize=7, rotation=20)
        ax.axhline(0, color='grey', lw=.6)
        ax.set_title(f'{ev} by {m}', fontsize=9)
    axs[i, 0].set_ylabel('mean cos(phase − group phase)\n(event level ± SEM)')
axs[0, 0].legend(fontsize=7)
fig.suptitle('Where is coupling concentrated? Clean NREM events, projection on leave-one-out group phase', fontsize=10)
fig.tight_layout()
fig.savefig(FIG / 'step3_explore_modulators.png', dpi=80)
