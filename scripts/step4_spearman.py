"""Step 4 follow-up: plain Spearman correlations (no covariates) between severity measures and respiration modulation.

Two per-subject modulation measures in clean N2 (primaries per family, all events):
- projected coupling: occupancy-normalised resultant projected on the leave-one-out group phase (polarity-unambiguous
  subjects; from Step 3b `subject_coupling.csv`);
- modulation strength: pairwise phase consistency (PPC) of event phases — direction-free, so independent of polarity
  calibration (all subjects) and, unlike MVL, not inflated in subjects with few events.
BH-FDR over all measure x marker x family correlations.
"""
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt.coupling import ppc  # noqa: E402
from snoredreamt.io import load_config, participants  # noqa: E402

cfg = load_config()
DER, FIG = Path(cfg['paths']['derivatives']), Path(cfg['paths']['figures'])
MIN = {'SO': 20, 'spindle': 20, 'complex': 10}
Sv = pd.DataFrame([json.load(open(f)) for f in sorted(DER.glob('step4/S*_severity.json'))]).set_index('sid')
Sv = Sv.join(participants(cfg))
markers = {'AHI': Sv.AHI, 'OAHI': Sv.OAHI, 'ODI3': Sv.odi3, 'T90': Sv.t90, 'hypoxic burden': Sv.hypoxic_burden,
           'mean event duration': Sv.mean_event_duration_s, 'snore index': Sv.snore_breath_frac,
           'mean SaO2': Sv.mean_sao2, 'arousal index (clinical)': Sv['Arousal Index'],
           'arousal index (detected)': Sv.arousal_index_detected, 'age': Sv.AGE, 'BMI': Sv.BMI}

U = pd.read_csv(DER / 'step3b/N2/subject_coupling.csv')
U = U[U.primary & (U.exclusion == 'none') & (U.condition == 'clean')]
prim = U.groupby('family')[['event_set', 'variant']].first()

E = pd.concat(pd.read_parquet(f) for f in sorted(DER.glob('step3b/S*/events.parquet')))
E = E[(E.stage == 2) & (E.label == 0)]
measures = {}
for fam, r in prim.iterrows():
    u = U[U.family == fam]
    measures[(fam, 'projected coupling')] = u[u.n >= MIN[fam]].set_index('sid').pcw
    g = E[E.event_set == r.event_set].dropna(subset=[f'phase_{r.variant}'])
    cnt = g.groupby('sid').size()
    p = g.groupby('sid')[f'phase_{r.variant}'].apply(lambda a: ppc(a.to_numpy()))
    measures[(fam, 'modulation strength (PPC)')] = p[cnt >= MIN[fam]]

rows = []
for (fam, mname), m in measures.items():
    for kname, k in markers.items():
        d = pd.concat([m.rename('m'), k.rename('k')], axis=1, join='inner').dropna()
        rho, p = stats.spearmanr(d.m, d.k)
        rows.append(dict(family=fam, measure=mname, marker=kname, n=len(d), rho=rho, p=p,
                         event_set=prim.loc[fam, 'event_set'], variant=prim.loc[fam, 'variant']))
R = pd.DataFrame(rows)
R['p_fdr'] = multipletests(R.p, method='fdr_bh')[1]
R.to_csv(DER / 'step4/spearman.csv', index=False)
with pd.option_context('display.width', 200):
    print(R.pivot_table(index='marker', columns=['measure', 'family'], values='rho').round(2).to_string())
    print('\nUncorrected p < .05:\n', R[R.p < .05][['family', 'measure', 'marker', 'n', 'rho', 'p', 'p_fdr']].round(4)
          .to_string(index=False))

fig, axs = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
for ax, mname in zip(axs, ['projected coupling', 'modulation strength (PPC)']):
    sub = R[R.measure == mname]
    piv = sub.pivot(index='marker', columns='family', values='rho').reindex(list(markers))
    pp = sub.pivot(index='marker', columns='family', values='p').reindex(list(markers))
    im = ax.imshow(piv.to_numpy(), cmap='RdBu_r', vmin=-.4, vmax=.4, aspect='auto')
    for (yi, xi), v in np.ndenumerate(piv.to_numpy()):
        star = '**' if pp.iloc[yi, xi] < .01 else ('*' if pp.iloc[yi, xi] < .05 else '')
        ax.text(xi, yi, f'{v:.2f}{star}', ha='center', va='center', fontsize=8)
    ns = sub.groupby('family').n.max()
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels([f'{c}\n(N≤{ns[c]})' for c in piv.columns], fontsize=8)
    ax.set_yticks(range(len(markers)))
    ax.set_yticklabels(list(markers), fontsize=8)
    ax.set_title(mname, fontsize=10)
fig.colorbar(im, ax=axs, label='Spearman ρ', shrink=.7)
fig.suptitle('Spearman correlations: severity measures vs respiration modulation in clean N2 (* p<.05, ** p<.01, '
             'uncorrected; none survives FDR)' if R.p_fdr.min() >= .05 else
             'Spearman correlations: severity measures vs respiration modulation in clean N2 (* p<.05, ** p<.01)',
             fontsize=10)
fig.savefig(FIG / 'step4_spearman.png', dpi=80, bbox_inches='tight')
print('\nmin FDR p:', R.p_fdr.min().round(4))
