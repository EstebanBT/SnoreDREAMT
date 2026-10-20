"""Explore SNORE channel: per-subject epoch distributions + raw examples of high vs low snore epochs."""
import glob, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import signal
D = '/home/esteban/Documents/SnoreDREAMT/derivatives/survey'
FIG = '/home/esteban/Documents/SnoreDREAMT/figures'
E = pd.concat([pd.read_parquet(f).assign(sid=f.split('/')[-1][:4]) for f in sorted(glob.glob(f'{D}/S*_epochs.parquet'))])
E['nrem'] = E.stage.isin(['N2', 'N3'])
E['lsn'] = np.log10(E.snore_hf_p90)
# per-subject distributions of log snore power in NREM vs wake
fig, axs = plt.subplots(10, 10, figsize=(25, 22), sharex=True)
for ax, (sid, g) in zip(axs.flat, E.groupby('sid')):
    bins = np.linspace(-7, -4, 60)
    for st, c in [('W', 'grey'), ('N2', 'C0'), ('N3', 'C2'), ('R', 'C3')]:
        ax.hist(g.lsn[g.stage == st], bins, histtype='step', color=c, density=True)
    ax.set_title(sid, fontsize=8); ax.tick_params(labelsize=6)
fig.suptitle('log10 SNORE 20-45 Hz RMS (90th pct of 1 s windows per 30 s epoch); grey W, blue N2, green N3, red REM')
fig.tight_layout(); fig.savefig(f'{FIG}/snore_hist_per_subject.png', dpi=70)
# relation snore vs PTAF flutter & resp events
N = E[E.nrem]
print('corr log snore vs log ptaf_hf (NREM, per subject median):',
      N.groupby('sid').apply(lambda g: np.corrcoef(g.lsn, np.log10(g.ptaf_hf_med + 1e-12))[0, 1]).median())
print('corr log snore vs resp_event frac:', N.groupby('sid').apply(lambda g: np.corrcoef(g.lsn, g.resp_event)[0, 1]).median())
S = N.groupby('sid').lsn.describe(percentiles=[.1, .5, .9])
S['spread90_10'] = S['90%'] - S['10%']
print(S.sort_values('spread90_10').round(2).to_string())
