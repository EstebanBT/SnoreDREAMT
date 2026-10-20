"""Raw-trace examples: quiet vs loud N2 epochs (by within-subject snore percentile)."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import signal
DATA = '/media/longterm1/esteban/DREAMT/data_100Hz'
D = '/home/esteban/Documents/SnoreDREAMT/derivatives/survey'
FIG = '/home/esteban/Documents/SnoreDREAMT/figures'
FS = 100
CH = ['F4-M1', 'C4-M1', 'SNORE', 'PTAF', 'FLOW', 'THORAX']
for sid in sys.argv[1:]:
    E = pd.read_parquet(f'{D}/{sid}_epochs.parquet')
    df = pd.read_csv(f'{DATA}/{sid}_PSG_df.csv', engine='pyarrow', usecols=CH)
    N = E[(E.stage == 'N2') & (E.resp_event == 0)]
    q = N.snore_hf_p90.rank(pct=True)
    picks = [('quiet (~10th pct)', N.epoch[(q - .1).abs().idxmin()]), ('loud (~90th pct)', N.epoch[(q - .9).abs().idxmin()])]
    fig, axs = plt.subplots(len(CH) + 1, 2, figsize=(18, 13), sharex=True)
    for j, (lab, e) in enumerate(picks):
        sl = slice(e * 30 * FS, (e + 1) * 30 * FS); t = np.arange(30 * FS) / FS
        for i, c in enumerate(CH):
            x = df[c].values[sl]
            if c in ['F4-M1', 'C4-M1', 'SNORE']: x = x * 1e6
            axs[i, j].plot(t, x, lw=.5, color='k'); axs[i, j].set_ylabel(c)
        f, tt, Sx = signal.spectrogram(df.SNORE.values[sl], FS, nperseg=50, noverlap=40)
        axs[-1, j].pcolormesh(tt, f, 10 * np.log10(Sx + 1e-20), shading='auto')
        axs[-1, j].set_ylabel('SNORE Hz')
        axs[0, j].set_title(f'{sid} N2 epoch {e} {lab}')
    for i in range(len(CH)):
        lo = min(axs[i, 0].get_ylim()[0], axs[i, 1].get_ylim()[0]); hi = max(axs[i, 0].get_ylim()[1], axs[i, 1].get_ylim()[1])
        axs[i, 0].set_ylim(lo, hi); axs[i, 1].set_ylim(lo, hi)
    fig.tight_layout(); fig.savefig(f'{FIG}/snore_example_{sid}.png', dpi=60)
