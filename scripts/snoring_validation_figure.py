"""Visual validation of breath-by-breath snoring labels: for 4 subjects, a 40 s N2 snippet with SNORE, ABDOMEN and
the breath labels, plus per-breath snore loudness vs PTAF flutter (event-free NREM breaths)."""
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt import resp  # noqa: E402
from snoredreamt.io import load_config, load_subject  # noqa: E402

cfg = load_config()
fs = cfg['fs']
DER = cfg['paths']['derivatives']
sids = ['S022', 'S094', 'S046', 'S070']
fig, axs = plt.subplots(3, len(sids), figsize=(20, 9))
for j, sid in enumerate(sids):
    d = load_subject(sid, cfg)
    br = pd.read_parquet(f'{DER}/step3/{sid}/snore_breaths_primary.parquet')
    st = d['stage']
    # snippet: a 40 s N2 window containing both snoring and non-snoring breaths if possible
    br['n2'] = st[br.start.to_numpy()] == 2
    cand = br[br.n2].reset_index(drop=True)
    sw = np.flatnonzero(cand.snore.to_numpy()[1:] != cand.snore.to_numpy()[:-1])
    c = cand.start.iloc[sw[len(sw) // 2]] if len(sw) else cand.start.iloc[len(cand) // 2]
    sl = slice(c - 20 * fs, c + 20 * fs)
    t = np.arange(40 * fs) / fs
    ab = resp.highpass(resp.lowpass(d['ABDOMEN'], fs), fs)[sl]
    axs[0, j].plot(t, d['SNORE'][sl], 'k', lw=.3)
    axs[0, j].set_title(f'{sid}: SNORE (shaded = snoring breaths)', fontsize=9)
    axs[1, j].plot(t, ab, 'k', lw=.8)
    axs[1, j].set_title('ABDOMEN', fontsize=9)
    for _, b in br[(br.start >= sl.start) & (br.end <= sl.stop)].iterrows():
        for ax in axs[:2, j]:
            if b.snore:
                ax.axvspan((b.start - sl.start) / fs, (b.end - sl.start) / fs, color='orange', alpha=.25)
    x = br[br.n2]
    ax = axs[2, j]
    ok = x.ptaf_flutter > 0
    ax.scatter(x.snore_db[ok], np.log10(x.ptaf_flutter[ok]), s=2, c=np.where(x.snore[ok], 'orange', 'grey'), alpha=.4)
    ax.set_xlabel('breath snore loudness (dB re quiet reference)')
    ax.set_ylabel('log10 PTAF 20-45 Hz flutter')
    ax.axvline(cfg['snoring']['min_db'], color='r', lw=.6)
    ax.set_title(f'N2 breaths: {x.snore.mean():.0%} snoring', fontsize=9)
fig.suptitle('Snoring detection check: orange = labelled snoring (>= 10 dB, phasic >= 2, >= 3 breaths)', fontsize=10)
fig.tight_layout()
fig.savefig(f"{cfg['paths']['figures']}/step3_snoring_validation.png", dpi=70)
