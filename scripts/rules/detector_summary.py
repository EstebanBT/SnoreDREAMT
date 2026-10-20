"""Snakemake rule `detector_summary`: event densities per stage/subject/channel and specificity ratios."""
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

cfg = snakemake.config  # noqa: F821
STAGES = {0: 'W', 1: 'N1', 2: 'N2', 3: 'N3', 5: 'REM'}
rows = []
for kind in ['sos', 'spindles']:
    for f in getattr(snakemake.input, kind):  # noqa: F821
        f = Path(f)
        sid, ch = f.parent.name, f.stem.split('_')[1]
        ev = pd.read_parquet(f)
        meta = json.load(open(f.with_name(f'{f.stem}_thresholds.json')))
        for s, name in STAGES.items():
            m = meta['valid_minutes'][str(s)]
            e = ev[ev.stage == s]
            rows.append(dict(sid=sid, ch=ch, kind=kind, stage=name, minutes=m, n=len(e),
                             density=len(e) / m if m >= 5 else np.nan,
                             duration=e.duration.median() if len(e) else np.nan,
                             freq=e.freq.median() if kind == 'spindles' and len(e) else np.nan,
                             amp=(e.pk2pk if kind == 'sos' else e.env_max).median() if len(e) else np.nan))
for f in snakemake.input.complexes:  # noqa: F821
    f = Path(f)
    cx = pd.read_parquet(f)
    sid, ch = f.parent.name, f.stem.split('_')[1]
    for s in [2, 3]:
        m = json.load(open(f.with_name(f'sos_{ch}_thresholds.json')))['valid_minutes'][str(s)]
        n = (cx.stage == s).sum()
        rows.append(dict(sid=sid, ch=ch, kind='complexes', stage=STAGES[s], minutes=m, n=n,
                         density=n / m if m >= 5 else np.nan))
T = pd.DataFrame(rows)
T.to_csv(snakemake.output.table, index=False)  # noqa: F821

# figure: density per stage, per event type and channel
kinds = ['sos', 'spindles', 'complexes']
chs = sorted(T.ch.unique())
fig, axs = plt.subplots(len(chs), 3, figsize=(13, 4 * len(chs)), squeeze=False)
for i, ch in enumerate(chs):
    for j, k in enumerate(kinds):
        sub = T[(T.ch == ch) & (T.kind == k)]
        order = [s for s in STAGES.values() if s in sub.stage.unique()]
        data = [sub[sub.stage == s].density.dropna() for s in order]
        ax = axs[i, j]
        ax.boxplot(data, labels=order, showfliers=False)
        for x, dd in enumerate(data, start=1):
            ax.scatter(np.random.default_rng(0).normal(x, 0.06, len(dd)), dd, s=6, alpha=0.5, color='k')
        ax.set_title(f'{k} — {ch}')
        ax.set_ylabel('events / min (valid data)')
fig.suptitle('Step 1 detector densities (W/N1/REM = detections with NREM thresholds; subjects with >= 5 min)')
fig.tight_layout()
fig.savefig(snakemake.output.fig, dpi=90)  # noqa: F821

# console summary: median densities and within-subject NREM/REM ratios
med = T.pivot_table(index=['kind', 'ch'], columns='stage', values='density', aggfunc='median')
print(med.round(2).to_string())
for k in ['sos', 'spindles']:
    for ch in chs:
        p = T[(T.kind == k) & (T.ch == ch)].pivot(index='sid', columns='stage', values='density')
        r = (p['N2'] / p['REM']).replace([np.inf], np.nan).dropna()
        print(f'{k} {ch}: N2/REM density ratio median {r.median():.2f} (IQR {r.quantile(.25):.2f}-{r.quantile(.75):.2f}), '
              f'subjects with ratio < 2: {(r < 2).sum()}/{len(r)}')
