"""Validate snoredreamt.arousals against the clinical DREAMT Arousal Index (per-subject ArI).

Usage: MPLCONFIGDIR=$TMPDIR/mpl .venv/bin/python derivatives/arousals/validate_subset.py [--all]
Writes derivatives/arousals/validation_subset.csv (or validation_all.csv), figures/arousal_validation_subset.png
and figures/arousal_examples.png (subset run only).
"""
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from snoredreamt.arousals import SLEEP_STAGES, detect_arousals  # noqa: E402
from snoredreamt.io import ROOT, load_subject, participants, subject_ids  # noqa: E402
from snoredreamt.quality import load_masks  # noqa: E402

SUBSET = ['S013', 'S049', 'S034', 'S024', 'S084', 'S022', 'S036', 'S070', 'S010', 'S003', 'S040', 'S009',
          'S029', 'S020', 'S008', 'S007', 'S094', 'S037', 'S095', 'S030', 'S055', 'S002', 'S046', 'S059']
OUT = ROOT / 'derivatives' / 'arousals'
FIG = ROOT / 'figures'
BLUE, ORANGE = '#2a78d6', '#eb6834'


def run(sids):
    P = participants()
    rows, events, data = [], {}, {}
    for sid in sids:
        d = load_subject(sid)
        fs = d['fs']
        valid = load_masks(OUT.parent / 'quality' / f'{sid}.npz')['C4-M1']
        t = time.time()
        df = detect_arousals(d['C4-M1'], d['CHIN'], d['ECG'], d['stage'], fs, valid=valid, return_all=True)
        rt = time.time() - t
        ok = df[df.ok]
        sleep_h = np.isin(d['stage'], SLEEP_STAGES).sum() / fs / 3600
        rows.append({'sid': sid, 'clinical_ArI': P.loc[sid, 'Arousal Index'], 'AHI': P.loc[sid, 'AHI'],
                     'sleep_h': sleep_h, 'n_arousals': len(ok), 'ArI': len(ok) / sleep_h,
                     'ArI_clean': (ok.valid_frac >= 0.99).sum() / sleep_h,
                     'n_alpha': (ok.kind == 'alpha').sum(), 'n_beta': (ok.kind == 'beta').sum(),
                     'n_rem': (ok.stage == 5).sum(), 'median_dur_s': ok.duration_s.median(),
                     'rej_awake': (df.reason == 'awake').sum(), 'rej_rem': (df.reason == 'rem').sum(),
                     'rej_close': (df.reason == 'within 10 s of previous').sum(),
                     'rej_noemg': (df.reason == 'no EMG').sum(), 'rej_spindles': (df.reason == 'spindles').sum(),
                     'runtime_s': rt})
        events[sid] = ok
        data[sid] = d
        print(rows[-1], flush=True)
    return pd.DataFrame(rows), events, data


def summarise(tab):
    x, y = tab.clinical_ArI.to_numpy(float), tab.ArI.to_numpy(float)
    r, pr = stats.pearsonr(x, y)
    rho, ps = stats.spearmanr(x, y)
    diff = y - x
    t, pw = stats.wilcoxon(diff)
    b = np.polyfit(x, y, 1)
    s = {'n': len(tab), 'pearson_r': r, 'pearson_p': pr, 'spearman_rho': rho, 'spearman_p': ps,
         'mean_diff_auto_minus_clin': diff.mean(), 'sd_diff': diff.std(ddof=1), 'wilcoxon_p': pw,
         'slope': b[0], 'intercept': b[1],
         'spearman_rho_clean': stats.spearmanr(x, tab.ArI_clean)[0]}
    return s


def scatter(tab, s, path):
    fig, ax = plt.subplots(figsize=(5.2, 5))
    lim = max(tab.clinical_ArI.max(), tab.ArI.max()) * 1.08
    ax.plot([0, lim], [0, lim], color='0.6', lw=1, ls='--', zorder=1)
    ax.text(lim * 0.97, lim * 0.93, 'identity', color='0.4', ha='right', fontsize=8)
    ax.scatter(tab.clinical_ArI, tab.ArI, s=40, color=BLUE, edgecolor='white', linewidth=1.5, zorder=3)
    for _, r in tab[tab.clinical_ArI >= (40 if len(tab) <= 30 else 80)].iterrows():
        ax.annotate(r.sid, (r.clinical_ArI, r.ArI), xytext=(4, 2), textcoords='offset points', fontsize=7,
                    color='0.35')
    xx = np.array([0, lim])
    ax.plot(xx, s['intercept'] + s['slope'] * xx, color=BLUE, lw=1.5, alpha=0.6, zorder=2)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel('Clinical Arousal Index (/h)')
    ax.set_ylabel('Detected arousals per hour of scored sleep')
    ax.set_title(f"Arousal detector vs clinical ArI (n = {s['n']})\n"
                 f"Spearman ρ = {s['spearman_rho']:.2f}, Pearson r = {s['pearson_r']:.2f}, "
                 f"mean diff = {s['mean_diff_auto_minus_clin']:+.1f}/h", fontsize=9)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def examples(events, data, path, n=6, seed=0, half=15):
    rng = np.random.default_rng(seed)
    pool = [(sid, i) for sid, ev in events.items() for i in range(len(ev))]
    pick = [pool[k] for k in rng.choice(len(pool), n, replace=False)]
    stage_names = {0: 'W', 1: 'N1', 2: 'N2', 3: 'N3', 5: 'R', -1: '?'}
    fig, axes = plt.subplots(n, 1, figsize=(11, 2.3 * n))
    for ax, (sid, i) in zip(axes, pick):
        d, ev = data[sid], events[sid].iloc[i]
        fs = d['fs']
        a, b = max(ev.start - half * fs, 0), min(ev.start + half * fs, len(d['stage']))
        t = (np.arange(a, b) - ev.start) / fs
        eeg = d['C4-M1'][a:b] - np.median(d['C4-M1'][a:b])
        emg = d['CHIN'][a:b] - np.median(d['CHIN'][a:b])
        sc_e = np.percentile(np.abs(eeg), 99) or 1
        sc_m = np.percentile(np.abs(emg), 99) or 1
        ax.axvspan(0, ev.duration_s, color=ORANGE, alpha=0.15, lw=0)
        ax.plot(t, 1.2 + np.clip(eeg / sc_e, -1.2, 1.2) * 0.5, color='0.15', lw=0.5)
        ax.plot(t, -0.2 + np.clip(emg / sc_m, -1.2, 1.2) * 0.5, color=BLUE, lw=0.5)
        ax.text(-half - 0.3, 1.2, f'C4-M1\n±{sc_e:.0f} µV', fontsize=7, ha='right', va='center')
        ax.text(-half - 0.3, -0.2, f'CHIN\n±{sc_m:.0f} µV', fontsize=7, ha='right', va='center', color=BLUE)
        ax.set_ylim(-0.9, 1.9)
        st = d['stage'][a:b]
        ch = np.r_[0, np.flatnonzero(np.diff(st)) + 1]
        lab = ' | '.join(f'{(c + a - ev.start) / fs:+.0f}s:{stage_names.get(int(st[c]), st[c])}' for c in ch)
        ax.set_title(f'{sid}  onset {ev.onset_s / 3600:.2f} h  {ev.kind}  dur {ev.duration_s:.1f} s  '
                     f'stage {lab}', fontsize=8, loc='left')
        ax.set_yticks([])
        ax.set_xlim(-half, half)
        for sp in ('top', 'right', 'left'):
            ax.spines[sp].set_visible(False)
    axes[-1].set_xlabel('time from detected onset (s); shaded = detected arousal')
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


if __name__ == '__main__':
    full = '--all' in sys.argv
    sids = subject_ids() if full else SUBSET
    tab, events, data = run(sids)
    s = summarise(tab)
    print(s)
    OUT.mkdir(parents=True, exist_ok=True)
    name = 'validation_all' if full else 'validation_subset'
    tab.to_csv(OUT / f'{name}.csv', index=False)
    pd.Series(s).to_csv(OUT / f'{name}_summary.csv', header=['value'])
    scatter(tab, s, FIG / f'arousal_{name}.png')
    if not full:
        examples(events, data, FIG / 'arousal_examples.png')
