"""Snakemake rule `arousal_agreement`: agreement of the detected arousal index with the clinical arousal index.

Replaces the earlier ratio-vs-clinical correlation, which is misleading: a ratio y/x of two noisy measures falls with x
and rises with y even without any proportional bias. Reported per detector variant (primary, liberal, extended):
  - ICC(2,1) absolute agreement and ICC(3,1) consistency (two-way, single measures; Shrout & Fleiss);
  - Bland-Altman bias (detected - clinical) with 95 % limits of agreement;
  - Spearman rho;
  - log-log slopes in both directions (log detected on log clinical, and the reverse). If under-detection were
    proportionally worse in high-arousal subjects, both slopes would be < 1; slopes bracketing 1 mean that the
    apparent trend is regression to the mean.
The detected index is per hour of scored sleep (N1, N2, N3, REM); the clinical index comes from participant_info.csv
with unknown scoring rules and denominator.
"""
import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt.arousals import SLEEP_STAGES  # noqa: E402
from snoredreamt.io import participants  # noqa: E402

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
fs = cfg['fs']
plt.rcParams.update({'font.size': 12, 'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 200, 'savefig.bbox': 'tight'})


def icc(x, y):
    """Shrout & Fleiss ICC(2,1) (absolute agreement) and ICC(3,1) (consistency) for two raters."""
    Y = np.column_stack([x, y]).astype(float)
    n, k = Y.shape
    gm = Y.mean()
    msr = k * ((Y.mean(1) - gm) ** 2).sum() / (n - 1)
    msc = n * ((Y.mean(0) - gm) ** 2).sum() / (k - 1)
    sse = ((Y - Y.mean(1, keepdims=True) - Y.mean(0, keepdims=True) + gm) ** 2).sum()
    mse = sse / ((n - 1) * (k - 1))
    icc2 = (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
    icc3 = (msr - mse) / (msr + (k - 1) * mse)
    return icc2, icc3


P = participants(cfg)
sleep_h = {}
for f in inp.cache:
    sid = f.split('/')[-1].split('.')[0]
    st = np.load(f)['stage']
    sleep_h[sid] = np.isin(st, SLEEP_STAGES).sum() / fs / 3600
rows, per_subject = [], []
variants = {'primary': inp.primary, 'liberal': inp.liberal, 'extended': inp.extended}
for variant, files in variants.items():
    for f in files:
        sid = f.split('/')[-1].split('_')[0]
        ev = pd.read_parquet(f)
        per_subject.append(dict(variant=variant, sid=sid, n=len(ev), sleep_h=sleep_h[sid],
                                detected=len(ev) / sleep_h[sid], clinical=P.loc[sid, 'Arousal Index']))
S = pd.DataFrame(per_subject).dropna(subset=['clinical'])
for variant, d in S.groupby('variant', sort=False):
    x, y = d.clinical.to_numpy(float), d.detected.to_numpy(float)
    i2, i3 = icc(x, y)
    diff = y - x
    lx, ly = np.log(x[(x > 0) & (y > 0)]), np.log(y[(x > 0) & (y > 0)])
    rows.append(dict(variant=variant, n=len(d), detected_median=np.median(y), clinical_median=np.median(x),
                     icc_agreement=i2, icc_consistency=i3, bias=diff.mean(), loa_lo=diff.mean() - 1.96 * diff.std(ddof=1),
                     loa_hi=diff.mean() + 1.96 * diff.std(ddof=1), spearman=stats.spearmanr(x, y)[0],
                     slope_log_detected_on_clinical=np.polyfit(lx, ly, 1)[0],
                     slope_log_clinical_on_detected=1 / np.polyfit(ly, lx, 1)[0],
                     ratio_vs_clinical_rho=stats.spearmanr(x, y / x)[0], ratio_vs_detected_rho=stats.spearmanr(y, y / x)[0]))
T = pd.DataFrame(rows)
T.to_csv(out.table, index=False)
S.to_csv(out.subjects, index=False)
with pd.option_context('display.width', 200):
    print(T.round(3).to_string(index=False))

COL = {'primary': '#0072B2', 'liberal': '#999999', 'extended': '#E69F00'}
LABEL = {'primary': 'original rules', 'liberal': 'liberal', 'extended': 'extended (used)'}
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
for variant in ['primary', 'extended']:
    d = S[S.variant == variant]
    x, y = d.clinical.to_numpy(float), d.detected.to_numpy(float)
    r = T[T.variant == variant].iloc[0]
    axs[0].scatter(x, y, s=16, alpha=.6, color=COL[variant],
                   label=f'{LABEL[variant]}: ICC = {r.icc_agreement:.2f}, ρ = {r.spearman:.2f}')
    m, dd = (x + y) / 2, y - x
    axs[1].scatter(m, dd, s=16, alpha=.6, color=COL[variant], label=f'{LABEL[variant]}: bias {r.bias:+.1f}/h')
    axs[1].axhline(r.bias, color=COL[variant], lw=1.5)
    for v in (r.loa_lo, r.loa_hi):
        axs[1].axhline(v, color=COL[variant], lw=1, ls='--')
lim = [1, S.clinical.max() * 1.2]
axs[0].plot(lim, lim, 'k--', lw=1, label='identity')
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_xlabel('clinical arousal index (/h)')
axs[0].set_ylabel('detected arousal index (/h)')
axs[0].legend(fontsize=10)
axs[1].axhline(0, color='grey', lw=.8)
axs[1].set_xlabel('mean of detected and clinical (/h)')
axs[1].set_ylabel('detected − clinical (/h)')
axs[1].set_title('Bland–Altman (solid = bias, dashed = 95 % limits)', fontsize=12)
axs[1].legend(fontsize=10, loc='lower left')
fig.savefig(out.fig)
