"""Snakemake rule `step2_group`: group-level replication of respiration-SO/spindle coupling (Schreiner 2023) and
arousal-breathing coupling, in clean NREM.

Group phase tests (Rayleigh on subject preferred phases, V-test vs 0 = inhalation peak) use subjects with an
unambiguous polarity for the phase channel (primary) and all subjects (sensitivity). Subjects need >= min_events.
"""
import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt import coupling  # noqa: E402

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
MIN = cfg['coupling']['min_events']
EVENTS = ['SO', 'spindle', 'complex', 'arousal', 'SO_REM_control']
LABEL = {'SO': 'SO trough (F4)', 'spindle': 'spindle onset (C4)', 'complex': 'SO-spindle complex (F4)',
         'arousal': 'arousal onset', 'SO_REM_control': 'REM pseudo-SO (control)'}
VARIANTS = ['FLOW_hilbert', 'FLOW_hilbertgated', 'FLOW_cycle', 'ABDOMEN_hilbert', 'ABDOMEN_hilbertgated', 'ABDOMEN_cycle']
COL = {'FLOW': '#1f77b4', 'ABDOMEN': '#d62728'}

S = pd.concat(pd.read_parquet(f) for f in inp.summary)
R = pd.concat(pd.read_parquet(f) for f in inp.rates)
Pt = pd.concat(pd.read_parquet(f) for f in inp.peth)
L = pd.concat(pd.read_parquet(f) for f in inp.locked)
Q = pd.read_csv(inp.resp_summary).set_index('sid')
amb = {'FLOW': set(Q.index[Q.flow_ambiguous]), 'ABDOMEN': set(Q.index[Q.abdomen_ambiguous])}


def subjects(ev, var, polarity_ok=True):
    s = S[(S.event == ev) & (S.variant == var) & (S.n >= MIN)]
    if polarity_ok:
        s = s[~s.sid.isin(amb[var.split('_')[0]])]
    return s


# ------------------------------------------------------------------ statistics table
rows = []
for ev in EVENTS:
    for var in VARIANTS:
        for pol in [True, False]:
            s = subjects(ev, var, pol)
            if len(s) < 3:
                continue
            a = s['mean'].to_numpy()
            pr, zr = coupling.rayleigh(a)
            pv, V = coupling.vtest(a, 0.0)
            zw = s.mvl_z.dropna()
            rows.append(dict(event=ev, variant=var, polarity_unambiguous_only=pol, n_subjects=len(s),
                             median_occupancy_R=s.occupancy_R.median(),
                             median_events=s.n.median(), group_mean_deg=np.degrees(coupling.circ_mean(a)),
                             group_R=coupling.circ_r(a), rayleigh_p=pr, rayleigh_z=zr, vtest_p=pv,
                             median_subject_mvl=s.mvl.median(), median_mvl_z=zw.median(),
                             wilcoxon_mvl_z_gt0_p=stats.wilcoxon(zw, alternative='greater').pvalue if len(zw) > 5 else np.nan,
                             frac_subjects_p_null_lt05=(s.mvl_p_null < 0.05).mean(),
                             frac_subjects_rayleigh_lt05=(s.rayleigh_p < 0.05).mean()))
T = pd.DataFrame(rows)
T.to_csv(out.table, index=False)


# ------------------------------------------------------------------ figures
def polar_panel(ax, s, title):
    a = s['mean'].to_numpy()
    ax.scatter(a, s.mvl.clip(upper=1), s=10, alpha=.6, color='k')
    if len(a) >= 3:
        mu, r = coupling.circ_mean(a), coupling.circ_r(a)
        ax.annotate('', xy=(mu, r), xytext=(0, 0), arrowprops=dict(color='crimson', width=2, headwidth=7))
        pr, _ = coupling.rayleigh(a)
        pv, _ = coupling.vtest(a, 0)
        title += f'\nN={len(a)} mean {np.degrees(mu):.0f}° R={r:.2f}\nRayleigh p={pr:.1g}  V(0°) p={pv:.1g}'
    ax.set_title(title, fontsize=7)
    ax.set_theta_zero_location('E')
    ax.set_ylim(0, max(0.6, s.mvl.max() if len(s) else 0.6))
    ax.set_yticks([])
    ax.set_xticks(np.radians([0, 90, 180, 270]))
    ax.set_xticklabels(['0° insp. peak', '90°', '±180°', '-90°'], fontsize=6)


fig = plt.figure(figsize=(21, 17))
for i, ev in enumerate(EVENTS):
    for j, var in enumerate(VARIANTS):
        ax = fig.add_subplot(len(EVENTS), len(VARIANTS), i * len(VARIANTS) + j + 1, projection='polar')
        polar_panel(ax, subjects(ev, var), f'{LABEL[ev]} | {var}')
fig.suptitle('Step 2: subject preferred respiratory phases in clean NREM (dot radius = subject MVL; red = group mean '
             'vector)\npolarity-unambiguous subjects, >= %d events; phases equalised (uniform occupancy), 0° = '
             'detected peak of the channel; hilbert = peak freq +- 0.05 Hz (not equalised), cycle = equalised' % MIN, fontsize=10)
fig.tight_layout()
fig.savefig(out.polar, dpi=80)


def rate_mod(ev, var, pol=True):
    sids = subjects(ev, var, pol).sid
    r = R[(R.event == ev) & (R.variant == var) & R.sid.isin(sids)].pivot(index='sid', columns='bin', values='rate')
    return 100 * (r.div(r.mean(1), axis=0) - 1)


fig, axs = plt.subplots(len(EVENTS), 2, figsize=(11, 15), sharex=True)
nb = cfg['coupling']['n_bins']
centers = np.degrees(-np.pi + (np.arange(nb) + .5) * 2 * np.pi / nb)
for i, ev in enumerate(EVENTS):
    for j, ch in enumerate(['FLOW', 'ABDOMEN']):
        ax = axs[i, j]
        for method, ls in [('hilbert', '-'), ('hilbertgated', ':'), ('cycle', '--')]:
            m = rate_mod(ev, f'{ch}_{method}')
            if len(m) < 3:
                continue
            mu, se = m.mean(), m.sem()
            ax.plot(centers, mu, ls, color=COL[ch], label=f'{method} (N={len(m)})')
            ax.fill_between(centers, mu - se, mu + se, color=COL[ch], alpha=.15)
        ax.axhline(0, color='grey', lw=.5)
        ax2 = ax.twinx()
        ax2.plot(centers, np.cos(np.radians(centers)), color='grey', lw=.6, alpha=.5)
        ax2.set_yticks([])
        ax.set_title(f'{LABEL[ev]} — {ch}', fontsize=8)
        ax.legend(fontsize=6)
        ax.set_ylabel('rate vs mean (%)')
axs[-1, 0].set_xlabel('respiratory phase (deg, 0 = peak)')
axs[-1, 1].set_xlabel('respiratory phase (deg, 0 = peak)')
fig.suptitle('Occupancy-normalised event rate across the respiratory cycle (group mean ± SEM; grey = cos(phase))',
             fontsize=10)
fig.tight_layout()
fig.savefig(out.rates_fig, dpi=80)

fig, axs = plt.subplots(len(EVENTS), 2, figsize=(11, 15), sharex=True)
for i, ev in enumerate(EVENTS):
    for j, ch in enumerate(['FLOW', 'ABDOMEN']):
        ax = axs[i, j]
        sids = subjects(ev, f'{ch}_hilbert').sid
        p = Pt[(Pt.event == ev) & (Pt.channel == ch) & Pt.sid.isin(sids)].pivot(index='sid', columns='bin_center',
                                                                                values='pct')
        if len(p) >= 3:
            ax.bar(p.columns, p.mean(), width=cfg['coupling']['peth_bin_s'], color=COL[ch], alpha=.7)
            ax.errorbar(p.columns, p.mean(), p.sem(), fmt='none', color='k', lw=.5)
            ax.axhline(100 / p.shape[1], color='grey', ls='--', lw=.7)
        ax.axvline(0, color='k', lw=.7)
        ax.set_title(f'{LABEL[ev]} around {ch} peaks (N={len(p)})', fontsize=8)
        ax.set_ylabel('% of events')
axs[-1, 0].set_xlabel('time from inhalation peak (s)')
axs[-1, 1].set_xlabel('time from inhalation peak (s)')
fig.suptitle('Peri-event time histograms around respiratory peaks (Schreiner 2023 Fig. 2d style)', fontsize=10)
fig.tight_layout()
fig.savefig(out.peth_fig, dpi=80)

fig, axs = plt.subplots(len(EVENTS), 1, figsize=(8, 15), sharex=True)
for i, ev in enumerate(EVENTS):
    ax = axs[i]
    for ch in ['FLOW', 'ABDOMEN']:
        sids = subjects(ev, f'{ch}_hilbert').sid
        l = L[(L.event == ev) & (L.channel == ch) & L.sid.isin(sids)].pivot(index='sid', columns='lag_s', values='value')
        if len(l) >= 3:
            ax.plot(l.columns, l.mean(), color=COL[ch], label=f'{ch} (N={len(l)})')
            ax.fill_between(l.columns, l.mean() - l.sem(), l.mean() + l.sem(), color=COL[ch], alpha=.2)
    ax.axvline(0, color='k', lw=.7)
    ax.set_title(LABEL[ev], fontsize=8)
    ax.set_ylabel('breathing (z, insp. up)')
    ax.legend(fontsize=6)
axs[-1].set_xlabel('time from event (s)')
fig.suptitle('Event-locked breathing (calibrated, z-scored; group mean ± SEM)', fontsize=10)
fig.tight_layout()
fig.savefig(out.locked_fig, dpi=80)

# ------------------------------------------------------------------ dedicated arousal figure
fig = plt.figure(figsize=(16, 8))
for j, var in enumerate(['FLOW_hilbert', 'FLOW_cycle', 'ABDOMEN_hilbert', 'ABDOMEN_cycle']):
    ax = fig.add_subplot(2, 4, j + 1, projection='polar')
    polar_panel(ax, subjects('arousal', var), f'arousal onset | {var}')
ax = fig.add_subplot(2, 4, 5)
for ch in ['FLOW', 'ABDOMEN']:
    m = rate_mod('arousal', f'{ch}_cycle')
    if len(m) >= 3:
        ax.plot(centers, m.mean(), color=COL[ch], label=f'{ch} cycle (N={len(m)})')
        ax.fill_between(centers, m.mean() - m.sem(), m.mean() + m.sem(), color=COL[ch], alpha=.2)
ax.axhline(0, color='grey', lw=.5)
ax.set_xlabel('phase (deg, 0 = peak)')
ax.set_ylabel('arousal rate vs mean (%)')
ax.legend(fontsize=7)
ax = fig.add_subplot(2, 4, 6)
for ch in ['FLOW', 'ABDOMEN']:
    sids = subjects('arousal', f'{ch}_hilbert').sid
    l = L[(L.event == 'arousal') & (L.channel == ch) & L.sid.isin(sids)].pivot(index='sid', columns='lag_s',
                                                                              values='value')
    if len(l) >= 3:
        ax.plot(l.columns, l.mean(), color=COL[ch], label=ch)
        ax.fill_between(l.columns, l.mean() - l.sem(), l.mean() + l.sem(), color=COL[ch], alpha=.2)
ax.axvline(0, color='k', lw=.7)
ax.set_xlabel('time from arousal onset (s)')
ax.set_ylabel('breathing (z, insp. up)')
ax.legend(fontsize=7)
for k, ch in enumerate(['FLOW', 'ABDOMEN']):
    ax = fig.add_subplot(2, 4, 7 + k)
    sids = subjects('arousal', f'{ch}_hilbert').sid
    p = Pt[(Pt.event == 'arousal') & (Pt.channel == ch) & Pt.sid.isin(sids)].pivot(index='sid', columns='bin_center',
                                                                                     values='pct')
    if len(p) >= 3:
        ax.bar(p.columns, p.mean(), width=cfg['coupling']['peth_bin_s'], color=COL[ch], alpha=.7)
        ax.axhline(100 / p.shape[1], color='grey', ls='--', lw=.7)
    ax.axvline(0, color='k', lw=.7)
    ax.set_xlabel(f'time from {ch} peak (s)')
    ax.set_ylabel('% of arousals')
s = S[(S.event == 'arousal') & (S.variant == 'FLOW_cycle')]
fig.suptitle(f'Arousal-breathing coupling in clean N1-N3 (outside apneas/hypopneas +-margins). Arousals per subject: '
             f'median {s.n.median():.0f} (IQR {s.n.quantile(.25):.0f}-{s.n.quantile(.75):.0f}); subjects with >= {MIN}: '
             f'{(s.n >= MIN).sum()}', fontsize=10)
fig.tight_layout()
fig.savefig(out.arousal_fig, dpi=80)

with pd.option_context('display.width', 250, 'display.max_columns', 30):
    print(T[T.polarity_unambiguous_only].round(4).to_string())
