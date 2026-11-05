"""Snakemake rule `report_figures`: compact, readable figures for the Quarto report.

The step rules write diagnostic figures with many panels; this rule redraws only the primary analyses from the saved
tables (no new statistics), with few panels and large fonts. Outputs go to figures/report/.
"""
import json

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt import resp  # noqa: E402
from snoredreamt.io import participants  # noqa: E402

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
fs = cfg['fs']

plt.rcParams.update({'font.size': 12, 'axes.titlesize': 13, 'axes.labelsize': 12, 'legend.fontsize': 11,
                     'xtick.labelsize': 11, 'ytick.labelsize': 11, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 200, 'savefig.bbox': 'tight'})
C_ALL, C_EXCL, C_CTRL = '#222222', '#7B3294', '#999999'          # all events / snore-evoked excluded / control
C_EV = {'SO': '#0072B2', 'spindle': '#D55E00', 'complex': '#009E73', 'arousal': '#CC79A7'}
FAMILY_LABEL = {'SO': 'Large SOs', 'spindle': 'Consensus spindles', 'complex': 'Pooled complexes'}
COND_LABEL = {'clean_odd': 'clean vs clean\n(control)', 'resp_event': 'apneas +\nhypopneas', 'hypopnea': 'hypopnea',
              'recovery': 'recovery\n(15 s after)', 'snoring': 'snoring'}
CONDS = ['resp_event', 'hypopnea', 'recovery', 'snoring']


def sem(x, axis=0):
    return np.nanstd(x, axis=axis, ddof=1) / np.sqrt(np.sum(np.isfinite(x), axis=axis))


# ------------------------------------------------------------------ Step 1: detector densities
D = pd.read_csv(inp.detector)
fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
stages = ['W', 'N1', 'N2', 'N3', 'REM']
for ax, (kind, ch, lab) in zip(axs, [('sos', 'F4', 'SOs (F4-M1)'), ('spindles', 'C4', 'Spindles (C4-M1)')]):
    d = D[(D.kind == kind) & (D.ch == ch) & (D.minutes >= 5)]
    vals = [d[d.stage == s].density.dropna().to_numpy() for s in stages]
    ax.boxplot(vals, showfliers=False, widths=.6, medianprops=dict(color='k', lw=2))
    for i, v in enumerate(vals):
        ax.scatter(np.random.default_rng(i).normal(i + 1, .07, len(v)), v, s=6, alpha=.35, color='#0072B2')
    ax.set_xticks(range(1, 6), [f'{s}\n(n={len(v)})' for s, v in zip(stages, vals)])
    ax.set_ylabel('events per minute')
    ax.set_title(lab)
fig.suptitle('Event density per sleep stage (stages with ≥ 5 min; REM and wake use N2 thresholds)', y=1.02)
fig.savefig(out.density)
plt.close(fig)

# ------------------------------------------------------------------ Step 2: phase example (obstructive apnea)
cache = np.load(inp.example_cache)
ph = np.load(inp.example_phase)
info = json.load(open(inp.example_info))
st, oa = cache['stage'], cache['ev_oa'].astype(bool)
dd = np.diff(np.r_[0, oa.astype(int), 0])
on, off = np.flatnonzero(dd == 1), np.flatnonzero(dd == -1)
fl_ph, ab_ph = ph['FLOW_hilbertgated'].astype(float), ph['ABDOMEN_hilbertgated'].astype(float)
runs = [(a, b) for a, b in zip(on, off) if 15 * fs <= b - a <= 30 * fs and a > 25 * fs
        and (st[a - 25 * fs:b + 15 * fs] == 2).all() and not oa[a - 25 * fs:a].any() and not oa[b:b + 15 * fs].any()
        and np.isfinite(fl_ph[a + 2 * fs:b]).mean() < .2 and np.isfinite(ab_ph[a + 2 * fs:b]).mean() > .8
        and np.isfinite(fl_ph[a - 25 * fs:a]).mean() > .9]
a, b = runs[0]
sl = slice(a - 25 * fs, b + 15 * fs)
t = (np.arange(sl.start, sl.stop) - a) / fs
fig, axs = plt.subplots(3, 1, figsize=(10, 6.5), sharex=True)
for ax, ch in zip(axs[:2], ['FLOW', 'ABDOMEN']):
    x = resp.highpass(resp.lowpass(cache[ch], fs), fs) * (info[ch].get('sign') or 1)
    ax.plot(t, x[sl], 'k', lw=1)
    ax.set_ylabel(f'{ch}\n(a.u.)')
    ax.set_yticks([])
for ch, col in [('FLOW', '#0072B2'), ('ABDOMEN', '#D55E00')]:
    axs[2].plot(t, np.degrees(ph[f'{ch}_hilbertgated'][sl].astype(float)), '.', ms=1.5, color=col, label=f'{ch} gated phase')
axs[2].set_ylabel('phase (°)\n0 = inspiration peak')
axs[2].set_yticks([-180, 0, 180])
axs[2].legend(loc='lower left', markerscale=6, ncol=2)
axs[2].set_xlabel('time from apnea onset (s)')
for ax in axs:
    ax.axvspan(0, (b - a) / fs, color='#E69F00', alpha=.2, lw=0)
axs[0].set_title('Obstructive apnea (shaded): airflow stops while the abdominal belt keeps moving', fontsize=13)
fig.savefig(out.phase_example)
plt.close(fig)

# ------------------------------------------------------------------ Step 2: preferred phases (polar)
S = pd.concat(pd.read_parquet(f) for f in inp.step2_summary)
Q = pd.read_csv(inp.resp_summary).set_index('sid')
G = pd.read_csv(inp.step2_group)
G = G[G.polarity_unambiguous_only]
MIN = cfg['coupling']['min_events']
panels = [('SO', 'SO trough'), ('spindle', 'spindle onset'), ('complex', 'SO–spindle complex'),
          ('arousal', 'arousal onset\n(control)'), ('SO_REM_control', 'REM pseudo-SO\n(control)')]
for key, var, sensor in [('polar_flow', 'FLOW_hilbert', 'FLOW'), ('polar_abd', 'ABDOMEN_hilbertgated', 'ABDOMEN')]:
    amb = set(Q.index[Q[f'{sensor.lower()}_ambiguous']])
    fig, axs = plt.subplots(1, 5, figsize=(16, 4.4), subplot_kw={'projection': 'polar'},
                            gridspec_kw={'wspace': .45})
    for ax, (ev, lab) in zip(axs, panels):
        d = S[(S.event == ev) & (S.variant == var) & (S.n >= MIN) & ~S.sid.isin(amb)]
        col = C_EV.get(ev.split('_')[0], '#777777')
        ax.scatter(d['mean'], d.mvl, s=14, color=col, alpha=.6)
        z = np.mean(np.exp(1j * d['mean']))
        rmax = max(.3, d.mvl.quantile(.95) * 1.1) if len(d) else .3
        ax.annotate('', xy=(np.angle(z), np.abs(z) * rmax), xytext=(0, 0),
                    arrowprops=dict(arrowstyle='-|>', color='k', lw=2.5))
        ax.set_ylim(0, rmax)
        ax.set_yticks([])
        ax.set_xticks(np.radians([0, 90, 180, 270]), ['0°', '90°', '180°', '−90°'])
        g = G[(G.event == ev) & (G.variant == var)]
        p = g.rayleigh_p.iloc[0] if len(g) else np.nan
        ax.set_title(f'{lab}\nN={len(d)}, R={np.abs(z):.2f}, p={p:.1g}', fontsize=12, pad=12)
    fig.suptitle(f'Preferred phase of each subject in clean NREM, {sensor} phase (0° = inspiration peak)\n'
                 f'dot = subject (distance from centre = coupling strength); arrow = group mean', y=1.12)
    fig.savefig(out[key])
    plt.close(fig)

# ------------------------------------------------------------------ Step 2: occupancy-normalised rate per phase bin
R = pd.concat(pd.read_parquet(f) for f in inp.step2_rates)
nb = R.bin.max() + 1
centers = np.degrees(-np.pi + (np.arange(nb) + .5) * 2 * np.pi / nb)
fig, axs = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
for ax, ev in zip(axs, ['SO', 'spindle', 'complex']):
    for var, sensor, ls in [('FLOW_hilbert', 'FLOW', '-'), ('ABDOMEN_hilbertgated', 'ABDOMEN', '--')]:
        amb = set(Q.index[Q[f'{sensor.lower()}_ambiguous']])
        n_ok = S[(S.event == ev) & (S.variant == var) & (S.n >= MIN) & ~S.sid.isin(amb)].sid
        d = R[(R.event == ev) & (R.variant == var) & R.sid.isin(n_ok)]
        m = d.pivot_table(index='sid', columns='bin', values='rate')
        m = 100 * (m.div(m.mean(axis=1), axis=0) - 1)
        ax.errorbar(centers, m.mean(), sem(m.to_numpy()), color=C_EV[ev], ls=ls, marker='o', ms=3, capsize=0,
                    label=f'{sensor} (N={len(m)})')
    ax.axhline(0, color='grey', lw=.8)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_xlabel('respiratory phase (°)')
    ax.set_title({'SO': 'SOs', 'spindle': 'Spindles', 'complex': 'Complexes'}[ev])
    ax.legend(loc='upper left', fontsize=10)
axs[0].set_ylabel('event rate vs subject mean (%)')
fig.suptitle('Event rate across the breathing cycle, clean NREM (0° = inspiration peak; mean ± SEM over subjects)',
             y=1.03)
fig.savefig(out.phase_rates)
plt.close(fig)

# ------------------------------------------------------------------ Step 3: event-locked rates
L = pd.read_parquet(inp.locked)
# SOs and spindles: rate relative to the subject's clean-NREM rate. Arousals: absolute rate, because many subjects
# have no arousal in clean NREM (relative rate undefined).
fig, axs = plt.subplots(2, 2, figsize=(11, 6.5), sharex=True, sharey='row', gridspec_kw={'height_ratios': [1.3, 1]})
for j, (anchor, lab) in enumerate([('onset', 'respiratory event onset'), ('end', 'respiratory event end')]):
    for i, (evs, col) in enumerate([(['SO', 'spindle'], 'rel'), (['arousal'], 'rate')]):
        ax = axs[i, j]
        for ev in evs:
            m = L[(L.anchor == anchor) & (L.event == ev)].pivot_table(index='sid', columns='lag_s', values=col)
            mu, se = m.mean(), sem(m.to_numpy())
            ax.plot(m.columns, mu, color=C_EV[ev], lw=2, label=ev)
            ax.fill_between(m.columns, mu - se, mu + se, color=C_EV[ev], alpha=.2)
        ax.axvline(0, color='k', lw=.8)
        if col == 'rel':
            ax.axhline(1, color='grey', lw=.8, ls='--')
        else:
            base = L[(L.event == 'arousal')].groupby('sid').baseline.first().mean()
            ax.axhline(base, color='grey', lw=.8, ls='--')
        if anchor == 'end':
            ax.axvspan(0, 15, color='#E69F00', alpha=.15, lw=0)
        if i == 0:
            ax.set_title(f'aligned to {lab}' + ('\n(shaded = 15 s recovery window)' if anchor == 'end' else ''))
        else:
            ax.set_xlabel(f'time from {lab} (s)')
axs[0, 0].legend(loc='upper left')
axs[1, 0].legend(loc='upper left')
axs[0, 0].set_ylabel('rate relative\nto clean NREM')
axs[1, 0].set_ylabel('arousals per min\n(dashed = clean NREM)')
fig.savefig(out.locked)
plt.close(fig)

# ------------------------------------------------------------------ Step 3b/3c: contrasts and GLMM
CS = pd.read_csv(inp.clean_signal)
PRIM = {r.family: (r.event_set, r.variant) for r in CS[CS.primary].itertuples()}
K = pd.read_csv(inp.contrasts)
K = K[K.primary]
M = pd.read_csv(inp.glmm)
fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
for ax, fam in zip(axs, ['SO', 'spindle', 'complex']):
    k = K[K.family == fam]
    for j, cond in enumerate(CONDS):
        for dy, excl, col, lab in [(.15, 'none', C_ALL, 'all events'), (-.15, 'evoked', C_EXCL, 'snore-evoked excluded')]:
            r = k[(k.condition == cond) & (k.exclusion == excl)]
            if len(r):
                r = r.iloc[0]
                ax.errorbar(r.delta_mean, j + dy, xerr=[[r.delta_mean - r.ci_lo], [r.ci_hi - r.delta_mean]], fmt='o',
                            color=col, ms=6, capsize=3, label=lab if j == 0 else None)
    clean = k.clean_level.median()
    ax.axvline(0, color='grey', lw=1)
    ax.axvline(-clean, color='k', ls=':', lw=1.2, label='complete loss' if fam == 'SO' else None)
    ax.set_yticks(range(len(CONDS)), [COND_LABEL[c] for c in CONDS] if fam == 'SO' else [''] * len(CONDS))
    ax.invert_yaxis()
    ax.set_title(f'{FAMILY_LABEL[fam]} (clean = {clean:.3f})')
    ax.set_xlabel('coupling: condition − clean')
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=3, bbox_to_anchor=(.5, -.08), frameon=False)
fig.suptitle('Coupling strength per breathing condition (N2, occupancy-normalised; mean and 95 % CI)', y=1.02)
fig.savefig(out.contrasts)
plt.close(fig)

fig, axs = plt.subplots(1, 3, figsize=(14, 4.2), sharex=True)
for ax, fam in zip(axs, ['SO', 'spindle', 'complex']):
    es, var = PRIM[fam]
    m = M[(M.event_set == es) & (M.variant == var)]
    for j, cond in enumerate(CONDS):
        for dy, excl, col, lab in [(.15, 'none', C_ALL, 'all events'), (-.15, 'evoked', C_EXCL, 'snore-evoked excluded')]:
            r = m[(m.condition == cond) & (m.exclusion == excl)]
            if len(r):
                r = r.iloc[0]
                v = 100 * np.array([r.rel_change_median, r.rel_change_lo, r.rel_change_hi])
                ax.errorbar(v[0], j + dy, xerr=[[v[0] - v[1]], [v[2] - v[0]]], fmt='o', color=col, ms=6, capsize=3,
                            label=lab if j == 0 else None)
    ax.axvline(0, color='grey', lw=1)
    ax.axvline(-100, color='k', ls=':', lw=1.2)
    ax.axvline(-50, color='k', ls='--', lw=.8)
    ax.set_xlim(-250, 300)
    ax.set_yticks(range(len(CONDS)), [COND_LABEL[c] for c in CONDS] if fam == 'SO' else [''] * len(CONDS))
    ax.invert_yaxis()
    ax.set_title(FAMILY_LABEL[fam])
    ax.set_xlabel('change of respiratory modulation (%)')
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=2, bbox_to_anchor=(.5, -.08), frameon=False)
fig.suptitle('Poisson GLMM: relative change of the modulation vs clean breathing '
             '(dotted −100 % = complete loss, dashed −50 %;\nintervals are cut at the axis limits)', y=1.06)
fig.savefig(out.glmm)
plt.close(fig)

# ------------------------------------------------------------------ Step 3c: snore-evoked events
E = pd.read_parquet(inp.snore_evoked)
fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
for ax, ev, lab, win in [(axs[0], 'SO_large trough', 'Large SOs (trough)', (.5, 1.5)),
                         (axs[1], 'spindle_cons onset', 'Consensus spindles (onset)', (.25, 1.0))]:
    for ref, col, rl in [('snore burst peak', '#E69F00', 'after snore burst'), ('quiet pseudo-burst', C_CTRL,
                                                                              'same moment of quiet breaths')]:
        m = E[(E.event == ev) & (E.ref == ref)].pivot_table(index='sid', columns='lag', values='rate')
        mu, se = m.mean(), sem(m.to_numpy())
        ax.plot(m.columns, mu, color=col, lw=2, label=f'{rl} (N={len(m)})')
        ax.fill_between(m.columns, mu - se, mu + se, color=col, alpha=.25)
    ax.axvspan(*win, color=C_EXCL, alpha=.12, lw=0, label='excluded window')
    ax.axvline(0, color='k', lw=.8)
    ax.set_title(lab)
    ax.set_xlabel('time from snore-burst peak (s)')
axs[0].set_ylabel('events per s')
axs[0].legend(loc='upper left', fontsize=10)
fig.savefig(out.snore_evoked)
plt.close(fig)

# ------------------------------------------------------------------ Step 4: coupling vs AHI
Pp = participants(cfg)
U = pd.read_csv(inp.subject_coupling)
U = U[U.primary & (U.exclusion == 'none') & (U.condition == 'clean')]
U = U[U.n >= U.family.map({'SO': 20, 'spindle': 20, 'complex': 10})]
T4 = pd.read_csv(inp.step4)
fig, axs = plt.subplots(1, 3, figsize=(14, 4))
for ax, fam in zip(axs, ['SO', 'spindle', 'complex']):
    d = U[U.family == fam].set_index('sid').join(Pp[['AHI']], how='inner')
    x = np.log1p(d.AHI)
    ax.scatter(x, d.pcw, s=18, color=C_EV[fam], alpha=.7)
    b = np.polyfit(x, d.pcw, 1)
    xx = np.linspace(x.min(), x.max(), 10)
    ax.plot(xx, np.polyval(b, xx), 'k', lw=1.5)
    ax.axhline(0, color='grey', lw=.8)
    r = T4[(T4.family == fam) & (T4.marker == 'log_ahi')].iloc[0]
    ax.set_title(f'{FAMILY_LABEL[fam]} (N={int(r.n)})\nβ = {r.beta_rlm:+.3f} per SD, p = {r.p_rlm:.2f}', fontsize=12)
    ticks = [0, 5, 15, 30, 60, 120]
    ax.set_xticks(np.log1p(ticks), ticks)
    ax.set_xlabel('AHI (events/h, log scale)')
axs[0].set_ylabel('clean-N2 coupling (projected)')
fig.suptitle('Coupling in clean breathing vs AHI (robust regression, adjusted for age, sex, BMI)', y=1.04)
fig.savefig(out.severity)
plt.close(fig)

# ------------------------------------------------------------------ Step 5: GLMM rotation and CircStat
GP = pd.read_csv(inp.step5)
rows5 = [('clean_odd', 'none')] + [(c, 'none') for c in CONDS] + [('snoring', 'evoked')]
lab5 = [COND_LABEL['clean_odd']] + [COND_LABEL[c] for c in CONDS] + ['snoring, evoked\nexcluded']
fig, axs = plt.subplots(1, 3, figsize=(14, 4.6), sharex=True, gridspec_kw={'wspace': .35})
for ax, fam in zip(axs, ['SO', 'spindle', 'complex']):
    g = GP[GP.family == fam]
    for j, (cond, excl) in enumerate(rows5):
        r = g[(g.condition == cond) & (g.exclusion == excl)]
        if not len(r) or not np.isfinite(r.glmm_rotation_deg.iloc[0]):
            continue
        r = r.iloc[0]
        sig = r.glmm_orthogonal_p < .05
        col = C_CTRL if cond == 'clean_odd' else (C_EXCL if excl == 'evoked' else C_ALL)
        ax.errorbar(r.glmm_rotation_deg, j, xerr=[[r.glmm_rotation_deg - r.glmm_rot_lo], [r.glmm_rot_hi - r.glmm_rotation_deg]],
                    fmt='o', mfc=col if sig else 'white', color=col, ms=7, capsize=3)
        ax.text(1.02, j, f'p={r.glmm_orthogonal_p:.3f}', va='center', fontsize=10, transform=ax.get_yaxis_transform())
    ax.axvline(0, color='grey', lw=1)
    ax.set_xlim(-180, 180)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_yticks(range(len(rows5)), lab5 if fam == 'SO' else [''] * len(rows5))
    ax.invert_yaxis()
    ax.set_title(FAMILY_LABEL[fam])
    ax.set_xlabel('shift of preferred phase (°)')
fig.suptitle('Preferred-phase shift vs clean breathing, pooled over events (GLMM; filled = orthogonal-shift p < .05)',
             y=1.02)
fig.savefig(out.rotation)
plt.close(fig)

CG = pd.read_csv(inp.circstat)
fig, ax = plt.subplots(figsize=(10, 3.8))
conds5 = ['clean_odd'] + CONDS
w = .25
for i, fam in enumerate(['SO', 'spindle', 'complex']):
    g = CG[(CG.family == fam) & (CG.exclusion == 'none')].set_index('condition')
    v = [100 * g.frac_sig_cm.get(c, np.nan) for c in conds5]
    ax.bar(np.arange(len(conds5)) + (i - 1) * w, v, w, color=C_EV[fam], label=FAMILY_LABEL[fam])
ax.axhline(5, color='k', ls='--', lw=1, label='5 % (chance)')
ax.set_xticks(range(len(conds5)), [COND_LABEL[c] for c in conds5])
ax.set_ylabel('subjects with a changed\nmedian phase (%)')
ax.set_title('Per participant: CircStat common median test, condition vs clean (N2, all events)')
ax.legend(ncol=4, loc='upper left', fontsize=10)
ax.set_ylim(0, 30)
fig.savefig(out.circstat)
plt.close(fig)

# ------------------------------------------------------------------ Step 6: clean preferred phase vs AHI
from snoredreamt import circstat as cs  # noqa: E402

codes = cfg['condition_codes']
stages6 = cfg['step3b']['stage_sets']['N2']
amb = {'FLOW': set(Q.index[Q.flow_ambiguous]), 'ABDOMEN': set(Q.index[Q.abdomen_ambiguous])}
amb['BEST'] = amb['FLOW']
EV = pd.concat(pd.read_parquet(f) for f in inp.events)
EV = EV[EV.stage.isin(stages6) & (EV.label == codes['clean'])]
B6 = pd.read_csv(inp.step6)
fig, axs = plt.subplots(1, 3, figsize=(14, 4), sharey=True)
for ax, fam in zip(axs, ['SO', 'spindle', 'complex']):
    es, var = PRIM[fam]
    g = EV[(EV.event_set == es) & EV[f'phase_{var}'].notna() & ~EV.sid.isin(amb[var.split('_')[0]])]
    cnt = g.groupby('sid').size()
    mu = g.groupby('sid')[f'phase_{var}'].apply(lambda a: cs.circ_mean(a.to_numpy()))[cnt >= {'SO': 20, 'spindle': 20,
                                                                                                 'complex': 10}[fam]]
    d = mu.rename('phase').to_frame().join(Pp[['AHI']], how='inner')
    ax.scatter(np.log1p(d.AHI), np.degrees(d.phase), s=18, color=C_EV[fam], alpha=.7)
    r = B6[(B6.family == fam) & (B6.marker == 'log AHI')].iloc[0]
    ax.set_title(f'{FAMILY_LABEL[fam]} (N={int(r.n)})\nρ = {r.stat:.2f}, p = {r.p:.2f}', fontsize=12)
    ticks = [0, 5, 15, 30, 60, 120]
    ax.set_xticks(np.log1p(ticks), ticks)
    ax.set_xlabel('AHI (events/h, log scale)')
    ax.set_yticks([-180, -90, 0, 90, 180])
axs[0].set_ylabel('clean preferred phase (°)')
fig.suptitle('Clean-breathing preferred phase vs AHI (circular–linear correlation)', y=1.06)
fig.savefig(out.step6)
plt.close(fig)
