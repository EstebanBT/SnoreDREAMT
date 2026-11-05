"""Snakemake rule `step3_group`: group statistics and figures for coupling per breathing condition.

Primary analysis: N2, matched clean controls, ABDOMEN_hilbertgated phase, main arousal exclusion (config).
Per (condition, event): subjects with >= min_events on both sides. Tests:
  - Wilcoxon signed-rank on projected coupling PC (condition vs matched control), median difference with bootstrap
    95 % CI, matched-pairs rank-biserial r; BH-FDR across conditions x events of the primary family;
  - phase shift: one-sample Hotelling T^2 on paired differences of subject mean vectors rotated by mu_ref;
  - coupling present at all: Wilcoxon PC > 0 for control and condition sides;
  - event density ratio (condition / control).
Sensitivity grid: all variants x stage sets x control types x arousal variants.
"""
import matplotlib
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
MIN = cfg['step3']['min_events']
COND = ['clean_split', 'snoring', 'hypopnea', 'recovery', 'obstructive', 'central', 'mixed']
EVENTS = ['SO', 'SO_large', 'spindle', 'complex']
rng = np.random.default_rng(0)

frames = []
MAIN = cfg['arousals']['exclusion']
for arv, files in [(MAIN, inp.main), ('liberal', inp.liberal)]:
    for f in files:
        frames.append(pd.read_parquet(f).assign(arousal_variant=arv))
D = pd.concat(frames)
key = ['arousal_variant', 'stage_set', 'control', 'variant', 'event', 'condition', 'sid']
W = D.pivot_table(index=key, columns='side', values=['pc', 'n', 'density', 'amp_median', 'vec_x', 'vec_y', 'mu_ref'],
                  aggfunc='first')
W.columns = [f'{a}_{b}' for a, b in W.columns]
W = W.reset_index()
W['mu_ref'] = W['mu_ref_condition']


def rotate(x, y, mu):
    c, s = np.cos(-mu), np.sin(-mu)
    return x * c - y * s, x * s + y * c


def hotelling_one_sample(X):
    X = np.asarray(X)
    n, p = X.shape
    if n <= p + 1:
        return np.nan, np.nan
    m = X.mean(0)
    S = np.cov(X, rowvar=False)
    t2 = n * m @ np.linalg.solve(S, m)
    F = (n - p) / (p * (n - 1)) * t2
    return t2, 1 - stats.f.cdf(F, p, n - p)


def rank_biserial(d):
    d = d[d != 0]
    if len(d) == 0:
        return np.nan
    r = stats.rankdata(np.abs(d))
    return (r[d > 0].sum() - r[d < 0].sum()) / r.sum()


def boot_median_ci(d, n_boot=2000):
    if len(d) < 3:
        return np.nan, np.nan
    b = np.median(rng.choice(d, (n_boot, len(d)), replace=True), axis=1)
    return np.percentile(b, [2.5, 97.5])


# ---------------------------------------------------------------- group reference direction (leave-one-subject-out)
# Subject-level preferred phases are unreliable (split-half: SO 30 % within 45 deg vs 25 % chance), so a second
# metric projects each subject's events onto the group preferred phase of Step 2 computed without that subject
# (V-test logic, Schreiner 2023); only subjects with unambiguous polarity for the phase channel enter it.
S2 = pd.concat(pd.read_parquet(f) for f in inp.step2)
Q = pd.read_csv(inp.resp_summary).set_index('sid')
amb = {'FLOW': set(Q.index[Q.flow_ambiguous]), 'ABDOMEN': set(Q.index[Q.abdomen_ambiguous])}
S2 = S2[S2.n >= cfg['coupling']['min_events']]
S2 = S2[[sid not in amb[v.split('_')[0]] for sid, v in zip(S2.sid, S2.variant)]]
vec = S2.assign(c=np.cos(S2['mean']), s=np.sin(S2['mean'])).groupby(['event', 'variant'])[['c', 's']].sum()
cnt = S2.groupby(['event', 'variant']).size()


def mu_loso(ev, var, sid):
    ev = {'SO_large': 'SO'}.get(ev, ev)          # large-SO subset uses the SO group phase
    if (ev, var) not in vec.index:
        return np.nan
    c, s_ = vec.loc[(ev, var)]
    me = S2[(S2.event == ev) & (S2.variant == var) & (S2.sid == sid)]
    if len(me):
        c, s_ = c - np.cos(me['mean'].iloc[0]), s_ - np.sin(me['mean'].iloc[0])
    return np.arctan2(s_, c)


W['mu_group'] = [mu_loso(e, v, sd) for e, v, sd in zip(W.event, W.variant, W.sid)]
W['polarity_ok'] = [sd not in amb[v.split('_')[0]] for sd, v in zip(W.sid, W.variant)]
for side in ['condition', 'control']:
    W[f'pcg_{side}'] = W[f'vec_x_{side}'] * np.cos(W.mu_group) + W[f'vec_y_{side}'] * np.sin(W.mu_group)

KEYS = ['arousal_variant', 'stage_set', 'control', 'variant', 'event', 'condition']
rows = []
for metric in ['pc', 'pcg']:
    for k, g in W.groupby(KEYS):
        g = g[(g.n_condition >= MIN) & (g.n_control >= MIN)].dropna(subset=[f'{metric}_condition', f'{metric}_control'])
        if metric == 'pcg':
            g = g[g.polarity_ok]
        base = dict(zip(KEYS, k), metric=metric, n_subjects=len(g))
        if len(g) < 3:
            rows.append(base)
            continue
        dpc = (g[f'{metric}_condition'] - g[f'{metric}_control']).to_numpy()
        lo, hi = boot_median_ci(dpc)
        mu = g.mu_ref if metric == 'pc' else g.mu_group
        xc, yc = rotate(g.vec_x_condition, g.vec_y_condition, mu)
        xm, ym = rotate(g.vec_x_control, g.vec_y_control, mu)
        t2, p_h = hotelling_one_sample(np.c_[xc - xm, yc - ym])
        shift = np.degrees(np.angle(np.mean(xc + 1j * yc)) - np.angle(np.mean(xm + 1j * ym)))
        rows.append(dict(base, median_events_condition=g.n_condition.median(), median_events_control=g.n_control.median(),
                         control_median=g[f'{metric}_control'].median(), condition_median=g[f'{metric}_condition'].median(),
                         delta_median=np.median(dpc), delta_ci_lo=lo, delta_ci_hi=hi,
                         wilcoxon_p=stats.wilcoxon(dpc).pvalue, rank_biserial=rank_biserial(dpc),
                         p_control_gt0=stats.wilcoxon(g[f'{metric}_control'], alternative='greater').pvalue,
                         p_condition_gt0=stats.wilcoxon(g[f'{metric}_condition'], alternative='greater').pvalue,
                         phase_shift_deg=(shift + 180) % 360 - 180, hotelling_p=p_h,
                         density_ratio_median=np.nanmedian(g.density_condition / g.density_control),
                         density_wilcoxon_p=stats.wilcoxon(g.density_condition - g.density_control).pvalue,
                         amp_ratio_median=np.nanmedian(g.amp_median_condition / g.amp_median_control)))
T = pd.DataFrame(rows)
prim = (T.arousal_variant == MAIN) & (T.stage_set == 'N2') & (T.control == 'matched') & \
       (T.variant == 'ABDOMEN_hilbertgated') & (T.condition != 'clean_split') & T.wilcoxon_p.notna()
T['primary_family'] = prim
T['wilcoxon_p_fdr'] = np.nan
for m in ['pc', 'pcg']:
    sel = prim & (T.metric == m)
    T.loc[sel, 'wilcoxon_p_fdr'] = multipletests(T.loc[sel, 'wilcoxon_p'], method='fdr_bh')[1]

# ---------------------------------------------------------------- event-level model (group reference)
import statsmodels.formula.api as smf  # noqa: E402
E = pd.concat(pd.read_parquet(f) for f in inp.events_main)
E = E[E.variant == 'ABDOMEN_hilbertgated']
E['polarity_ok'] = ~E.sid.isin(amb['ABDOMEN'])
E = E[E.polarity_ok]
E['mu_group'] = [mu_loso(e, 'ABDOMEN_hilbertgated', sd) for e, sd in zip(E.event, E.sid)]
E['y'] = np.cos(E.phase - E.mu_group)
E['is_condition'] = (E.side == 'condition').astype(int)
mm = []
for (c, ev), g in E.groupby(['condition', 'event']):
    g = g.dropna(subset=['y'])
    if g.sid.nunique() < 5:
        continue
    try:
        # event-level OLS with subject-clustered SEs (MixedLM collapsed the intercept to 0 with p = 1 on these data)
        r = smf.ols('y ~ is_condition', g).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(g['sid'])[0]})
        mm.append(dict(condition=c, event=ev, n_subjects=g.sid.nunique(), n_events=len(g),
                       intercept_control=r.params['Intercept'], beta_condition=r.params['is_condition'],
                       ci_lo=r.conf_int().loc['is_condition', 0], ci_hi=r.conf_int().loc['is_condition', 1],
                       p=r.pvalues['is_condition'], p_control_intercept=r.pvalues['Intercept']))
    except Exception as e:  # report and keep going
        mm.append(dict(condition=c, event=ev, error=repr(e)))
MM = pd.DataFrame(mm)
MM.to_csv(out.mixed, index=False)
T.to_csv(out.table, index=False)
W.to_parquet(out.wide)

# ---------------------------------------------------------------- figures
P = T[(T.arousal_variant == MAIN) & (T.stage_set == 'N2') & (T.control == 'matched') &
      (T.variant == 'ABDOMEN_hilbertgated')]
MLAB = {'pc': 'subject reference (mu_ref, own clean data)', 'pcg': 'group reference (LOSO Step 2 phase, polarity-ok)'}

fig, axs = plt.subplots(2, 4, figsize=(21, 9), sharey=True)
for i, metric in enumerate(['pcg', 'pc']):
    for j, ev in enumerate(EVENTS):
        ax = axs[i, j]
        p = P[(P.event == ev) & (P.metric == metric)].set_index('condition').reindex(COND)
        y = np.arange(len(COND))
        ok = p.delta_median.notna()
        ax.errorbar(p.delta_median[ok], y[ok.to_numpy()],
                    xerr=[(p.delta_median - p.delta_ci_lo)[ok], (p.delta_ci_hi - p.delta_median)[ok]],
                    fmt='o', color='k', capsize=3)
        ax.axvline(0, color='grey', lw=.7)
        for yi, (c, r) in zip(y, p.iterrows()):
            nn = int(r.n_subjects) if np.isfinite(r.n_subjects) else 0
            txt = f'N={nn} p={r.wilcoxon_p:.2g}' if np.isfinite(r.get('wilcoxon_p', np.nan)) else f'N={nn} (too few)'
            ax.text(0.02, yi + 0.3, txt, fontsize=7, transform=ax.get_yaxis_transform(),
                    color='k' if np.isfinite(r.get('wilcoxon_p', np.nan)) else 'grey')
        ax.set_yticks(y)
        ax.set_yticklabels(COND)
        ax.set_title(f'{ev} — {MLAB[metric]}', fontsize=9)
        ax.set_xlabel('Δ projected coupling (condition − matched clean), median ± boot 95% CI', fontsize=8)
axs[0, 0].invert_yaxis()
fig.suptitle('Step 3: N2, ABDOMEN gated Hilbert (±0.05 Hz), matched clean controls (±10 min, same stage); '
             'clean_split = clean-vs-clean control', fontsize=10)
fig.tight_layout()
fig.savefig(out.forest, dpi=80)

G = W[(W.arousal_variant == MAIN) & (W.stage_set == 'N2') & (W.control == 'matched') &
      (W.variant == 'ABDOMEN_hilbertgated') & (W.n_condition >= MIN) & (W.n_control >= MIN) & W.polarity_ok]
fig, axs = plt.subplots(len(EVENTS), len(COND), figsize=(20, 9), sharey='row')
for i, ev in enumerate(EVENTS):
    for j, c in enumerate(COND):
        g = G[(G.event == ev) & (G.condition == c)]
        ax = axs[i, j]
        for _, r in g.iterrows():
            ax.plot([0, 1], [r.pcg_control, r.pcg_condition], color='grey', alpha=.4, lw=.7)
        if len(g):
            ax.plot([0, 1], [g.pcg_control.median(), g.pcg_condition.median()], 'r-o', lw=2)
        ax.axhline(0, color='k', lw=.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['matched\nclean', c], fontsize=7)
        ax.set_title(f'{ev} | {c} (N={len(g)})', fontsize=8)
    axs[i, 0].set_ylabel('coupling projected on group phase')
fig.suptitle('Per-subject coupling (projection on LOSO group phase): matched clean vs condition (N2, ABDOMEN gated); '
             'red = median', fontsize=10)
fig.tight_layout()
fig.savefig(out.paired, dpi=80)

S = T[T.event.isin(['SO', 'SO_large', 'spindle']) & T.wilcoxon_p.notna()].copy()
S['label'] = S.metric + '|' + S.arousal_variant.str[:3] + '|' + S.stage_set + '|' + S.control + '|' + S.variant
fig, axs = plt.subplots(1, 3, figsize=(22, 14), sharey=True)
for ax, ev in zip(axs, ['SO', 'SO_large', 'spindle']):
    m = S[S.event == ev].pivot(index='label', columns='condition', values='delta_median').reindex(columns=COND)
    pv = S[S.event == ev].pivot(index='label', columns='condition', values='wilcoxon_p').reindex(columns=COND)
    im = ax.imshow(m.to_numpy(), cmap='RdBu_r', vmin=-0.15, vmax=0.15, aspect='auto')
    for (yi, xi), v in np.ndenumerate(pv.to_numpy()):
        if np.isfinite(v) and v < 0.05:
            ax.text(xi, yi, '*' if v >= 0.01 else '**', ha='center', va='center', fontsize=9)
    ax.set_xticks(range(len(COND)))
    ax.set_xticklabels(COND, rotation=45, fontsize=8)
    ax.set_yticks(range(len(m)))
    ax.set_yticklabels(m.index, fontsize=6)
    ax.set_title(ev)
fig.colorbar(im, ax=axs, label='median Δ coupling (condition − control)')
fig.suptitle('Sensitivity grid: metric | arousal detector | stages | control | phase variant '
             '(* p<.05, ** p<.01, uncorrected)', fontsize=10)
fig.savefig(out.grid, dpi=80, bbox_inches='tight')

cols = ['metric', 'event', 'condition', 'n_subjects', 'control_median', 'condition_median', 'delta_median',
        'delta_ci_lo', 'delta_ci_hi', 'wilcoxon_p', 'wilcoxon_p_fdr', 'p_control_gt0', 'phase_shift_deg', 'hotelling_p',
        'density_ratio_median', 'density_wilcoxon_p']
with pd.option_context('display.width', 280, 'display.max_columns', 30):
    print(P[P.event != 'complex'][cols].round(3).to_string())
    print(MM.round(4).to_string())
