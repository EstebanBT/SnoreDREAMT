"""Snakemake rule `step5_circstat`: Step 5 with CircStat tests on event phases (Esteban's request, 2026-10-26).

5a: per subject, compare the respiratory phases of events during each condition with those during clean N2 using
    circ_wwtest (Watson-Williams, equal mean directions; parametric, assumes von Mises with sufficient concentration
    -> CircStat's applicability rule is recorded) and circ_cmtest (common median test; non-parametric). Count subjects
    with p < .05 and compare with the 5 % expected under H0 (binomial) and with the clean-vs-clean control
    (alternating 10-min blocks of clean N2: clean_odd vs clean_even).
5b: per subject mean phase difference d_mu = circ_dist(circ_mean(cond), circ_mean(clean)) (and median difference);
    population: does the distribution of d_mu differ from 0?  circ_mtest(d_mu, 0) (0 outside the 95 % CI of the mean),
    circ_medtest(d_mu, 0) (sign test for median 0), circ_rtest for non-uniformity (context only).
Event phases are used directly (Hilbert gated phase occupancy is ~uniform, R ~ 0.01). With snore-evoked exclusion,
the same breath-phase window is removed in both groups (phase-matched).
"""
import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt import circstat as cs  # noqa: E402
from snoredreamt.coupling import rayleigh  # noqa: E402

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
fs = cfg['fs']
P = cfg['step3b']
SS = snakemake.wildcards.ss  # noqa: F821
STAGES = P['stage_sets'][SS]
MIN_N = 10                                              # CircStat circ_cmtest requirement per group
codes = cfg['condition_codes']
COND = {'resp_event': P['pooled_event_codes'], 'snoring': [codes['snoring']], 'hypopnea': [codes['hypopnea']],
        'recovery': [codes['recovery']], 'obstructive': [codes['obstructive']], 'central': [codes['central']],
        'mixed': [codes['mixed']]}

A = pd.read_csv(inp.clean_signal)
PRIM = {r.family: (r.event_set, r.variant) for r in A[A.primary].itertuples()}
E = pd.concat(pd.read_parquet(f) for f in inp.events)
E = E[E.stage.isin(STAGES) & E.event_set.isin([p[0] for p in PRIM.values()])]
E['block'] = (E['sample'] // (600 * fs)) % 2

rows = []
for fam, (es, var) in PRIM.items():
    col = f'phase_{var}'
    for excl in ['none', 'evoked']:
        g0 = E[(E.event_set == es) & E[col].notna()]
        if excl == 'evoked':
            g0 = g0[~g0.evoked_window]
        for sid, gs in g0.groupby('sid'):
            clean = gs[gs.label == codes['clean']]
            pairs = {c: (gs[gs.label.isin(cc)][col].to_numpy(), clean[col].to_numpy()) for c, cc in COND.items()}
            pairs['clean_odd'] = (clean[clean.block == 1][col].to_numpy(), clean[clean.block == 0][col].to_numpy())
            for cond, (a1, a0) in pairs.items():
                # snore-evoked exclusion removes a breath-phase window from snoring and quiet-clean breaths only, so
                # it is phase-matched (valid) for snoring vs clean and clean vs clean, but would bias other contrasts
                if excl == 'evoked' and cond not in ('snoring', 'clean_odd'):
                    continue
                if len(a1) < MIN_N or len(a0) < MIN_N:
                    continue
                p_ww, F, ww_ok = cs.circ_wwtest(a1, a0)
                p_cm = cs.circ_cmtest(a1, a0)[0]
                rows.append(dict(family=fam, event_set=es, variant=var, exclusion=excl, condition=cond, sid=sid,
                                 n_cond=len(a1), n_ref=len(a0), r_cond=cs.circ_r(a1), r_ref=cs.circ_r(a0),
                                 p_ww=p_ww, ww_applicable=ww_ok, p_cm=p_cm,
                                 d_mean=cs.circ_dist(cs.circ_mean(a1), cs.circ_mean(a0)),
                                 d_median=cs.circ_dist(cs.circ_median(a1), cs.circ_median(a0))))
S = pd.DataFrame(rows)
S.to_csv(out.subjects, index=False)


def summarise(g):
    n = len(g)
    res = {'n_subjects': n, 'median_n_cond': g.n_cond.median(), 'frac_ww_applicable': g.ww_applicable.mean()}
    for t in ['ww', 'cm']:
        k = int((g[f'p_{t}'] < .05).sum())
        res[f'n_sig_{t}'] = k
        res[f'frac_sig_{t}'] = k / n
        res[f'p_binom_{t}'] = stats.binomtest(k, n, .05, alternative='greater').pvalue
    for d in ['d_mean', 'd_median']:
        a = g[d].to_numpy()
        for k in ['pop_mean_deg', 'ci_lo_deg', 'ci_hi_deg', 'mtest_h_differs_from_0', 'medtest_p', 'pop_median_deg',
                  'rayleigh_p', 'R']:
            res[f'{d}_{k}'] = np.nan
        if n >= 5:
            h, mu, ul, ll = cs.circ_mtest(a, 0.0)
            res.update({f'{d}_pop_mean_deg': np.degrees(mu), f'{d}_ci_lo_deg': np.degrees(ll),
                        f'{d}_ci_hi_deg': np.degrees(ul), f'{d}_mtest_h_differs_from_0': h,
                        f'{d}_medtest_p': cs.circ_medtest(a, 0.0), f'{d}_pop_median_deg': np.degrees(cs.circ_median(a)),
                        f'{d}_rayleigh_p': rayleigh(a)[0], f'{d}_R': cs.circ_r(a)})
    return pd.Series(res)


import warnings  # noqa: E402
warnings.filterwarnings('ignore', message='Ties detected')
T = S.groupby(['family', 'event_set', 'variant', 'exclusion', 'condition'])[S.columns.drop(
    ['family', 'event_set', 'variant', 'exclusion', 'condition'])].apply(summarise).reset_index()
# compare the fraction of significant subjects with the clean-vs-clean control (Fisher exact, one-sided)
for t in ['ww', 'cm']:
    ctrl = T[T.condition == 'clean_odd'].set_index(['family', 'exclusion'])
    pv = []
    for r in T.itertuples():
        c = ctrl.loc[(r.family, r.exclusion)]
        tab = [[getattr(r, f'n_sig_{t}'), r.n_subjects - getattr(r, f'n_sig_{t}')],
               [c[f'n_sig_{t}'], c.n_subjects - c[f'n_sig_{t}']]]
        pv.append(stats.fisher_exact(tab, alternative='greater')[1] if r.condition != 'clean_odd' else np.nan)
    T[f'p_vs_control_{t}'] = pv
T.to_csv(out.table, index=False)
order = ['clean_odd', 'resp_event', 'snoring', 'hypopnea', 'recovery', 'obstructive', 'central', 'mixed']
T['o'] = T.condition.map({c: i for i, c in enumerate(order)})
with pd.option_context('display.width', 280, 'display.max_columns', 40):
    print(T.sort_values(['family', 'exclusion', 'o'])[
        ['family', 'exclusion', 'condition', 'n_subjects', 'median_n_cond', 'frac_ww_applicable', 'n_sig_ww', 'frac_sig_ww',
         'p_binom_ww', 'p_vs_control_ww', 'n_sig_cm', 'frac_sig_cm', 'p_binom_cm', 'p_vs_control_cm',
         'd_mean_pop_mean_deg', 'd_mean_ci_lo_deg', 'd_mean_ci_hi_deg', 'd_mean_mtest_h_differs_from_0',
         'd_mean_medtest_p', 'd_mean_R']].round(3).to_string(index=False))

# ------------------------------------------------------------------ figure
fams = list(PRIM)
fig, axs = plt.subplots(2, len(fams), figsize=(6 * len(fams), 9))
for j, fam in enumerate(fams):
    ax = axs[0, j]
    t = T[(T.family == fam) & (T.exclusion == 'none')].set_index('condition').reindex(order).dropna(subset=['n_subjects'])
    x = np.arange(len(t))
    ax.bar(x - .2, t.frac_sig_ww, .4, label='Watson-Williams')
    ax.bar(x + .2, t.frac_sig_cm, .4, label='common median')
    ax.axhline(.05, color='r', ls='--', lw=.8, label='5 % (H0)')
    ax.set_xticks(x)
    ax.set_xticklabels([f'{c}\nN={int(n)}' for c, n in zip(t.index, t.n_subjects)], fontsize=7, rotation=30)
    ax.set_ylabel('fraction of subjects with p < .05')
    ax.set_title(f'{fam}: {PRIM[fam][0]} | {PRIM[fam][1]} (all events)', fontsize=9)
    ax.legend(fontsize=7)
    ax = fig.add_subplot(2, len(fams), len(fams) + j + 1, projection='polar')
    axs[1, j].set_visible(False)
    for cond, c in [('resp_event', 'tab:blue'), ('snoring', 'tab:orange'), ('clean_odd', 'grey')]:
        g = S[(S.family == fam) & (S.exclusion == 'none') & (S.condition == cond)]
        h, e = np.histogram(g.d_mean, np.linspace(-np.pi, np.pi, 19))
        ax.plot(np.r_[e[:-1], e[0]] + np.pi / 18, np.r_[h, h[0]] / max(len(g), 1), color=c, label=f'{cond} (N={len(g)})')
    ax.set_title('per-subject mean phase difference (condition − clean)', fontsize=8)
    ax.legend(fontsize=6, loc='lower left')
fig.tight_layout()
fig.savefig(out.fig, dpi=80)
