"""Snakemake rule `step3b_group`: power-improved Step 3 (ideas 2 + 3) and data export for the GLMM (idea 4).

A. Clean-N2 signal per event set x phase variant (selection uses clean data only, not the condition contrasts):
   per-event y = cos(phase - mu_LOSO), mu_LOSO = leave-one-subject-out group phase of subjects' clean-N2 mean phases
   (polarity-unambiguous subjects for the phase channel; BEST follows FLOW polarity). Group signal z = mean/SEM of
   subject means. The primary (event set, variant) maximises z.
B. Contrasts (N2): pooled respiratory events and single conditions vs whole-night clean N2; per-subject mean y,
   one-sided Wilcoxon for loss (condition < clean) plus two-sided; event-level OLS with subject-clustered SEs.
C. Dose-response (N2 events in clean, snoring, events, recovery): y ~ snore dB + airflow log2 + effort log2 of the
   breath containing the event, all centred within subject and scaled by the pooled within-subject SD;
   subject-clustered OLS; loss hypotheses: louder snoring -> lower y, lower airflow -> lower y.
"""
import matplotlib
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt.coupling import circ_mean  # noqa: E402

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
codes = cfg['condition_codes']
P = cfg['step3b']
SETS = ['SO', 'SO_large', 'SO_cons', 'SO_large_cons', 'spindle', 'spindle_cons']
VARS = P['phase_variants']
Q = pd.read_csv(inp.resp_summary).set_index('sid')
amb = {'FLOW': set(Q.index[Q.flow_ambiguous]), 'ABDOMEN': set(Q.index[Q.abdomen_ambiguous])}
amb['BEST'] = amb['FLOW']

E = pd.concat(pd.read_parquet(f) for f in inp.events)
E = E[E.stage == 2]                                    # N2 primary
inv = {v: k for k, v in codes.items()}
E['condition'] = E.label.map(inv)
E.loc[E.label.isin(P['pooled_event_codes']), 'pooled'] = 'resp_event'

# ------------------------------------------------------------------ A. clean signal and LOSO group phase
mu = {}                                                 # (set, var) -> Series sid -> mu_LOSO
rowsA = []
for s in SETS:
    for v in VARS:
        col = f'phase_{v}'
        c = E[(E.event_set == s) & (E.label == codes['clean'])].dropna(subset=[col])
        c = c[~c.sid.isin(amb[v.split('_')[0]])]
        cnt = c.groupby('sid').size()
        m = c.groupby('sid')[col].apply(lambda a: circ_mean(a.to_numpy()))[cnt >= 20]
        cs, sn = np.cos(m).sum(), np.sin(m).sum()
        mu[(s, v)] = pd.Series({sid: np.arctan2(sn - np.sin(m[sid]), cs - np.cos(m[sid])) for sid in m.index})
        c = c[c.sid.isin(m.index)]
        y = np.cos(c[col] - c.sid.map(mu[(s, v)]))
        sm = y.groupby(c.sid).mean()
        rowsA.append(dict(event_set=s, variant=v, n_subjects=len(sm), median_events=cnt[m.index].median(),
                          group_phase_deg=np.degrees(np.arctan2(sn, cs)), mean_y=sm.mean(), sem_y=sm.sem(),
                          z=sm.mean() / sm.sem(), p_one_sided=stats.wilcoxon(sm, alternative='greater').pvalue,
                          event_level_mean=y.mean()))
A = pd.DataFrame(rowsA).sort_values('z', ascending=False)
A['family'] = np.where(A.event_set.str.startswith('SO'), 'SO', 'spindle')
# primary per family: max clean-N2 z among phase sources that also cover apneas (FLOW gated is undefined there)
cover = A[A.variant.isin(['BEST_hilbertgated', 'ABDOMEN_hilbertgated'])]
PRIMARIES = {fam: tuple(g.iloc[0][['event_set', 'variant']]) for fam, g in cover.groupby('family')}
A['primary'] = [(s_, v_) in PRIMARIES.values() for s_, v_ in zip(A.event_set, A.variant)]
A.to_csv(out.clean_signal, index=False)
print('Clean-N2 coupling signal per event set x phase variant:\n', A.round(4).to_string(index=False))
print('PRIMARIES (max clean z per family, apnea-covering phase):', PRIMARIES)
PRIMARY = PRIMARIES['SO']

# per-event y for every set/variant
for v in VARS:
    col = f'phase_{v}'
    E[f'y_{v}'] = np.nan
    for s in SETS:
        sel = (E.event_set == s) & ~E.sid.isin(amb[v.split('_')[0]])
        E.loc[sel, f'y_{v}'] = np.cos(E.loc[sel, col] - E.loc[sel, 'sid'].map(mu[(s, v)]))

# ------------------------------------------------------------------ B. contrasts
CONDS = ['resp_event', 'snoring', 'hypopnea', 'recovery', 'obstructive', 'central', 'mixed']
rowsB = []
for s in SETS:
    for v in VARS:
        yc = f'y_{v}'
        g0 = E[(E.event_set == s) & E[yc].notna()]
        clean = g0[g0.label == codes['clean']]
        cm = clean.groupby('sid')[yc].agg(['mean', 'size'])
        cm = cm[cm['size'] >= 20]
        for cond in CONDS:
            gc = g0[g0.pooled == 'resp_event'] if cond == 'resp_event' else g0[g0.condition == cond]
            km = gc.groupby('sid')[yc].agg(['mean', 'size'])
            km = km[km['size'] >= P['min_events']]
            sids = km.index.intersection(cm.index)
            base = dict(event_set=s, variant=v, condition=cond, n_subjects=len(sids))
            if len(sids) < 5:
                rowsB.append(base)
                continue
            dlt = km.loc[sids, 'mean'] - cm.loc[sids, 'mean']
            ev = pd.concat([clean[clean.sid.isin(sids)].assign(x=0), gc[gc.sid.isin(sids)].assign(x=1)])
            r = smf.ols(f'{yc} ~ x', ev).fit(cov_type='cluster', cov_kwds={'groups': pd.factorize(ev.sid)[0]})
            b, se = r.params['x'], r.bse['x']
            rowsB.append(dict(base, median_events_condition=km.loc[sids, 'size'].median(),
                              clean_mean=cm.loc[sids, 'mean'].mean(), condition_mean=km.loc[sids, 'mean'].mean(),
                              delta_mean=dlt.mean(), delta_sem=dlt.sem(),
                              p_loss_wilcoxon=stats.wilcoxon(dlt, alternative='less').pvalue,
                              p_two_sided_wilcoxon=stats.wilcoxon(dlt).pvalue,
                              ev_beta=b, ev_ci_lo=b - 1.96 * se, ev_ci_hi=b + 1.96 * se,
                              p_loss_event=stats.norm.cdf(b / se), p_two_sided_event=r.pvalues['x'],
                              n_events_condition=int((ev.x == 1).sum()), n_events_clean=int((ev.x == 0).sum())))
B = pd.DataFrame(rowsB)
B['primary'] = [(s_, v_) in PRIMARIES.values() for s_, v_ in zip(B.event_set, B.variant)]
B.to_csv(out.contrasts, index=False)
for fam, pr in PRIMARIES.items():
    print(f'\nContrasts, primary {fam}: {pr}\n',
          B[(B.event_set == pr[0]) & (B.variant == pr[1])].drop(columns=['primary', 'event_set', 'variant'])
          .round(4).to_string(index=False))

# ------------------------------------------------------------------ C. dose-response
rowsC, figdata = [], {}
for s in SETS:
    for v in VARS:
        yc = f'y_{v}'
        g = E[(E.event_set == s) & E[yc].notna() & E.label.isin([0, 1, 2, 3, 4, 5, 6])].copy()
        for pr in ['snore_db', 'flow_log2', 'effort_log2']:
            g[pr + '_c'] = g[pr] - g.groupby('sid')[pr].transform('mean')
            g[pr + '_c'] /= g[pr + '_c'].std()
        g = g.dropna(subset=['snore_db_c', 'flow_log2_c', 'effort_log2_c'])
        g = g[np.isfinite(g[['snore_db_c', 'flow_log2_c', 'effort_log2_c']]).all(1)]
        if g.sid.nunique() < 10:
            continue
        r = smf.ols(f'{yc} ~ snore_db_c + flow_log2_c + effort_log2_c', g).fit(
            cov_type='cluster', cov_kwds={'groups': pd.factorize(g.sid)[0]})
        row = dict(event_set=s, variant=v, n_subjects=g.sid.nunique(), n_events=len(g), intercept=r.params['Intercept'])
        for pr, loss_sign in [('snore_db_c', -1), ('flow_log2_c', +1), ('effort_log2_c', +1)]:
            b, se = r.params[pr], r.bse[pr]
            # one-sided p for the loss hypothesis: snoring louder -> beta < 0; airflow/effort lower -> beta > 0
            p_loss = stats.norm.cdf(b / se) if loss_sign < 0 else 1 - stats.norm.cdf(b / se)
            row.update({f'{pr}_beta': b, f'{pr}_se': se, f'{pr}_p_two_sided': r.pvalues[pr], f'{pr}_p_loss': p_loss})
        rowsC.append(row)
        if (s, v) == PRIMARY:
            figdata = g
C = pd.DataFrame(rowsC)
C['primary'] = [(s_, v_) in PRIMARIES.values() for s_, v_ in zip(C.event_set, C.variant)]
C.to_csv(out.dose, index=False)
print('\nDose-response (per within-subject SD):\n', C.round(4).to_string(index=False))

# ------------------------------------------------------------------ export phase-bin data for the GLMM (idea 4)
Bn = pd.concat(pd.read_parquet(f) for f in inp.bins)
Bn = Bn[(Bn.stage_set == 'N2')]
Bn['polarity_ok'] = [sd not in amb[v.split('_')[0]] for sd, v in zip(Bn.sid, Bn.variant)]
gph = {(r.event_set, r.variant): np.radians(r.group_phase_deg) for r in A.itertuples()}
Bn['mu_group'] = [gph.get((s, v), np.nan) for s, v in zip(Bn.event_set, Bn.variant)]
Bn['primary'] = [(s_, v_) in PRIMARIES.values() for s_, v_ in zip(Bn.event_set, Bn.variant)]
Bn.to_csv(out.bins_export, index=False)

# ------------------------------------------------------------------ figures
fig, axs = plt.subplots(1, 3, figsize=(17, 5))
ax = axs[0]
A2 = A.set_index(['event_set', 'variant'])
for k, v in enumerate(VARS):
    sub = A2.xs(v, level='variant').reindex(SETS)
    ax.errorbar(np.arange(len(SETS)) + (k - 1) * .2, sub.mean_y, 1.96 * sub.sem_y, fmt='o', capsize=3, label=v)
ax.set_xticks(range(len(SETS)))
ax.set_xticklabels(SETS, rotation=30, fontsize=8)
ax.axhline(0, color='grey', lw=.6)
ax.set_ylabel('clean-N2 coupling: mean cos(phase − group phase)\nsubject means ± 95% CI')
ax.set_title(f'A. signal per event set / phase source (primary: {PRIMARY[0]}, {PRIMARY[1]})', fontsize=9)
ax.legend(fontsize=7)
ax = axs[1]
y = np.arange(len(CONDS))
for k, (fam, pr) in enumerate(PRIMARIES.items()):
    pb = B[(B.event_set == pr[0]) & (B.variant == pr[1])].set_index('condition').reindex(CONDS)
    yy = y + (k - .5) * .3
    col = 'k' if fam == 'SO' else 'tab:green'
    ax.errorbar(pb.ev_beta, yy, xerr=[pb.ev_beta - pb.ev_ci_lo, pb.ev_ci_hi - pb.ev_beta], fmt='o', color=col,
                capsize=3, label=f'{pr[0]} | {pr[1]}')
    for yi, (cn, r) in zip(yy, pb.iterrows()):
        if np.isfinite(r.get('ev_beta', np.nan)):
            ax.text(0.01 + .5 * k, yi + .12, f'N={int(r.n_subjects)} p_loss={r.p_loss_event:.2g}', fontsize=6,
                    color=col, transform=ax.get_yaxis_transform())
    ax.axvline(-A2.loc[pr, 'mean_y'], color=col, ls='--', lw=.8)
ax.axvline(0, color='grey', lw=.7)
ax.set_yticks(y)
ax.set_yticklabels(CONDS)
ax.invert_yaxis()
ax.legend(fontsize=7)
ax.set_xlabel('condition − clean (event-level, subject-clustered 95% CI)')
ax.set_title('B. condition − clean (N2); dashed = complete loss of clean coupling', fontsize=9)
ax = axs[2]
if len(figdata):
    yc = f'y_{PRIMARY[1]}'
    for pr, c in [('snore_db_c', 'tab:orange'), ('flow_log2_c', 'tab:blue'), ('effort_log2_c', 'tab:red')]:
        q = pd.qcut(figdata[pr], 8, duplicates='drop')
        m = figdata.groupby(q, observed=True)[yc].agg(['mean', 'sem'])
        xc = figdata.groupby(q, observed=True)[pr].mean()
        ax.errorbar(xc, m['mean'], m['sem'], fmt='-o', color=c, capsize=2, label=pr.replace('_c', ''))
    ax.axhline(0, color='grey', lw=.6)
    ax.set_xlabel('breath predictor (within-subject z)')
    ax.set_ylabel('mean cos(phase − group phase)')
    ax.legend(fontsize=7)
    ax.set_title('C. dose-response (primary set, octiles)', fontsize=9)
fig.tight_layout()
fig.savefig(out.fig, dpi=80)
