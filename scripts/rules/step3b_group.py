"""Snakemake rule `step3b_group`: power-improved Step 3 (ideas 2 + 3), snore-evoked exclusion, data export for GLMM.

A. Clean-N2 signal per event set x phase variant (selection on clean data only): per-event y = cos(phase - mu_LOSO),
   mu_LOSO = leave-one-subject-out group phase of subjects' clean-N2 mean phases (polarity-unambiguous subjects; BEST
   follows FLOW polarity). Group z = mean/SEM of subject means. Primary per family (SO, spindle, complex) = max z
   among phase sources that cover apneas (BEST, ABDOMEN).
B. Contrasts (N2): pooled respiratory events and single conditions vs whole-night clean N2; per-subject mean y,
   one-sided Wilcoxon for loss plus two-sided; event-level OLS with subject-clustered SEs. Run without and with the
   phase-matched exclusion of snore-evoked events (events in the family window after snore-burst peaks in snoring
   breaths and after pseudo-bursts at the same breath position in quiet clean breaths).
C. Dose-response (N2): y ~ snore dB + airflow log2 + effort log2 of the containing breath, within-subject centred and
   scaled; subject-clustered OLS; both exclusion levels.
D. Snoring diagnostics for the primaries: per-subject vectors rotated by mu_LOSO (phase shift), event densities.
"""
import matplotlib
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt.coupling import circ_mean, rayleigh, vtest  # noqa: E402

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
codes = cfg['condition_codes']
P = cfg['step3b']
SETS = ['SO', 'SO_large', 'SO_cons', 'SO_large_cons', 'spindle', 'spindle_cons', 'complex', 'complex_large',
        'complex_cons', 'complex_pooled', 'complex_pooled_large']
SS = snakemake.wildcards.ss  # noqa: F821  stage set: N2 or N2N3
STAGES = P['stage_sets'][SS]
VARS = P['phase_variants']
EXCL = ['none', 'evoked']
CONDS = ['resp_event', 'snoring', 'hypopnea', 'recovery', 'obstructive', 'central', 'mixed']


def family(s):
    return 'complex' if s.startswith('complex') else ('SO' if s.startswith('SO') else 'spindle')


MIN_MU = {'SO': 20, 'spindle': 20, 'complex': 10}     # clean events needed for a subject's mean phase
MIN_CLEAN = {'SO': 20, 'spindle': 20, 'complex': 10}  # clean events needed to enter a contrast
Q = pd.read_csv(inp.resp_summary).set_index('sid')
amb = {'FLOW': set(Q.index[Q.flow_ambiguous]), 'ABDOMEN': set(Q.index[Q.abdomen_ambiguous])}
amb['BEST'] = amb['FLOW']

E = pd.concat(pd.read_parquet(f) for f in inp.events)
E = E[E.stage.isin(STAGES)]
inv = {v: k for k, v in codes.items()}
E['condition'] = E.label.map(inv)
E['pooled'] = np.where(E.label.isin(P['pooled_event_codes']), 'resp_event', None)

# ------------------------------------------------------------------ A. clean signal and LOSO group phase
mu, rowsA = {}, []
for s in SETS:
    for v in VARS:
        col = f'phase_{v}'
        c = E[(E.event_set == s) & (E.label == codes['clean'])].dropna(subset=[col])
        c = c[~c.sid.isin(amb[v.split('_')[0]])]
        cnt = c.groupby('sid').size()
        m = c.groupby('sid')[col].apply(lambda a: circ_mean(a.to_numpy()))[cnt >= MIN_MU[family(s)]]
        if len(m) < 5:
            continue
        cs, sn = np.cos(m).sum(), np.sin(m).sum()
        mu[(s, v)] = pd.Series({sid: np.arctan2(sn - np.sin(m[sid]), cs - np.cos(m[sid])) for sid in m.index})
        c = c[c.sid.isin(m.index)]
        y = np.cos(c[col] - c.sid.map(mu[(s, v)]))
        sm = y.groupby(c.sid).mean()
        rowsA.append(dict(event_set=s, family=family(s), variant=v, n_subjects=len(sm),
                          median_events=cnt[m.index].median(), group_phase_deg=np.degrees(np.arctan2(sn, cs)),
                          group_R=np.hypot(cs, sn) / len(m), rayleigh_p=rayleigh(m.to_numpy())[0],
                          vtest0_p=vtest(m.to_numpy(), 0.0)[0], mean_y=sm.mean(), sem_y=sm.sem(),
                          z=sm.mean() / sm.sem(), p_one_sided=stats.wilcoxon(sm, alternative='greater').pvalue,
                          event_level_mean=y.mean()))
A = pd.DataFrame(rowsA).sort_values('z', ascending=False)
cover = A[A.variant.isin(['BEST_hilbertgated', 'ABDOMEN_hilbertgated'])]
PRIMARIES = {fam: tuple(g.iloc[0][['event_set', 'variant']]) for fam, g in cover.groupby('family')}
A['primary'] = [(s_, v_) in PRIMARIES.values() for s_, v_ in zip(A.event_set, A.variant)]
A.to_csv(out.clean_signal, index=False)
print('Clean-N2 coupling signal:\n', A.round(4).to_string(index=False))
print('PRIMARIES:', PRIMARIES)

for v in VARS:
    E[f'y_{v}'] = np.nan
    for s in SETS:
        if (s, v) not in mu:
            continue
        sel = (E.event_set == s) & ~E.sid.isin(amb[v.split('_')[0]])
        E.loc[sel, f'y_{v}'] = np.cos(E.loc[sel, f'phase_{v}'] - E.loc[sel, 'sid'].map(mu[(s, v)]))

# ------------------------------------------------------------------ B. contrasts (both exclusion levels)
rowsB = []
for excl in ['none']:       # per-event mean cos assumes uniform phase occupancy -> not valid after phase-window exclusion
    Ex = E
    for s in SETS:
        for v in VARS:
            yc = f'y_{v}'
            g0 = Ex[(Ex.event_set == s) & Ex[yc].notna()]
            clean = g0[g0.label == codes['clean']]
            cm = clean.groupby('sid')[yc].agg(['mean', 'size'])
            cm = cm[cm['size'] >= MIN_CLEAN[family(s)]]
            for cond in CONDS:
                gc = g0[g0.pooled == 'resp_event'] if cond == 'resp_event' else g0[g0.condition == cond]
                km = gc.groupby('sid')[yc].agg(['mean', 'size'])
                km = km[km['size'] >= P['min_events']]
                sids = km.index.intersection(cm.index)
                base = dict(exclusion=excl, event_set=s, family=family(s), variant=v, condition=cond,
                            n_subjects=len(sids))
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
                                  clean_event_level=ev[ev.x == 0][yc].mean(),
                                  n_events_condition=int((ev.x == 1).sum()), n_events_clean=int((ev.x == 0).sum())))
B = pd.DataFrame(rowsB)
B['primary'] = [(s_, v_) in PRIMARIES.values() for s_, v_ in zip(B.event_set, B.variant)]
B['loss_fraction_excluded'] = -B.ev_ci_lo / B.clean_event_level    # loss larger than this fraction is outside the CI
B.to_csv(out.contrasts, index=False)
cols = ['exclusion', 'condition', 'n_subjects', 'n_events_condition', 'clean_event_level', 'ev_beta', 'ev_ci_lo',
        'ev_ci_hi', 'p_loss_event', 'p_two_sided_event', 'loss_fraction_excluded']
for fam, pr in PRIMARIES.items():
    print(f'\nContrasts, primary {fam}: {pr}\n',
          B[(B.event_set == pr[0]) & (B.variant == pr[1]) & B.condition.isin(['resp_event', 'snoring', 'hypopnea',
                                                                                'recovery'])][cols]
          .round(4).to_string(index=False))

# ------------------------------------------------------------------ C. dose-response
rowsC = []
for excl in ['none']:       # see above; exclusion is handled by the occupancy-normalised analyses (B2, GLMM)
    Ex = E
    for s in SETS:
        for v in VARS:
            yc = f'y_{v}'
            g = Ex[(Ex.event_set == s) & Ex[yc].notna() & Ex.label.isin([0, 1, 2, 3, 4, 5, 6])].copy()
            for pr in ['snore_db', 'flow_log2', 'effort_log2']:
                g[pr + '_c'] = g[pr] - g.groupby('sid')[pr].transform('mean')
                g[pr + '_c'] /= g[pr + '_c'].std()
            g = g[np.isfinite(g[['snore_db_c', 'flow_log2_c', 'effort_log2_c']]).all(1)]
            if g.sid.nunique() < 10:
                continue
            r = smf.ols(f'{yc} ~ snore_db_c + flow_log2_c + effort_log2_c', g).fit(
                cov_type='cluster', cov_kwds={'groups': pd.factorize(g.sid)[0]})
            row = dict(exclusion=excl, event_set=s, variant=v, n_subjects=g.sid.nunique(), n_events=len(g),
                       intercept=r.params['Intercept'])
            for pr, loss_sign in [('snore_db_c', -1), ('flow_log2_c', +1), ('effort_log2_c', +1)]:
                b, se = r.params[pr], r.bse[pr]
                # one-sided p for loss: louder snoring -> beta < 0; lower airflow / effort -> beta > 0
                p_loss = stats.norm.cdf(b / se) if loss_sign < 0 else 1 - stats.norm.cdf(b / se)
                row.update({f'{pr}_beta': b, f'{pr}_se': se, f'{pr}_p_two_sided': r.pvalues[pr], f'{pr}_p_loss': p_loss})
            rowsC.append(row)
C = pd.DataFrame(rowsC)
C['primary'] = [(s_, v_) in PRIMARIES.values() for s_, v_ in zip(C.event_set, C.variant)]
C.to_csv(out.dose, index=False)
print('\nDose-response, primaries:\n', C[C.primary][['exclusion', 'event_set', 'variant', 'n_subjects', 'n_events',
                                                      'snore_db_c_beta', 'snore_db_c_p_two_sided', 'flow_log2_c_beta',
                                                      'flow_log2_c_p_two_sided', 'effort_log2_c_beta',
                                                      'effort_log2_c_p_two_sided']].round(4).to_string(index=False))

# ------------------------------------------------------------------ D. snoring diagnostics (phase shift, density)
rowsD = []
for fam, (s, v) in PRIMARIES.items():
    for excl in ['none']:
        Ex = E
        g = Ex[(Ex.event_set == s) & Ex[f'phase_{v}'].notna() & Ex.label.isin([0, 1]) & Ex.sid.isin(mu[(s, v)].index)]
        rel = g[f'phase_{v}'] - g.sid.map(mu[(s, v)])                 # phase relative to the LOSO group direction
        g = g.assign(cx=np.cos(rel), sx=np.sin(rel))
        V = g.groupby(['sid', 'label'])[['cx', 'sx']].mean().unstack('label')
        cnt = g.groupby(['sid', 'label']).size().unstack('label')
        ok = (cnt[0] >= MIN_CLEAN[fam]) & (cnt[1] >= P['min_events'])
        V = V[ok]
        dx, dy = V['cx'][1] - V['cx'][0], V['sx'][1] - V['sx'][0]
        Xd = np.c_[dx, dy]
        n_ = len(Xd)
        if n_ > 3:
            mvec = Xd.mean(0)
            S = np.cov(Xd, rowvar=False)
            t2 = n_ * mvec @ np.linalg.solve(S, mvec)
            p_h = 1 - stats.f.cdf((n_ - 2) / (2 * (n_ - 1)) * t2, 2, n_ - 2)
        else:
            p_h = np.nan
        rowsD.append(dict(family=fam, event_set=s, variant=v, exclusion=excl, n_subjects=n_,
                          clean_angle_deg=np.degrees(np.arctan2(V['sx'][0].mean(), V['cx'][0].mean())),
                          snore_angle_deg=np.degrees(np.arctan2(V['sx'][1].mean(), V['cx'][1].mean())),
                          clean_length=np.hypot(V['sx'][0].mean(), V['cx'][0].mean()),
                          snore_length=np.hypot(V['sx'][1].mean(), V['cx'][1].mean()),
                          orthogonal_shift=dy.mean(), orthogonal_shift_p=stats.wilcoxon(dy).pvalue if n_ > 5 else np.nan,
                          hotelling_p=p_h))
D = pd.DataFrame(rowsD)
D.to_csv(out.snoring_diag, index=False)
print('\nSnoring vs clean (vectors relative to group phase):\n', D.round(4).to_string(index=False))

# ------------------------------------------------------------------ export phase-bin data for the GLMM (idea 4)
Bn = pd.concat(pd.read_parquet(f) for f in inp.bins)
Bn = Bn[Bn.stage_set == SS]

# ------------------------------------------------------------------ B2. occupancy-normalised contrasts (valid with exclusion)
# subject coupling = projection on mu_LOSO of the rate-weighted resultant: M = sum_b r_b e^{i phi_b} / sum_b r_b,
# r_b = count_b / occupancy_b (time in bin outside excluded windows); immune to non-uniform phase occupancy.
rowsB2 = []
Bx = Bn[Bn.condition.isin(['clean'] + CONDS)].copy()
Bx['rate'] = Bx['count'] / Bx.occupancy_min.where(Bx.occupancy_min > 0)
for (s, v, excl), g in Bx.groupby(['event_set', 'variant', 'exclusion']):
    if (s, v) not in mu:
        continue
    g = g[g.sid.isin(mu[(s, v)].index)]
    agg = g.groupby(['sid', 'condition']).apply(
        lambda x: pd.Series({'n': x['count'].sum(),
                             'pcw': np.real(np.nansum(x.rate * np.exp(1j * x.bin_center)) / np.nansum(x.rate)
                                            * np.exp(-1j * mu[(s, v)][x.sid.iloc[0]])) if np.nansum(x.rate) > 0 else np.nan}))
    agg = agg.reset_index()
    clean = agg[(agg.condition == 'clean') & (agg.n >= MIN_CLEAN[family(s)])].set_index('sid')
    for cond in CONDS:
        cc = agg[(agg.condition == cond) & (agg.n >= P['min_events'])].set_index('sid')
        sids = cc.index.intersection(clean.index)
        base = dict(exclusion=excl, event_set=s, family=family(s), variant=v, condition=cond, n_subjects=len(sids))
        if len(sids) < 5:
            rowsB2.append(base)
            continue
        dlt = (cc.loc[sids, 'pcw'] - clean.loc[sids, 'pcw']).dropna()
        boot = [np.mean(np.random.default_rng(i).choice(dlt, len(dlt))) for i in range(2000)]
        rowsB2.append(dict(base, clean_level=clean.loc[sids, 'pcw'].mean(), condition_level=cc.loc[sids, 'pcw'].mean(),
                           delta_mean=dlt.mean(), ci_lo=np.percentile(boot, 2.5), ci_hi=np.percentile(boot, 97.5),
                           p_loss=stats.wilcoxon(dlt, alternative='less').pvalue,
                           p_two_sided=stats.wilcoxon(dlt).pvalue, n_events_condition=int(cc.loc[sids, 'n'].sum())))
B2 = pd.DataFrame(rowsB2)
B2['primary'] = [(s_, v_) in PRIMARIES.values() for s_, v_ in zip(B2.event_set, B2.variant)]
B2.to_csv(out.contrasts_occ, index=False)
for fam, pr in PRIMARIES.items():
    print(f'\nOccupancy-normalised contrasts, primary {fam}: {pr}\n',
          B2[(B2.event_set == pr[0]) & (B2.variant == pr[1]) & B2.condition.isin(['resp_event', 'snoring', 'hypopnea',
                                                                                  'recovery'])]
          [['exclusion', 'condition', 'n_subjects', 'n_events_condition', 'clean_level', 'condition_level', 'delta_mean',
            'ci_lo', 'ci_hi', 'p_loss', 'p_two_sided']].round(4).to_string(index=False))
Bn['polarity_ok'] = [sd not in amb[v.split('_')[0]] for sd, v in zip(Bn.sid, Bn.variant)]
gph = {(r.event_set, r.variant): np.radians(r.group_phase_deg) for r in A.itertuples()}
Bn['mu_group'] = [gph.get((s, v), np.nan) for s, v in zip(Bn.event_set, Bn.variant)]
Bn['primary'] = [(s_, v_) in PRIMARIES.values() for s_, v_ in zip(Bn.event_set, Bn.variant)]
Bn['family'] = Bn.event_set.map(family)
Bn[Bn.primary].to_csv(out.bins_export, index=False)

# ------------------------------------------------------------------ figure: contrasts with / without exclusion
fams = list(PRIMARIES)
fig, axs = plt.subplots(1, len(fams), figsize=(6 * len(fams), 5), sharey=True)
show = ['resp_event', 'snoring', 'hypopnea', 'recovery']
for ax, fam in zip(np.atleast_1d(axs), fams):
    s, v = PRIMARIES[fam]
    for k, (excl, colr) in enumerate([('none', 'k'), ('evoked', 'tab:purple')]):
        pb = B2[(B2.event_set == s) & (B2.variant == v) & (B2.exclusion == excl)].set_index('condition').reindex(show)
        yy = np.arange(len(show)) + (k - .5) * .3
        ax.errorbar(pb.delta_mean, yy, xerr=[pb.delta_mean - pb.ci_lo, pb.ci_hi - pb.delta_mean], fmt='o', color=colr,
                    capsize=3, label='all events' if excl == 'none' else 'snore-evoked excluded (phase-matched)')
        for yi, (_, r) in zip(yy, pb.iterrows()):
            if np.isfinite(r.get('delta_mean', np.nan)):
                ax.text(0.01, yi + .12, f'N={int(r.n_subjects)} n={int(r.n_events_condition)} p={r.p_two_sided:.2g}',
                        fontsize=6, color=colr, transform=ax.get_yaxis_transform())
        lvl = pb.clean_level.mean()
        ax.axvline(-lvl, color=colr, ls='--', lw=.7)
    ax.axvline(0, color='grey', lw=.7)
    ax.set_yticks(range(len(show)))
    ax.set_yticklabels(show)
    ax.set_title(f'{fam}: {s} | {v}', fontsize=9)
    ax.set_xlabel('condition − clean (N2), occupancy-normalised coupling,\nmean ± bootstrap 95% CI (dashed = complete loss)')
np.atleast_1d(axs)[0].invert_yaxis()
np.atleast_1d(axs)[0].legend(fontsize=7, loc='lower right')
fig.tight_layout()
fig.savefig(out.fig, dpi=80)
