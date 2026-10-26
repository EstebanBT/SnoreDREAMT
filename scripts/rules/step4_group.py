"""Snakemake rule `step4_group`: does sleep-disordered-breathing severity relate to respiration-oscillation coupling?

Coupling per subject (Step 3b/3c, clean N2, all events, primaries per family): occupancy-normalised coupling
projected on the leave-one-subject-out group phase (polarity-unambiguous subjects only).
1. Primary (one test per family, BH-FDR over the 3 families): robust regression (Huber, statsmodels RLM)
   coupling ~ z(log1p AHI) + z(age) + sex + z(BMI). Arousal index is a mediator and is not adjusted for.
   Sensitivity: OLS (HC3), WLS weighted by event count, + depression (antidepressant proxy), Spearman.
2. Secondary (BH-FDR over markers x families): same model with each data-derived marker instead of AHI.
3. Exploratory: within-subject change (snoring / pooled events - clean; all events and snore-evoked excluded) vs AHI.
Writes the clean-N2 phase-bin data with subject covariates for the GLMM moderator model (step4_glmm.R).
"""
import json

import matplotlib
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt.io import participants  # noqa: E402

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
MIN_CLEAN = {'SO': 20, 'spindle': 20, 'complex': 10}

Sv = pd.DataFrame([json.load(open(f)) for f in inp.severity]).set_index('sid')
Pp = participants(cfg)
Sv = Sv.join(Pp)
Sv['male'] = (Sv.GENDER == 'M').astype(int)
Sv['depression'] = Sv.MEDICAL_HISTORY.fillna('').str.contains('Depression', case=False).astype(int)
Sv['log_ahi'] = np.log1p(Sv.AHI)
Sv['log_odi3'] = np.log1p(Sv.odi3)
Sv['log_hb'] = np.log1p(Sv.hypoxic_burden)
Sv['log_t90'] = np.log1p(Sv.t90)
Sv['log_oahi'] = np.log1p(Sv.OAHI)
Sv['log_arousal_clin'] = np.log1p(Sv['Arousal Index'])
Sv['log_arousal_det'] = np.log1p(Sv.arousal_index_detected)
Sv['snore_frac'] = Sv.snore_breath_frac
Sv['nadir_proxy'] = -Sv.mean_sao2                     # higher = worse
for c in ['AGE', 'BMI']:
    Sv[f'z_{c}'] = (Sv[c] - Sv[c].mean()) / Sv[c].std()

U = pd.read_csv(inp.subjects)
U = U[U.primary & (U.exclusion == 'none')]
clean = U[U.condition == 'clean'].copy()
clean = clean[clean.n >= clean.family.map(MIN_CLEAN)]
FAM = {f: g.set_index('sid') for f, g in clean.groupby('family')}


def zs(x):
    return (x - x.mean()) / x.std()


def fit_all(fam, marker, data):
    d = data.dropna(subset=['coupling', marker, 'z_AGE', 'z_BMI']).copy()
    d['zm'] = zs(d[marker])
    X = sm.add_constant(d[['zm', 'z_AGE', 'male', 'z_BMI']])
    rlm = sm.RLM(d.coupling, X, M=sm.robust.norms.HuberT()).fit()
    ols = smf.ols('coupling ~ zm + z_AGE + male + z_BMI', d).fit(cov_type='HC3')
    wls = smf.wls('coupling ~ zm + z_AGE + male + z_BMI', d, weights=d.n).fit(cov_type='HC3')
    dep = smf.ols('coupling ~ zm + z_AGE + male + z_BMI + depression', d).fit(cov_type='HC3')
    rho, prho = stats.spearmanr(d[marker], d.coupling)
    return dict(family=fam, marker=marker, n=len(d), beta_rlm=rlm.params['zm'], se_rlm=rlm.bse['zm'],
                p_rlm=rlm.pvalues['zm'], beta_ols=ols.params['zm'], p_ols=ols.params['zm'] * 0 + ols.pvalues['zm'],
                beta_wls=wls.params['zm'], p_wls=wls.pvalues['zm'], beta_dep=dep.params['zm'], p_dep=dep.pvalues['zm'],
                spearman_rho=rho, spearman_p=prho, coupling_mean=d.coupling.mean())


rows = []
MARKERS = ['log_ahi', 'log_oahi', 'log_odi3', 'log_t90', 'log_hb', 'mean_event_duration_s', 'snore_frac',
           'log_arousal_clin', 'log_arousal_det', 'nadir_proxy']
for fam, g in FAM.items():
    data = g[['pcw', 'n']].rename(columns={'pcw': 'coupling'}).join(Sv, how='inner')
    for m in MARKERS:
        rows.append(fit_all(fam, m, data))
R = pd.DataFrame(rows)
prim = R.marker == 'log_ahi'
R['primary'] = prim
R['p_fdr'] = np.nan
R.loc[prim, 'p_fdr'] = multipletests(R.loc[prim, 'p_rlm'], method='fdr_bh')[1]
sec = ~prim & ~R.marker.isin(['log_arousal_clin', 'log_arousal_det'])   # arousal indices: mediators, reported only
R.loc[sec, 'p_fdr'] = multipletests(R.loc[sec, 'p_rlm'], method='fdr_bh')[1]
R.to_csv(out.table, index=False)
with pd.option_context('display.width', 250):
    print(R.round(4).to_string(index=False))

# ------------------------------------------------------------------ within-subject change vs AHI (exploratory)
U2 = pd.read_csv(inp.subjects)
U2 = U2[U2.primary]
rows = []
for (fam, excl), g in U2.groupby(['family', 'exclusion']):
    piv = g.pivot_table(index='sid', columns='condition', values='pcw')
    cnt = g.pivot_table(index='sid', columns='condition', values='n')
    for cond in ['snoring', 'resp_event', 'hypopnea']:
        if cond not in piv:
            continue
        ok = (cnt['clean'] >= MIN_CLEAN[fam]) & (cnt[cond] >= cfg['step3b']['min_events'])
        dlt = (piv[cond] - piv['clean'])[ok].rename('coupling').to_frame().join(Sv, how='inner')
        dlt['n'] = 1
        if len(dlt) >= 10:
            rows.append(dict(fit_all(fam, 'log_ahi', dlt), condition=cond, exclusion=excl))
W = pd.DataFrame(rows)
W.to_csv(out.change, index=False)
print('\nWithin-subject change vs AHI:\n', W[['family', 'condition', 'exclusion', 'n', 'beta_rlm', 'p_rlm',
                                              'spearman_rho', 'spearman_p']].round(4).to_string(index=False))

# ------------------------------------------------------------------ GLMM export (clean N2 bins + covariates)
Bn = pd.read_csv(inp.bins)
Bn = Bn[(Bn.condition == 'clean') & (Bn.exclusion == 'none')]
cov = Sv[['log_ahi', 'z_AGE', 'male', 'z_BMI']].copy()
cov['z_log_ahi'] = zs(cov.log_ahi)
Bn.merge(cov, left_on='sid', right_index=True).to_csv(out.glmm_data, index=False)

# ------------------------------------------------------------------ figures
fams = list(FAM)
fig, axs = plt.subplots(1, len(fams) + 1, figsize=(5 * (len(fams) + 1), 4.5))
for ax, fam in zip(axs, fams):
    d = FAM[fam][['pcw']].join(Sv, how='inner').dropna(subset=['pcw', 'log_ahi'])
    ax.scatter(d.AHI + 1, d.pcw, s=12, alpha=.7)
    ax.set_xscale('log')
    b = np.polyfit(d.log_ahi, d.pcw, 1)
    xs = np.linspace(d.log_ahi.min(), d.log_ahi.max(), 50)
    ax.plot(np.expm1(xs) + 1, np.polyval(b, xs), 'r')
    r = R[(R.family == fam) & prim].iloc[0]
    ev = FAM[fam].event_set.iloc[0]
    ax.set_title(f'{fam} ({ev}): β={r.beta_rlm:.3f}/SD, p={r.p_rlm:.2g}, FDR={r.p_fdr:.2g}\nρ={r.spearman_rho:.2f}, N={r.n}',
                 fontsize=8)
    ax.axhline(0, color='grey', lw=.6)
    ax.set_xlabel('AHI + 1 (log scale)')
    ax.set_ylabel('clean-N2 coupling (projection on group phase)')
ax = axs[-1]
Rs = R[~R.marker.isin(['log_oahi'])]
markers = [m for m in MARKERS if m != 'log_oahi']
for k, fam in enumerate(fams):
    rr = Rs[Rs.family == fam].set_index('marker').reindex(markers)
    y = np.arange(len(markers)) + (k - 1) * .25
    ax.errorbar(rr.beta_rlm, y, xerr=1.96 * rr.se_rlm, fmt='o', capsize=2, label=fam)
ax.axvline(0, color='grey', lw=.6)
ax.set_yticks(range(len(markers)))
ax.set_yticklabels(markers, fontsize=7)
ax.invert_yaxis()
ax.set_xlabel('β per SD of marker (robust, adj. age/sex/BMI) ± 95% CI')
ax.legend(fontsize=7)
ax.set_title('Severity markers vs clean coupling', fontsize=9)
fig.tight_layout()
fig.savefig(out.fig, dpi=80)
