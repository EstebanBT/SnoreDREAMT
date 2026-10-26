"""Snakemake rule `step6_phase_severity`: do respiratory-phase effects relate to sleep-disordered-breathing severity?

6a (within-subject phase change vs severity), per event family (Step 3b primaries) x condition (pooled events,
   snoring, hypopnea, recovery; snoring also with snore-evoked events excluded):
   - signed per-subject mean phase difference d_mu (Step 5 CircStat analysis) vs each marker: circ_corrcl;
   - size of the change: |d_mu| and the subject's common-median evidence (-log10 p_cm) vs each marker: Spearman, and
     rank regression adjusted for age, sex and BMI (robust to outliers in outcome and marker).
   Within-subject differences are polarity-invariant (all subjects).
6b (clean preferred phase vs severity): subject clean-N2 preferred phase (event circular mean; polarity-unambiguous
   subjects only) vs each marker: circ_corrcl; AHI groups (< 5, 5-15, 15-30, >= 30) compared with the k-group common
   median test and Watson-Williams (+ CircStat applicability flag).
BH-FDR within 6a and within 6b. Also writes the phase-bin data with covariates for the GLMM (step6_glmm.R).
"""
import json
import warnings

import matplotlib
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt import circstat as cs  # noqa: E402
from snoredreamt.io import participants  # noqa: E402

warnings.filterwarnings('ignore', message='Ties detected')
cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
SS = snakemake.wildcards.ss  # noqa: F821
STAGES = cfg['step3b']['stage_sets'][SS]
codes = cfg['condition_codes']

Sv = pd.DataFrame([json.load(open(f)) for f in inp.severity]).set_index('sid').join(participants(cfg))
Sv['male'] = (Sv.GENDER == 'M').astype(int)
MARK = {'log AHI': np.log1p(Sv.AHI), 'log OAHI': np.log1p(Sv.OAHI), 'log ODI3': np.log1p(Sv.odi3),
        'log T90': np.log1p(Sv.t90), 'log hypoxic burden': np.log1p(Sv.hypoxic_burden),
        'mean event duration': Sv.mean_event_duration_s, 'snore index': Sv.snore_breath_frac}
for k, v in MARK.items():
    Sv[k] = v
for c in ['AGE', 'BMI']:
    Sv[f'z_{c}'] = (Sv[c] - Sv[c].mean()) / Sv[c].std()


def robust_beta(y, x, d):
    """Rank regression (outcome, marker and covariates rank-transformed, HC3 SEs): robust to outliers in both the
    outcome and the marker (Huber RLM is not robust to leverage: one subject with mean event duration z = 4.3 produced
    p ~ 0). Returns the standardised slope (partial rank association) and its p value."""
    dd = d[[y, x, 'AGE', 'male', 'BMI']].dropna().rank()
    dd = (dd - dd.mean()) / dd.std()
    r = sm.OLS(dd[y], sm.add_constant(dd[[x, 'AGE', 'male', 'BMI']])).fit(cov_type='HC3')
    return r.params[x], r.pvalues[x]


# ------------------------------------------------------------------ 6a
C = pd.read_csv(inp.circstat_subjects)
C = C[((C.exclusion == 'none') & C.condition.isin(['resp_event', 'snoring', 'hypopnea', 'recovery'])) |
      ((C.exclusion == 'evoked') & (C.condition == 'snoring'))].copy()
C['abs_dmu'] = np.abs(C.d_mean)
C['cm_evidence'] = -np.log10(C.p_cm)
rows = []
for (fam, excl, cond), g in C.groupby(['family', 'exclusion', 'condition']):
    d = g.set_index('sid').join(Sv, how='inner')
    if len(d) < 15:
        continue
    for m in MARK:
        dm = d.dropna(subset=[m])
        rho_c, p_c = cs.circ_corrcl(dm.d_mean.to_numpy(), dm[m].to_numpy())
        rho_a, p_a = stats.spearmanr(dm.abs_dmu, dm[m])
        rho_e, p_e = stats.spearmanr(dm.cm_evidence, dm[m])
        b_a, pb_a = robust_beta('abs_dmu', m, dm)
        b_e, pb_e = robust_beta('cm_evidence', m, dm)
        rows.append(dict(family=fam, exclusion=excl, condition=cond, marker=m, n=len(dm),
                         circlin_rho=rho_c, circlin_p=p_c, absdmu_spearman=rho_a, absdmu_spearman_p=p_a,
                         absdmu_rankadj_beta=b_a, absdmu_rankadj_p=pb_a, cm_spearman=rho_e, cm_spearman_p=p_e,
                         cm_rankadj_beta=b_e, cm_rankadj_p=pb_e))
A6 = pd.DataFrame(rows)
allp = pd.concat([A6.circlin_p, A6.absdmu_rankadj_p, A6.cm_rankadj_p])
fdr = multipletests(allp, method='fdr_bh')[1]
n6 = len(A6)
A6['circlin_p_fdr'], A6['absdmu_rankadj_p_fdr'], A6['cm_rankadj_p_fdr'] = fdr[:n6], fdr[n6:2 * n6], fdr[2 * n6:]
A6.to_csv(out.table_6a, index=False)
print('6a: smallest p values (FDR over all 6a tests):')
for col in ['circlin_p', 'absdmu_rankadj_p', 'cm_rankadj_p']:
    print(A6.sort_values(col).head(5)[['family', 'exclusion', 'condition', 'marker', 'n', col, f'{col}_fdr']].round(4)
          .to_string(index=False))

# ------------------------------------------------------------------ 6b
Q = pd.read_csv(inp.resp_summary).set_index('sid')
amb = {'FLOW': set(Q.index[Q.flow_ambiguous]), 'ABDOMEN': set(Q.index[Q.abdomen_ambiguous])}
amb['BEST'] = amb['FLOW']
Acs = pd.read_csv(inp.clean_signal)
PRIM = {r.family: (r.event_set, r.variant) for r in Acs[Acs.primary].itertuples()}
MIN = {'SO': 20, 'spindle': 20, 'complex': 10}
E = pd.concat(pd.read_parquet(f) for f in inp.events)
E = E[E.stage.isin(STAGES) & (E.label == codes['clean'])]
rows, phase_rows = [], []
for fam, (es, var) in PRIM.items():
    g = E[(E.event_set == es) & E[f'phase_{var}'].notna() & ~E.sid.isin(amb[var.split('_')[0]])]
    cnt = g.groupby('sid').size()
    mu = g.groupby('sid')[f'phase_{var}'].apply(lambda a: cs.circ_mean(a.to_numpy()))[cnt >= MIN[fam]]
    d = mu.rename('phase').to_frame().join(Sv, how='inner')
    phase_rows.append(d.assign(family=fam)[['family', 'phase', 'AHI']])
    for m in MARK:
        dm = d.dropna(subset=[m])
        rho, p = cs.circ_corrcl(dm.phase.to_numpy(), dm[m].to_numpy())
        rows.append(dict(family=fam, analysis='circular-linear', marker=m, n=len(dm), stat=rho, p=p))
    grp = pd.cut(d.AHI, [-1, 5, 15, 30, 1000], labels=['<5', '5-15', '15-30', '>=30'])
    groups = [d.phase[grp == k].to_numpy() for k in grp.cat.categories]
    p_cm, _, P = cs.circ_cmtest_k(groups)
    p_ww, F, ok = cs.circ_wwtest_k(groups)
    gm = {k: np.degrees(cs.circ_mean(a)) for k, a in zip(grp.cat.categories, groups)}
    rows.append(dict(family=fam, analysis='AHI groups: common median', marker='AHI group', n=len(d), stat=P, p=p_cm,
                     **{f'mean_deg_{k}': v for k, v in gm.items()},
                     **{f'n_{k}': len(a) for k, a in zip(grp.cat.categories, groups)}))
    rows.append(dict(family=fam, analysis=f'AHI groups: Watson-Williams (applicable={ok})', marker='AHI group',
                     n=len(d), stat=F, p=p_ww))
B6 = pd.DataFrame(rows)
B6['p_fdr'] = multipletests(B6.p, method='fdr_bh')[1]
B6.to_csv(out.table_6b, index=False)
with pd.option_context('display.width', 250, 'display.max_columns', 20):
    print('\n6b:\n', B6.round(4).to_string(index=False))

# ------------------------------------------------------------------ GLMM export
Bn = pd.read_csv(inp.bins)
Bn = Bn[Bn.primary & Bn.polarity_ok]
Bn = Bn[(Bn.exclusion == 'none') | (Bn.condition.isin(['clean', 'snoring']))]
cov = pd.DataFrame({'z_log_ahi': (np.log1p(Sv.AHI) - np.log1p(Sv.AHI).mean()) / np.log1p(Sv.AHI).std()})
Bn.merge(cov, left_on='sid', right_index=True).to_csv(out.glmm_data, index=False)

# ------------------------------------------------------------------ figure
P6 = pd.concat(phase_rows)
fams = list(PRIM)
fig = plt.figure(figsize=(6 * len(fams), 10))
for j, fam in enumerate(fams):
    ax = fig.add_subplot(2, len(fams), j + 1)
    d = P6[P6.family == fam]
    ax.scatter(np.log1p(d.AHI), np.degrees(d.phase), s=12)
    r = B6[(B6.family == fam) & (B6.marker == 'log AHI')].iloc[0]
    ax.set_title(f'6b {fam}: clean preferred phase vs AHI\ncirc-lin ρ={r.stat:.2f}, p={r.p:.2g} (N={r.n})', fontsize=9)
    ax.set_xlabel('log(1 + AHI)')
    ax.set_ylabel('clean preferred phase (deg)')
    ax.set_ylim(-180, 180)
    ax = fig.add_subplot(2, len(fams), len(fams) + j + 1)
    a6 = A6[(A6.family == fam) & (A6.exclusion == 'none')]
    piv = a6.pivot(index='marker', columns='condition', values='circlin_p')
    im = ax.imshow(-np.log10(piv.to_numpy()), cmap='viridis', vmin=0, vmax=3, aspect='auto')
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels(piv.columns, fontsize=7, rotation=30)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels(piv.index, fontsize=7)
    ax.set_title(f'6a {fam}: −log10 p, circ-lin corr. of phase shift with markers', fontsize=9)
    fig.colorbar(im, ax=ax, shrink=.7)
fig.tight_layout()
fig.savefig(out.fig, dpi=80)
