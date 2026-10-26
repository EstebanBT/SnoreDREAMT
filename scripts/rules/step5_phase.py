"""Snakemake rule `step5_phase`: does the preferred respiratory phase of SOs / spindles / complexes change during
breathing conditions (vs whole-night clean breathing)?

Data: Step 3b phase-bin counts and occupancy per subject x condition x exclusion (primaries per family).
Per subject and condition: occupancy-normalised resultant M = sum_b r_b e^{i phi_b} / sum_b r_b, r_b = count/occupancy
(bins where both conditions have occupancy). Signed phase difference dphi = arg(M_cond / M_ref).
5a (any change, direction-free): |dphi| against a parametric null in which both conditions share the pooled rate
   profile rho_b = (k_c + k_ref) / (o_c + o_ref) and counts are drawn Multinomial(n, rho_b * o_b) for each condition
   (keeps each condition's occupancy and event count) -> per-subject z, p. Group: one-sided Wilcoxon of z > 0,
   Stouffer combination of p, binomial test of the fraction p < .05.
5b (consistent rotation): Rayleigh test on the subjects' dphi (all; and reliability-filtered = subjects whose clean
   resultant exceeds a uniform-rate null at p < .1), circular mean rotation with bootstrap CI; GLMM rotation estimate
   (step3b_glmm.R) for comparison.
Within-subject dphi is polarity-invariant, so all subjects are used. Control: clean_odd vs clean_even (10-min blocks).
"""
import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from snoredreamt.coupling import circ_mean, circ_r, rayleigh  # noqa: E402

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
P = cfg['step3b']
rng = np.random.default_rng(5)
N_NULL = 2000
MIN_REF = {'SO': 20, 'spindle': 20, 'complex': 10}
CONDS = ['resp_event', 'snoring', 'hypopnea', 'recovery', 'obstructive', 'central', 'mixed', 'clean_odd']

B = pd.read_csv(inp.bins)
B = B[B.primary]
G = pd.read_csv(inp.glmm)


def resultant(k, o, phi):
    r = np.where(o > 0, k / np.where(o > 0, o, 1), 0.0)
    s = r.sum(-1, keepdims=True)
    return (r * np.exp(1j * phi)).sum(-1) / np.where(s[..., 0] > 0, s[..., 0], np.nan)


def subject_stats(kc, oc, k0, o0, phi):
    ok = (oc > 0) & (o0 > 0)
    kc, oc, k0, o0, phi = kc[ok], oc[ok], k0[ok], o0[ok], phi[ok]
    nc, n0 = int(kc.sum()), int(k0.sum())
    Mc, M0 = resultant(kc, oc, phi), resultant(k0, o0, phi)
    dphi = np.angle(Mc / M0)
    rho = (kc + k0) / (oc + o0)
    pc, p0 = rho * oc, rho * o0
    if pc.sum() == 0 or p0.sum() == 0:
        return None
    sc = rng.multinomial(nc, pc / pc.sum(), size=N_NULL)
    s0 = rng.multinomial(n0, p0 / p0.sum(), size=N_NULL)
    nd = np.abs(np.angle(resultant(sc, oc, phi) / resultant(s0, o0, phi)))
    nd = nd[np.isfinite(nd)]
    # reliability of the clean (reference) coupling: |M0| vs a uniform-rate null (counts proportional to occupancy)
    u0 = rng.multinomial(n0, o0 / o0.sum(), size=N_NULL)
    p_rel = (np.sum(np.abs(resultant(u0, o0, phi)) >= np.abs(M0)) + 1) / (N_NULL + 1)
    return dict(n_cond=nc, n_ref=n0, len_cond=np.abs(Mc), len_ref=np.abs(M0), dphi=dphi, abs_dphi=np.abs(dphi),
                z=(np.abs(dphi) - nd.mean()) / nd.std(), p=(np.sum(nd >= np.abs(dphi)) + 1) / (len(nd) + 1),
                p_ref_reliable=p_rel)


rows = []
for (es, var, excl), g in B.groupby(['event_set', 'variant', 'exclusion']):
    fam = g.family.iloc[0]
    for cond in CONDS:
        ref = 'clean_even' if cond == 'clean_odd' else 'clean'
        for sid, gs in g[g.condition.isin([cond, ref])].groupby('sid'):
            piv_k = gs.pivot_table(index='condition', columns='bin', values='count', aggfunc='sum')
            piv_o = gs.pivot_table(index='condition', columns='bin', values='occupancy_min', aggfunc='sum')
            if cond not in piv_k.index or ref not in piv_k.index:
                continue
            kc, k0 = piv_k.loc[cond].to_numpy(float), piv_k.loc[ref].to_numpy(float)
            if kc.sum() < P['min_events'] or k0.sum() < MIN_REF[fam]:
                continue
            phi = gs.drop_duplicates('bin').set_index('bin').loc[piv_k.columns, 'bin_center'].to_numpy()
            r = subject_stats(kc, piv_o.loc[cond].to_numpy(float), k0, piv_o.loc[ref].to_numpy(float), phi)
            if r is not None:
                rows.append(dict(r, sid=sid, family=fam, event_set=es, variant=var, exclusion=excl, condition=cond))
S = pd.DataFrame(rows)
S.to_csv(out.subjects, index=False)


KEYS = ['p_5b_rotation_ne0', 'p_5b_rotation_ne0_reliable', 'median_z', 'p_5a_wilcoxon', 'stouffer_z', 'p_5a_stouffer', 'frac_p05', 'p_5a_binomial', 'mean_rotation_deg',
        'rot_ci_lo', 'rot_ci_hi', 'rotation_R', 'p_5b_rayleigh', 'n_reliable', 'p_5b_rayleigh_reliable',
        'mean_rotation_reliable_deg']


def group(g):
    res = dict(n_subjects=len(g), median_events_cond=g.n_cond.median(), **{k: np.nan for k in KEYS})
    if len(g) < 5:
        return pd.Series(res)
    zs = stats.norm.isf(np.clip(g.p, 1e-6, 1 - 1e-6))
    pr, zr = rayleigh(g.dphi.to_numpy())
    rel = g[g.p_ref_reliable < 0.1]
    boot = [circ_mean(rng.choice(g.dphi.to_numpy(), len(g))) for _ in range(2000)]
    mu = circ_mean(g.dphi.to_numpy())
    bc = np.degrees(np.angle(np.exp(1j * (np.array(boot) - mu))))
    # 5b proper: is the mean rotation different from 0? (Rayleigh only says dphi is non-uniform, which is also true
    # when there is NO change and both conditions are coupled: dphi then clusters at 0, see clean_odd control)
    res['p_5b_rotation_ne0'] = stats.wilcoxon(np.sin(g.dphi)).pvalue
    res['p_5b_rotation_ne0_reliable'] = stats.wilcoxon(np.sin(rel.dphi)).pvalue if len(rel) >= 6 else np.nan
    res.update(median_z=g.z.median(), p_5a_wilcoxon=stats.wilcoxon(g.z, alternative='greater').pvalue,
               stouffer_z=zs.sum() / np.sqrt(len(zs)), p_5a_stouffer=stats.norm.sf(zs.sum() / np.sqrt(len(zs))),
               frac_p05=(g.p < .05).mean(), p_5a_binomial=stats.binomtest(int((g.p < .05).sum()), len(g), .05,
                                                                          alternative='greater').pvalue,
               mean_rotation_deg=np.degrees(mu), rot_ci_lo=np.degrees(mu) + np.percentile(bc, 2.5),
               rot_ci_hi=np.degrees(mu) + np.percentile(bc, 97.5), rotation_R=circ_r(g.dphi), p_5b_rayleigh=pr,
               n_reliable=len(rel),
               p_5b_rayleigh_reliable=rayleigh(rel.dphi.to_numpy())[0] if len(rel) >= 5 else np.nan,
               mean_rotation_reliable_deg=np.degrees(circ_mean(rel.dphi)) if len(rel) >= 5 else np.nan)
    return pd.Series(res)


T = S.groupby(['family', 'event_set', 'variant', 'exclusion', 'condition']).apply(group).reset_index()

# 5a calibrated: the parametric null ignores temporal clustering of events (spindle trains, merged SOs), so the
# clean-vs-clean control is itself 'significant' for spindles/complexes; test each condition's z against the same
# subject's control z (one-sided paired Wilcoxon: phase difference larger during the condition than between two
# halves of clean breathing).
ctrl = S[S.condition == 'clean_odd'][['sid', 'family', 'exclusion', 'z']].rename(columns={'z': 'z_ctrl'})
SC = S.merge(ctrl, on=['sid', 'family', 'exclusion'], how='inner')
cal = []
for key, g in SC[SC.condition != 'clean_odd'].groupby(['family', 'exclusion', 'condition']):
    d = (g.z - g.z_ctrl).dropna()
    cal.append(dict(zip(['family', 'exclusion', 'condition'], key), n_paired=len(d),
                    median_z_minus_ctrl=d.median() if len(d) else np.nan,
                    p_5a_calibrated=stats.wilcoxon(d, alternative='greater').pvalue if len(d) >= 6 else np.nan))
T = T.merge(pd.DataFrame(cal), on=['family', 'exclusion', 'condition'], how='left')
Gs = G[['event_set', 'condition', 'exclusion', 'rotation_deg', 'rotation_lo', 'rotation_hi', 'gamma_sg', 'gamma_sg_p',
         'length_ratio', 'length_ratio_lo', 'length_ratio_hi']].rename(
    columns={'rotation_deg': 'glmm_rotation_deg', 'rotation_lo': 'glmm_rot_lo', 'rotation_hi': 'glmm_rot_hi',
             'gamma_sg': 'glmm_orthogonal_shift', 'gamma_sg_p': 'glmm_orthogonal_p', 'length_ratio': 'glmm_length_ratio',
             'length_ratio_lo': 'glmm_len_lo', 'length_ratio_hi': 'glmm_len_hi'}) \
    if 'rotation_deg' in G else pd.DataFrame(columns=['event_set', 'condition', 'exclusion'])
T = T.merge(Gs, on=['event_set', 'condition', 'exclusion'], how='left')
T.to_csv(out.table, index=False)
with pd.option_context('display.width', 280, 'display.max_columns', 40):
    print(T.round(4).to_string(index=False))

# ------------------------------------------------------------------ figures
fams = sorted(S.family.unique())
show = ['resp_event', 'snoring', 'hypopnea', 'recovery', 'clean_odd']
fig = plt.figure(figsize=(4 * len(show), 4 * len(fams) * 2))
k = 0
for fam in fams:
    for excl in ['none', 'evoked']:
        for cond in show:
            k += 1
            ax = fig.add_subplot(len(fams) * 2, len(show), k, projection='polar')
            g = S[(S.family == fam) & (S.exclusion == excl) & (S.condition == cond)]
            if len(g):
                ax.scatter(g.dphi, np.clip(g.z, -2, 6) + 2.5, s=8, alpha=.6, c=np.where(g.p < .05, 'crimson', 'k'))
                t = T[(T.family == fam) & (T.exclusion == excl) & (T.condition == cond)]
                if len(t) and np.isfinite(t.iloc[0].get('mean_rotation_deg', np.nan)):
                    r = t.iloc[0]
                    ax.annotate('', xy=(np.radians(r.mean_rotation_deg), 8 * r.rotation_R), xytext=(0, 0),
                                arrowprops=dict(color='tab:blue', width=2, headwidth=7))
                    p5a = r.p_5a_calibrated if cond != 'clean_odd' else r.p_5a_wilcoxon
                    ax.set_title(f'{fam} | {cond} | excl={excl}\nN={int(r.n_subjects)} 5a(cal) p={p5a:.2g} | rot '
                                 f'{r.mean_rotation_deg:.0f}° p(≠0)={r.p_5b_rotation_ne0:.2g}', fontsize=7)
            ax.set_yticks([])
            ax.set_ylim(0, 8.5)
            ax.set_theta_zero_location('E')
fig.suptitle('Step 5: per-subject phase shift (condition − clean). Angle = Δφ; radius = 5a z (red: subject p < .05); '
             'blue = mean rotation (length = consistency R)', fontsize=10)
fig.tight_layout()
fig.savefig(out.fig, dpi=70)
