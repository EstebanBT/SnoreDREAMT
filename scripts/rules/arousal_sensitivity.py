"""Snakemake rule `arousal_sensitivity`: do the Step 3b results depend on missed arousals?

Compares the primary Step 3b/3c results with the same analysis run with the extended arousal exclusion (primary
arousals + recovered no-EMG NREM alpha candidates + candidates that start in sleep and end in wake), for the primary
event families: occupancy-normalised contrasts (condition - clean) and the GLMM relative change of the modulation.
Also reports how much extra time the extended exclusion removes per condition.
"""
import json

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

inp, out = snakemake.input, snakemake.output  # noqa: F821
plt.rcParams.update({'font.size': 12, 'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 200, 'savefig.bbox': 'tight'})
FAM = {'SO': 'Large SOs', 'spindle': 'Consensus spindles', 'complex': 'Pooled complexes'}
CONDS = ['resp_event', 'hypopnea', 'recovery', 'snoring']
LAB = {'resp_event': 'apneas + hypopneas', 'hypopnea': 'hypopnea', 'recovery': 'recovery', 'snoring': 'snoring'}

# extra excluded time (N2 minutes per condition, summed over subjects)
mins = []
for variant, files in [('primary', inp.info_primary), ('extended', inp.info_extended)]:
    for f in files:
        j = json.load(open(f))
        for key, m in j['minutes'].items():
            if key.endswith('_N2'):
                mins.append(dict(variant=variant, sid=j['sid'], condition=key[:-3], minutes=m))
M = pd.DataFrame(mins).groupby(['variant', 'condition']).minutes.sum().unstack(0)
M['extended_minus_primary_pct'] = 100 * (M.extended - M.primary) / M.primary

# The event set / phase source of each family is held fixed at the primary analysis's choice. (step3b_group.py
# re-selects the 'primary' set from the clean data; with the extended exclusion it would pick a different spindle set,
# which would make the comparison not like-for-like.)
CS = pd.read_csv(inp.clean_primary)
prim = {r.family: (r.event_set, r.variant) for r in CS[CS.primary].itertuples()}
print('primary sets (held fixed):', prim)
rows = []
for variant, k_file, g_file in [('primary', inp.contrasts_primary, inp.glmm_primary),
                                ('extended', inp.contrasts_extended, inp.glmm_extended)]:
    K = pd.read_csv(k_file)
    G = pd.read_csv(g_file)
    for fam, (es, var) in prim.items():
        for cond in CONDS:
            for excl in ['none', 'evoked']:
                k = K[(K.event_set == es) & (K.variant == var) & (K.condition == cond) & (K.exclusion == excl)]
                g = G[(G.event_set == es) & (G.variant == var) & (G.condition == cond) & (G.exclusion == excl)]
                if not len(k):
                    continue
                k, g = k.iloc[0], (g.iloc[0] if len(g) else None)
                rows.append(dict(arousal_exclusion=variant, family=fam, event_set=es, phase=var, condition=cond,
                                 snore_evoked_excluded=excl == 'evoked', n_subjects=k.n_subjects,
                                 clean_level=k.clean_level, delta=k.delta_mean, ci_lo=k.ci_lo, ci_hi=k.ci_hi,
                                 p_loss=k.p_loss,
                                 glmm_rel_change=np.nan if g is None else g.rel_change_median,
                                 glmm_rel_lo=np.nan if g is None else g.rel_change_lo,
                                 glmm_rel_hi=np.nan if g is None else g.rel_change_hi,
                                 glmm_p_loss=np.nan if g is None else g.gamma_cg_p_loss,
                                 glmm_p_complete_loss=np.nan if g is None else g.p_complete_loss))
R = pd.DataFrame(rows)
W = R.pivot_table(index=['family', 'condition', 'snore_evoked_excluded'], columns='arousal_exclusion',
                  values=['delta', 'p_loss', 'glmm_rel_change', 'glmm_p_loss'])
W.columns = [f'{a}_{b}' for a, b in W.columns]
W['delta_shift'] = W.delta_extended - W.delta_primary
R.to_csv(out.table, index=False)
W.reset_index().to_csv(out.comparison, index=False)
M.reset_index().to_csv(out.minutes, index=False)
print(M.round(1).to_string())
with pd.option_context('display.width', 200):
    print(W.round(3).to_string())

fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))
for ax, fam in zip(axs, FAM):
    for j, cond in enumerate(CONDS):
        for k_, (variant, col, dy) in enumerate([('primary', '#222222', .2), ('extended', '#E69F00', -.2)]):
            for excl, mk in [(False, 'o'), (True, 's')]:
                r = R[(R.arousal_exclusion == variant) & (R.family == fam) & (R.condition == cond)
                      & (R.snore_evoked_excluded == excl)]
                if not len(r) or (excl and cond != 'snoring'):
                    continue
                r = r.iloc[0]
                y = j + dy + (-.1 if excl else 0)
                ax.errorbar(r.delta, y, xerr=[[r.delta - r.ci_lo], [r.ci_hi - r.delta]], fmt=mk, color=col, ms=6,
                            mfc=col if not excl else 'white', capsize=3,
                            label=f'{variant} arousal exclusion' + (', snore-evoked excluded' if excl else '')
                            if (fam == 'SO' and j == 3) else None)
    clean = R[(R.family == fam) & (R.arousal_exclusion == 'primary')].clean_level.median()
    ax.axvline(0, color='grey', lw=1)
    ax.axvline(-clean, color='k', ls=':', lw=1.2)
    ax.set_yticks(range(len(CONDS)), [LAB[c] for c in CONDS] if fam == 'SO' else [''] * len(CONDS))
    ax.invert_yaxis()
    ax.set_title(FAM[fam])
    ax.set_xlabel('coupling: condition − clean')
h, lab = axs[0].get_legend_handles_labels()
fig.legend(h, lab, loc='lower center', ncol=2, bbox_to_anchor=(.5, -.16), frameon=False)
fig.suptitle('Arousal sensitivity: Step 3b/3c contrasts with the primary vs the extended arousal exclusion (N2)',
             y=1.02)
fig.savefig(out.fig)
