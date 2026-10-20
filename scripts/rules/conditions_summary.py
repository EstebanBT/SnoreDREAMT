"""Snakemake rule `conditions_summary`: feasibility table (minutes per condition and stage per subject) + figure."""
import json

import matplotlib
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

rows = []
for f in snakemake.input.info:  # noqa: F821
    i = json.load(open(f))
    rows.append({'sid': i['sid'], 'snore_breath_frac': i['snore_breath_frac_eventfree_nrem'], **i['minutes'],
                 **{f'n_{k[3:]}': v for k, v in i['n_events'].items()}})
T = pd.DataFrame(rows).sort_values('sid')
T.to_csv(snakemake.output.table, index=False)  # noqa: F821
conds = ['clean', 'snoring', 'hypopnea', 'obstructive', 'central', 'mixed', 'recovery']
feas = pd.DataFrame({c: {f'>= {m} min N2': int((T[f'{c}_N2'] >= m).sum()) for m in (1, 5, 10, 20, 60)} for c in conds})
print(feas.to_string())
fig, ax = plt.subplots(figsize=(10, 5))
ax.boxplot([T[f'{c}_N2'] for c in conds], labels=conds, showfliers=True)
ax.set_yscale('symlog', linthresh=1)
ax.set_ylabel('minutes in N2 (per subject)')
ax.set_title(f'Time per breathing condition in N2 (n = {len(T)} subjects)')
fig.tight_layout()
fig.savefig(snakemake.output.fig, dpi=80)  # noqa: F821
