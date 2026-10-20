"""Snakemake rule `resp_summary`: polarity decisions, conflicts, breathing frequency and phase coverage per subject."""
import json

import pandas as pd

rows = []
for f in snakemake.input.info:  # noqa: F821
    i = json.load(open(f))
    r = {'sid': i['sid']}
    for ch in ['ABDOMEN', 'THORAX']:
        a = i[ch]
        r.update({f'{ch}_sign': a['sign'], f'{ch}_conflict': a['conflict'], f'{ch}_snore_z': a.get('snore_z'),
                  f'{ch}_rsa_z': a.get('rsa_z'), f'{ch}_ie_ratio': a.get('ie_ratio')})
    r.update(pooled_snore_z=i['ABDOMEN'].get('pooled_snore_z'), pooled_rsa_z=i['ABDOMEN'].get('pooled_rsa_z'),
             r_abd_tho=i['ABDOMEN'].get('r_abd_tho'))
    for ch in ['FLOW', 'PTAF']:
        r.update({f'{ch}_sign': i[ch]['sign'], f'{ch}_lag_s': i[ch]['lag_s'], f'{ch}_r': i[ch]['r_with_dabdomen'],
                  f'{ch}_snore_z': i[ch].get('snore_z'), f'{ch}_conflict': i[ch].get('conflict')})
    for ch in ['FLOW', 'ABDOMEN']:
        r.update({f'{ch}_f0': i[ch]['f0'], f'{ch}_cycle_defined': i[ch]['cycle_defined_frac_nrem'],
                  f'{ch}_n_breaths': i[ch]['n_breaths']})
    rows.append(r)
T = pd.DataFrame(rows).sort_values('sid')
# belts: conflicting or undecided anchors; flow: neither a strong snore anchor nor a clear belt relation
T['abdomen_ambiguous'] = T.ABDOMEN_conflict | (T.ABDOMEN_sign == 0)
T['flow_ambiguous'] = (T.FLOW_sign == 0) | ((T.FLOW_snore_z.abs() < 5) & (T.FLOW_r.abs() < 0.3))
T['needs_visual_check'] = T.abdomen_ambiguous | T.flow_ambiguous
T.to_csv(snakemake.output[0], index=False)  # noqa: F821
print(T[['ABDOMEN_sign', 'ABDOMEN_conflict', 'FLOW_sign', 'needs_visual_check']].apply(pd.Series.value_counts).to_string())
print(T[['FLOW_f0', 'ABDOMEN_f0', 'FLOW_cycle_defined', 'ABDOMEN_cycle_defined', 'FLOW_lag_s', 'FLOW_r']]
      .describe().round(3).to_string())
print('needs visual check:', T.sid[T.needs_visual_check].tolist())
