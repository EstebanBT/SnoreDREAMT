"""Snakemake rule `conditions`: per-sample breathing-condition labels and per-breath snoring labels for one subject.

Label codes (int8), assigned with this priority (first match wins):
  3 obstructive apnea, 4 central apnea, 5 mixed (Multiple_Events), 2 hypopnea,
  6 recovery (event end -> + post_event_s), 7 pre-event margin (pre_event_s before onset; excluded),
  8 arousal window (detected arousal +- pad; excluded), 1 snoring (>= min_run consecutive snoring breaths),
  0 clean (quiet, event-free breathing), -1 not analysed (stage not N2/N3 or invalid ABDOMEN signal).
Stage is kept separately so N2-only and N2+N3 analyses can be run.
"""
import json

import numpy as np
import pandas as pd

from snoredreamt import conditions, snoring
from snoredreamt.io import load_subject
from snoredreamt.quality import load_masks

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
sid = snakemake.wildcards.sid  # noqa: F821
fs, c, sc = cfg['fs'], cfg['clean'], cfg['snoring']
variant = snakemake.wildcards.arv  # noqa: F821  arousal detector variant used for exclusion

d = load_subject(sid, cfg)
st = d['stage']
n = len(st)
masks = load_masks(inp.quality)
ar = pd.read_parquet(inp.arousals)
b = pd.read_parquet(inp.breaths)
b = b[b.channel == 'ABDOMEN']

# snoring per breath
ev_any = conditions.respiratory_events(d)
br = snoring.breath_features(d['SNORE'], b.peak.to_numpy(), b.trough.to_numpy(), fs, band=sc['band'], ptaf=d['PTAF'])
br = snoring.label_snoring(br, st, np.isin(st, c['stages']) & ~ev_any, st == 0, sc['min_db'], sc['min_phasic'],
                           sc['min_run'])
snore = np.zeros(n, bool)
sb = br[br.snore]
snore = conditions.paint(n, sb.start.to_numpy(), sb.end.to_numpy() - 1, 0, 0)

lab = np.zeros(n, np.int8)                       # clean by default
lab[snore] = 1
lab[conditions.arousal_exclusion(n, ar, fs, cfg['arousals']['exclude_pad_s'])] = 8
s, e = conditions.runs(ev_any)
lab[conditions.paint(n, s, s, int(c['pre_event_s'] * fs), -1)] = 7          # [onset - pre, onset - 1]
lab[conditions.paint(n, e, e, -1, int(c['post_event_s'] * fs))] = 6          # [end + 1, end + post]
for code, key in [(2, 'ev_hyp'), (5, 'ev_mixed'), (4, 'ev_ca'), (3, 'ev_oa')]:
    lab[d[key]] = code
lab[~np.isin(st, [2, 3]) | ~masks['ABDOMEN']] = -1

np.savez_compressed(out.labels, labels=lab, stage=st)
br.to_parquet(out.breaths)
names = {0: 'clean', 1: 'snoring', 2: 'hypopnea', 3: 'obstructive', 4: 'central', 5: 'mixed', 6: 'recovery',
         7: 'pre_event', 8: 'arousal_window'}
minutes = {f'{nm}_{stg}': float(((lab == k) & (st == sv)).sum() / fs / 60)
           for k, nm in names.items() for stg, sv in [('N2', 2), ('N3', 3)]}
nrem_free = np.isin(st[br.start.to_numpy()], [2, 3]) & ~ev_any[br.start.to_numpy()]
with open(out.info, 'w') as f:
    json.dump({'sid': sid, 'arousal_variant': variant, 'minutes': minutes,
               'n_events': {k: int(len(conditions.runs(d[k])[0])) for k in ['ev_oa', 'ev_ca', 'ev_hyp', 'ev_mixed']},
               'snore_breath_frac_eventfree_nrem': float(br.snore[nrem_free].mean()) if nrem_free.any() else None,
               'quiet_ref': br.attrs['quiet_ref']}, f, indent=1)
