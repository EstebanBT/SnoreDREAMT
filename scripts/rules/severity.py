"""Snakemake rule `severity`: sleep-disordered-breathing severity markers for one subject, computed from the data.

SaO2 (finger probe; DREAMT 'SAO2') cleaned: values < 50 % or deviating > 10 % from a 5 s running median are set to
NaN, then averaged to 1 Hz. Sleep = scored N1/N2/N3/REM.
- ODI3: desaturations >= 3 % below the running maximum of the preceding 120 s, lasting >= 10 s (gaps < 10 s merged),
  per hour of sleep.
- T90: % of sleep time with SaO2 < 90 %.
- Hypoxic burden (simplified after Azarbarzin 2019): for each scored apnea/hypopnea, area of (baseline - SaO2)^+
  from event onset to 60 s after event end, baseline = max SaO2 in the 100 s before onset; summed (%·min) per hour
  of sleep.
- event index (scored events per hour of sleep, check vs clinical AHI), mean event duration, detected arousal index,
  snore index (% of event-free NREM breaths snoring; % of N2 time labelled snoring).
"""
import json

import numpy as np
import pandas as pd
from scipy import ndimage

from snoredreamt import conditions
from snoredreamt.io import load_subject

cfg = snakemake.config  # noqa: F821
inp, out = snakemake.input, snakemake.output  # noqa: F821
sid = snakemake.wildcards.sid  # noqa: F821
fs = cfg['fs']
d = load_subject(sid, cfg)
st = d['stage']
n_s = len(st) // fs
sleep_s = np.isin(st[:n_s * fs:fs], [1, 2, 3, 5])
sleep_h = sleep_s.sum() / 3600

x = d['SAO2'].astype(float)
x[x < 50] = np.nan
med = ndimage.median_filter(np.nan_to_num(x, nan=np.nanmedian(x)), size=5 * fs)
x[np.abs(x - med) > 10] = np.nan
s1 = np.nanmean(x[:n_s * fs].reshape(n_s, fs), 1)                       # 1 Hz
valid = np.isfinite(s1)

# ODI3
base = pd.Series(s1).rolling(120, min_periods=30).max().shift(1).to_numpy()
desat = (base - s1 >= 3) & valid & sleep_s
r_s, r_e = conditions.runs(desat)
if len(r_s):
    keep_s, keep_e = [r_s[0]], [r_e[0]]
    for a, b in zip(r_s[1:], r_e[1:]):
        if a - keep_e[-1] < 10:
            keep_e[-1] = b
        else:
            keep_s.append(a)
            keep_e.append(b)
    dur = np.array(keep_e) - np.array(keep_s) + 1
    n_desat = int((dur >= 10).sum())
else:
    n_desat = 0
odi3 = n_desat / sleep_h if sleep_h > 0 else np.nan
t90 = 100 * np.mean(s1[sleep_s & valid] < 90) if (sleep_s & valid).any() else np.nan

# hypoxic burden
ev = conditions.respiratory_events(d)
es, ee = conditions.runs(ev)
es, ee = es // fs, ee // fs
area = 0.0
for a, b in zip(es, ee):
    if not sleep_s[min(a, n_s - 1)]:
        continue
    pre = s1[max(a - 100, 0):a]
    if not np.isfinite(pre).any():
        continue
    bl = np.nanmax(pre)
    seg = s1[a:min(b + 60, n_s)]
    area += np.nansum(np.clip(bl - seg, 0, None)) / 60                   # %·min
hb = area / sleep_h if sleep_h > 0 else np.nan

ev_sleep = [(a, b) for a, b in zip(es, ee) if sleep_s[min(a, n_s - 1)]]
ar = pd.read_parquet(inp.arousals)
cinfo = json.load(open(inp.cond_info))
lab = np.load(inp.labels)['labels']
n2 = st == 2
out_d = dict(sid=sid, sleep_h=sleep_h, odi3=odi3, t90=t90, hypoxic_burden=hb,
             event_index=len(ev_sleep) / sleep_h if sleep_h else np.nan,
             mean_event_duration_s=float(np.mean([b - a + 1 for a, b in ev_sleep])) if ev_sleep else np.nan,
             arousal_index_detected=len(ar) / sleep_h if sleep_h else np.nan,
             snore_breath_frac=cinfo['snore_breath_frac_eventfree_nrem'],
             snore_time_frac_n2=float(((lab == 1) & n2).sum() / max(((lab >= 0) & n2).sum(), 1)),
             mean_sao2=float(np.nanmean(s1[sleep_s])), sao2_valid_frac=float(valid[sleep_s].mean()))
with open(out[0], 'w') as f:
    json.dump(out_d, f, indent=1, default=float)
