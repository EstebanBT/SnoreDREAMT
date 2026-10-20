"""Snakemake rule `detect_sos`: slythm SO detection on one channel of one subject.

NREM events use per-stage thresholds (N2/N3); stages in `check_stages` (W, N1, REM) are re-run with the N2
threshold so their 'events' measure the detector's false-positive floor.
"""
import json

import pandas as pd

from snoredreamt import detectors as D
from snoredreamt.io import load_subject
from snoredreamt.quality import load_masks

cfg = snakemake.config  # noqa: F821
p = cfg['detect']['so']
ch = cfg['channel_alias'][snakemake.wildcards.ch]  # noqa: F821
d = load_subject(snakemake.wildcards.sid, cfg)  # noqa: F821
valid = load_masks(snakemake.input.quality)[ch]  # noqa: F821

ev, thr = D.detect_sos(d[ch], d['stage'], valid, cfg['fs'], p)
n2 = thr.get(p['stages'][0])
if n2 is not None:
    chk, _ = D.detect_sos(d[ch], d['stage'], valid, cfg['fs'], p, detect_stages=cfg['check_stages'],
                          thresholds={s: n2 for s in cfg['check_stages']})
    ev = pd.concat([ev.assign(nrem_event=True), chk.assign(nrem_event=False)]).sort_values('start')
ev.reset_index(drop=True).to_parquet(snakemake.output.events)  # noqa: F821

minutes = {int(s): float(((d['stage'] == s) & valid).sum() / cfg['fs'] / 60) for s in [0, 1, 2, 3, 5]}
with open(snakemake.output.thresholds, 'w') as f:  # noqa: F821
    json.dump({'channel': ch, 'thresholds': {int(k): list(map(float, v)) for k, v in thr.items()},
               'valid_minutes': minutes}, f, indent=1)
