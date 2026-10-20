"""Snakemake rule `detect_spindles`: slythm spindle detection on one channel of one subject.

NREM events use the threshold from pooled N2/N3; stages in `check_stages` are re-run with that same threshold.
"""
import json

import pandas as pd

from snoredreamt import detectors as D
from snoredreamt.io import load_subject
from snoredreamt.quality import load_masks

cfg = snakemake.config  # noqa: F821
p = cfg['detect']['spindle']
ch = cfg['channel_alias'][snakemake.wildcards.ch]  # noqa: F821
d = load_subject(snakemake.wildcards.sid, cfg)  # noqa: F821
valid = load_masks(snakemake.input.quality)[ch]  # noqa: F821

ev, thr = D.detect_spindles(d[ch], d['stage'], valid, cfg['fs'], p)
chk, _ = D.detect_spindles(d[ch], d['stage'], valid, cfg['fs'], p, detect_stages=cfg['check_stages'], threshold=thr)
ev = pd.concat([ev.assign(nrem_event=True), chk.assign(nrem_event=False)]).sort_values('start')
ev.reset_index(drop=True).to_parquet(snakemake.output.events)  # noqa: F821

minutes = {int(s): float(((d['stage'] == s) & valid).sum() / cfg['fs'] / 60) for s in [0, 1, 2, 3, 5]}
with open(snakemake.output.thresholds, 'w') as f:  # noqa: F821
    json.dump({'channel': ch, 'threshold': float(thr), 'valid_minutes': minutes}, f, indent=1)
