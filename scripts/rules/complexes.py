"""Snakemake rule `complexes`: SO-spindle complexes from the NREM SOs and spindles of one channel."""

import pandas as pd

from snoredreamt.detectors import find_complexes

cfg = snakemake.config  # noqa: F821
sos = pd.read_parquet(snakemake.input.sos)  # noqa: F821
sp = pd.read_parquet(snakemake.input.spindles)  # noqa: F821


def nrem(events):
    # channels without any detected event (e.g. an invalid derivation) are saved without the nrem_event column
    return events[events['nrem_event']] if 'nrem_event' in events.columns else events.iloc[:0]


cx = find_complexes(nrem(sos), nrem(sp), cfg['fs'], cfg['detect']['complex']['max_lag_s'], spindle_time='trough')
cx.to_parquet(snakemake.output[0])  # noqa: F821
