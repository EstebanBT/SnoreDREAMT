"""Snakemake rule `complexes`: SO-spindle complexes from the NREM SOs and spindles of one channel."""

import pandas as pd

from snoredreamt.detectors import find_complexes

cfg = snakemake.config  # noqa: F821
sos = pd.read_parquet(snakemake.input.sos)  # noqa: F821
sp = pd.read_parquet(snakemake.input.spindles)  # noqa: F821
cx = find_complexes(sos[sos.nrem_event], sp[sp.nrem_event], cfg['fs'], cfg['detect']['complex']['max_lag_s'],
                    spindle_time='trough')
cx.to_parquet(snakemake.output[0])  # noqa: F821
