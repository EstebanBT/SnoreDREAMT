"""Snakemake rule `quality`: per-channel validity masks for one subject."""

from snoredreamt import quality
from snoredreamt.io import load_subject

cfg = snakemake.config  # noqa: F821
d = load_subject(snakemake.wildcards.sid, cfg)  # noqa: F821
masks, tables = quality.compute_masks(d, cfg)
quality.save_masks(snakemake.output[0], masks, tables, len(d['stage']))  # noqa: F821
