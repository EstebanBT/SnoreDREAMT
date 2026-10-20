"""Snakemake rule `arousals`: arousal detection for one subject and one parameter variant."""
from snoredreamt.arousals import detect_arousals
from snoredreamt.ecg import detect_rpeaks
from snoredreamt.io import load_subject
from snoredreamt.quality import load_masks

cfg = snakemake.config  # noqa: F821
a = cfg['arousals']
d = load_subject(snakemake.wildcards.sid, cfg)  # noqa: F821
valid = load_masks(snakemake.input.quality)[a['eeg']]  # noqa: F821
rp, _ = detect_rpeaks(d['ECG'], cfg['fs'], band=cfg['ecg']['band'], min_rr_s=cfg['ecg']['min_rr_s'])
ev = detect_arousals(d[a['eeg']], d[a['emg']], d['ECG'], d['stage'], cfg['fs'], valid=valid, rpeaks=rp,
                     params=a['variants'][snakemake.wildcards.variant] or None)  # noqa: F821
ev.to_parquet(snakemake.output[0])  # noqa: F821
