"""Snakemake rule `arousals_extended`: sensitivity variant of the arousal exclusion.

Runs the primary detector with all rejected candidates returned and keeps, in addition to the primary arousals,
the candidates most likely to be missed arousals (adversarial review of the detector, 2026-11-05):
  - alpha candidates rejected only for lacking chin-EMG activity, with onset in N2/N3 (AASM requires an EMG
    increase only for REM arousals);
  - candidates rejected by the awake rule although their onset is in sleep (they end in a wake-scored epoch).
Column `source` says why each row is included. All rows are used as exclusion windows downstream.
"""
from snoredreamt.arousals import detect_arousals
from snoredreamt.ecg import detect_rpeaks
from snoredreamt.io import load_subject
from snoredreamt.quality import load_masks

cfg = snakemake.config  # noqa: F821
a = cfg['arousals']
x = a['extended']
d = load_subject(snakemake.wildcards.sid, cfg)  # noqa: F821
valid = load_masks(snakemake.input.quality)[a['eeg']]  # noqa: F821
rp, _ = detect_rpeaks(d['ECG'], cfg['fs'], band=cfg['ecg']['band'], min_rr_s=cfg['ecg']['min_rr_s'])
ev = detect_arousals(d[a['eeg']], d[a['emg']], d['ECG'], d['stage'], cfg['fs'], valid=valid, rpeaks=rp,
                     params=a['variants'][x['base']] or None, return_all=True)
no_emg = (ev.reason == 'no EMG') & ev.stage.isin(x['no_emg_stages'])
awake = (ev.reason == 'awake') & ev.stage.isin(x['awake_onset_stages'])
ev['source'] = 'primary'
ev.loc[no_emg, 'source'] = 'recovered: no EMG'
ev.loc[awake, 'source'] = 'recovered: ends in wake'
ev = ev[ev.ok | no_emg | awake].sort_values('start').reset_index(drop=True)
ev.to_parquet(snakemake.output[0])  # noqa: F821
