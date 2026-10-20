"""Load DREAMT 100 Hz PSG files into tidy arrays and cache them per subject.

Units: EEG/EOG/EMG/ECG/SNORE and respiratory channels are stored in volts in the CSVs
(Metadata.txt says µV, which is wrong); we multiply by 1e6. SAO2 is a fraction -> %.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
EVENT_COLS = {'oa': 'Obstructive_Apnea', 'ca': 'Central_Apnea', 'hyp': 'Hypopnea', 'mixed': 'Multiple_Events'}
VOLT_CHANNELS = ['C4-M1', 'F4-M1', 'O2-M1', 'Fp1-O2', 'T3 - CZ', 'CZ - T4', 'CHIN', 'E1', 'E2', 'ECG',
                 'LAT', 'RAT', 'SNORE', 'PTAF', 'FLOW', 'THORAX', 'ABDOMEN']


def load_config(path=ROOT / 'config.yaml'):
    with open(path) as f:
        return yaml.safe_load(f)


def subject_ids(cfg=None, include_excluded=False):
    cfg = cfg or load_config()
    sids = sorted(p.name[:4] for p in Path(cfg['paths']['data']).glob('S*_PSG_df.csv'))
    return sids if include_excluded else [s for s in sids if s not in cfg['exclude_subjects']]


def _cache_path(sid, cfg):
    return Path(cfg['paths']['derivatives']) / 'cache' / f'{sid}.npz'


def read_csv(sid, cfg):
    """Parse the raw CSV into a dict of numpy arrays (slow path, ~1 GB file)."""
    df = pd.read_csv(Path(cfg['paths']['data']) / f'{sid}_PSG_df.csv', engine='pyarrow',
                     usecols=['TIMESTAMP', 'SAO2', 'Sleep_Stage'] + VOLT_CHANNELS + list(EVENT_COLS.values()))
    t = df.TIMESTAMP.to_numpy()
    dt = np.diff(t)
    if not np.allclose(dt, 1 / cfg['fs'], atol=1e-6):
        raise ValueError(f'{sid}: irregular sampling (dt range {dt.min()}-{dt.max()})')
    codes = cfg['stages']
    stage = df.Sleep_Stage.map(codes).fillna(-1).to_numpy(np.int8)
    out = {'t0': np.float64(t[0]), 'stage': stage, 'SAO2': (df.SAO2.to_numpy() * 100).astype(np.float32)}
    for ch in VOLT_CHANNELS:
        out[ch] = (df[ch].to_numpy() * 1e6).astype(np.float32)
    for key, col in EVENT_COLS.items():
        out['ev_' + key] = df[col].fillna(0).to_numpy() > 0
    return out


def load_subject(sid, cfg=None, use_cache=True):
    """Return a dict: channel name -> float32 array, 'stage' (int8 codes), 'ev_*' bool masks, 't0', 'fs'."""
    cfg = cfg or load_config()
    cp = _cache_path(sid, cfg)
    if use_cache and cp.exists():
        with np.load(cp) as z:
            d = {k: z[k] for k in z.files}
    else:
        d = read_csv(sid, cfg)
        cp.parent.mkdir(parents=True, exist_ok=True)
        np.savez(cp, **d)
    d['fs'] = cfg['fs']
    d['sid'] = sid
    return d


def participants(cfg=None):
    cfg = cfg or load_config()
    p = pd.read_csv(cfg['paths']['participants'])
    p['Mean_SaO2'] = p['Mean_SaO2'].astype(str).str.rstrip('%').astype(float)
    return p.set_index('SID')
