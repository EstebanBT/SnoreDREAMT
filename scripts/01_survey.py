"""Survey all DREAMT 100 Hz nights: durations, staging, events, signal quality, snore-channel envelope per 30 s epoch."""
import sys, glob, os
import numpy as np, pandas as pd
from scipy import signal
from multiprocessing import Pool

DATA = '/media/longterm1/esteban/DREAMT/data_100Hz'
OUT = '/home/esteban/Documents/SnoreDREAMT/derivatives/survey'
FS = 100
COLS = ['TIMESTAMP', 'C4-M1', 'F4-M1', 'SNORE', 'PTAF', 'FLOW', 'THORAX', 'ABDOMEN', 'Sleep_Stage',
        'Obstructive_Apnea', 'Central_Apnea', 'Hypopnea', 'Multiple_Events']


def band_rms(x, lo, hi, win):
    sos = signal.butter(4, [lo, hi], btype='band', fs=FS, output='sos')
    y = signal.sosfiltfilt(sos, x)
    n = len(y) // win
    return np.sqrt((y[:n * win].reshape(n, win) ** 2).mean(1))


def run(f):
    sid = os.path.basename(f)[:4]
    with open(f) as fh:
        hdr = fh.readline().strip().split(',')
    df = pd.read_csv(f, engine='pyarrow', usecols=[c for c in COLS if c in hdr])
    n_ep = len(df) // (30 * FS)
    df = df.iloc[:n_ep * 30 * FS]
    ep = lambda a: a.reshape(n_ep, 30 * FS)
    stage = ep(df.Sleep_Stage.astype(str).values)[:, 0]
    ev = np.zeros(len(df), bool)
    for c in ['Obstructive_Apnea', 'Central_Apnea', 'Hypopnea', 'Multiple_Events']:
        if c in df:
            ev |= df[c].fillna(0).values > 0
    rows = dict(epoch=np.arange(n_ep), t0=ep(df.TIMESTAMP.values)[:, 0], stage=stage, resp_event=ep(ev).mean(1))
    # snore: 1 s band RMS in 20-45 Hz, summarised per epoch
    for name, (lo, hi) in {'snore_hf': (20, 45), 'snore_lf': (5, 20)}.items():
        r = band_rms(np.nan_to_num(df.SNORE.values), lo, hi, FS)
        r = r[:n_ep * 30].reshape(n_ep, 30)
        rows[name + '_med'] = np.median(r, 1)
        rows[name + '_p90'] = np.percentile(r, 90, 1)
    # PTAF high-frequency flutter (classic nasal-pressure snoring marker)
    r = band_rms(np.nan_to_num(df.PTAF.values), 20, 45, FS)[:n_ep * 30].reshape(n_ep, 30)
    rows['ptaf_hf_med'] = np.median(r, 1)
    for c in ['PTAF', 'THORAX', 'ABDOMEN', 'FLOW', 'C4-M1', 'F4-M1']:
        x = df[c].values
        rows[c + '_clip'] = ep(np.abs(x) >= 0.99 * np.nanmax(np.abs(x))).mean(1)
        rows[c + '_sd'] = np.nanstd(ep(x), 1)
    E = pd.DataFrame(rows)
    E.to_parquet(f'{OUT}/{sid}_epochs.parquet')
    # resp peak frequency during N2/N3 from THORAX and FLOW
    nrem = np.repeat(np.isin(stage, ['N2', 'N3']), 30 * FS)
    summ = dict(sid=sid, hours=len(df) / FS / 3600, n_ep=n_ep)
    for s in ['W', 'N1', 'N2', 'N3', 'R', 'Missing', 'P']:
        summ['min_' + s] = (stage == s).sum() / 2
    summ['resp_event_frac_nrem'] = E.resp_event[np.isin(stage, ['N2', 'N3'])].mean()
    for c in ['THORAX', 'FLOW', 'PTAF']:
        x = df[c].values[nrem]
        if len(x) > 600 * FS:
            fr, p = signal.welch(np.nan_to_num(x), FS, nperseg=60 * FS)
            m = (fr > 0.1) & (fr < 0.6)
            summ[c + '_peakHz'] = fr[m][np.argmax(p[m])]
            summ[c + '_peakness'] = p[m].max() / np.median(p[(fr > 0.05) & (fr < 2)])
    for c in ['PTAF', 'THORAX', 'ABDOMEN', 'FLOW']:
        summ[c + '_clipfrac'] = E[c + '_clip'].mean()
    return summ


if __name__ == '__main__':
    files = sorted(glob.glob(f'{DATA}/S*_PSG_df.csv'))
    with Pool(10) as p:
        res = p.map(run, files)
    S = pd.DataFrame(res)
    S.to_csv(f'{OUT}/summary.csv', index=False)
    print(S.round(3).to_string())
