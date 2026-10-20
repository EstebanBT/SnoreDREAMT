"""Compare the Python slythm port with the original MATLAB detectors on the same 2-h segment (S003).

1. In MATLAB (needs the license server):  matlab -batch "run('validation/matlab_parity/run_parity.m')"
2. Then:                                    .venv/bin/python validation/matlab_parity/compare.py
Pass criterion (plan): >= 95 % of events matched within +-2 samples.
"""
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from snoredreamt import detectors as D  # noqa: E402
from snoredreamt.io import load_config  # noqa: E402

HERE = Path(__file__).parent
cfg = load_config()
seg = sio.loadmat(HERE / 'seg.mat')
ml = sio.loadmat(HERE / 'matlab_out.mat')
X, scoring = seg['X'], seg['scoring'].ravel().astype(int)
valid = np.ones(X.shape[1], bool)
p_so = dict(cfg['detect']['so'], min_stage_minutes=0)   # MATLAB splitStages: every present stage gets its own threshold


def match(a, b, tol=2):
    """Fraction of a with a partner in b within tol samples (and vice versa)."""
    a, b = np.sort(a), np.sort(b)
    if len(a) == 0 or len(b) == 0:
        return 0.0, 0.0
    ia = np.clip(np.searchsorted(b, a), 1, len(b) - 1)
    da = np.minimum(np.abs(b[ia] - a), np.abs(b[ia - 1] - a))
    ib = np.clip(np.searchsorted(a, b), 1, len(a) - 1)
    db = np.minimum(np.abs(a[ib] - b), np.abs(a[ib - 1] - b))
    return (da <= tol).mean(), (db <= tol).mean()


ok = True
for c, ch in enumerate(['F4-M1', 'C4-M1'], start=1):
    so, thr = D.detect_sos(X[c - 1], scoring, valid, cfg['fs'], p_so)
    sp, sthr = D.detect_spindles(X[c - 1], scoring, valid, cfg['fs'], cfg['detect']['spindle_slythm'])
    m_so = ml[f'so_tr_{c}'].ravel() - 1
    m_sp = ml[f'sp_sta_{c}'].ravel() - 1
    f1, f2 = match(so.trough.to_numpy(), m_so)
    g1, g2 = match(sp.start.to_numpy(), m_sp)
    print(f'{ch}: SO python {len(so)} / matlab {len(m_so)}; matched {f1:.3f} / {f2:.3f}')
    print(f'      thresholds python {thr} | matlab trough {ml["so_thr"][:, c - 1]}, pk2pk {ml["so_thr2"][:, c - 1]}')
    print(f'{ch}: spindles python {len(sp)} / matlab {len(m_sp)}; matched {g1:.3f} / {g2:.3f}; '
          f'threshold python {sthr:.3f} | matlab {ml["sp_thr"].ravel()[c - 1]:.3f}')
    ok &= min(f1, f2, g1, g2) >= 0.95
print('PASS' if ok else 'FAIL')
