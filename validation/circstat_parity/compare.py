"""Compare snoredreamt.circstat with CircStat 2012a on fixed inputs. 1) python compare.py (writes inputs.mat if
missing); 2) matlab -batch "run('validation/circstat_parity/run_parity.m')"; 3) python compare.py again."""
from pathlib import Path

import numpy as np
import scipy.io as sio

from snoredreamt import circstat as cs

HERE = Path(__file__).parent
if not (HERE / 'inputs.mat').exists():
    rng = np.random.default_rng(42)
    sio.savemat(HERE / 'inputs.mat', {'a1': rng.vonmises(0.3, 1.5, 61)[:, None], 'a2': rng.vonmises(1.0, 1.5, 48)[:, None]})
z = sio.loadmat(HERE / 'inputs.mat')
a1, a2 = z['a1'].ravel(), z['a2'].ravel()
py = {'p_ww': cs.circ_wwtest(a1, a2)[0], 'F_ww': cs.circ_wwtest(a1, a2)[1], 'p_cm': cs.circ_cmtest(a1, a2)[0],
      'P_cm': cs.circ_cmtest(a1, a2)[2], 'med1': cs.circ_median(a1), 'h_mt': cs.circ_mtest(a1, 0)[0],
      'mu_mt': cs.circ_mtest(a1, 0)[1], 'ul_mt': cs.circ_mtest(a1, 0)[2], 'll_mt': cs.circ_mtest(a1, 0)[3],
      'p_med': cs.circ_medtest(a1, 0)}
if not (HERE / 'matlab_out.mat').exists():
    print('python values:', {k: float(v) for k, v in py.items()}, '\nRun run_parity.m in MATLAB, then compare again.')
else:
    ml = sio.loadmat(HERE / 'matlab_out.mat')
    ok = True
    for k, v in py.items():
        m = float(np.squeeze(ml[k]))
        same = np.isclose(float(v), m, atol=1e-8) or (k == 'med1' and np.isclose(np.angle(np.exp(1j * (v - m))), 0, atol=1e-8))
        ok &= same
        print(f'{k:6s} python {float(v): .10f}  matlab {m: .10f}  {"OK" if same else "DIFF"}')
    print('PASS' if ok else 'FAIL')
