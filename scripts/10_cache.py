"""Step 0: parse every raw CSV once into derivatives/cache/SXXX.npz (float32, µV)."""
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from snoredreamt.io import load_config, load_subject, subject_ids  # noqa: E402


def run(sid):
    try:
        d = load_subject(sid)
        return sid, len(d['stage']), None
    except Exception as e:  # report and continue; one broken file must not stop the batch
        return sid, 0, repr(e)


if __name__ == '__main__':
    cfg = load_config()
    sids = subject_ids(cfg, include_excluded=True)
    with Pool(8) as p:
        for sid, n, err in p.imap_unordered(run, sids):
            print(sid, f'{n / cfg["fs"] / 3600:.2f} h' if n else '', err or 'ok', flush=True)
