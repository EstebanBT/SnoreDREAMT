# SnoreDREAMT

Does the coupling of slow oscillations, sleep spindles and SO–spindle complexes to the breathing cycle change during
apneas, hypopneas and snoring, and does it scale with sleep-disordered-breathing severity?
Analysis of the DREAMT sleep-clinic cohort (100 overnight PSGs) with a Snakemake pipeline (Steps 0–6).

**Repository:** <https://github.com/EstebanBT/SnoreDREAMT>
**Report:** `reports/SnoreDREAMT_report.html` (source `reports/SnoreDREAMT_report.qmd`).
**Progress / decisions:** `notes/PLAN.md` and `notes/STEP*_REPORT.md`.

## Setup

Requires [pixi](https://pixi.sh) (Linux x86-64). The lock file pins Python, Snakemake 7.32, R + lme4, graphviz and Quarto.

```bash
pixi install
pixi run test          # unit tests (detectors, circular statistics, arousal detector, coupling)
```

### Data (not included)
DREAMT is a credentialed PhysioNet dataset (PhysioNet Restricted Health Data License 1.5.0). Request access on
PhysioNet, download it, and set `paths.data` and `paths.participants` in `config.yaml` to the `data_100Hz` folder and
`participant_info.csv`. Verify the download with the dataset's `SHA256SUMS.txt`.

### Third-party arousal detector
The arousal detector is a Python port of the GPL-3 MATLAB toolbox by Fernández-Varela / Alvarez-Estevez. The MATLAB
original is only needed for comparison:

```bash
git clone https://github.com/bigsasi/arousals-detection third_party/arousals-detection
git -C third_party/arousals-detection checkout 7034ef739e0d34e5a8b1e0bed0f6ec8d2ab7624b
```

## Reproduce

```bash
pixi run snakemake -n      # dry run: what would be built and why
pixi run pipeline          # snakemake -j 8: cache → quality → detection → phase → Steps 2–6
pixi run dag               # figures/pipeline_dag.png
pixi run report            # renders reports/SnoreDREAMT_report.html
```

A run from scratch is long (98 subjects; per-subject detection and permutation nulls); `-j` sets the number of cores.
Step outputs land in `derivatives/` (large intermediate files are git-ignored), figures in `figures/`.

## Layout

| Path | Content |
|---|---|
| `Snakefile`, `config.yaml` | pipeline (one rule per step) and all parameters (one primary analysis per step) |
| `snoredreamt/` | library: I/O, quality masks, detectors (slythm port), respiration phase, conditions, snoring, arousals, coupling, CircStat port |
| `scripts/rules/` | one script per Snakemake rule (Python and R/lme4) |
| `tests/` | unit tests |
| `validation/` | MATLAB parity scripts (slythm detectors, CircStat) |
| `notes/` | plan, progress and per-step reports |
| `literature/` | structured notes of the key papers (PDFs not included) |
| `reports/` | Quarto report |
| `CLAUDE.md`, `.claude/skills/` | instructions and Skills for coding agents |

## Validation
- Detectors: MATLAB parity with the original `+slythm` code (`validation/matlab_parity/run_parity.m`, then `compare.py`).
- CircStat port: `validation/circstat_parity/run_parity.m`, then `compare.py`.
