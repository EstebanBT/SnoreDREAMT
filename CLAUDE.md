# SnoreDREAMT — instructions for coding agents

Respiration coupling of slow oscillations (SOs), spindles and SO–spindle complexes during apnea, hypopnea and
snoring in the DREAMT sleep-clinic cohort. Owner: Esteban Bullón Tarrasó (sleep/memory neuroscientist, MATLAB and
FieldTrip background; the project is in Python + R).

## Start here
- **Progress file:** `notes/PLAN.md` (step table and status). Per-step results: `notes/STEP3_REPORT.md` …
  `notes/STEP6_REPORT.md`; the decision log is in `notes/STEP3_GOAL.md`. Read these before changing anything and
  update them when a step finishes.
- **Report:** `reports/SnoreDREAMT_report.qmd` (render with `pixi run report`).

## Environment and commands
- `pixi install` once; then `pixi run test`, `pixi run pipeline` (= `snakemake -j 8`), `pixi run dag`,
  `pixi run report`. Use `pixi run snakemake -n` before any real run and read the reasons it gives.
- In the Claude Code sandbox, caches must point to `$TMPDIR`: `MPLCONFIGDIR`, `XDG_CACHE_HOME`, `PIXI_CACHE_DIR`.
- MATLAB has no license inside the sandbox. Parity checks against MATLAB (`validation/*/run_parity.m`) are run by the
  owner; prepare the script and the comparison (`compare.py`), do not claim parity without his output.

## Data
- DREAMT (PhysioNet, credentialed, **read-only**) at the path in `config.yaml` → `paths.data`. Never copy raw data
  into the repository; integrity is given by the dataset's `SHA256SUMS.txt`.
- EEG is stored in volts (converted to µV in the cache). S097 and S077 are excluded (broken recordings).
- Arousals and snoring are not labelled in DREAMT; they are detected (`snoredreamt/arousals.py`, `snoredreamt/snoring.py`).

## Conventions
- **One Snakemake rule per analysis step** (`Snakefile`, scripts in `scripts/rules/`). No analysis results from ad hoc
  scripts end up in a report without a rule.
- **All parameters live in `config.yaml`.** Each step names one primary analysis; everything else is labelled
  sensitivity. Do not change a primary after seeing the condition contrasts.
- Derived files go to `derivatives/`, figures to `figures/`; never edit them by hand.
- Library code in `snoredreamt/` gets unit tests in `tests/` (`pixi run test` must pass before committing).
- Commit at milestones only when the owner asks; end commit messages with the co-author line he uses.

## Scientific guardrails (learned the hard way in this project)
- **Never invent reference values.** Test expectations come from formulas, simulations with known ground truth, or
  real MATLAB/CircStat output — not from memory. Same for citations: no guessed DOIs, volumes or pages.
- **Every new test needs a negative control that shows it is calibrated** (clean-vs-clean split, `clean_odd` vs
  `clean_even`, REM pseudo-SOs, random-time surrogates). If the control is "significant", the test is wrong.
- **Phase occupancy:** compare event phases only with occupancy-normalised measures (rate per phase bin, Poisson GLMM
  with occupancy offset). Excluding phase windows (e.g. snore-evoked events) makes per-event mean-cos invalid.
- **Respiratory polarity differs between subjects:** group phase analyses only on polarity-unambiguous subjects;
  within-subject contrasts are polarity-free.
- **No circularity:** a reference phase must not be estimated from the same events it is compared with (odd/even split).
- Rayleigh on phase *differences* does not test rotation; Watson–Williams is invalid at the low concentrations seen
  here (use the common median test / GLMM orthogonal component). Use rank regression rather than Huber when leverage
  points are possible.
- Report N per comparison and correct for multiple comparisons (BH-FDR within a family).

## Skills
- `literature-prep` (`.claude/skills/literature-prep/`): rename PDFs to FirstAuthor_Year, convert with docling, write
  structured notes to `literature/notes/` and update `literature/index.json`. Use it before relying on a paper's methods.
