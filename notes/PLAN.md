# Project plan and status (SnoreDREAMT)
Full approved plan: `~/.claude/plans/fuzzy-shimmying-breeze.md`. Reports: `notes/STEP3_REPORT.md`, `notes/STEP4_REPORT.md`.

| Step | Question | Status |
|---|---|---|
| 0 | Data cache, quality masks, literature | done |
| 1 | SO / spindle / complex detection (slythm port, MATLAB parity passed) | done |
| 2 | Respiratory phase + replication of Schreiner 2023 coupling; arousal–breathing coupling | done |
| 3 (3b–3d) | Coupling **strength** per breathing condition (apneas, hypopneas, snoring) | done |
| 4 | Coupling strength vs AHI and data-derived severity markers | done |
| 5 | Coupling **preferred phase** per breathing condition: any change? consistent rotation? | done (`notes/STEP5_REPORT.md`) |
| 6 | Phase changes vs AHI and severity markers | planned |

## Step 5 — Does the preferred phase change?
Same conditions, event families (large SOs × BEST phase, consensus spindles × ABDOMEN gated, pooled complexes ×
ABDOMEN gated), stage sets and snore-evoked exclusion as Step 3b/3c; whole-night clean N2 as reference. Within-subject
phase differences are polarity-invariant, so all subjects can be used.
- **5a — any phase difference (direction-free):** per subject, angle between the occupancy-normalised mean directions of
  condition and clean events; within-subject label-permutation null (1000×) → z/p per subject; group: one-sided
  Wilcoxon on z, Stouffer combination.
- **5b — consistent rotation across participants:** signed Δφ per subject → Rayleigh (any common rotation), weighted and
  reliability-filtered versions; group estimate of the rotation from the Poisson GLMM (condition × cos/sin terms) with
  simulation CI.
- Controls: clean-vs-clean pseudo-condition, REM pseudo-SOs.

## Step 6 — Phase changes vs severity
- **6a:** Δφ (circular–linear correlation) and |Δφ| / 5a-z (Spearman, robust regression adjusted for age, sex, BMI) vs
  log1p AHI, ODI3, T90, hypoxic burden, snore index, mean event duration; BH-FDR.
- **6b:** clean preferred phase vs severity (follow-up of the Step 4 PPC finding): circular–linear correlations, AHI
  groups compared (Watson–Williams / Hotelling), GLMM AHI × sin term; polarity-unambiguous subjects only.
