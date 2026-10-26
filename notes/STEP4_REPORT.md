# Step 4 report — sleep-disordered-breathing severity and respiration–oscillation coupling (2026-10-26)
Rules `severity`, `step4_group`, `step4_glmm` (R/lme4); numbers in `derivatives/step4/`; figure `step4_severity.png`.

## Design
- **Coupling per subject** (from Step 3b/3c, clean N2, all events): occupancy-normalised coupling projected on the
  leave-one-subject-out group phase, for the three primaries selected on clean data — large SOs × BEST phase,
  consensus spindles × ABDOMEN, pooled large complexes × ABDOMEN. Polarity-unambiguous subjects only; N = 89 / 72 / 70.
- **Primary test** (one per family, BH-FDR over the 3): robust regression (Huber) coupling ~ z(log1p AHI) + z(age) +
  sex + z(BMI); arousal index not adjusted (mediator). Sensitivity: OLS (HC3), event-count weighted WLS,
  + depression (antidepressant proxy), Spearman.
- **Secondary markers** computed from the data (`severity` rule): OAHI, ODI3, T90, hypoxic burden (simplified
  Azarbarzin), mean event duration, snore index, mean SaO2; arousal indices reported but not FDR-tested (mediators).
  Checks: scored-event index vs clinical AHI ρ = 1.00; mean SaO2 vs clinical ρ = 0.95; ODI3 vs AHI ρ = 0.71; detected
  vs clinical arousal index ρ = 0.65.
- **Hierarchical**: Poisson GLMM of clean-N2 phase-bin counts, (cos + sin) × (z log1p AHI + age + sex + BMI),
  random intercept and random cos/sin slopes per subject — AHI as moderator of the phase modulation.
- **Exploratory**: within-subject change (snoring / pooled events / hypopnea − clean) vs AHI.

## Results
**Primary: no association between AHI and clean coupling.**

| Family | N | β per SD log1p(AHI) (robust) | p | FDR | Spearman ρ |
|---|---|---|---|---|---|
| Large SOs | 89 | −0.003 ± 0.008 (clean mean 0.041) | 0.72 | 0.77 | −0.005 |
| Pooled complexes | 70 | +0.004 ± 0.012 (clean mean 0.051) | 0.77 | 0.77 | +0.02 |
| Consensus spindles | 72 | −0.005 ± 0.006 (clean mean 0.027) | 0.38 | 0.77 | −0.15 |

All sensitivity models agree (p ≥ 0.24).

**GLMM (AHI × phase modulation):** large SOs −2 % per SD [−218, +140] (p = 0.94; SO modulation itself weak,
p = 0.047, so the ratio is unstable); complexes −8 % [−70, +35] (p = 0.73); **spindles −12 % per SD [−51, +16]
(p = 0.41)** — i.e. a severity-related loss of more than ~half the spindle modulation per SD of log AHI is unlikely.

**Secondary markers:** nothing survives FDR (all FDR ≥ 0.62). Uncorrected: complexes vs ODI3 (+0.027/SD, p = 0.048)
and T90 (p = 0.059) — positive direction, i.e. *more* coupling with more desaturation; detected arousal index vs SO
coupling +0.020/SD (p = 0.011; clinical arousal index p = 0.65) — reported only.

**Within-subject change vs AHI (exploratory, uncorrected):** with all events, the change in large-SO coupling during
pooled respiratory events (β = −0.058/SD, p = 0.045) and hypopneas (−0.059/SD, p = 0.035) becomes more negative with
higher AHI; with snore-evoked events excluded these vanish (p ≥ 0.48). Spindles and complexes: nothing (p ≥ 0.20).

## Interpretation
Sleep-disordered-breathing severity — AHI or any oxygen/event marker — is **not related to how strongly SOs,
spindles or complexes lock to breathing during undisturbed (clean) N2 breathing**. With ~70–90 subjects and
β ≈ 0 ± 0.01 (≈ 20–30 % of the mean coupling per SD), moderate trait effects are unlikely; the GLMM bounds a
severity-related spindle loss to < ~50 % per SD. Together with Step 3, the picture is: coupling is a **state**
feature of ongoing breathing — present in clean breathing irrespective of severity, preserved through apneas and
hypopneas, and reduced (trend) only during snoring once snore-evoked events are removed.

## Caveats
- Coupling per subject is noisy (Step 3 split-half reliability), which attenuates any correlation; the GLMM is the
  better-powered estimate.
- AHI range is wide but skewed; severe OSA subjects contribute little clean N2 (and therefore few events).
- Hypoxic burden is a simplified version (fixed windows); the hypopnea scoring rule of DREAMT is unknown.
- The within-subject × AHI trend depends on including snore-evoked events and is exploratory.

## Follow-up: plain Spearman correlations (`scripts/step4_spearman.py`, `derivatives/step4/spearman.csv`, `figures/step4_spearman.png`)
Two per-subject modulation measures in clean N2: projected coupling (direction on the group phase, polarity-ok
subjects) and **direction-free modulation strength (PPC)** (all subjects; unbiased by event count).
- **Projected coupling:** no correlation with any severity marker (|ρ| ≤ 0.22; only detected arousal index vs SOs
  ρ = 0.22, p = 0.04).
- **Modulation strength (PPC):** SO modulation is *stronger* in more severe patients — ODI3 ρ = 0.30 (p = 0.003),
  snore index ρ = 0.26 (p = 0.013), detected arousal index ρ = 0.25 (p = 0.017), hypoxic burden ρ = 0.20 (p = 0.06),
  AHI ρ = 0.19 (p = 0.07); complexes similar for ODI3 (ρ = 0.24, p = 0.02) and hypoxic burden (ρ = 0.23, p = 0.03);
  spindles nothing. **Nothing survives FDR** over 72 correlations (min FDR p = 0.25).
- Checks for SO PPC vs ODI3: partial ρ controlling for event count and the share of clean events < 60 s from a
  respiratory event = 0.25 (p = 0.016); using only SOs > 60 s from any event ρ = 0.14 (p = 0.20); within subjects,
  near-event vs far SOs do not differ in PPC (p = 0.58).
- Reading: in more severe patients the SO–breathing modulation tends to be stronger but **not at the common preferred
  phase** (the projected measure shows nothing), so it is not "more of the same coupling"; it may reflect residual
  event-related responses (arousal/K-complex) in fragmented clean sleep, but this is not settled and is uncorrected.

## Detector validation
MATLAB parity of the slythm port (run by Esteban 2026-10-26, `validation/matlab_parity/`): S003 2-h segment, F4-M1 and
C4-M1 — SOs 285/285 and 255/255, spindles 902/902 and 1165/1165, 100 % matched within ±2 samples, thresholds identical.
