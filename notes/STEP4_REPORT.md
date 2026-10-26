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
