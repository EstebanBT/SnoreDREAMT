# Step 3 report — respiration coupling during snoring, hypopneas and apneas (DREAMT)

*Prepared 2026-10-20 for review on Monday 2026-10-26. Pipeline: `Snakefile` (rules `event_locked_rates`,
`conditions`, `conditions_summary`, `coupling_step3`, `step3_group`); numbers in `derivatives/step3/*.csv`.*

## Bottom line
1. **No breathing condition changed SO, spindle or complex coupling detectably** after FDR correction
   (primary family: N2, gated ABDOMEN Hilbert ±0.05 Hz, matched clean controls). Nothing was lost or gained beyond
   chance.
2. **The design cannot exclude loss of coupling.** Within-subject coupling in this cohort is tiny (clean NREM,
   projection on the group phase: SO 0.009 ± 0.0035 and spindle 0.020 ± 0.003 per event), while the smallest
   detectable change with ~40 subjects is ≈ 0.065 per subject (clean-vs-clean null SD ≈ 0.15). Even a complete loss
   of coupling during hypopneas would fall inside the confidence intervals. The result is therefore
   **inconclusive, not a demonstration that coupling survives**.
3. **One consistent but uncorrected trend goes the opposite way to the hypothesis:** SO coupling is slightly
   *higher* during snoring than in matched clean breathing (Δ ≈ +0.04; positive in 31/32 sensitivity cells, p < .05
   in ~31 %; primary p = 0.019 / 0.070 for the two metrics, FDR 0.23 / 0.43). It is not explained by more
   breathing-locked EEG during snoring (bias-controlled EEG–ABDOMEN coherence 0.15 vs 0.13, p = 0.14); the effect size
   is similar in the half of subjects with the lowest EEG–breathing coherence (Δ +0.04 / +0.055, N = 16 / 14, n.s.)
   and does not scale with the coherence increase during snoring (ρ = −0.10, p = 0.59).
   A weaker opposite trend: **large SOs lose coupling during hypopneas** (group-reference Δ = −0.06,
   95 % CI −0.16 to 0.01, p = 0.08; negative in 31/32 cells).
4. **Obstructive, central and mixed apneas cannot be analysed as group tests** (≥ 10 events per side in only 0–5
   subjects); they stay descriptive, as anticipated in the plan.

## What was done
| Item | Result |
|---|---|
| Post-event window from data | SO rate 1.6× and arousal rate 13× clean within 0–4 s after event end, back to baseline by ~8–10 s → kept **15 s** recovery window (`figures/step3_event_locked_rates.png`). At event onsets SO/spindle rates drop to ~0.5×. |
| Snoring per breath | SNORE 20–45 Hz RMS ≥ quiet reference + 10 dB (reference = min(10th pct of event-free NREM breaths, median of wake breaths)), phasic (louder/quieter half-cycle ≥ 2), ≥ 3 consecutive breaths. In clear snorers snoring breaths carry 1.6–4× more PTAF flutter (r 0.6–0.96 where PTAF is not clipped); tonic noise (S052) is rejected by the phasic criterion. Visual check (`figures/step3_snoring_validation.png`): labelled breaths show clear inspiratory bursts; the detector is conservative — quieter burst trains below +10 dB stay unlabelled (they end up in "clean"). |
| Condition labels | priority OA > CA > mixed > hypopnea > recovery (15 s) > pre-event (5 s, excluded) > arousal (±5 s, excluded) > snoring > clean; N2 primary, N2+N3 secondary; ABDOMEN must be valid. |
| Feasibility (≥ 5 min in N2) | clean 96, snoring 52, hypopnea 58, recovery 52, obstructive 12, central 7, mixed 2 subjects (`derivatives/step3/conditions_summary.csv`). |
| Matched controls | clean samples of the same stage within ±10 min of each condition run, equal duration, no reuse; runs with < 50 % match dropped. |
| Metrics | (a) projected coupling on the subject's own clean phase μ_ref (planned primary); (b) projection on the leave-one-subject-out **group** phase from Step 2 (added: subject μ is unreliable — split-half agreement within 45° only 30 % for SOs vs 25 % chance; 47 % spindles); PPC, densities, amplitudes. |
| Negative control | `clean_split`: random clean runs vs matched clean → null behaviour after fixing a circularity bias (below). |
| Event-level model | OLS on cos(phase − group phase), condition vs control, subject-clustered SEs (MixedLM collapsed the intercept, replaced). |
| Sensitivity grid | 2 metrics × 2 arousal detectors × 2 stage sets × 2 control types (matched / whole-night) × 4 phase variants (`figures/step3_sensitivity_grid.png`). |
| Exploratory | coupling modulators in clean NREM (`figures/step3_explore_modulators.png`): carried by **larger SOs** (FLOW: Q3–Q4 ≈ 0.03–0.05 vs Q1 ≈ 0); no AHI gradient; no benefit of distance from events; N3 too scarce. Termination SOs (0–5 s after event end) are not phase-locked (Rayleigh p = 0.95). |

## Primary results (N2, ABDOMEN gated, matched controls, primary arousal detector)
Group-reference metric (subject-reference metric in `derivatives/step3/group_stats.csv`):

| Event | Condition | N | control | condition | Δ median [95 % CI] | p (FDR) | density ratio |
|---|---|---|---|---|---|---|---|
| SO | snoring | 35 | −0.001 | 0.040 | +0.047 [−0.009, 0.082] | 0.070 (0.43) | 0.97 |
| SO | hypopnea | 36 | −0.022 | 0.009 | −0.005 [−0.086, 0.049] | 0.80 (0.95) | 0.85* |
| SO | recovery | 28 | 0.000 | 0.006 | +0.011 [−0.090, 0.112] | 0.87 (0.95) | 1.60** |
| SO large | hypopnea | 24 | −0.001 | −0.104 | −0.059 [−0.156, 0.008] | 0.079 (0.43) | 0.81* |
| spindle | snoring | 37 | 0.000 | 0.007 | +0.014 [−0.052, 0.056] | 0.36 (0.90) | 1.04 |
| spindle | hypopnea | 37 | 0.030 | 0.040 | −0.002 [−0.063, 0.105] | 0.79 (0.95) | 0.97 |
| spindle | recovery | 33 | 0.023 | 0.087 | +0.040 [−0.038, 0.116] | 0.20 (0.71) | 0.85* |
| complex | snoring | 10 | 0.017 | 0.044 | +0.019 [−0.083, 0.158] | 0.49 (0.90) | 0.90 |

\* density difference p < .05; \*\* p < .01. Event-level model (all events, subject-clustered SE): snoring SO
β = +0.031 [−0.004, 0.066], p = 0.08; hypopnea large SO β = −0.042 [−0.100, 0.015], p = 0.15; no other condition
effect.

## Methodological lessons (now fixed in the pipeline)
- **Phase occupancy:** interpolated cycle phase over-represents expiration (occupancy R ≈ 0.09–0.16); the ±0.05 Hz
  Hilbert phase is uniform (R ≈ 0.002). Cycle phase is equalised; Hilbert is used raw (Esteban's recommendation).
- **Filter ringing in apneas:** the narrow-band Hilbert phase keeps "breathing" through apneas; the gated variant
  (phase only where a breath is detected in the raw belt) is the Step 3 primary.
- **Circularity in whole-night controls:** estimating μ_ref from the same clean events used as controls inflated the
  control side and made every condition (including the clean-vs-clean null) look like a loss of coupling; fixed
  by an even/odd split of clean events.
- **Polarity:** respiratory polarity differs between subjects; within-subject contrasts are polarity-free, group-phase
  analyses use only the 79/95 polarity-unambiguous subjects.

## Caveats
- Arousal detector finds ~0.38× the clinical Arousal Index; residual arousals remain in "clean" and "snoring".
- Snoring is detected without ground truth; thresholds (10 dB, phasic 2, 3 breaths) were not tuned.
- EEG–breathing coherence in clean N2 is clearly above estimation bias in many subjects (median 0.13, up to 0.6–0.8):
  some breathing-locked EEG (artifact or physiology) is present throughout and should be examined subject by subject.
- Hypopnea scoring rule (3 % vs 4 %) unknown; apnea labels were added retrospectively by the DREAMT team.

## Suggested next steps
1. **Power:** the question needs more events per condition or a stronger coupling estimate. Options: pool N2+N3 and
   all EEG channels per event (multichannel consensus), restrict to larger SOs throughout, or a hierarchical
   Bayesian model to quantify evidence *for* no change (equivalence) instead of failing to reject.
2. **Snoring SO effect:** test whether it survives restriction to subjects with low EEG–breathing coherence and with
   the FLOW phase; check SO waveform/topography during snoring.
3. **Complexes (open concern):** complex definition (spindle onset vs minTime, F4 vs C4, 1.5 s window), complexes
   per subject are few (~60), and Step 2 coupling was only weak (V-test p = 0.03).
4. Step 4 (AHI) can proceed on the Step 2 per-subject coupling (trait effect), with the caveat that subject-level
   coupling estimates are noisy.
