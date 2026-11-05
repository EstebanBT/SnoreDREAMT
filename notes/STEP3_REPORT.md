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

---

# Step 3b — more statistical power (2026-10-26)
Rules `coupling_step3b`, `step3b_group`, `step3b_glmm` (R/lme4); scripts `snore_evoked_so.py`;
numbers in `derivatives/step3b/`; figures `step3b_power.png`, `step3b_snore_evoked_so.png`.

## What changed
- **More signal per event (idea 2)**: large SOs (≥ subject median |trough|), multichannel consensus events
  (SO: F4 + C4/O2 trough within ±0.3 s; spindle: C4 overlapping another channel), and a combined phase
  (`BEST`: FLOW gated Hilbert where a breath is detected in airflow, else ABDOMEN gated rotated by the subject's
  FLOW–ABDOMEN offset; offset consistency median R = 0.86; FLOW supplies 96 % of N2 samples).
  Primary per family chosen on **clean data only** (max group z among apnea-covering phase sources):
  **SO_large × BEST** (clean coupling 0.041 per event, z = 3.4; all SOs on ABDOMEN were 0.010, z = 1.3) and
  **spindle_cons × ABDOMEN** (0.024, z = 3.7).
- **Fewer, bigger questions (idea 3)**: pooled respiratory events (OA + CA + mixed + hypopnea) vs whole-night clean
  N2 (group-phase projection, no circularity), one-sided tests for loss; per-breath dose–response.
- **Hierarchical model (idea 4)**: Poisson GLMM on event counts per phase bin with occupancy offset,
  (1 | subject:condition) + (0 + cos + sin | subject); relative change of the phase modulation with simulation
  intervals.

## Results
**Contrasts (event-level, subject-clustered 95 % CI; condition − clean):**

| | Large SOs (clean ≈ 0.05) | Consensus spindles (clean ≈ 0.02) |
|---|---|---|
| pooled respiratory events | +0.020 [−0.033, 0.072], N = 54, p_loss = 0.77 | +0.017 [−0.010, 0.045], N = 48, p_loss = 0.89 |
| snoring | +0.006 [−0.030, 0.042], N = 45 | −0.002 [−0.031, 0.026], N = 46 |
| hypopnea | +0.023 [−0.032, 0.078], N = 51 | +0.016 [−0.014, 0.046], N = 46 |
| recovery | +0.023 [−0.025, 0.071], N = 47 | +0.037 [−0.022, 0.095], N = 43 |
| obstructive | −0.022 [−0.181, 0.137], N = 7 | +0.043 [−0.043, 0.129], N = 10 |

**GLMM (relative change of phase modulation, 95 % interval; P(loss > 50 %)):**
large SOs — events +27 % [−63 %, +182 %], P = 0.045; snoring +34 % [−54 %, +282 %], P = 0.028;
hypopnea +37 % [−63 %, +226 %], P = 0.044. Consensus spindles — events +59 % [−44 %, +216 %], P = 0.019;
snoring −27 % [−89 %, +48 %], P = 0.24; hypopnea +58 % [−58 %, +235 %], P = 0.036.
P(complete loss) ≤ 1 % everywhere. Clean modulation itself is significant in every model (p ≤ 0.014).

**Interpretation:** with the stronger event sets and phase source, **a complete loss of respiration coupling
during respiratory events can now be excluded**, and a loss of more than ~50–60 % is unlikely (P ≈ 2–5 %), for
both large SOs and consensus spindles. The data are compatible with no change. The one exception with real
uncertainty is **spindles during snoring** (point estimate −27 %, interval to −89 %).

**Dose–response (per within-subject SD):** louder snoring → *higher* large-SO coupling (β = +0.022, p = 1e-4;
consistent in every SO set and phase source); airflow and effort reduction: no linear effect (the octile plot
suggests a U shape — higher coupling at both extremes — not modelled). Spindles: no reliable loudness effect on
ABDOMEN (β = +0.003, p = 0.53).

**Confound for the snoring effect — snore-evoked SOs:** large-SO troughs are more frequent ~0.8–1.1 s after
snore-burst peaks than at the equivalent moment of quiet breaths (excess 0.3–1.0 s, p = 0.02, N = 38), and the
probability of a large SO per breath rises with snore loudness (logit +0.15 per 10 dB, p = 1e-5). Because bursts are
inspiration-locked, **the snoring-related increase in SO "coupling" is at least partly auditory/arousal-evoked
K-complexes timed by the snore**, not intrinsic respiratory pacing. This matches the concern raised in
Bullón Tarrasó et al. (sounds evoke K-complex-like SOs).

## Caveats
- The primary event sets and phase source were chosen on clean data (not on the contrasts), but after Step 3 had
  been seen; they should be treated as a second, confirmatory-style pass, not as pre-registered.
- Whole-night clean controls (instead of ±10 min matched) gain events but are exposed to time-of-night / sleep-depth
  differences; the matched analysis from Step 3 remains as a sensitivity check.
- Some GLMM fits are singular (random-slope variance ≈ 0 for spindles), i.e. effectively fixed modulation.

---

# Step 3c — snore-evoked events excluded; complexes rebuilt (2026-10-26)
Scripts/rules: `scripts/snore_evoked_events.py`, `snoring.burst_references`, `coupling_step3b` (complex sets,
evoked flags and masks), `step3b_group` (occupancy-normalised contrasts `contrasts_occupancy.csv`), `step3b_glmm`
(both exclusion levels). Figures: `step3b_snore_evoked_events.png`, `step3b_power.png`.

## Snore-evoked events
Event rates around snore-burst peaks vs the same breath position in quiet clean breaths (N = 38, 0.5 s windows):
**large SOs** excess 0.5–1.5 s after the burst (peak 0.75–1.25 s, p = 1e-4); **consensus spindles** onsets excess
0.25–1.0 s (p = 0.01–0.05; envelope peaks 0.5–1.25 s). Snore bursts therefore evoke K-complex/spindle responses.
Exclusion windows: SO and complex troughs 0.5–1.5 s, spindle onsets 0.25–1.0 s after each burst peak in snoring
breaths **and** after phase-matched pseudo-bursts in quiet clean breaths; excluded time is removed from phase
occupancy too.

**Method note:** excluding a fixed breath-phase window makes phase occupancy non-uniform, which biases any per-event
mean-cos metric (clean coupling turned negative). With exclusion, only occupancy-normalised measures are valid:
subject coupling from rate-weighted resultants (`contrasts_occupancy.csv`) and the Poisson GLMM.

## Results with snore-evoked events excluded (occupancy-normalised; N2; Δ = condition − clean)
| | all events | evoked excluded | GLMM rel. change (excluded) |
|---|---|---|---|
| Large SOs, snoring | +0.022 [−0.023, 0.070] | **−0.044 [−0.108, 0.021], p_loss = 0.050** | **−68 % [−170, +5], p_loss = 0.031**, P(loss > 50 %) = 0.70 |
| Consensus spindles, snoring | −0.006 [−0.058, 0.044] | −0.029 [−0.089, 0.027], p_loss = 0.14 | **−71 % [−149, +17], p_loss = 0.047**, P(loss > 50 %) = 0.71 |
| Large SOs, respiratory events | +0.031 | +0.029 [−0.025, 0.080] | +34 % [−53, +195], P(complete loss) = 0.003 |
| Spindles, respiratory events | +0.031 | +0.035 [−0.011, 0.080] | +80 % [−42, +317], P(complete loss) = 0.001 |
| Hypopnea / recovery | no loss | no loss | P(complete loss) ≤ 0.003 |

(clean levels ≈ 0.05 for large SOs, ≈ 0.02 for spindles; two-sided GLMM p for snoring 0.063 (SO) and 0.094 (spindle).)

**Interpretation:** the apparent *increase* of SO coupling during snoring was produced by snore-evoked
K-complexes timed by the inspiratory burst. Once those are removed, **both SO and spindle coupling to respiration
are reduced during snoring by roughly two thirds** (one-sided p ≈ 0.03–0.05; two-sided trend), whereas
**apneas/hypopneas show no loss** (complete loss excluded). This is the most direct evidence so far for the
original hypothesis, but it is trend-level and depends on an exclusion window chosen from the same data — it needs
confirmation (e.g. pre-specified windows, other datasets, sensitivity to window bounds).

## Complexes rebuilt (F4 SO trough → consensus C4 spindle envelope peak within 1.5 s)
Median 38 complexes per subject over the night (IQR 14–72); consensus/large variants fewer. Clean-N2 group phase
tests (`derivatives/step3b/clean_group_phase_tests.csv`): best variant (all SOs, BEST phase, ≥ 10 complexes/subject,
N = 66) group mean +18°, R = 0.21, **Rayleigh p = 0.056, V-test(0°) p = 0.011**; large/consensus complexes are
weaker (R ≤ 0.18, p ≥ 0.2) because they leave even fewer events (median 28–34). For comparison: large SOs
R = 0.43 (p < 1e-6), consensus spindles R = 0.34 (p = 2e-4). Condition contrasts for complexes include only 9–26
subjects and are uninformative (snoring −0.04 to −0.06, n.s.). **The complex problem is mainly event count:**
coupled SO–spindle complexes are rare in this cohort's N2-dominated sleep, so subject preferred phases are too
noisy; the group V-test is significant but the Rayleigh is not.

---

# Step 3d — complexes pooled over channels and N2+N3 (2026-10-26)
**Definition:** SO troughs from F4, C4, O2 (mastoid-referenced) merged within ±0.3 s (earliest trough, largest
amplitude), followed within (0, 1.5] s by a spindle envelope peak on any of F4, C4, O2, T3-CZ, CZ-T4
(`complex_pooled`; `complex_pooled_large` = pooled SO ≥ subject median amplitude). Group analysis now runs per stage
set (`derivatives/step3b/{N2,N2N3}/`, figures `step3b_power_{N2,N2N3}.png`).

**Counts:** ~10× more complexes (S094: 479 vs 49; median ~145 per subject in clean N2 vs ~45).

**Clean coupling (group phase tests on subject mean phases):**

| | N2 | N2+N3 |
|---|---|---|
| complex_pooled_large × ABDOMEN (primary) | N = 70, mean −16°, R = 0.28, Rayleigh p = 0.004, V(0°) p = 0.0008 | R = 0.27, p = 0.007, V p = 0.001 |
| complex_pooled × ABDOMEN | N = 74, R = 0.30, p = 0.001, V p = 0.0009 | R = 0.29, p = 0.002 |
| complex_pooled_large × BEST | N = 86, mean +47°, R = 0.27, p = 0.002 | R = 0.25, p = 0.004 |
| old single-channel complex × BEST | N = 66, R = 0.21, p = 0.056 | R = 0.20, p = 0.077 |

Complexes therefore **do** show significant respiration coupling once pooled. Adding N3 changes little (N3 is scarce).

**But most complexes are coincidental:** shifting spindle times randomly (±5–15 s) still yields 16 % of pooled SOs
followed by a spindle within 1.5 s vs 22 % observed (single-channel C4: 6 % vs 7 %), i.e. only ~25–30 % of
"complexes" exceed chance pairing. Their preferred phase equals that of large SOs (ABDOMEN −16° vs −21°), so complex
coupling to breathing largely inherits the SO coupling; genuine SO–spindle nesting is weak in this cohort.

**Conditions (pooled large complexes, ~25–32 subjects):** no change with all events. With snore-evoked events
excluded, **snoring −147 % [−539, +62] (N2, one-sided p = 0.032) / −144 % [−361, 0] (N2+N3, one-sided p = 0.023,
two-sided p = 0.047)** — same direction as SOs and spindles but very imprecise. Respiratory events, hypopnea,
recovery: no change. N2+N3 for spindles during snoring (evoked excluded): −83 % [−167, +3], one-sided p = 0.027.

**Summary across event types (snore-evoked excluded, GLMM, N2+N3):** snoring reduces respiration coupling of large
SOs (−65 %, p_loss = 0.033), consensus spindles (−83 %, p_loss = 0.027) and pooled complexes (−144 %, p_loss = 0.023);
apneas/hypopneas do not. Consistent direction, trend-to-nominal significance, wide intervals.

---

# Step 3e — arousal-exclusion sensitivity (2026-11-05)
Rules `arousals_extended`, `coupling_step3b_arv`, `step3b_group_arv` (primary event sets held fixed), `step3b_glmm_arv`,
`arousal_sensitivity`, `arousal_agreement`; numbers in `derivatives/arousal_sensitivity/` and
`derivatives/arousals/agreement.csv`.

- Agreement of the primary detector with the clinical index: median 11 vs 28/h, bias −24/h (LoA −61 to +13), ICC(A,1)
  0.32, ICC(C,1) 0.55, ρ 0.65; log–log slopes 0.59 / 1.36 bracket 1 (no proportional undercount; the earlier
  "ratio falls with clinical index" argument was regression to the mean — adversarial review).
- Extended exclusion = primary + 3735 N2/N3 alpha candidates rejected only for missing EMG + 664 candidates starting in
  sleep and ending in wake (detected median 20/h, ICC 0.47). Removes 4.7 % of clean N2 and 4.0 % of snoring N2.
- Step 3b/3c with the extended exclusion: apneas/hypopneas unchanged (complete loss still excluded, P ≤ 0.012); snoring
  with snore-evoked events excluded: SOs Δ −0.044 → −0.050 (GLMM −68 % → −60 %), spindles GLMM −68 % → −73 %.
  Conclusions do not depend on missed arousals.
- Pitfall: step3b_group.py re-selects the primary event set from the clean data; with the extended exclusion it would
  pick all spindles × BEST instead of consensus spindles × ABDOMEN. The sensitivity run holds the primary sets fixed.
- **Adopted 2026-11-05:** the extended detection replaced the original rules in all analyses (`arousals.exclusion:
  extended` in config.yaml; Step 3b primary event sets fixed in `step3b.primary_sets`). Full rerun. No conclusion
  changed; three nominal side results disappeared (Step 4 detected-arousal-index correlation, Step 5 CircStat SO
  recovery, Step 6 SO Watson–Williams across AHI groups). Snoring loss (evoked excluded, N2 GLMM): SOs −60 %
  (p_loss 0.073), spindles −73 % (0.052). The original-detection run is kept as the comparison
  (`derivatives/step3b_arousal_primary/`, rule `arousal_sensitivity`).
