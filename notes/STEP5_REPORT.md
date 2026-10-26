# Step 5 report — does the preferred respiratory phase change during breathing conditions? (2026-10-26)
Rule `step5_phase` (+ updated `step3b_glmm.R`, `coupling_step3b` control condition); numbers in
`derivatives/step5/{N2,N2N3}/` and `derivatives/step3b/{N2,N2N3}/glmm.csv`; figures `step5_phase_{N2,N2N3}.png`.

## Design
Event families / phase sources as Step 3b–3d (large SOs × BEST, consensus spindles × ABDOMEN gated, pooled large
complexes × ABDOMEN gated); conditions vs whole-night clean N2; snoring with and without phase-matched exclusion of
snore-evoked events. Control: clean breathing split into alternating 10-min blocks (clean_odd vs clean_even).
Within-subject phase differences are polarity-invariant → all subjects used.
- **5a (any change, any direction):** per subject |Δφ| between occupancy-normalised mean directions, against a
  parametric null (shared pooled rate profile, each condition's own occupancy and event count).
- **5b (consistent rotation):** per-subject signed Δφ; test of mean rotation ≠ 0 (Wilcoxon on sin Δφ), bootstrap CI;
  population estimate from the Poisson GLMM: orthogonal (cond × sin) component — a linear, unbiased test of rotation —
  plus rotation angle and modulation-length ratio with simulation CIs.

## Calibration lessons (fixed)
1. **Rayleigh on Δφ is the wrong 5b test:** with no change, Δφ clusters at 0 when both conditions are coupled, so the
   clean-vs-clean control was "significant" (SO R ≈ 0.4, p ≈ 0). Replaced by a test of mean rotation ≠ 0.
2. **The parametric 5a null is anti-conservative for spindles/complexes** (events cluster in time): the control gave
   p = 0.02–0.05 (spindles) and p < 0.001 (complexes); SOs were calibrated (p = 0.3–0.5). 5a is therefore tested
   against each subject's own clean-vs-clean z (paired, "calibrated").
3. **GLMM random slopes:** the condition × modulation terms now have random subject slopes; fits are singular
   (no detectable heterogeneity) and results are unchanged versus the previous model — Step 3c conclusions stand.
4. **Vector-length ratios are biased upward** in sparse conditions (noise adds length) — interpret
   `glmm_length_ratio` > 1 cautiously; rotation is judged from the linear orthogonal component.

## Results (N2; N2+N3 equivalent)
**5a — any phase change per subject (calibrated):** essentially none. SOs during snoring (snore-evoked excluded)
p = 0.038 and recovery p = 0.087 (N2; N2+N3 recovery p = 0.040); all spindle/complex conditions p ≥ 0.4.

**5b — subject-level consistent rotation:** none (all p ≥ 0.17). Individual subjects have only ~20–50 condition
events, so single-subject preferred phases are too noisy to show a shared rotation.

**5b — population (GLMM) rotation:** controls flat (clean_odd: rotation −5° to +3°, orthogonal p ≥ 0.76).

| Family | Condition | Orthogonal shift p | Rotation [95 % CI] |
|---|---|---|---|
| **Spindles** | hypopnea | **0.002** (evoked-excl. 0.003) | **+46° [+4, +93]** |
| **Spindles** | pooled respiratory events | **0.001** (0.001) | **+43° [+6, +83]** |
| **Spindles** | snoring | **0.014** (evoked-excl. 0.032) | **−42° [−85, −6]** (excl. −72° [−134, −15]) |
| Spindles | recovery | 0.92 | −5° |
| Complexes | hypopnea / events | 0.038 / 0.032 (excl. 0.004 / 0.003) | +54° / +58° (excl. +110° / +114°) |
| Large SOs | events / hypopnea / recovery | 0.19–0.31 | +20° to +29° (n.s.) |
| Large SOs | snoring | 0.27–0.71 | n.s. |

Robustness of the spindle shifts: dropping the 5 or 10 subjects with most condition events keeps the orthogonal
shift significant (hypopnea p = 0.02, events p = 0.01/0.002, snoring p = 0.04/0.016) with similar rotation
(CIs widen). Phase occupancy is uniform in every condition (R ≈ 0.01), so the shift is not an occupancy artefact.
The FLOW–ABDOMEN phase offset changes by only −9° (hypopnea) / −6° (snoring), too little to explain ±40–70°.
Event-pooled rate profiles reproduce the GLMM (hypopnea +25° vs clean +7°; snoring −51°).

## Interpretation
- **No evidence that the preferred phase changes consistently at the level of individual participants** — but the
  per-subject tests are underpowered (few events per condition).
- **At the population level, spindles shift their timing within the breathing cycle without losing modulation:**
  ~+45° (≈ 0.5 s later at 0.25 Hz breathing) during hypopneas/respiratory events and ~−40° to −70° (earlier) during
  snoring, relative to clean breathing. Complexes (spindle-defined) shift the same way during events. SOs show no
  significant rotation (same positive direction, n.s.).
- **This reframes the Step 3c spindle result:** the apparent ~70 % "loss" of spindle coupling during snoring in the
  projected measure is largely a **phase shift** (length ratio ≈ 1.1–1.4) rather than a weakening. The large-SO loss
  during snoring (evoked excluded) is not a rotation (orthogonal p = 0.27) — it remains a candidate loss.
- Caveat: these are population (event-weighted) effects across many models (3 families × 5 conditions × 2 exclusion
  levels × 2 stage sets); the spindle hypopnea/event shifts (p ≈ 0.001–0.003) are the most robust, snoring
  (p ≈ 0.01–0.04) less so. Phase reference is the abdominal belt.

---

## Step 5 with CircStat tests (Esteban's suggestion; rule `step5_circstat`, `snoredreamt/circstat.py`)
Faithful Python ports of CircStat 2012a `circ_wwtest`, `circ_cmtest`, `circ_median`, `circ_mtest`, `circ_confmean`,
`circ_medtest`, `circ_kappa` (tests in `tests/test_circstat.py`; MATLAB parity script in
`validation/circstat_parity/` — run `run_parity.m`, then `compare.py`). One deviation: `circ_medtest` can return
p > 1 in MATLAB (overlapping tails when n1 = n2); capped at 1. Event phases used directly (occupancy uniform).
Per subject: condition events vs clean N2 events, ≥ 10 per group (CircStat's requirement for the median test).
Snore-evoked exclusion only for snoring and the control (it is phase-matched only there).

**5a — how many subjects have different circular means / medians?**
- **Watson–Williams is not usable here:** CircStat's applicability rule (average resultant ≥ 0.45 for n ≥ 11) fails
  in 92–100 % of comparisons (within-subject coupling R ≈ 0.05–0.15), and the test is badly anti-conservative — the
  clean-vs-clean control flags 28–67 % of subjects (the small-κ correction β = 1 + 3/(8κ) explodes as κ → 0).
- **Common median test is calibrated** (control: 3–11 % of subjects at p < .05). Fraction of subjects with
  different medians (N2, all events; control in brackets): SOs — events 6 %, snoring 12.5 % (p vs control 0.058),
  hypopnea 10 %, recovery 14 % (p = 0.03) [3 %]; spindles — events 4 %, snoring 10 % (p = 0.09), hypopnea 0 %,
  recovery 8 % [3 %]; complexes — events 13 %, snoring 0 %, hypopnea 14 %, recovery 12 % [4 %]. Snoring with
  snore-evoked events excluded: SOs 11 % [11 %], spindles 15 % [7 %] (p = 0.11), complexes 3 % [6 %].
  N2+N3 is equivalent. → no reliable excess of subjects with a changed median over the clean-vs-clean control.

**5b — does the population distribution of per-subject mean phase differences differ from 0?**
`circ_mtest(Δμ, 0)` never rejects (0 lies inside the 95 % CI of the population mean in every case; the CI is undefined
for some spindle/complex conditions because the Δμ distribution is too dispersed); `circ_medtest(Δμ, 0)` p ≥ 0.17
(median differences: p ≥ 0.07). Population mean differences: SOs +15° to +40° (CIs ≈ ±30–60°), spindles +0° to +20°,
complexes inconsistent.

**Conclusion of the CircStat analysis:** at the level of individual participants (each subject weighted equally) there
is no evidence for a change of preferred phase in any breathing condition — consistent with the per-subject analyses
above. The spindle rotations reported earlier come only from the event-weighted population models (GLMM / pooled
rates), which treat each event, not each subject, as the unit of information; the two views answer different
questions and should be reported as such.
