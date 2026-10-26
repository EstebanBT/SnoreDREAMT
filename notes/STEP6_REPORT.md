# Step 6 report — do respiratory-phase effects relate to sleep-disordered-breathing severity? (2026-10-26)
Rules `step6_phase_severity`, `step6_glmm` (R/lme4); numbers in `derivatives/step6/{N2,N2N3}/`; figure
`step6_{N2,N2N3}.png`. CircStat ports extended with `circ_corrcl`, k-group `circ_wwtest` / `circ_cmtest` (tested).

## Design
Event families as Step 3b–5 (large SOs × BEST, consensus spindles × ABDOMEN gated, pooled large complexes ×
ABDOMEN gated). Severity markers: log1p AHI, OAHI, ODI3, T90, hypoxic burden; mean event duration; snore index.
- **6a — within-subject phase change vs severity** (pooled events, snoring, hypopnea, recovery; snoring also with
  snore-evoked events excluded; all subjects, polarity-invariant): signed per-subject mean phase difference
  (Step 5 CircStat) vs markers with `circ_corrcl`; size of change (|Δμ| and −log10 common-median p) vs markers with
  Spearman and rank regression adjusted for age, sex, BMI. BH-FDR over all 6a tests.
- **6b — clean preferred phase vs severity** (polarity-unambiguous subjects): `circ_corrcl` with markers; AHI groups
  (< 5, 5–15, 15–30, ≥ 30) with k-group common median and Watson–Williams (+ CircStat applicability). BH-FDR over 6b.
- **Population level:** Poisson GLMM, (clean modulation + condition change) × z(log1p AHI), random subject slopes:
  AHI × clean cos/sin (6b) and AHI × condition cos/sin (6a).

Method note: Huber robust regression was replaced by rank regression after one subject with extreme mean event
duration (z = 4.3) produced a spurious p ≈ 0 (Huber protects against outcome outliers, not leverage; OLS-HC3 p = 0.17).

## Results
**6a — no association between severity and within-subject phase changes.** Smallest FDR-corrected p = 0.79 (N2) /
0.99 (N2+N3) across 315 tests. Uncorrected hits (p ≈ 0.007–0.04: SOs during recovery |Δμ| vs ODI3/T90; complexes
during snoring vs ODI3/AHI) are scattered and not consistent between N2 and N2+N3.

**6b — no association between severity and the clean preferred phase.**
| | log AHI circ-lin ρ (p) | other markers | AHI groups, common median |
|---|---|---|---|
| Large SOs (N = 89) | 0.18 (0.23) | ρ ≤ 0.23, p ≥ 0.10 | p = 0.13 (WW p = 0.0045, applicable; FDR 0.12) |
| Consensus spindles (N = 72) | 0.23 (0.15) | p ≥ 0.12 | p = 0.82–0.95 |
| Pooled complexes (N = 70) | 0.11 (0.66) | p ≥ 0.31 | p = 0.55–0.60 |
The only nominal result (SO Watson–Williams across AHI groups) is non-monotonic — group means < 5: +33°, 5–15: +95°,
15–30: −13°, ≥ 30: +32° — and not supported by the common median test.

**GLMM (population level):** AHI does not moderate the clean modulation (strength p ≥ 0.18; direction p ≥ 0.07 —
spindle clean phase −0.04 per SD, p = 0.07–0.09) nor the condition-related rotation (p ≥ 0.12) or strength change
(p ≥ 0.10) in any family/condition. (The spindle snoring rotation itself remains, p = 0.011–0.02.)

## Interpretation
Neither the preferred respiratory phase of SOs, spindles or complexes in clean breathing nor its change during
apneas/hypopneas/snoring depends on sleep-disordered-breathing severity (AHI, oxygen markers, event duration, snore
index), at the individual or population level. This parallels Step 4 (no severity effect on coupling strength):
respiration–oscillation coupling behaves as a state property of ongoing breathing, not a trait that scales with
disease severity. The Step 4 hint (stronger but not phase-aligned SO modulation with ODI3) is not explained by a
severity-dependent shift of the clean preferred phase.

## Caveats
Per-subject phase estimates are noisy (Step 5), which limits 6a/6b power; the GLMM is the better-powered test and is
also null. Severe-OSA subjects contribute little clean N2.
