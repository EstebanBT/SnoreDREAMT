# Goal: complete Step 3 (coupling per breathing condition) — autonomous run until Mon 2026-10-26

Set by Esteban on 2026-10-20 ("jump to Step 3 … feel free to try multiple things until step 3 is over").
Approved plan: `~/.claude/plans/fuzzy-shimmying-breeze.md` (Step 3 section). Decisions since then are listed below.

## Decisions in force
- Phase source for condition contrasts: **ABDOMEN_hilbertgated** (peak ± 0.05 Hz Hilbert, NaN where no breath is
  detected in the raw belt signal). Sensitivity: ABDOMEN_cycle (equalised), FLOW_hilbertgated (conditions with airflow).
- All breathing conditions analysed as group tests, Ns reported, low-N conditions flagged exploratory
  (obstructive ~13 subjects, central/mixed ≤ 6).
- Arousal exclusion layered: detected arousals ±5 s + data-driven post-event window; liberal detector as sensitivity.
- Group phase tests: polarity-unambiguous subjects primary, all subjects sensitivity. Within-subject contrasts are
  polarity-invariant (projected coupling relative to the subject's own clean reference).
- Spindles: 2.0 × median RMS + spectral-peak criterion (slythm 1.5 × as sensitivity).

## Open issues to come back to (not blocking Step 3)
- **SO–spindle complexes do not show significant group coupling** in Step 2 (FLOW Hilbert ±0.05: +4°, Rayleigh
  p = 0.15, V p = 0.03; ~60 complexes/subject, ~50 subjects with ≥ 20). Esteban is troubled by this — revisit:
  complex definition (spindle minTime vs onset, same vs different channel, 1.5 s window), stricter spindle definition
  reducing complexes, N3 scarcity, SO threshold choice (2.5× median), polarity subset, per-subject counts.
- MATLAB parity of the slythm port not run yet (needs the MATLAB license server; Esteban to run
  `validation/matlab_parity/run_parity.m`).
- Arousal detector undercounts (~0.38× clinical ArI).

## Step 3 checklist
- [x] Post-event window from data (event-locked SO / spindle / arousal rates around respiratory-event ends)
- [x] Snoring detection per breath (+ validation figure)
- [x] Condition labels per sample (OA, CA, mixed, hypopnea, recovery, snoring, clean) + feasibility table
- [x] Matched clean controls (same subject & stage, ±10 min, equal duration)
- [x] Coupling per condition: projected coupling (primary), PPC, preferred phase, densities, amplitudes
- [x] Group statistics + figures; sensitivity analyses
- [x] Event-locked analysis around event onset/termination
- [x] Summary report for Monday

## Log
- 2026-10-20: Steps 0–2 committed (79953e2). Step 3 started.
- 2026-10-20: Step 3 complete — see `notes/STEP3_REPORT.md`. Result: no detectable change (underpowered: baseline
  coupling per event ≈ 0.01–0.02 vs detectable change ≈ 0.065); uncorrected trend of higher SO coupling during
  snoring; fixed whole-night circularity bias; added group-phase metric (subject phases unreliable).

# Step 3b: more statistical power (requested 2026-10-26, ideas 2 + 3, then 4)
- [x] 2a large SOs (top half amplitude within subject)
- [x] 2b multichannel consensus events (SO: F4 + C4/O2 within ±0.3 s; spindle: C4 overlapping another channel)
- [x] 2c best respiratory phase: FLOW gated where defined, else ABDOMEN gated rotated by subject FLOW–ABD offset
- [x] 3a pooled "respiratory event" condition (OA+CA+mixed+hypopnea) vs clean N2, one-sided
- [x] 3b per-breath dose–response (snore dB, airflow reduction, within-subject centred)
- [x] 4  hierarchical Poisson model on phase-bin counts with occupancy offset + ROPE bounds
- 2026-10-26: Step 3b done (see STEP3_REPORT.md, section Step 3b): complete loss excluded for pooled events
  (large SOs, consensus spindles); snoring SO increase partly snore-evoked K-complexes.
