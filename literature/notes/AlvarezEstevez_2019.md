---
key: AlvarezEstevez_2019
title: "Large-scale validation of an automatic EEG arousal detection algorithm using different heterogeneous databases"
authors: [Diego Alvarez-Estevez, Isaac Fernández-Varela]
year: 2019
journal: Sleep Medicine 57:6-14
doi: 10.1016/j.sleep.2019.01.025
type: methods
tags: [EEG arousal, automatic scoring, polysomnography, validation, SHHS, chin EMG, open-source]
fulltext: AlvarezEstevez_2019.fulltext.md
---

# Alvarez-Estevez & Fernández-Varela (2019) — Large-scale validation of an automatic EEG arousal detector

## TL;DR
Validates an updated version of the authors' single-EEG + chin-EMG arousal detector on 2768 full-night clinical/cohort PSGs (two HMC clinical sets, SHHS-2). Epoch-based agreement with clinical scorers (κ 0.56–0.60) matches estimated human inter-scorer agreement (κ 0.54–0.59), and per-recording Arousal Index agrees with clinical ArI (r/ICC ≈ 0.65–0.79). **The algorithm itself is not specified here** — details are in ref. 16 (Alvarez-Estevez et al. 2017, *Comput Biol Med* 87:77–86) and in open-source MATLAB code (ref. 29: github.com/bigsasi/arousals-detection).

## Research question & hypotheses
Does the detector generalise to large, heterogeneous, non-curated clinical data and a different external database, at a level comparable to human inter-scorer variability? (Introduction)

## Sample & design
- HMC (Haaglanden Medisch Centrum, The Hague): 472 clinical recordings (Apr–Oct 2017), SOMNOscreen plus, EDF+, AASM scoring, TST ≥ 4 h, no exclusion for signal quality. Split into HMC-S (n = 220, scored semi-automatically with the previous algorithm version, then reviewed) and HMC-M (n = 252, purely manual). (Methods › 2.1.1) Note: APSG/HPSG counts given (176+45, 176+98, 352+143) do not sum consistently with 220/252/472.
- SHHS-2: 2296 of 2561 recordings (montage v3, arousals scored, no total technical failure, TST ≥ 4 h); R&K staging, ASDA1992 arousal rules; C3/A2 EEG + default bipolar submental EMG; HMC signals at 256 Hz; SHHS EEG/EMG sampling rate not stated in the paper (ECG 125 Hz SHHS-1, 250 Hz SHHS-2) (Methods › 2.1.2, 2.3).
- Demographics (Table 1): HMC-S age 53.0 ± 14.3, ArI 12.9 ± 7.9, AHI 13.3 ± 15.3; HMC-M 51.6 ± 16.2, ArI 12.5 ± 10.5, AHI 14.8 ± 21.0; SHHS2 67.4 ± 10.0, ArI 12.9 ± 7.0, AHI 16.3 ± 15.6.

## Methods (replication-level detail)
### Arousal definition used (Introduction; AASM)
Abrupt EEG frequency shift (alpha, theta and/or > 16 Hz, not spindles), ≥ 3 s, preceded by ≥ 10 s of stable sleep; in REM additionally a concurrent chin-EMG increase ≥ 1 s.

### Algorithm overview (Methods › 2.2) — only a summary is given
| Stage | Content |
|---|---|
| Inputs | 1 EEG derivation (HMC: C4/M1, or Cz/O2 in APSG; SHHS: C3/A2) + 1 chin EMG; optional ECG (modified lead II) for artifact removal |
| Preprocessing | Notch filter on EEG and EMG (60 Hz US / 50 Hz Europe); optional adaptive-filter ECG artifact removal (refs 16, 28); high-pass filtering of EMG |
| Candidate detection | Power content in alpha (8–12 Hz) and beta (> 16 Hz) bands of the EEG |
| Pattern recognition | (i) EEG power-based, (ii) EEG amplitude-based, (iii) EMG amplitude-based subroutines |
| Region adjustment | Candidates merged or offset extended if arousal activity follows |
| False-positive rejection | (i) duration criteria, (ii) not sleep-spindle activity, (iii) no EMG activity during REM → reject, (iv) ≥ 10 s stable sleep before onset |
| Changes vs ref. 16 | Thresholds tuned for sensitivity to alpha-prevalent arousals and better spindle discrimination; skip first false-detection check of EEG power pattern for alpha arousals; dropped the concurrent-EMG procedure in some steps |
- Parameters were set on two development sets (HMC-22, SHHS1-26) independent of test data. Analysis restricted to lights-out → lights-on (TIB). Sampling 256 Hz for HMC. (Methods › 2.3)
- Runtime ≈ 30 s per full PSG on a laptop. (Discussion)

### Validation (Methods › 2.3)
- Event-by-event on 30 s epochs (arousal assigned to the epoch of its midpoint): sensitivity, specificity, precision, F1, Cohen's κ.
- Arousal Index (ArI) per recording vs clinical ArI: Anderson-Darling/Lilliefors normality, Wilcoxon signed-rank, Pearson r, ICC (two-way absolute, single measures), GLMM-based repeatability (rptR, log link).
- Inter-scorer variability: 5 recordings per dataset at κ percentiles 12.5/37.5/50/62.5/87.5 re-scored blind by an independent expert.

## Key results
- Pooled epoch agreement (Table 2): HMC-S κ = 0.600 (sens 0.580, spec 0.972, prec 0.707, F1 0.637); HMC-M κ = 0.559 (0.563/0.953/0.641/0.600); SHHS2 κ = 0.573 (0.517/0.979/0.743/0.610).
- Per-recording medians (Table 3): κ 0.609 / 0.529 / 0.552 (HMC-S / HMC-M / SHHS2).
- ArI (Table 4): mean ref − auto differences 0.84 ± 5.41 (p = 0.023), −0.52 ± 6.68 (p = 0.224), 0.35 ± 4.89 (p < 0.001) — automatic ArI close to clinical.
- Repeatability (Table 5): r = 0.771 / 0.791 / 0.780; ICC = 0.770 / 0.790 / 0.780 (95 % CIs ≈ ±0.06).
- Human inter-scorer κ (Table 6): 0.594 / 0.561 / 0.543 → algorithm behaves like "one expert more".
- SHHS2: performance robust across EEG signal-quality grades; only the worst EMG-quality group differed (Fig. 1).
- Outperforms earlier approaches on SHHS1-26 and HMC-22 (Table 7: current κ 0.597 and 0.773).

## Figures
- Fig. 1: κ distribution by SHHS signal-quality grade (EEG and chin EMG), with recording counts per grade.

## Limitations
- Stated by authors: inter-scorer estimates from only 5 recordings per dataset; algorithm details deferred to ref. 16 and code; only one EEG channel; reference scorings partly semi-automatic (HMC-S, potential bias toward the algorithm); development sets partly overlap the original algorithm's development data.
- Not stated by authors: sensitivity is only ~0.52–0.58 (many arousals missed; specificity carries κ); validation populations have lower ArI (~13/h) than DREAMT (mean 37.7/h); validated at 256 Hz (HMC) and SHHS rates (not stated here; SHHS EEG is usually 125 Hz) — beta (> 16 Hz) and EMG analysis at DREAMT's 100 Hz are untested; 30-s epoch agreement says little about onset timing precision, which matters for phase analyses.

## Relevance to this project
- We need arousals (unannotated in DREAMT) to (a) exclude them from clean/snore/event conditions, since arousals are breathing-locked and evoke K-complex-like SOs and spindle-band bursts, and (b) test arousal–breathing coupling as a supplementary analysis.
- Inputs map directly onto DREAMT: C4-M1 (EEG), CHIN (EMG), ECG. Port should come from the public MATLAB code (github.com/bigsasi/arousals-detection) plus ref. 16; this paper alone is insufficient to implement it.
- 100 Hz adaptations: beta band must stay below Nyquist (50 Hz) and avoid a possible 40 Hz alias of 60 Hz mains; EMG high-pass + envelope limited to < 50 Hz; notch at 60 Hz is moot (alias instead).
- Validation strategy to copy: per-subject ArI vs clinical ArI (DREAMT `participant_info.csv` "Arousal Index"), Spearman/ICC; expect r ≈ 0.7 as benchmark, lower is plausible given DREAMT's high ArI.
- Onset timing matters for phase analysis; check onset precision visually since the paper only validates at 30-s epoch level.

## Follow-up references
- Ref. 16: Alvarez-Estevez D, Fernández-Varela I et al. 2017, *Comput Biol Med* 87:77–86 — the actual algorithm description; needed for the port.
- Ref. 28: Alvarez-Estevez D et al. 2017, *Biomed Signal Process Control* 35:87–99 — adaptive ECG artifact removal details.
- Ref. 29: github.com/bigsasi/arousals-detection — MATLAB source code.
- Refs. 38–43: inter-scorer reliability of arousal scoring (κ ≈ 0.47–0.57, ICC 0.54–0.98).
