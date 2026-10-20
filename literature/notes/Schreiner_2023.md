---
key: Schreiner_2023
title: "Respiration modulates sleep oscillations and memory reactivation in humans"
authors: [Thomas Schreiner, Marit Petzka, Tobias Staudigl, Bernhard P. Staresina]
year: 2023
journal: Nature Communications
doi: 10.1038/s41467-023-43450-5
type: empirical
tags: [respiration, slow oscillations, spindles, SO-spindle coupling, memory reactivation, NREM sleep, EEG, nap, phase-amplitude coupling, BreathMetrics]
fulltext: Schreiner_2023.fulltext.md
---

# Schreiner et al. (2023) — Respiration modulates sleep oscillations and memory reactivation

## TL;DR
Re-analysis of a nap EEG dataset (N = 20 healthy young adults, thermistor airflow) asking whether breathing paces NREM sleep rhythms. SOs, fast spindles (12–18 Hz) and coupled SO_spindle complexes are non-uniformly distributed over the respiratory cycle, clustering around the inhalation peak (SOs/SO_spindles just before, spindles just after), and the strength of respiration–SO_spindle coupling predicts memory-reactivation decoding accuracy across participants. Correlational only; healthy sleepers not screened for sleep-disordered breathing.

## Research question & hypotheses
- Is respiration an underlying pacemaker for NREM sleep oscillations (SOs, spindles, SO-spindle coupling)? (Intro)
- Hypothesis: respiration–sleep-oscillation association should appear as respiration-locked EEG power changes and as non-uniform distribution of SO/spindle events across the respiratory phase (Results › Respiration modulates NREM sleep rhythms).
- Further question: is respiration–SO_spindle coupling linked to memory reactivation during SO_spindles (decoded in ref. 13, same dataset)? (Results › Respiration modulates SO_spindle complexes)
- Motivation includes parallels between age-related breathing changes / sleep-disordered breathing and age-related changes in SO-spindle coupling precision (Intro).

## Sample & design
- Re-analysis of the dataset from Schreiner et al. 2021 (ref. 13, Nat Commun 12:3112). (Methods)
- N = 20 healthy participants; mean age 20.75 ± 0.35 years; 17 female. 5 additional participants excluded for insufficient sleep (< 30 min sleep in one session). No a-priori sample-size calculation. (Methods)
- Screening: PSQI, morningness-eveningness questionnaire, self-developed health/stimulant questionnaire; no medication, no neurological/psychiatric disorder; good self-reported sleep quality. **Not screened for sleep-disordered breathing**; no pulse oximetry or plethysmography; no clinician review. (Methods; Discussion)
- Ethics: University of Birmingham Research Ethics Committee. (Methods)
- Design: two sessions ≥ 1 week apart (mean 8.5 ± 0.85 days), order counterbalanced, investigators not blinded. Adaptation nap beforehand. Arrival 11 a.m.; ~12 p.m. PVT + encoding of 120 verb–object or verb–scene pairs (category per session counterbalanced); memory test before and after a 120 min nap opportunity starting ~1 p.m.; post-nap PVT and memory test; localizer task (new object and scene images) at the end of each session. (Methods › Experimental overview; Fig. 1)
- Total sleep time 101.63 ± 2.23 min (Methods) / 101.63 ± 2.3 min (Results) — minor inconsistency in the text.
- Recording: Brain Products 64-channel EEG, impedances < 10 kΩ, online reference FCz, 1000 Hz sampling; EMG and ECG for PSG; respiration via **Embla thermistor airflow sensor**. Sleep scored offline by two independent raters (AASM 2007 manual, ref. 91). (Methods › EEG and respiration)
- Software: MATLAB 2018b, FieldTrip (v.09/01/2020), Psychophysics Toolbox 3, BreathMetrics, CircStat v1, MVPA-Light. (Methods)
- Data/code: processed data and code on OSF, https://doi.org/10.17605/OSF.IO/D6YHB; raw data on request only. (Data/Code availability)

## Methods (replication-level detail)

### Preprocessing (Methods › Data analysis)
- Downsampled to 200 Hz. Noisy channels visually identified, removed, interpolated (weighted neighbour average).
- Sleep data re-referenced to **linked mastoids**.
- Segmented into 4 s epochs time-locked to inhalation peaks [−2, +2 s].
- Representative electrodes: **F2** for SO and SO_spindle analyses, **P1** for spindles (Supplementary Figs 10–11 show averages of significant electrodes / strongest electrodes; not in the full text).
- Only movement-free / artifact-free data from **NREM stages 2 and 3** used for event detection (N2 + SWS pooled).
- Data from both sessions collapsed per participant → 20 datasets for statistics (Methods › Statistics).

### Event detection (Methods › Event detection)

| Parameter | Slow oscillations | Fast spindles |
|---|---|---|
| Band-pass | 0.3–1.25 Hz | 12–18 Hz |
| Filter | two-pass FIR, order = 3 cycles of low cutoff | two-pass FIR, order = 3 cycles of low cutoff |
| Stages | N2 + N3, movement-free | N2 + N3, artifact-free |
| Criterion | consecutive positive-to-negative zero-crossings; SO candidate = downstate followed by upstate | RMS with 200 ms moving average; threshold = 75th percentile of RMS values |
| Duration | 0.8–2 s | > 0.5 s and < 3 s above threshold |
| Amplitude criterion | **none reported** | 75th percentile RMS (reference distribution not specified) |
| Phase reference time | SO downstate (negative peak) | spindle onset |

- 12–18 Hz (rather than 12–16 Hz) chosen for consistency with ref. 13 and because fast spindles can exceed 16 Hz.
- Slow spindles (control analysis): 9–12 Hz (ref. 34); detection presumably same RMS algorithm (not explicitly stated). (Results; Supplementary Fig. 3)
- **SO_spindle complex**: SO with a spindle detected within SO downstate + 1.5 s. (Methods; Results)
- Event-free control segments: 4 s NREM intervals with no SO, spindle or SO-spindle event, drawn from ±5 min around the corresponding event.
- **Inhalation peaks: BreathMetrics** (Noto et al. 2018, ref. 26) applied to the thermistor signal.

### Respiration-locked time–frequency analysis (Methods › Time-frequency analysis)
- Inhalation-peak-centred segments [−4, +4 s] (padding), window of interest [−2, +2 s].
- Sliding-window Fourier, 50 ms steps, window = 5 cycles per frequency, 1–25 Hz in 1 Hz steps, Hanning taper; power z-scored across time [−2, 2 s].
- Contrast: inhalation-locked vs randomly selected segments ("randomly centred data segments"); cluster-based permutation across time × frequency × electrode.
- Source analysis: LCMV beamformer (FieldTrip), 10 mm grid, standard BEM, AAL-derived gray-matter grid → 1457 virtual sensors; TFR 1–25 Hz in 2 Hz steps, 100 ms steps, 5 cycles, Hanning; z-scored [−2, 2 s]. (Methods › Source analysis)

### Respiratory phase extraction (Methods › Respiration-sleep oscillation coupling)
- Per participant, respiratory peak frequency determined (mean **0.25 ± 0.005 Hz**).
- Respiration band-pass filtered at **peak frequency ± 0.05 Hz**, two-pass Butterworth, "order = three cycles of the low-frequency cutoff" (unusual specification for an IIR Butterworth — see ambiguities).
- Hilbert transform → instantaneous phase. Phase convention: **0° = inhalation peak**.
- Filtering described as applied to respiratory data "locked to inhalation peaks" — unclear whether filtering was done on continuous data or on epochs (a 4 s epoch is ~1 breath at 0.25 Hz).
- Alternative phase estimate via double interpolation (Kluger et al. 2021, ref. 21) in Supplementary Fig. 12 (not in full text).
- Phase taken at SO downstates (SO and SO_spindle coupling) or spindle onsets (spindle coupling). Per-participant preferred phase = circular mean over events; group test across participants' preferred phases.

### Modulation index (Methods › Modulation Index)
- Tort et al. 2010 MI (KL distance from uniform), 20 phase bins.
- Phase: continuous respiration filtered 0.25 ± 0.05 Hz (two-pass Butterworth) — here a fixed band, not individual peak.
- Amplitude: electrode F2, two-pass FIR (order = 3 × lower frequency bound), 1 Hz bandwidth bins, centre frequencies 0.5–24.5 Hz in 1 Hz steps, Hilbert envelope. (Text says "20 equally spaced frequency bins", which is inconsistent with 0.5–24.5 Hz in 1 Hz steps = 25 bins.)
- Surrogates: random shifts of phase time series, 200 repetitions; MI z-scored against surrogate mean/SD; significance |z| > 1.96.

### Phase-slope index (Methods › Phase-slope index)
- Cross-frequency PSI (Jiang et al. 2015) between respiration (peak ± 0.05 Hz) and F2 (0.3–1.25 Hz) for SOs, P1 (12–18 Hz) for spindles. Positive = respiration leads. Tested against zero (t test).

### Peri-event time histograms (Methods)
- Seed = inhalation peak (t = 0), targets = SO downstates / spindle onsets, 50 ms bins, window ±1.5 s; normalised to % of events in the window, per participant and electrode; compared with event-free segments (cluster permutation).

### SO_spindle rate across respiration cycle (Methods)
- At F2: count SO_spindles and all SOs in 20 phase bins (−π to π, 0.31 rad steps); rate = SO_spindles/SOs × 100 per bin; per-participant Kolmogorov–Smirnov test of uniformity, FDR across participants.

### Memory reactivation decoding (Methods › Multivariate analysis; details in ref. 13)
- MVPA-Light, LDA; CAR re-reference; temporal generalisation trained on localizer (objects vs scenes), tested on SO_spindle segments [−1.5, 1.5 s]; baseline over whole trial; z-scored across trials; sessions pooled; PCA, 30 components; 150 ms running-average smoothing; AUC; surrogate = 250 label shuffles.

### Statistics (Methods › Statistics)
- Cluster-based permutation (FieldTrip): dependent-samples t, cluster-forming P = 0.05, 1000 randomisations, sum-of-t cluster statistic, Monte Carlo, two-tailed α = 0.05.
- Circular: V-test against 0° (inhalation peak), CircStat; unimodality checked with Watson's test vs von Mises (all P > 0.1).
- Robust regression (per electrode) of respiration–SO_spindle vector length vs memory reactivation (decoding accuracy averaged over the significant cluster from ref. 13).
- FDR (Benjamini–Hochberg) across electrodes (circular stats, regression) or participants (KS tests).

## Key results
- **Respiration-locked TFR**: increased power around inhalation peak in SO/theta range followed by fast-spindle burst (12–18 Hz) vs random segments; P < 0.001, cluster-corrected over time/frequency/electrodes. Source: frontoparietal + right MTL, P = 0.037. (Results; Fig. 2a–b)
- **MI** at F2: respiration phase modulates EEG amplitude with local peaks at 0.5 Hz (SO) and 15.5 Hz (spindle); z > 1.96, corrected across frequencies. (Fig. 2c)
- **SOs**: non-uniform respiratory phase at SO downstates over frontal, central, parietal electrodes; V-test vs 0°, Vmean = 7.42 ± 0.28, all P < 0.05 FDR. Mean angle of significant electrodes −9.4° ± 0.22, mean vector length (MVL) 0.51; F2: −21.36° ± 0.20, MVL 0.58. Significant individually in 15/20 participants. PSI (F2) = 0.0003 ± 0.0001, t(19) = 2.17, P = 0.042. PETH vs event-free: P = 0.004. (Fig. 2d–e)
- **Fast spindles**: non-uniform phase at spindle onsets over parietal electrodes; Vmean = 6.84 ± 0.24, all P < 0.05 FDR. Mean angle +16.95° ± 0.23, MVL 0.52; P1: +14.06° ± 0.23, MVL 0.44. Individually significant in 18/20. PSI (P1) = 0.0024 ± 0.0007, t(19) = 3.1, P = 0.005. PETH vs event-free: P = 0.002. (Fig. 2d–e)
- SOs group at earlier respiratory phases than spindles (phase shift; Supplementary Fig. 2).
- **Slow spindles (9–12 Hz)**: no significant non-uniformity (all P > 0.05) → effect specific to SOs and fast spindles. (Supplementary Fig. 3)
- **SO_spindle complexes**: non-uniform phase (locked to SO downstate) over frontal, central, parietal sites; Vmean = 9.46 ± 0.23, P < 0.05 corrected; mean angle −4.5° ± 0.16, MVL 0.73; F2: −24.62° ± 0.20, MVL 0.57. Individually significant in 18/20. (Fig. 3a)
- **SO_spindle rate (SO_spindles / SOs) across 20 respiratory phase bins at F2**: non-uniform in all 20 participants (KS, all P < 0.001, FDR). (Fig. 3a right; Supplementary Figs 5–6)
- **Memory link**: robust regression, respiration–SO_spindle vector length vs memory reactivation, significant positive relationship in a frontal cluster (P < 0.05 FDR); F2: R² = 0.27, F(1,19) = 6.69, P = 0.018. Respiration–SO or respiration–spindle coupling alone: no significant relation (Supplementary Fig. 8). Behavioural consolidation relation in Supplementary Fig. 7 (not described in main text). (Fig. 3b)
- Results unchanged with female-only subsample (Supplementary Fig. 9).
- My conversion (not in paper): at 0.25 Hz (4 s cycle), 1° ≈ 11 ms, so −9.4° ≈ −100 ms (SOs), +16.95° ≈ +190 ms (spindles), −4.5° ≈ −50 ms (SO_spindles) relative to inhalation peak.

## Figures
- Fig. 1: Paradigm (encoding 120 verb–object/scene pairs, pre/post-nap tests, 120 min nap, localizer); example hypnogram; 30 s NREM segment at Cz with respiration; IRASA 1/f-corrected spectra of EEG and respiration.
- Fig. 2: (a) inhalation-locked TFR vs random segments with topography; (b) source localisation; (c) MI spectrum at F2; (d) V-test topographies, circular plots (SO at F2, spindle at P1) and PETHs ±1.5 s; (e) PSI for SOs and spindles.
- Fig. 3: (a) respiration–SO_spindle phase topography, circular plot at F2, SO_spindle rate across respiratory phase; (b) robust regression of coupling strength vs memory reactivation, topography + F2 scatter.

## Limitations
- Stated by authors:
  - Correlational; causal manipulation of breathing needed (Discussion).
  - No screening for sleep-disordered breathing; no pulse oximetry/plethysmography; no clinician evaluation. Argued against major SDB influence: no fragmentation, no observed breathing cessation, spindle frequency/topography in normal range (Discussion).
  - 17/20 female; sex effects on sleep/oscillations possible (female-only reanalysis unchanged) (Discussion).
  - Scalp EEG cannot address hippocampal ripples (Discussion).
  - Sensory (olfactory/nasal) vs brainstem routes cannot be dissociated; breathing route (nose vs mouth) not controlled (Discussion).
- Not stated by authors:
  - Daytime nap (~100 min) in young healthy adults; generalisation to full nights, older or clinical populations untested.
  - N2 and N3 pooled; no stage-specific breakdown, so it is unknown whether the effect depends on N3/SWS density.
  - SO detection uses only zero-crossing duration (0.8–2 s) with no amplitude threshold reported, so "SOs" include many low-amplitude slow waves; results may differ with amplitude-thresholded detectors.
  - Narrow-band (peak ± 0.05 Hz) Hilbert phase assumes regular, quasi-sinusoidal breathing; ill-suited to irregular breathing, and filter-order specification for the Butterworth is unclear.
  - Group V-test is computed on 20 per-participant preferred phases; individual-level non-uniformity only in 15/20 (SOs). Vector lengths reported are across-participant phase consistency, not within-participant event-level locking strength.
  - Thermistor signal (temperature-based, slower and phase-lagged relative to true airflow) used for inhalation peaks; the exact phase angles therefore depend on the sensor.
  - Memory regression: N = 20, one predictor, R² = 0.27 at a single representative electrode; susceptible to small-sample instability.
  - No control for respiratory-related EEG artifacts (movement/pulse) beyond the slow-spindle negative control.

## Relevance to this project
Our question: is SO/spindle–respiration coupling lost when people snore (DREAMT, 100 sleep-clinic patients)?

- **Comparison point / hypothesis anchor**: healthy-sleeper reference values — SOs at about −9° (just before inhalation peak), fast spindles at about +17°, SO_spindles at about −4.5° relative to inhalation peak, across-participant MVL ~0.5–0.7, individual-level significance in 15–18/20. The authors themselves propose OSA/arrhythmic breathing as a quasi-causal test (Discussion), which directly motivates the snoring contrast.
- **Methods to borrow**: SO detector (0.3–1.25 Hz, two-pass FIR, 0.8–2 s zero-crossing duration), spindle detector (12–18 Hz, 200 ms RMS, 75th percentile, 0.5–3 s), SO_spindle = spindle within SO downstate + 1.5 s, BreathMetrics inhalation peaks, Hilbert phase with 0° = inhalation peak, V-test against 0°, per-participant circular mean, SO_spindle/SO rate across 20 phase bins with KS test, event-free/random-segment controls.
- **Channel mapping**: Schreiner used F2 (SO, SO_spindle) and P1 (spindles), linked mastoids. DREAMT offers F4-M1 and C4-M1 (single contralateral mastoid). F4-M1 is the closest analogue for SO/SO_spindle; C4-M1 is the best available proxy for fast spindles (no parietal channel equivalent to P1). Expect somewhat different amplitudes/phases due to reference and location.
- **Respiration sensor**: Schreiner used a thermistor → DREAMT's thermistor FLOW channel is the like-for-like replication; PTAF (nasal pressure) can serve as a sensitivity analysis (more direct airflow, less phase lag). Thermistor vs pressure phase lag should be characterised before comparing angles with this paper.
- **Sampling**: they analysed at 200 Hz; DREAMT 100 Hz (Nyquist 50 Hz) is sufficient for 0.3–18 Hz bands.
- **Little N3 in DREAMT**: Schreiner pooled N2 + N3 and still found effects in a ~100 min nap, so N2-dominant data are in principle adequate; but SO counts will be lower and per-participant phase estimates noisier — consider minimum-event thresholds and reporting N2-only results.
- **Phase-estimation caveat for snoring/apnea**: the ±0.05 Hz narrow-band Hilbert approach will produce meaningless phase during apneas/hypopneas, flow-limited or irregular breathing. For the patient data prefer peak-based/interpolated phase (Kluger 2021 double interpolation, their Supplementary Fig. 12) or BreathMetrics-derived cycle phase, and exclude scored apnea/hypopnea segments (or analyse them separately) so that "loss of coupling" is not an artifact of phase estimation.
- **Confounds this paper raises for us**:
  - Population: their sample was young (≈21 y), mostly female, healthy; DREAMT patients differ in age, sex mix, BMI and SDB severity — all plausibly affect SO-spindle coupling (Intro cites age effects on coupling precision). Match or covary.
  - Snoring-related/respiratory arousals and K-complexes at event termination are respiration-locked and SO-like; they could inflate apparent respiration–SO coupling in snorers rather than abolish it. Use arousal annotations and consider amplitude/shape criteria.
  - Snoring vibration/movement artifacts in EEG near the respiratory frequency; check SNORE-locked EEG averages and use event-free/random-segment controls as Schreiner did.
  - Breathing route (mouth vs nose) — snoring often co-occurs with mouth breathing; the authors note nasal breathing may be key (Discussion).
- Statistical design: Schreiner's tests are within-group (vs uniform/vs random); our snoring vs non-snoring comparison needs a between-condition test (e.g., within-subject snore vs non-snore epochs, comparing MVL or SO_spindle rate modulation), with matched event counts since MVL is biased by N.

## Follow-up references
- Schreiner et al. 2021, Nat Commun 12:3112 (ref. 13) — original dataset; full details of detection and decoding.
- Bullón Tarrasó, Staudigl & Schreiner 2026 (manuscript, in this library) — follow-up on breathing-cycle gating of memory reactivation.
- Noto et al. 2018, Chem Senses (ref. 26) — BreathMetrics, inhalation-peak detection.
- Kluger et al. 2021, eLife (ref. 21) — double-interpolation respiratory phase, more robust to irregular breathing.
- Kluger & Gross 2021, PLOS Biol (ref. 22); Zelano et al. 2016, J Neurosci (ref. 17) — respiration modulation of brain oscillations in wake.
- Karalis & Sirota 2022, Nat Commun (ref. 46) — breathing coordinates cortico-hippocampal dynamics in sleeping mice (mechanism).
- Ngo et al. 2013, Neuron (ref. 31); Staresina et al. 2015 (ref. 93) — SO/spindle detection algorithms used.
- Mölle et al. 2011, Sleep (ref. 34) — slow vs fast spindle definitions.
- Tort et al. 2010 (ref. 30); Jiang et al. 2015 (ref. 33); Hülsemann et al. 2019 (ref. 97) — MI, PSI, PAC method comparison.
- Weiner & Dang-Vu 2016, Neural Plast (ref. 70) — spindle abnormalities in sleep disorders incl. OSA.
- Cellini 2017 (ref. 69); Ahuja et al. 2018 (ref. 68); Kimoff 1996 (ref. 71) — OSA, memory and sleep fragmentation.
- Arzi et al. 2010, Chem Senses (ref. 67) — odors alter sleep respiration (manipulation approach).
- Cox et al. 2018, Sleep (ref. 14); Helfrich et al. 2018, Neuron (ref. 10); Muehlroth et al. 2019 (ref. 16) — inter-individual and age variability in SO-spindle coupling (covariates for a clinical sample).
