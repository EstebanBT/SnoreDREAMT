---
key: FernandezVarela_2017b
title: "A simple and robust method for the automatic scoring of EEG arousals in polysomnographic recordings"
authors: [Isaac Fernandez-Varela, Diego Alvarez-Estevez, Elena Hernandez-Pereira, Vicente Moret-Bonillo]
year: 2017
journal: Computers in Biology and Medicine 87:77-86
doi: not printed in the text (Elsevier PII S0010-4825(17)30133-6)
type: methods
tags: [EEG arousal, automatic scoring, rule-based, adaptive ECG artifact removal, chin EMG, sleep spindles]
fulltext: FernandezVarela_2017b.fulltext.md
---

# Fernández-Varela et al. (2017b) — A simple and robust rule-based EEG arousal detector

## TL;DR
Rule-based (no classifier training) arousal detector using one EEG derivation (C4/A1), chin EMG and ECG: adaptive ECG-artifact removal → candidate events from abrupt alpha/beta power increases (STFT, 3 s windows, 0.2 s step, vs 10 s baseline) → event refinement by EEG power, EEG amplitude and EMG patterns → removal of false positives by AASM rules (3–15 s, spindle activity, REM EMG, 10 s stable sleep). 22 validation PSGs: κ = 0.775 (0.654–0.886), precision 0.86, F1 0.79. **This is the algorithm (ref. 16) that Alvarez-Estevez 2019 updated and validated at scale**; open-source code at github.com/bigsasi/arousals-detection.

## Research question & hypotheses
Can a simple, computationally light, artifact-tolerant signal-processing pipeline detect EEG arousals at expert level? (Introduction)

## Sample & design
- 30 in-lab clinical PSGs (Haaglanden Medisch Centrum; METC 16-027): 16 M / 14 F, age 50.4 (20.8), AHI 7.8 (8.0), ArI 13.9 (11.2); 16 with major OSAS component, others complex/central SAS, hypersomnia/narcolepsy, RBD, NREM parasomnia, 6 without major disorder. No rejection for artifacts. (§3)
- SOMNOscreen plus, EDF+, **all signals 256 Hz**; C4/A1, submental EMG, modified lead II ECG. AASM scoring + consensus review by a second expert. (§3)
- Split: 6 recordings for development/parameter setting, 22 for independent validation (6 + 22 = 28 ≠ 30 — inconsistency in the text). (§3)

## Methods (replication-level detail)
### 1. Signal conditioning (§2.1)
- Notch filter all signals at mains (50 Hz in their data); **EMG high-pass 15 Hz**. Filter implementations in ref. 20 (Alvarez-Estevez thesis).
- QRS detection: standard algorithm (ref. 21). **ECG reliability** per 30 s window: Hamilton-style "lump" signal (notch + 8–16 Hz band-pass, amplitude normalised; ref. 23); peaks = maxima of regions above the mean of 10 s non-overlapping windows; r = mean(A)/std(A) of peak amplitudes; reliable if r ≥ 0.05 (the comparison symbol was lost in conversion — "r  0:05"; ≥ assumed).
- **Adaptive ECG artifact removal** (§2.1.2): for each reliable R-peak, 0.5 s window W(t) centred on the R-peak in the EEG/EMG; template update T ← (1 − α)T + αW with α = 0.01 (formula not decoded — standard exponential-memory form assumed); template Hann-tapered and subtracted. Unreliable beats: W = 0 and no subtraction.

### 2. Candidate events (§2.2)
| Parameter | Value |
|---|---|
| Power estimate | STFT, 3 s Hamming window, 0.2 s step; band power = mean squared periodogram in band |
| Bands | alpha 8–12 Hz, beta > 16 Hz (upper edge not stated; Nyquist 128 Hz) |
| Baseline | mean band power of the previous 10 s of windows |
| Trigger | alpha power > 2.5 × baseline; beta power > 2 × baseline |
| Event | onset t1 / offset t2 = threshold crossing period |
Theta is deliberately not used (low specificity, §4).

### 3. Pattern recognition / event refinement (§2.3)
**EEG power-based (§2.3.1)**: V_ref = mean power over the t_ref pre-event windows (t_ref = 10 s beta, 3 s alpha). V_f(k) = mean power in the k-th 1 s block after t2. While V_f(k) > th·V_ref extend the end; stop after n consecutive failures (beta th = 3, n = 3; alpha th = 4, n = 1). Validation of the extension: mean event power or ≥ 1 future value > th·V_ref (th = 5.5 beta, 6.5 alpha), and EEG peak-to-peak in the event ≤ 8 × peak-to-peak of the previous 10 s (artifact guard).
**EEG amplitude-based (§2.3.2)**: if event p2p / mean p2p of the previous 5 non-overlapping 1 s windows > 4 (alpha and beta), append the next 1 s window; repeat while true.
**EMG (§2.3.3)**:
- Overlapping activity: V_ref = mean p2p of 0.1 s windows in [t1 − 15 s, t2 + 15 s] excluding the event; if V_ev > 1.4·V_ref, append consecutive 0.1 s windows while the condition holds.
- Non-continuous activity: rectified EMG, 3 s windows, 0.5 s step over [t1 − 3.5 s, t2]; s_i = mean(w_i − w_0); mark regions with s_i > 2·s_1 (1 < i ≤ n); μ_diffs, sd_diffs of s_i; append 0.5 s windows while s_n > μ_diffs − sd_diffs (sign ambiguous in conversion: "μdiffs sddiffs").
- Delayed (non-overlapping) activity: reference = mean p2p of 0.125 s windows in the 10 s before the event; if two consecutive future 0.5 s windows exceed 2.5 × reference, extend the end; stop when no update within the next 4 s.

### 4. False-positive removal (§2.4)
- Duration: remove < 3 s; remove > 15 s within one epoch (events spanning two epochs allowed up to 30 s); remove short **alpha** events of 3–4.5 s.
- Spindles: EEG band-passed 1–32 Hz; 0.5 s window, 0.125 s step; a window is spindle if its dominant frequency (formula not decoded) is 12–15 Hz; extend while dominant frequency stays in range; keep the arousal only if it contains ≥ 3 s spindle-free.
- REM without EMG: REM estimated per 30 s epoch as low chin tone (rectified EMG p2p < 2.5 × minimum epoch p2p); events in such epochs removed if event EMG p2p < 1.1 × epoch EMG p2p. (With real hypnograms, the scored REM can be used instead.)
- 10 s stable sleep: if an event starts < 10 s after the previous event's end, remove the later one; remove events in epochs scored W (needs a W/sleep hypnogram input).

### Validation (§3)
Per 30 s epoch (arousal assigned to the epoch of its midpoint): sensitivity, specificity, precision, F1, κ; onset/offset/duration differences of matched events (Fig. 5, Table 3).

## Key results
- 22 validation PSGs (Table 1): expert 1698 vs system 1496 arousals; TP 1278, FP 218, FN 420; sens 0.748, spec 0.988, precision 0.855, F1 0.793, **κ 0.775** (per recording 0.654–0.886; 10 recordings > 0.80).
- Timing (Table 3, seconds; minus signs probably lost in conversion for lower adjacent / 25th percentiles): onset median difference 0.38 s, IQR roughly −0.23 to 0.93 s; offset median 0.39 s; duration median 0.14 s. Onsets more precise than offsets (Wilcoxon p < 0.001).
- Only 8 % of epochs contain an arousal — error/specificity are uninformative (Table 2 trivial classifiers).

## Figures
- Fig. 1 pipeline; Fig. 2 ECG artifact removal example; Fig. 3 power-based extension (beta, alpha); Fig. 4 EMG-based extension cases (overlapping, non-continuous, delayed); Fig. 5 expert–method timing differences; Fig. 6 overlap agreement.

## Limitations
- Stated by authors: theta excluded; thresholds tuned on 6 recordings; requires W/sleep hypnogram for the stable-sleep rule.
- Not stated by authors: small development set (n = 6) and single-centre data at 256 Hz; parameter interplay (several hand-tuned thresholds) untested for different sampling rates or populations with high arousal indices; several formulas (template update, dominant frequency) and two inequality signs were not recoverable from the text conversion → cross-check with the GitHub code.

## Relevance to this project
- **This is the detector to port** (user decision: Alvarez-Estevez method). It needs no training, uses exactly the DREAMT channels C4-M1, CHIN and ECG, and can use DREAMT's scored hypnogram for the W and REM rules.
- 100 Hz adaptation: beta band must be capped < 50 Hz (e.g. 16–35/45 Hz) and avoid a possible 40 Hz alias of 60 Hz mains (check spectra); EMG high-pass 15 Hz leaves only 15–50 Hz of EMG; ECG template window 0.5 s = 50 samples (fine).
- Onset precision (~0.4 s median offset vs experts) is adequate for respiratory-phase analysis at ~0.25 Hz breathing (≈ 36° per 0.4 s) but adds phase noise — report it.
- Validation plan: per-subject Arousal Index vs DREAMT clinical ArI (Spearman/ICC), as in Alvarez-Estevez 2019; DREAMT ArI (mean 37.7/h) is far above this development set (13.9/h), so the "10 s stable sleep" rule will remove more events in fragmented nights.
- Where details are ambiguous here (template update, μ − sd extension, dominant-frequency formula), follow the GitHub implementation (bigsasi/arousals-detection) and the 2019 update (tuned thresholds; some EMG procedures dropped).

## Follow-up references
- Ref. 20: Alvarez-Estevez thesis — filter implementations.
- Ref. 21: QRS detection algorithm; ref. 22 Hamilton (lump signal); ref. 23 Alvarez-Estevez et al. — ECG reliability and lump-signal details.
- Alvarez-Estevez & Fernández-Varela 2019 (`AlvarezEstevez_2019`) — updated version and large-scale validation.
