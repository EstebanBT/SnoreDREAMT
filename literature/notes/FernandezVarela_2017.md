---
key: FernandezVarela_2017
title: "Combining machine learning models for the automatic detection of EEG arousals"
authors: [Isaac Fernández-Varela, Elena Hernández-Pereira, Diego Álvarez-Estévez, Vicente Moret-Bonillo]
year: 2017
journal: Neurocomputing 268:100-108
doi: not printed in the text (Elsevier PII S0925-2312(17)30750-6)
type: methods
tags: [EEG arousal, automatic scoring, machine learning, SVM, Hjorth parameters, SHHS, chin EMG]
fulltext: FernandezVarela_2017.fulltext.md
---

# Fernández-Varela et al. (2017) — Combining machine learning models for EEG arousal detection

## TL;DR
Two-stage arousal detector on SHHS data: (1) signal processing finds candidate intervals — abrupt power increases in θ/α/>16 Hz on two EEG derivations (3 s Hamming windows vs a 10 s baseline) and chin-EMG amplitude increases (3 s windows vs a 30 s baseline) — grouped per 30-s epoch; (2) **supervised** classifiers (SVM, ANN, CT, k-NN, …) trained on expert-labelled epochs decide arousal vs non-arousal, and their outputs are combined. Best combined model: sensitivity 0.78–0.81, specificity 0.88–0.89 per epoch. **This is ref. 15 of Alvarez-Estevez 2019, not the rule-based algorithm validated there** (that one is ref. 16, Fernández-Varela et al. 2017, *Comput Biol Med* 87:77–86).

## Research question & hypotheses
Can combining several ML classifiers on signal-processing-derived features improve epoch-level EEG arousal detection, and do Hjorth parameters and sleep stage help as features? (Introduction, §2)

## Sample & design
- SHHS recordings (C3/A2, C4/A1, submental EMG, ECG at 125 Hz), ASDA/AASM expert arousal scoring. (§3)
- Dataset 1: 20 recordings, 2981 arousals / 23,972 epochs → 15 training (2202 arousals, 18,094 epochs) + 5 test (779 arousals, 5878 epochs).
- Dataset 2 (independent, from ref. 21): 26 recordings, 4860 arousals / 31,070 epochs (the "SHHS1-26" set of the 2019 paper).

## Methods (replication-level detail)
### Signal conditioning (§2.1.1)
- ECG artifact removal: locate QRS start/end from the largest values of the ECG first derivative; where an artifact is detected in the EEG, **linearly interpolate** the EEG between QRS start and end (only affected segments). Detection criterion for "artifact present" is not specified.

### Candidate intervals (§2.1.2; parameters from §3)
| Signal | Window | Step | Feature | Baseline | Threshold |
|---|---|---|---|---|---|
| EEG (×2) | 3 s, Hamming, FFT power (formula not decoded in conversion) | 0.2 s | band power in θ 4–7, α 8–12, > 16 Hz (each band band-pass filtered on the whole signal first) | mean of previous 10 s | power / baseline > **1.5** |
| Chin EMG | 3 s | 0.4 s | peak-to-peak amplitude | mean of previous 30 s | amplitude / baseline > **2** |
- Intervals < 3 s discarded (per band / signal).
- Merging rule: consecutive intervals of durations t1, t2 separated by gap t3 are merged into one interval (t1 + t2 + t3) if t1 + t2 > t3, t1 < t3 and t2 < t3 — as printed; note the condition t1 < t3 and t2 < t3 seems at odds with merging short gaps, possibly a typo.
- Upper limit of the "> 16 Hz" band not stated (Nyquist 62.5 Hz at 125 Hz).

### Grouping (§2.1.3)
- Intervals assigned to the 30-s epoch containing their midpoint; a group requires one interval from **each** of the three signals (both EEGs + EMG); if several per signal, take the one with the highest θ/α/>16 Hz power. No group → no arousal possible in that epoch.

### Features & classifiers (§2.2, §3)
- Per EEG interval: total/max/min power in δ 0.5–4, θ, α, σ 12–15, > 16 Hz; Hjorth activity, mobility, complexity; duration (19 features each). EMG: total/max/min amplitude, duration. Context: sleep stage (automatic, ref. 16 of this paper), overlap time of the two EEG intervals.
- Models: FLD, SVM (RBF, C = 2^11, γ "S" = 2^-3), ANN (1 hidden layer, H = 40, conjugate-gradient backprop), CT (prune 15), k-NN (k = 5, Spearman distance), NB (kernel); training set under-sampled to balance classes, 10-fold CV.
- Combinations: Shortliffe–Buchanan certainty factors (categorical outputs mapped to ±0.7; threshold 0.7) and linear combination (equal weights chosen; threshold 0.3). Compared with RF (100 trees) and k-NN ensemble.

## Key results
- Test set (Table 2, all features): SVM error 0.134, sens 0.815, spec 0.874, AUC 0.845; ANN AUC 0.852; Hjorth + stage features improve AUC modestly (e.g. SVM 0.823 without both).
- Independent set (Table 3): S&B combination error 0.124, sens 0.781, spec 0.893, AUC 0.837; linear combination error 0.133, sens 0.810, spec 0.878, AUC 0.844; RF AUC 0.842.
- Upper bounds of the candidate stage (max sensitivity / min specificity) were used to tune thresholds but the values are not reported.

## Figures
- Fig. 1: pipeline overview. Fig. 2: ECG-artifact correction example. Fig. 3: sliding-window band-power detection. Fig. 4: grouping intervals per epoch.

## Limitations
- Stated by authors: ad-hoc artifact conditioning tuned to the available recordings; limited feature-selection study; only SHHS data.
- Not stated by authors: validation is per 30-s epoch (no onset timing accuracy); requires both EEG derivations and EMG to fire in the same epoch; the classifier stage needs **expert-labelled training data** and trained models are not provided; the 2019 paper (Table 7) shows this approach generalises less well than the later rule-based algorithm.

## Relevance to this project
- **Cannot be reproduced as a whole on DREAMT**: DREAMT has no arousal annotations, so the ML stage cannot be trained (only participant-level Arousal Index is available, usable for validation, not for training).
- The **signal-processing front-end is directly usable**, has explicit parameters, and maps onto DREAMT channels: C4-M1 and a second derivation (O2-M1 or CZ-T4; no C3), CHIN, ECG. At 100 Hz the > 16 Hz band must stop below 50 Hz and should avoid a possible 40 Hz mains alias.
- Recommended use: as a fallback/rule-based detector if the code of the validated 2019 algorithm (Comput Biol Med 2017, GitHub bigsasi/arousals-detection) cannot be obtained — candidate intervals + AASM rules (≥ 3 s, ≥ 10 s stable sleep before, EMG increase in REM), validated against each subject's clinical Arousal Index.
- For respiratory-phase analysis, interval onset (power ratio crossing 1.5) gives sub-second timing, unlike the epoch-level outputs of the classifiers.

## Follow-up references
- Fernández-Varela I, Alvarez-Estevez D, Hernández-Pereira E, et al. 2017, *Comput Biol Med* 87:77–86 — "A simple and robust method for the automatic scoring of EEG arousals" (ref. 16 in Alvarez-Estevez 2019): the rule-based algorithm actually validated in 2019; needed for the port.
- Ref. 11 (Shmiel et al. 2009) and ref. 21 (Álvarez-Estévez & Moret-Bonillo, source of the 26-recording set) — comparison baselines.
