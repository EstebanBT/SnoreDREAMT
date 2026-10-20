---
key: BullonTarraso_2026
title: "The breathing cycle gates memory reactivation during human NREM sleep"
authors: [Esteban Bullón Tarrasó, Tobias Staudigl, Thomas Schreiner]
year: 2026
journal: unpublished manuscript (LMU München; no DOI)
doi: null
type: empirical
tags: [respiration, slow oscillations, spindles, SO-spindle coupling, respiratory phase, closed-loop TMR, memory reactivation, MVPA decoding, NREM sleep, thermistor]
fulltext: BullonTarraso_2026.fulltext.md
---

# Bullón Tarrasó, Staudigl & Schreiner (2026): Breathing cycle gates memory reactivation in NREM sleep

> Unpublished manuscript (2026). Data and code are "available upon publication" (Data/Code availability). Numbers below come from the docling conversion; garbled parts are flagged at the end.

## TL;DR
The authors used respiration-locked closed-loop targeted memory reactivation (TMR) in 25 healthy young adults. Cues were delivered either toward each person's SO-spindle-favoured respiratory phase ("in-phase") or toward the opposite phase ("antiphase"). In-phase cueing improved associative memory retention relative to antiphase cueing (post/pre recall 102.12 % vs. 95.84 %). It also tightened spindle nesting in SO up-states and produced category-specific EEG reactivation locked to the evoked SO-spindle complex. Antiphase cues produced reactivation about 2 s later, at a similar respiratory phase and during a later SO-spindle complex. The authors conclude that respiration acts as a causal pacemaker of SO-spindle windows rather than a correlate. The paper also replicates the endogenous respiration–SO/spindle coupling of Schreiner et al. (2023) in a new sample (see `Schreiner_2023`).

## Research question & hypotheses
- Does respiratory phase *causally* shape SO-spindle coupling, memory reactivation and sleep-dependent memory retention, or does it only track oscillatory dynamics? (Introduction)
- Prediction: cues timed to the respiratory phase where SO-spindle complexes cluster lead to better consolidation, tighter SO-spindle coupling and better-timed reactivation than antiphase cues. (Introduction; Results)

## Sample & design
- **Recruitment and exclusions:** 36 volunteers were recruited. 2 withdrew after the nap. 1 was not invited because the nap yielded too little sleep. 5 were excluded during collection: technical errors n = 3, medical n = 1, and **predominant mouth breathing during sleep n = 1**. A further 3 were excluded after inspection: poor data quality n = 1, and no discernible cue-evoked ERP n = 2. (Methods › Participants)
- **Final sample:** N = 25, age 24.00 ± 3.37 y, 17 female. All were German native speakers with no respiratory or neurological disorders, including sleep apnea. Screening used PSQI, the Morningness-Eveningness Questionnaire and a custom health/stimulant questionnaire. None was on medication and none had done night-shift work in the previous 8 weeks. They abstained from alcohol the evening before and from caffeine on the test day. (Methods › Participants)
- **Sleep-disordered breathing screen:** everyone scored **STOP-BANG < 3**, i.e. low OSA risk. The authors also report no increase in awakenings and no breathing cessations on visual inspection, so this is explicitly a **non-SDB sample**. (Methods › Participants)
- **Adaptation ("habituation") nap:** held about 1 week before the main session at a self-chosen time, after the participant cut about 1 h of sleep the night before. It gave a 60-min opportunity with mean NREM = 31.480 ± 2.372 min. Respiration and 16-channel EEG were recorded (F3, F1, Fz, F2, F4, FC1, FCz, FC2, C3, Cz, C4, O1, Oz, O2, M1, M2). The nap served to estimate the preferred phase. (Methods › Habituation nap)
- **Main night:**
  - Set-up started about 20:00. Participants encoded 90 verb–image associations (objects vs. scenes), then took the pre-sleep test.
  - Lights out was about 23:30 with a 7–8 h sleep opportunity.
  - TMR was given during the first hour of stable N2/N3.
  - A 90-min post-TMR NREM period with no TMR served as the "no-TMR" reference.
  - The post-sleep test and an object/scene localizer followed in the morning. (Methods › Experimental task; Closed-loop TMR; Results)
- **Recording:** eego 65-channel EEG (ANT Neuro), impedances < 20 kΩ, online reference CPz, **512 Hz**. Respiration came from a **Braebon cTherm nasal cannula thermistor** under the nostrils, taped to the cheeks and recorded on the same amplifier, so it is synchronised with the EEG. The text spells the device both "Braedon" and "Braebon". (Methods › EEG and respiratory recordings)
- **Sleep architecture (Supp. Table 1, mean ± SD):**

  | Measure | Duration (min) | Share of night |
  |---|---|---|
  | Wake | 41.82 ± 40.77 | 9.13 % |
  | N1 | 15.30 ± 13.41 | 3.32 % |
  | N2 | 238.82 ± 41.29 | 53.22 % |
  | N3 | 67.56 ± 30.70 | 15.03 % |
  | REM | 77.70 ± 18.88 | 17.49 % |
  | Total sleep time | 399.38 ± 46.18 | 89.06 % |
  | Wake after sleep onset | 24.46 ± 37.91 | 5.47 % |

## Methods (replication-level detail)

### Behavioural task (Methods › Experimental task)
- **Stimuli:** 300 German verbs and 240 images (objects/scenes from Konkle et al. 2010), used across sessions. 90 associations were encoded.
- **Software:** Psychophysics Toolbox 3, MATLAB R2018b.
- **Encoding:**
  - A fixation cross (1.5 ± 0.1 s) was followed by the spoken verb and then the image for 4 s.
  - Participants rated their mental image as realistic or bizarre.
  - The encoding block was run twice in randomised order.
- **Memory tests (pre and post):**
  - All encoded verbs plus 50 foils.
  - Each verb was followed by a 3 s delay, then an old/new response.
  - After "old", participants typed a description of the image. A response counted as correct if it matched the familiarization caption or unambiguously described the image.
- **Localizer:** 180 novel images (90/90), each shown 2.5 ± 0.1 s and presented twice, with an old/new task.

### Preferred-phase estimation (adaptation nap) (Methods › Preferred phase computation)
1. **Extract phases:** for each detected SO and SO-spindle complex in N2 and N3, take the respiratory phase at each frontal and fronto-central electrode. A complex is an SO followed by a spindle within 1.5 s of the SO trough.
2. **Rayleigh test** per event type and per electrode.
3. **Preferred phase:** the circular mean of SO phases at the electrode with a significant distribution and the highest z.
   - If no electrode is significant for SOs, use SO-spindle complexes instead.
   - If neither is significant, use the event type and electrode with the highest z overall.
4. **Cue target phase:** the preferred phase plus a phase advance equivalent to **500 ms**, computed from the participant's mean respiratory frequency in N2/N3. This compensates for stimulus delivery and neural processing delays.

### Closed-loop TMR (Methods › Closed-loop TMR; EEG and respiration pre-processing)
- **Conditions:**
  - In-phase cues at the target phase.
  - Antiphase cues at the target phase ± π.
  - Uncued items.
  - 30 task-irrelevant control verbs (15 in-phase, 15 antiphase).
- **Item assignment:** items were assigned so that correctly recalled pre-sleep items were balanced across conditions. Per-condition item counts are not stated explicitly; 30/30/30 is presumably intended.
- **Software and timing:** OpenViBE. A cue fired when the online phase reached the target **and** at least **5 s** had passed since the previous cue. In-phase and antiphase cues were pseudorandomised.
- **Volume and arousals:** 40–65 dB, starting at 50 dB. Cueing was paused on signs of arousal and stopped after a cumulative 60 min of TMR.
- **Online phase estimation:**
  1. One-pass 2nd-order Butterworth band-pass, **0.05–0.5 Hz**.
  2. Take the last **16 s** of signal.
  3. Compute the median inhalation-peak→exhalation-trough (or trough→peak) half-cycle duration of the recent cycles.
  4. Interpolate the current phase linearly from the time elapsed since the last peak or trough.

### Offline pre-processing (Methods › EEG and respiration pre-processing)
- **Software:** MATLAB 2020a and FieldTrip (v.16/07/2023).
- **EEG:**
  - Downsampled to **256 Hz**.
  - Line-noise filter around 50 Hz (4th-order, two-pass Butterworth; stop-band width not given).
  - Re-referenced to the **average of the mastoids**.
  - Noisy channels were rejected visually and interpolated from neighbours.
- **TMR epochs:** −2 to 5 s around cue onset.
- **Sleep staging:** AASM (Iber 2007), scored by 2 independent raters.
- **Offline respiratory phase:**
  - Band-pass **±0.1 Hz around the individual peak respiratory frequency**. The filter is a two-pass Butterworth; its order is not given.
  - Phase is the Hilbert angle.
  - **Convention: 0° = inhalation peak (inhalation→exhalation transition); ±180° = exhalation trough (exhalation→inhalation transition).**

### Event detection (Methods › Event detection)

| Parameter | Slow oscillations | Spindles |
|---|---|---|
| Band-pass | 0.3–1.5 Hz | 12–16 Hz |
| Criterion | Successive positive→negative zero crossings 0.8–2 s apart | RMS envelope (200 ms sliding window) > mean envelope + 50 % |
| Amplitude / duration | Trough-to-peak amplitude > 2 × SD of the signal **within the extracted epoch** | Duration 0.5–3 s; threshold relative to "the analyzed epoch" |
| Cue association | SO trough 0.5–1.5 s after cue, in ≥ half of the "most prominent SO channels" (chosen from no-TMR detection, Supp. Fig. 2) | Spindle onset ≥ 0.5 s after cue and amplitude maximum within 3 s, in ≥ 1 prominent spindle channel |
| SO-spindle complex | SO followed by a spindle within **1.5 s of the SO trough** | |

- Thresholds are **epoch-relative**, not whole-night or stage-relative. That suits cue-locked epochs, but the text does not say how the "epoch" is defined for the continuous nap and no-TMR detection.

### Preferred-phase correction (post hoc) (Methods › Preferred phase correction)
- The preferred phase was re-estimated from SO-spindle complexes in the 90-min no-TMR NREM period after TMR.
- If this overnight estimate was closer to the antiphase target than to the in-phase target, **all cue labels for that participant were swapped. This happened for 8 of 25 participants.**

### Phase-amplitude coupling (Methods › Phase-amplitude coupling)
- SO phase came from the 0.3–1.5 Hz Hilbert transform and spindle amplitude from the 12–16 Hz Hilbert envelope.
- Estimator: Özkurt & Schnitzler (2011).
- Time-resolved: 1 s windows with a 50 ms step, baseline-corrected by subtracting the mean of the 2 s before the cue. Computed per channel and per condition.

### Peri-event histograms (Methods › Peri-event histograms)
- **Peri-event time histograms (PETHs):**
  - Targets: SO troughs, spindle "time points of minimum amplitude", and complex troughs.
  - 50 ms bins over −0.2 to 2.5 s (cue-locked) or −1 to 2.5 s (SO-trough-locked).
  - Normalised by total events in the window, smoothed with a 250 ms moving average, and z-scored to a 500 ms pre-cue baseline. Done per participant and electrode.
- **Phase histograms:**
  - SO phase from the Hilbert transform, with **0° = SO up-state**.
  - The SO phase at each "spindle minimum value" was binned at **π/60** (3°).
  - Rates were normalised to the mean spindle rate in the 2 s before the cue and smoothed over 9 bins (3π/20).
  - Restricted to the electrodes and time window of the significant PAC cluster.

### Multivariate decoding (Methods › Multivariate analysis)
- **Classifier:** MVPA-Light with linear discriminant analysis (LDA).
- **Localizer decoding:**
  - z-score across trials at each time point, then PCA keeping components that explain 95 % of variance.
  - 200 ms running-average smoothing. Channels are the features.
  - One classifier per time point, scored by AUC.
  - 5-fold cross-validation repeated 5 times; 250 label permutations for the null distribution.
  - Localizer epochs run −0.2 to 2.5 s.
- **Sleep decoding (temporal generalisation):**
  - Train on the localizer, then test on sleep epochs of −1 to 2.5 s around the trough of the cue-elicited SO-spindle complex.
  - Only cues for associations remembered before sleep are used, and no cross-validation is needed because train and test sets are independent.
  - Significance is tested against condition-specific label-shuffled null distributions.
- **Fig. 3b:** per-trial decision values averaged over significant localizer time points, re-aligned to *later* complexes 1.5–3.5 s post-cue, and compared with trials that had no late complex.

### GLMM (Methods › Generalized Linear Mixed-Effects Model)
- Binomial model with logit link. The outcome is post-sleep recall per verb, using items remembered before sleep only.
- The predictor is the interval from the late SO-spindle trough to the next cue onset.
- Random intercepts for participant and verb.
- The equation is not rendered in the conversion. From the text it is: logit P(R_ij = 1) = β0 + β1·Interval_ij + p_i + v_j, with p_i ~ N(0, σ_i²) and v_j ~ N(0, σ_j²).

### Statistics (Methods › Statistical analyses)
- Behaviour: repeated-measures ANOVA, then paired t-tests.
- Spindle rate vs. memory: ANCOVA with spindle rate as covariate.
- Circular statistics: CircStat toolbox (Berens).
  - Rayleigh test for non-uniformity.
  - V-test with an a priori mean of 0 rad, applied to within-participant phase *differences*, to test whether two conditions share the same phase.
- Time, phase and time-frequency data: FieldTrip cluster-based permutation tests.
  - Dependent-samples t at the sample level with a cluster-forming threshold of p = 0.05.
  - "maxsum" cluster statistic; Monte Carlo test at α = 0.05, two-tailed.
  - The number of permutations is not stated.

## Key results
- **Endogenous coupling replicates Schreiner 2023 (no-TMR 90 min; Results › Validation; Fig. 1c; Supp. Fig. 2):**
  - **Modulation index:** respiratory phase modulates delta (< 4 Hz, frontal) and sigma (13–17 Hz, parietal) amplitude, z > 2.58, p < 0.01. The main text gives the bands as SO < 2 Hz and sigma 12–16 Hz.
  - **Inhalation-peak vs. exhalation-trough-locked time-frequency maps:** more SO-range power (p = 0.029) and sigma power (p < 0.001) at the inhalation peak.
  - **Event phases** (circular mean ± error, Rayleigh test across participants; 0° = inhalation peak):

    | Event | Electrode | Preferred phase | Rayleigh z | p |
    |---|---|---|---|---|
    | SOs | AF8 | α = −81.978 ± 14.559° | 6.660 | < 0.001 |
    | SO-spindle complexes | FPz | α = −22.478 ± 15.944° | 5.688 | 0.003 |
    | Spindles | POz | α = +47.584 ± 12.770° | 8.216 | < 0.001 |

    SOs and complexes fall before the inhalation peak; spindles fall after it.
- **Delivered cue phases (Fig. 1d):** in-phase α = −144.976 ± 17.059°, antiphase α = 35.025 ± 17.059°. Rayleigh z = 5.046, p = 0.005. The authors say the mean offset from the previously reported preferred SO-spindle phase is "approximately π/2", which they consider consistent with the 500 ms advance.
- **Associative memory, post/pre ratio (Results › Respiratory phase during TMR; Fig. 1e; Supp. Table 3):**
  - In-phase 102.12 ± 1.45 %, antiphase 95.84 ± 1.74 %, uncued 100.05 ± 1.47 % (mean ± SEM).
  - Repeated-measures ANOVA: F(2,48) = 4.222, p = 0.018.
  - In-phase > antiphase: t(24) = 2.775, p = 0.008. In-phase vs. uncued: t = 1.003, p = 0.321. Antiphase vs. uncued: t = −1.850, p = 0.070.
  - **Recognition memory:** no effect, F(2,48) = 1.63, p = 0.203.
- **Event induction:** cue-locked rates of SOs, spindles and complexes did not differ between in-phase and antiphase (all p > 0.05; Supp. Fig. 3). Respiratory phase therefore does not change *whether* events are evoked.
  - Supp. Table 2 (reconstructed): about 40 % of cues had an SO, about 57 % a spindle and about 24 % a complex, in both conditions.
- **SO-spindle PAC:** in-phase > antiphase, p = 0.031, **1.273–2.070 s** post-cue. The cluster is frontal and centred on SO up-states. (Fig. 2a)
- **Spindle phase histogram:** more spindles at the SO up-state for in-phase, **−16.5° to 16.5°**, p = 0.009, corrected across phase. (Fig. 2b)
- **Spindle rate predicts memory:** the spindle rate in the significant cluster predicted memory, β = 2.602, p = 0.021, regardless of cue class. (Supp. Fig. 4, ANCOVA)
- **Decoding locked to cue-elicited complexes:** the analysis used 17.47 % of trials, a mean of 55.52 ± 6.51 trials per participant.
  - In-phase: above chance at **0.160–0.691 s** from the SO trough (p = 0.042).
  - Antiphase: above chance only at **1.965–2.461 s** (p = 0.032).
  - Direct in-phase > antiphase contrast: p = 0.029. (Fig. 2c; Supp. Fig. 5b)
- **Common respiratory phase for reactivation:** the within-participant phase difference between the in-phase and antiphase decoding windows clusters around 0 (V-test, V = 17.701, p < 0.001). (Fig. 3a)
- **Late complexes in antiphase trials:** 23.99 ± 1.61 % of remembered trials had a complex 1.5–3.5 s post-cue.
  - Decoding locked to these late complexes: 0.328–0.648 s, p = 0.040. (Fig. 3b)
  - Late complexes show tighter up-state coupling than early antiphase complexes: p = 0.012, −85.5° to −37.5°. The Fig. 3c caption lists two clusters: down-state p = 0.001 and up-state p = 0.012. (Fig. 3c)
  - Their coupling phase matches in-phase complexes (V = 15.275, p < 0.001). (Supp. Fig. 6b)
- **Interference from the next cue (antiphase only; GLMM):** a longer interval from the late complex to the next cue predicted recall, β = 0.230, p = 0.041. Binned recall probability rose with the interval (Pearson r = 0.619, p = 0.006). (Fig. 3e–f)

## Figures
- **Fig. 1:** (a) Design. (b) Cue-locked K-complex-like ERP and SO/sigma power increase. (c) SO (AF8) and spindle (POz) respiratory-phase distributions, no-TMR. (d) Delivered in-phase vs. antiphase cue phases. (e) Relative associative memory by condition.
- **Fig. 2:** (a) Time-resolved SO-spindle PAC, in-phase > antiphase (frontal cluster). (b) Spindle-in-SO-phase histogram. (c) Temporal-generalisation decoding locked to complexes, early for in-phase and about 2 s later for antiphase. (d) Respiratory phase trajectories relative to the in-phase cue angle. (e) SO-trough-locked PETH of complexes.
- **Fig. 3:** (a) Respiratory phase of the decoding windows is matched within participants. (b) Decoding locked to late complexes in antiphase trials. (c) Late vs. early antiphase complex coupling. (d) Late antiphase vs. in-phase coupling. (e–f) Interval to the next cue predicts recall.
- **Supplementary figures:**
  - S1: ERPs and TFRs per cue class; no condition differences.
  - S2: Replication of Schreiner 2023 (modulation index, complex phase at FPz, inhalation- vs. exhalation-locked TFR).
  - S3: Cue-locked PETHs, no differences.
  - S4: Spindle rate vs. memory.
  - S5: Localizer decoding from about 100 ms, plus the in-phase vs. antiphase decoding contrast.
  - S6: Coupling comparisons.

## Limitations
- **Stated by authors:**
  - Scalp EEG only, so the hippocampal ripple contribution is untested; they propose intracranial follow-up.
  - Delayed antiphase reactivation could also reflect weak early reactivation that scalp EEG cannot detect.
  - Data quality concerns: one participant was excluded for mouth breathing, and undiagnosed SDB was ruled out with STOP-BANG, scoring and visual inspection.
  - (Discussion; Methods › Participants)
- **Not stated by authors (my assessment):**
  - **Post hoc relabelling:** in-phase and antiphase labels were swapped for **8 of 25 participants (32 %)**, based on overnight no-TMR data. Phase estimates from about 30 min of nap NREM were therefore often unreliable. The behavioural contrast depends partly on a data-driven relabelling, and the paper reports no analysis without the swap.
  - **Small effect:** about 6 percentage points between in-phase and antiphase. Neither condition differs from uncued, so "in-phase helps" and "antiphase hurts" cannot be separated.
  - **"π/2 offset" claim:** at typical sleep breathing rates (about 0.2–0.3 Hz), a 500 ms advance is about 36–54° (about π/4). The observed in-phase mean (−145°) is about 63° from the SO preferred phase (−82°) and about 122° from the complex phase (−22°). The "approximately π/2" statement is unclear.
  - **Detection thresholds** are relative to cue-locked epochs, so they are not directly transferable to continuous-night detection.
  - **Thermistor caveats:** the signal lags true airflow. "Inhalation peak" in thermistor polarity is device-dependent and not spelled out.
  - **Missing details:** the offline respiratory filter order and the number of cluster permutations are not given. Sigma band definitions differ across analyses (12–16 vs. 13–17 Hz).
  - **Spindle timing:** spindle time is defined as the "minimum amplitude" point, which is ambiguous (see below).
  - **Population:** healthy young adults, so generalisation to clinical populations is untested, which the authors acknowledge only speculatively.

## Relevance to this project
- **Direct motivation.** The Discussion speculates that conditions with disrupted breathing such as **obstructive sleep apnea could impair memory by disturbing the temporal coordination of sleep oscillations and reactivation**. That is the hypothesis our snoring analysis tests. The authors deliberately excluded SDB (STOP-BANG < 3) and a mouth breather, so DREAMT's clinical sample covers exactly the population they left out.
- **Comparison values for the "intact coupling" baseline** (healthy, nasal thermistor, N2+N3, 0° = inhalation peak, ±180° = exhalation trough):

  | Event | Preferred phase | Rayleigh z |
  |---|---|---|
  | SOs | about −82° (before the inhalation peak) | 6.66 |
  | SO-spindle complexes | about −22° | 5.69 |
  | Spindles | about +48° (after the peak) | 8.22 |

  These are group-level tests across 25 participants of participant-level mean phases. The same direction should appear in DREAMT non-snoring N2/N3 if the effect generalises. Use the same convention, or convert explicitly.
- **Respiratory channel choice:**
  - Their sensor is a **nasal thermistor**, so DREAMT's **FLOW (thermistor)** is the like-for-like channel.
  - PTAF (nasal pressure) leads the thermistor in time and differs in shape, so phase angles would shift systematically. Use PTAF only as a sensitivity check, or report both.
  - Their mouth-breather exclusion signals that nasal sensors fail under oral breathing, which is likely common during snoring. Check for flat or low-amplitude FLOW and PTAF during snore segments and consider THORAX/ABDOMEN belts as a phase source robust to the breathing route.
- **Phase extraction recipe to borrow:** band-pass **±0.1 Hz around the individual peak respiratory frequency** (two-pass Butterworth), then the Hilbert angle. Exclude apnea/hypopnea segments, where the phase is meaningless. Narrow-band Hilbert phase will impose a sinusoid on irregular, flow-limited snoring breaths, so also consider a peak/trough-based (interpolated) phase like their online method (median half-cycle duration over a 16 s window).
- **EEG channels:** they find SOs strongest frontally (AF8, frontal PAC cluster), complexes at FPz and spindles parietally (POz). With DREAMT, use **F4-M1 for SOs and complexes** and **C4-M1 for spindles**. Their mastoid re-reference matches DREAMT's M1 derivations. Our 100 Hz sampling is sufficient for the 12–16 Hz spindle band and 0.3–1.5 Hz SOs.
- **Detection parameters to borrow:**
  - SO: 0.3–1.5 Hz, 0.8–2 s between zero crossings.
  - Spindle: 12–16 Hz, 200 ms RMS window, 0.5–3 s duration.
  - Complex: spindle within 1.5 s of the SO trough.
  - Their amplitude thresholds are per cue epoch, so for continuous data use stage-wide or per-recording thresholds (e.g. as in Schreiner 2023) and say so.
- **Little N3:** they pooled N2+N3 for phase estimation; we should do the same. Their sample still had about 68 min of N3.
- **Amount of data matters:** 30 min of nap NREM gave unstable individual preferred phases (labels swapped in 8 of 25 participants once 90 min was available).
  - For within-night snoring vs. non-snoring contrasts, segments may be short, so expect noisy per-participant phase estimates.
  - Prefer group-level tests on participant mean phases, and match event counts between conditions, because Rayleigh z scales with n.
  - Within-participant phase-difference V-tests (expected 0 rad), as in Fig. 3a and Supp. Fig. 6b, are a good template: snore vs. non-snore preferred phase (does the preferred phase stay put?) alongside coupling strength (does concentration drop?).
- **Confound raised by this paper:** auditory cues reliably evoke K-complex-like SOs and spindles *regardless of respiratory phase* (Fig. 1b; Supp. Fig. 3).
  - Snoring is itself a loud, inspiration-locked sound, and respiratory events end in inspiration-locked micro-arousals. Both could evoke K-complexes at a fixed respiratory phase.
  - That would create *artificial* respiration-SO locking during snoring, masking or even reversing a "loss of coupling".
  - Mitigations: separate K-complexes and arousal-related events (e.g. exclude events within a few seconds of arousals or respiratory events, or test SO amplitude and shape distributions), and check whether the SO phase in snoring shifts toward inspiration onset.
- **Interpretive frame:** SO-spindle *precision* (up-state nesting), not event rate, was the respiration-sensitive variable. Report SO-spindle coupling precision (phase-histogram concentration at the up-state, Özkurt PAC) as well as respiration-event phase locking when testing snoring.

## Follow-up references
- **Schreiner et al. 2023, Nat Commun 14, 8351** (`Schreiner_2023`): original respiration–SO/spindle coupling and reactivation paper. Its methods are the baseline this paper replicates.
- **Ghibaudo et al. 2024, Clin Neurophysiol 166, 252–261:** respiratory timing of spindles in humans, an independent replication with a clinical neurophysiology angle.
- **Sheriff et al. 2024, PNAS 121, e2405395121:** breathing orchestrates hippocampal sleep oscillations (intracranial); gives the mechanistic context.
- **Schwimmbeck, Bullón Tarrasó & Schreiner 2025, Trends Neurosci:** review of respiration in sleep-oscillation coordination; likely discusses SDB.
- **Corzo et al. 2026, Prog Neurobiol 256, 102857:** respiratory coordination of excitability across wake and sleep.
- **Casali et al. 2026, Nat Commun 17, 6620:** respiratory pauses and sleep architecture in mice, relevant to apneas and pauses.
- **Karalis & Sirota 2022, Nat Commun 13, 467:** breathing coordinates cortico-hippocampal dynamics in offline states (rodent).
- **Özkurt & Schnitzler 2011, J Neurosci Methods 201, 438–443:** the PAC estimator used here.
- **Helfrich et al. 2018, Neuron 97:221; Muehlroth et al. 2019, Sci Rep 9:1940:** SO-spindle coupling precision and its memory relevance and age effects; relevant because DREAMT is older and clinical.
- **Chung et al. 2016, Chest 149:631:** STOP-BANG screening, if we want a comparable risk score.

---

### Conversion issues / ambiguities in the full text
- **Supplementary Table 2** (cues associated with events) is column-misaligned in the conversion. My reconstruction, verified by encoded + control = condition and in-phase + antiphase ≈ all:

  | Cues | SO (n, %) | Spindle (n, %) | Complex (n, %) |
  |---|---|---|---|
  | All | 187 ± 10.18, 39.59 ± 1.78 % | 271 ± 10.66, 57.33 ± 1.16 % | 115 ± 8.59, 24.38 ± 1.48 % |
  | In-phase | 92.3, 39.15 % | 135, 57.07 % | 56.6, 23.88 % |
  | In-phase encoded | 63.6, 40.31 % | 91.2, 57.61 % | 39.7, 24.98 % |
  | In-phase control | 28.7, 36.80 % | 44.1, 55.98 % | 16.9, 21.65 % |
  | Antiphase | 94.2, 40.01 % | 136, 57.56 % | 58.9, 24.86 % |
  | Antiphase encoded | 63.8, 40.69 % | 91, 57.86 % | 39.8, 25.23 % |
  | Antiphase control | 30.4, 38.68 % | 45.1, 56.97 % | 19.1, 24.11 % |

  This implies about 470 cues per participant in total.
- **GLMM equation** is missing (`<!-- formula-not-decoded -->`); it is reconstructed above from the variable definitions.
- **Typos in numbers:**
  - "t24 = 2,775" in the text means 2.775. Supp. Table 3 gives t = 2.77, p = 0.01, while the text gives p = 0.008.
  - Supp. Table 3 "102.432.62" presumably means 102.43 ± 2.62 (antiphase recognition post/pre).
- **Spindle timing point:** spindles are timed at the "time points of minimum amplitude" and "spindle minimum value". This most likely means the most negative trough of the filtered spindle, but it is not defined.
- **Fig. 3c:** the text gives p = 0.012 for −85.5° to −37.5° and calls it up-state coupling. The caption lists two clusters (down-state p = 0.001, up-state p = 0.012). The window lies on the rising flank before the 0° up-state, so it is unclear which cluster is meant.
- **Supp. Fig. 6b caption:** says "no significant differences (V = 15.275, p < 0.001)". The significant V-test is read as phases being concentrated around a 0 difference, i.e. similar, which is consistent with the text.
- **Reference list:** PDF line numbers (717–838) are interleaved into the references, producing spurious list numbers, and refs 59–65 were converted as a table. The content is readable.
- **Konkle et al. (2010)** is cited by name and corresponds to ref. 60.
- **Figures:** images were not extracted; the Figures section relies on captions.
- **Minor naming inconsistencies:** "Habituation nap" vs. "adaptation nap" (same session). "Preferred-phase" vs. "in-phase" cueing (same condition). "Braedon" vs. "Braebon" (same device).
- **Nap EEG system:** the 16-channel nap montage is listed, but its recording system is not stated. The main night used 65 channels, and the electrodes reported in the results (AF8, POz, FPz) come from that montage.
