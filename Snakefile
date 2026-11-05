# SnoreDREAMT pipeline — respiration–SO/spindle coupling during apnea, hypopnea and snoring.
# Run:  .venv/bin/snakemake -j 8            (dry run: -n; one subject: derivatives/events/S002/complexes_F4.parquet)
from snoredreamt.io import subject_ids

configfile: "config.yaml"

DER = config["paths"]["derivatives"]
SUBJECTS = subject_ids(config)
CHANNELS = config["event_channels"]

wildcard_constraints:
    sid = r"S\d{3}",
    ch = "|".join(config["channel_alias"])


rule all:
    input:
        f"{DER}/step1/detector_summary.csv",
        f"{DER}/step2/resp_summary.csv",
        f"{DER}/step2/group_stats.csv",
        f"{DER}/step3/group_stats.csv",
        f"{DER}/step3/conditions_summary.csv",
        f"{DER}/step3/event_locked_rates.parquet",
        expand(f"{DER}/step3b/{{ss}}/glmm.csv", ss=["N2", "N2N3"]),
        f"{DER}/step4/glmm_ahi.csv",
        expand(f"{DER}/step5/{{ss}}/group_phase_shifts.csv", ss=["N2", "N2N3"]),
        expand(f"{DER}/step5/{{ss}}/circstat_group.csv", ss=["N2", "N2N3"]),
        expand(f"{DER}/step6/{{ss}}/glmm_ahi.csv", ss=["N2", "N2N3"]),
        expand(f"{config['paths']['figures']}/qc_polarity/{{sid}}.png", sid=SUBJECTS),
        f"{config['paths']['figures']}/report/fig01_density.png",
        f"{config['paths']['figures']}/step3_snoring_validation.png",


# ---------------------------------------------------------------- Step 0: data
rule cache:
    """Parse the raw 1 GB CSV once into float32 µV arrays."""
    output: f"{DER}/cache/{{sid}}.npz"
    resources: io = 1
    shell: "python -c \"from snoredreamt.io import load_subject; load_subject('{wildcards.sid}', use_cache=False)\""


rule quality:
    """Per-channel validity masks (EEG artifacts; flat/saturated/low-amplitude respiration)."""
    input: f"{DER}/cache/{{sid}}.npz"
    output: f"{DER}/quality/{{sid}}.npz"
    script: "scripts/rules/quality.py"


# ---------------------------------------------------------------- Step 1: sleep events
rule detect_sos:
    """slythm SO detector (N2/N3 thresholds), plus REM/W/N1 detections with N2 thresholds for specificity."""
    input:
        cache = f"{DER}/cache/{{sid}}.npz",
        quality = f"{DER}/quality/{{sid}}.npz",
    output:
        events = f"{DER}/events/{{sid}}/sos_{{ch}}.parquet",
        thresholds = f"{DER}/events/{{sid}}/sos_{{ch}}_thresholds.json",
    script: "scripts/rules/detect_sos.py"


rule detect_spindles:
    """slythm spindle detector (RMS envelope), plus REM/W/N1 detections with the NREM threshold."""
    input:
        cache = f"{DER}/cache/{{sid}}.npz",
        quality = f"{DER}/quality/{{sid}}.npz",
    output:
        events = f"{DER}/events/{{sid}}/spindles_{{ch}}.parquet",
        thresholds = f"{DER}/events/{{sid}}/spindles_{{ch}}_thresholds.json",
    script: "scripts/rules/detect_spindles.py"


rule complexes:
    """SO-spindle complexes: SO trough followed by a spindle within 1.5 s (same channel)."""
    input:
        sos = f"{DER}/events/{{sid}}/sos_{{ch}}.parquet",
        spindles = f"{DER}/events/{{sid}}/spindles_{{ch}}.parquet",
    output: f"{DER}/events/{{sid}}/complexes_{{ch}}.parquet"
    script: "scripts/rules/complexes.py"


# ---------------------------------------------------------------- Step 2: respiration
rule resp_phase:
    """Polarity calibration, breath detection, Hilbert and cycle-based phase for FLOW and ABDOMEN."""
    input:
        cache = f"{DER}/cache/{{sid}}.npz",
        quality = f"{DER}/quality/{{sid}}.npz",
    output:
        phase = f"{DER}/resp/{{sid}}_phase.npz",
        breaths = f"{DER}/resp/{{sid}}_breaths.parquet",
        info = f"{DER}/resp/{{sid}}_resp_info.json",
    script: "scripts/rules/resp_phase.py"


# ---------------------------------------------------------------- Step 3: breathing conditions
rule event_locked_rates:
    """SO / spindle / arousal rates around apnea/hypopnea onsets and ends -> post-event exclusion window."""
    input:
        sos = expand(f"{DER}/events/{{sid}}/sos_F4.parquet", sid=SUBJECTS),
        spindles = expand(f"{DER}/events/{{sid}}/spindles_C4.parquet", sid=SUBJECTS),
        arousals = expand(f"{DER}/arousals/{{sid}}_primary.parquet", sid=SUBJECTS),
    params: subjects = SUBJECTS
    output:
        table = f"{DER}/step3/event_locked_rates.parquet",
        fig = f"{config['paths']['figures']}/step3_event_locked_rates.png",
    script: "scripts/rules/event_locked_rates.py"


rule conditions:
    """Per-sample breathing-condition labels (+ per-breath snoring) using arousal variant {arv} for exclusion."""
    input:
        cache = f"{DER}/cache/{{sid}}.npz",
        quality = f"{DER}/quality/{{sid}}.npz",
        breaths = f"{DER}/resp/{{sid}}_breaths.parquet",
        arousals = f"{DER}/arousals/{{sid}}_{{arv}}.parquet",
    output:
        labels = f"{DER}/step3/{{sid}}/conditions_{{arv}}.npz",
        breaths = f"{DER}/step3/{{sid}}/snore_breaths_{{arv}}.parquet",
        info = f"{DER}/step3/{{sid}}/conditions_{{arv}}.json",
    script: "scripts/rules/conditions.py"


rule coupling_step3:
    """Coupling per breathing condition vs matched clean controls (projected coupling on the subject's clean mu)."""
    input:
        labels = f"{DER}/step3/{{sid}}/conditions_{{arv}}.npz",
        phase = f"{DER}/resp/{{sid}}_phase.npz",
        sos = f"{DER}/events/{{sid}}/sos_F4.parquet",
        spindles = f"{DER}/events/{{sid}}/spindles_C4.parquet",
        complexes = f"{DER}/events/{{sid}}/complexes_F4.parquet",
    output:
        table = f"{DER}/step3/{{sid}}/coupling_{{arv}}.parquet",
        events = f"{DER}/step3/{{sid}}/event_phases_{{arv}}.parquet",
    script: "scripts/rules/coupling_step3.py"


rule step3_group:
    """Group statistics + figures for coupling per breathing condition (primary + sensitivity grid)."""
    input:
        primary = expand(f"{DER}/step3/{{sid}}/coupling_primary.parquet", sid=SUBJECTS),
        liberal = expand(f"{DER}/step3/{{sid}}/coupling_liberal.parquet", sid=SUBJECTS),
        events_primary = expand(f"{DER}/step3/{{sid}}/event_phases_primary.parquet", sid=SUBJECTS),
        step2 = expand(f"{DER}/step2/{{sid}}/summary.parquet", sid=SUBJECTS),
        resp_summary = f"{DER}/step2/resp_summary.csv",
    output:
        table = f"{DER}/step3/group_stats.csv",
        mixed = f"{DER}/step3/mixed_model.csv",
        wide = f"{DER}/step3/subject_level.parquet",
        forest = f"{config['paths']['figures']}/step3_forest.png",
        paired = f"{config['paths']['figures']}/step3_paired.png",
        grid = f"{config['paths']['figures']}/step3_sensitivity_grid.png",
    script: "scripts/rules/step3_group.py"


rule coupling_step3b:
    """Step 3b data: event sets (large / multichannel consensus), best phase, breath predictors, phase-bin counts."""
    input:
        labels = f"{DER}/step3/{{sid}}/conditions_primary.npz",
        breaths = f"{DER}/step3/{{sid}}/snore_breaths_primary.parquet",
        phase = f"{DER}/resp/{{sid}}_phase.npz",
        cache = f"{DER}/cache/{{sid}}.npz",
        events = expand(f"{DER}/events/{{{{sid}}}}/{{kind}}_{{ch}}.parquet", kind=["sos", "spindles"],
                        ch=config["event_channels"]),
        snore_breaths = f"{DER}/step3/{{sid}}/snore_breaths_primary.parquet",
    output:
        events = f"{DER}/step3b/{{sid}}/events.parquet",
        bins = f"{DER}/step3b/{{sid}}/phase_bins.parquet",
        info = f"{DER}/step3b/{{sid}}/info.json",
    script: "scripts/rules/coupling_step3b.py"


rule step3b_group:
    """Step 3b/3c: clean-signal selection, contrasts (event-level + occupancy-normalised), dose-response; per stage set."""
    input:
        events = expand(f"{DER}/step3b/{{sid}}/events.parquet", sid=SUBJECTS),
        bins = expand(f"{DER}/step3b/{{sid}}/phase_bins.parquet", sid=SUBJECTS),
        resp_summary = f"{DER}/step2/resp_summary.csv",
    output:
        clean_signal = f"{DER}/step3b/{{ss}}/clean_signal.csv",
        contrasts = f"{DER}/step3b/{{ss}}/contrasts.csv",
        dose = f"{DER}/step3b/{{ss}}/dose_response.csv",
        snoring_diag = f"{DER}/step3b/{{ss}}/snoring_diagnostics.csv",
        contrasts_occ = f"{DER}/step3b/{{ss}}/contrasts_occupancy.csv",
        subjects = f"{DER}/step3b/{{ss}}/subject_coupling.csv",
        bins_export = f"{DER}/step3b/{{ss}}/phase_bins.csv",
        fig = f"{config['paths']['figures']}/step3b_power_{{ss}}.png",
    wildcard_constraints: ss = "N2|N2N3"
    script: "scripts/rules/step3b_group.py"


rule step3b_glmm:
    """Idea 4: hierarchical Poisson GLMM of event counts per phase bin (lme4), relative change in modulation."""
    input: f"{DER}/step3b/{{ss}}/phase_bins.csv"
    output: f"{DER}/step3b/{{ss}}/glmm.csv"
    wildcard_constraints: ss = "N2|N2N3"
    shell: "Rscript scripts/rules/step3b_glmm.R {input} {output}"


# ---------------------------------------------------------------- Step 5: preferred-phase changes
rule step5_phase:
    """Step 5: any phase change (5a, parametric occupancy-aware null) and consistent rotation (5b) per condition."""
    input:
        bins = f"{DER}/step3b/{{ss}}/phase_bins.csv",
        glmm = f"{DER}/step3b/{{ss}}/glmm.csv",
    output:
        subjects = f"{DER}/step5/{{ss}}/subject_phase_shifts.csv",
        table = f"{DER}/step5/{{ss}}/group_phase_shifts.csv",
        fig = f"{config['paths']['figures']}/step5_phase_{{ss}}.png",
    wildcard_constraints: ss = "N2|N2N3"
    script: "scripts/rules/step5_phase.py"


rule step5_circstat:
    """Step 5 with CircStat tests on event phases: per-subject Watson-Williams / common-median tests (5a) and the
    population distribution of per-subject mean phase differences vs 0 (5b)."""
    input:
        events = expand(f"{DER}/step3b/{{sid}}/events.parquet", sid=SUBJECTS),
        clean_signal = f"{DER}/step3b/{{ss}}/clean_signal.csv",
    output:
        subjects = f"{DER}/step5/{{ss}}/circstat_subjects.csv",
        table = f"{DER}/step5/{{ss}}/circstat_group.csv",
        fig = f"{config['paths']['figures']}/step5_circstat_{{ss}}.png",
    wildcard_constraints: ss = "N2|N2N3"
    script: "scripts/rules/step5_circstat.py"


# ---------------------------------------------------------------- Step 4: severity
rule severity:
    """Severity markers from the data: ODI3, T90, hypoxic burden, event index/duration, arousal and snore indices."""
    input:
        cache = f"{DER}/cache/{{sid}}.npz",
        arousals = f"{DER}/arousals/{{sid}}_primary.parquet",
        cond_info = f"{DER}/step3/{{sid}}/conditions_primary.json",
        labels = f"{DER}/step3/{{sid}}/conditions_primary.npz",
    output: f"{DER}/step4/{{sid}}_severity.json"
    script: "scripts/rules/severity.py"


rule step4_group:
    """Step 4: severity (AHI + data-derived markers) vs clean coupling; within-subject change vs AHI."""
    input:
        severity = expand(f"{DER}/step4/{{sid}}_severity.json", sid=SUBJECTS),
        subjects = f"{DER}/step3b/N2/subject_coupling.csv",
        bins = f"{DER}/step3b/N2/phase_bins.csv",
    output:
        table = f"{DER}/step4/severity_vs_coupling.csv",
        change = f"{DER}/step4/change_vs_ahi.csv",
        glmm_data = f"{DER}/step4/glmm_data.csv",
        fig = f"{config['paths']['figures']}/step4_severity.png",
    script: "scripts/rules/step4_group.py"


rule step4_glmm:
    """Step 4 hierarchical: AHI as moderator of the respiratory phase modulation (lme4)."""
    input: f"{DER}/step4/glmm_data.csv"
    output: f"{DER}/step4/glmm_ahi.csv"
    shell: "Rscript scripts/rules/step4_glmm.R {input} {output}"


# ---------------------------------------------------------------- Step 6: phase effects vs severity
rule step6_phase_severity:
    """Step 6: within-subject phase changes (6a) and clean preferred phase (6b) vs AHI and severity markers."""
    input:
        severity = expand(f"{DER}/step4/{{sid}}_severity.json", sid=SUBJECTS),
        circstat_subjects = f"{DER}/step5/{{ss}}/circstat_subjects.csv",
        events = expand(f"{DER}/step3b/{{sid}}/events.parquet", sid=SUBJECTS),
        clean_signal = f"{DER}/step3b/{{ss}}/clean_signal.csv",
        resp_summary = f"{DER}/step2/resp_summary.csv",
        bins = f"{DER}/step3b/{{ss}}/phase_bins.csv",
    output:
        table_6a = f"{DER}/step6/{{ss}}/6a_phase_change_vs_severity.csv",
        table_6b = f"{DER}/step6/{{ss}}/6b_clean_phase_vs_severity.csv",
        glmm_data = f"{DER}/step6/{{ss}}/glmm_data.csv",
        fig = f"{config['paths']['figures']}/step6_{{ss}}.png",
    wildcard_constraints: ss = "N2|N2N3"
    script: "scripts/rules/step6_phase_severity.py"


rule step6_glmm:
    """Step 6 population level: AHI as moderator of clean modulation and of condition-related rotation/change."""
    input: f"{DER}/step6/{{ss}}/glmm_data.csv"
    output: f"{DER}/step6/{{ss}}/glmm_ahi.csv"
    wildcard_constraints: ss = "N2|N2N3"
    shell: "Rscript scripts/rules/step6_glmm.R {input} {output}"


rule conditions_summary:
    """Feasibility: minutes per condition and stage per subject (primary arousal variant)."""
    input: info = expand(f"{DER}/step3/{{sid}}/conditions_primary.json", sid=SUBJECTS)
    output:
        table = f"{DER}/step3/conditions_summary.csv",
        fig = f"{config['paths']['figures']}/step3_condition_minutes.png",
    script: "scripts/rules/conditions_summary.py"


# ---------------------------------------------------------------- Step 2b: arousals
rule arousals:
    """Rule-based arousal detector (Fernandez-Varela 2017 / Alvarez-Estevez 2019 port); variant from config."""
    input:
        cache = f"{DER}/cache/{{sid}}.npz",
        quality = f"{DER}/quality/{{sid}}.npz",
    output: f"{DER}/arousals/{{sid}}_{{variant}}.parquet"
    wildcard_constraints: variant = "|".join(config["arousals"]["variants"])
    script: "scripts/rules/arousals.py"


rule coupling_step2:
    """Respiration coupling of SOs, spindles, complexes and arousals in clean NREM (four phase variants)."""
    input:
        cache = f"{DER}/cache/{{sid}}.npz",
        phase = f"{DER}/resp/{{sid}}_phase.npz",
        breaths = f"{DER}/resp/{{sid}}_breaths.parquet",
        info = f"{DER}/resp/{{sid}}_resp_info.json",
        sos = f"{DER}/events/{{sid}}/sos_F4.parquet",
        spindles = f"{DER}/events/{{sid}}/spindles_C4.parquet",
        complexes = f"{DER}/events/{{sid}}/complexes_F4.parquet",
        arousals = f"{DER}/arousals/{{sid}}_primary.parquet",
    output:
        summary = f"{DER}/step2/{{sid}}/summary.parquet",
        rates = f"{DER}/step2/{{sid}}/phase_rates.parquet",
        peth = f"{DER}/step2/{{sid}}/peth.parquet",
        phases = f"{DER}/step2/{{sid}}/event_phases.parquet",
        locked = f"{DER}/step2/{{sid}}/event_locked_breathing.parquet",
    script: "scripts/rules/coupling_step2.py"


rule step2_group:
    """Group replication figures + statistics (SOs, spindles, complexes, arousals, REM control)."""
    input:
        summary = expand(f"{DER}/step2/{{sid}}/summary.parquet", sid=SUBJECTS),
        rates = expand(f"{DER}/step2/{{sid}}/phase_rates.parquet", sid=SUBJECTS),
        peth = expand(f"{DER}/step2/{{sid}}/peth.parquet", sid=SUBJECTS),
        locked = expand(f"{DER}/step2/{{sid}}/event_locked_breathing.parquet", sid=SUBJECTS),
        resp_summary = f"{DER}/step2/resp_summary.csv",
    output:
        table = f"{DER}/step2/group_stats.csv",
        polar = f"{config['paths']['figures']}/step2_polar.png",
        rates_fig = f"{config['paths']['figures']}/step2_phase_rates.png",
        peth_fig = f"{config['paths']['figures']}/step2_peth.png",
        locked_fig = f"{config['paths']['figures']}/step2_event_locked_breathing.png",
        arousal_fig = f"{config['paths']['figures']}/step2_arousal_coupling.png",
    script: "scripts/rules/step2_group.py"


rule resp_qc_figure:
    """Visual polarity check figure for one subject (made for all; review the flagged ones)."""
    input:
        cache = f"{DER}/cache/{{sid}}.npz",
        quality = f"{DER}/quality/{{sid}}.npz",
        info = f"{DER}/resp/{{sid}}_resp_info.json",
    output: f"{config['paths']['figures']}/qc_polarity/{{sid}}.png"
    script: "scripts/rules/resp_qc_figure.py"


rule resp_summary:
    """Across-subject QC of polarity calibration and phase coverage."""
    input: info = expand(f"{DER}/resp/{{sid}}_resp_info.json", sid=SUBJECTS)
    output: f"{DER}/step2/resp_summary.csv"
    script: "scripts/rules/resp_summary.py"


rule detector_summary:
    """Event densities per stage and subject, NREM/REM specificity ratios, waveform/frequency checks."""
    input:
        sos = expand(f"{DER}/events/{{sid}}/sos_{{ch}}.parquet", sid=SUBJECTS, ch=CHANNELS),
        spindles = expand(f"{DER}/events/{{sid}}/spindles_{{ch}}.parquet", sid=SUBJECTS, ch=CHANNELS),
        complexes = expand(f"{DER}/events/{{sid}}/complexes_{{ch}}.parquet", sid=SUBJECTS, ch=CHANNELS),
        quality = expand(f"{DER}/quality/{{sid}}.npz", sid=SUBJECTS),
    output:
        table = f"{DER}/step1/detector_summary.csv",
        fig = f"{config['paths']['figures']}/step1_detector_densities.png",
    script: "scripts/rules/detector_summary.py"


# ---------------------------------------------------------------- scripts that were run by hand, now rules
rule snore_evoked_events:
    """Step 3c: event rates around snore-burst peaks vs the same moment of quiet breaths (chooses exclusion windows)."""
    input:
        summary = f"{DER}/step3/conditions_summary.csv",
        events = expand(f"{DER}/step3b/{{sid}}/events.parquet", sid=SUBJECTS),
        spindles = expand(f"{DER}/events/{{sid}}/spindles_C4.parquet", sid=SUBJECTS),
    output:
        table = f"{DER}/step3b/snore_evoked_events.parquet",
        fig = f"{config['paths']['figures']}/step3b_snore_evoked_events.png",
    shell: "python scripts/snore_evoked_events.py"


rule snoring_validation_figure:
    """Step 3: visual check of the breath-by-breath snoring labels (two subjects)."""
    input:
        cache = expand(f"{DER}/cache/{{sid}}.npz", sid=["S022", "S046"]),
        breaths = expand(f"{DER}/step3/{{sid}}/snore_breaths_primary.parquet", sid=["S022", "S046"]),
    output: f"{config['paths']['figures']}/step3_snoring_validation.png"
    shell: "python scripts/snoring_validation_figure.py"


# ---------------------------------------------------------------- report figures
REPORT_FIGS = dict(density="fig01_density", phase_example="fig02_phase_example", polar_flow="fig03_polar_flow",
                   polar_abd="fig04_polar_abdomen", phase_rates="fig05_phase_rates", locked="fig06_event_locked",
                   contrasts="fig07_contrasts", glmm="fig08_glmm_change", snore_evoked="fig09_snore_evoked",
                   severity="fig10_severity", rotation="fig11_rotation", circstat="fig12_circstat",
                   step6="fig13_phase_vs_ahi")

rule report_figures:
    """Compact, readable figures of the primary analyses for the Quarto report (redrawn from saved tables)."""
    input:
        detector = f"{DER}/step1/detector_summary.csv",
        example_cache = f"{DER}/cache/S047.npz",
        example_phase = f"{DER}/resp/S047_phase.npz",
        example_info = f"{DER}/resp/S047_resp_info.json",
        step2_summary = expand(f"{DER}/step2/{{sid}}/summary.parquet", sid=SUBJECTS),
        step2_rates = expand(f"{DER}/step2/{{sid}}/phase_rates.parquet", sid=SUBJECTS),
        step2_group = f"{DER}/step2/group_stats.csv",
        resp_summary = f"{DER}/step2/resp_summary.csv",
        locked = f"{DER}/step3/event_locked_rates.parquet",
        clean_signal = f"{DER}/step3b/N2/clean_signal.csv",
        contrasts = f"{DER}/step3b/N2/contrasts_occupancy.csv",
        glmm = f"{DER}/step3b/N2/glmm.csv",
        snore_evoked = f"{DER}/step3b/snore_evoked_events.parquet",
        subject_coupling = f"{DER}/step3b/N2/subject_coupling.csv",
        step4 = f"{DER}/step4/severity_vs_coupling.csv",
        step5 = f"{DER}/step5/N2/group_phase_shifts.csv",
        circstat = f"{DER}/step5/N2/circstat_group.csv",
        events = expand(f"{DER}/step3b/{{sid}}/events.parquet", sid=SUBJECTS),
        step6 = f"{DER}/step6/N2/6b_clean_phase_vs_severity.csv",
    output: **{k: f"{config['paths']['figures']}/report/{v}.png" for k, v in REPORT_FIGS.items()}
    script: "scripts/rules/report_figures.py"
