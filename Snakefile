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
        expand(f"{config['paths']['figures']}/qc_polarity/{{sid}}.png", sid=SUBJECTS),


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
