"""
config.py -- every path, protocol dimension and processing parameter used by
the pipeline. Nothing else in the repository hard-codes a path or a constant.

The dataset (Zenodo) is NOT part of this repository. Point the code at it
either by placing it in ./data/ or by setting the environment variable
PPG_PHANTOM_DATA to the folder that contains PPG/, Pressure/ and Flow/.
"""

import os

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

DATA_DIR = os.environ.get("PPG_PHANTOM_DATA", os.path.join(REPO_ROOT, "data"))
PPG_DIR = os.path.join(DATA_DIR, "PPG")
PRESSURE_DIR = os.path.join(DATA_DIR, "Pressure")
FLOW_DIR = os.path.join(DATA_DIR, "Flow")

RESULTS_DIR = os.path.join(REPO_ROOT, "results")
TABLES_DIR = os.path.join(RESULTS_DIR, "tables")      # step 01 output (tracked in git)
MANUSCRIPT_TABLES_DIR = os.path.join(RESULTS_DIR, "tables_manuscript")
FIGURES_DIR = os.path.join(RESULTS_DIR, "figures")
HAEMO_DIR = os.path.join(RESULTS_DIR, "haemodynamics")
BEAT_LEVEL_DIR = os.path.join(RESULTS_DIR, "beat_level")   # large, git-ignored

RECORDING_SQI_CSV = os.path.join(TABLES_DIR, "recording_sqi_table.csv")
RECORDING_MORPHOLOGY_CSV = os.path.join(TABLES_DIR, "recording_morphology_table.csv")
TEMPLATES_CSV = os.path.join(TABLES_DIR, "pulse_templates.csv")
BEAT_LEVEL_SQI_CSV = os.path.join(BEAT_LEVEL_DIR, "beat_level_sqi_table.csv")
BEAT_LEVEL_MORPHOLOGY_CSV = os.path.join(BEAT_LEVEL_DIR, "beat_level_morphology_table.csv")

# ---------------------------------------------------------------------------
# Acquisition
# ---------------------------------------------------------------------------

FS = 200                 # Hz -- PPG (AFE4420) and decimated pressure
EDGE_CROP_SECONDS = 5    # removed from each end after zero-phase filtering

# ---------------------------------------------------------------------------
# Filtering (zero-phase 4th-order Butterworth)
# ---------------------------------------------------------------------------

BP_LOW_HZ = 0.5          # AC band-pass lower cutoff
BP_HIGH_HZ = 8.0         # AC band-pass upper cutoff
DC_LP_HZ = 0.5           # DC low-pass cutoff
HP_CUTOFF_HZ = 8.0       # noise high-pass cutoff
FILTER_ORDER = 4

# ---------------------------------------------------------------------------
# Protocol: 3 skin tones x 3 replicates x 3 heart rates x 3 target flows = 81
# ---------------------------------------------------------------------------

SKIN_TONES = ["P", "M", "D"]                 # ordered pale -> dark (JT order)
SKIN_LABELS = {"P": "Pale", "M": "Medium", "D": "Dark"}
REPLICATES = [1, 2, 3]
HEART_RATES = [60, 90, 120]                  # bpm (pump rate)
TARGET_FLOWS = [5, 6, 7]                     # L/min (pump target flow)

# Transmittance GREEN is recorded in the dataset but excluded from all
# analysis: penetration depth leaves no usable pulsatile signal.
CHANNELS = [
    "Reflectance RED",
    "Reflectance IR",
    "Reflectance GREEN",
    "Transmittance RED",
    "Transmittance IR",
]


def recording_id(skin: str, rep: int, hr: int, flow: int) -> str:
    """File stem used throughout the dataset, e.g. 'DR1_60HR_5CO'
    (Dark, replicate 1, 60 bpm, 5 L/min target flow)."""
    return f"{skin}R{rep}_{hr}HR_{flow}CO"


def iter_recordings():
    """Yield (skin, rep, hr, flow) for all 81 recordings in a fixed order."""
    for skin in SKIN_TONES:
        for rep in REPLICATES:
            for hr in HEART_RATES:
                for flow in TARGET_FLOWS:
                    yield skin, rep, hr, flow
