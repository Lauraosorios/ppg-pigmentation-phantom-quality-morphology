"""
segmentation.py -- pulse container and pooled pulse template.
"""

from dataclasses import dataclass

import numpy as np

TEMPLATE_POINTS = 100


@dataclass
class PulseSegment:
    start: int             # onset sample (in the cropped signal)
    end: int               # next onset sample
    ac: np.ndarray         # band-passed AC, non-normalised
    dc: np.ndarray         # low-passed DC
    noise: np.ndarray      # high-passed noise
    ac_norm: np.ndarray    # AC resampled to the recording's median pulse length, scaled to [0, 1]
    ibi_s: float
    ibi_samples: int


def resample_pulse(ac_norm: np.ndarray, n_points: int = TEMPLATE_POINTS) -> np.ndarray:
    """Linear resampling of one normalised pulse to a common length, so pulses
    from different recordings (different median pulse lengths) can be pooled."""
    return np.interp(np.linspace(0, 1, n_points), np.linspace(0, 1, len(ac_norm)), ac_norm)


def pooled_template(pulses: np.ndarray) -> tuple:
    """Median pulse and 5th / 95th percentile envelope of a (n_pulses, n_points) stack."""
    return (np.median(pulses, axis=0),
            np.percentile(pulses, 5, axis=0),
            np.percentile(pulses, 95, axis=0))
