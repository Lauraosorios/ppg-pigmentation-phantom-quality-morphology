"""
filtering.py -- AC / DC / noise separation (zero-phase 4th-order Butterworth).

  AC    : mean-subtract -> band-pass 0.5-8 Hz -> polarity flip -> crop edges
  DC    : low-pass 0.5 Hz -> crop edges (absolute level preserved)
  noise : high-pass 8 Hz -> crop edges

The AFE decodes photocurrent with the opposite sign to the usual PPG
convention, so the AC component is negated once here (and only here) to give
a fast systolic upstroke and slower diastolic decay. All three outputs are
cropped by the same number of samples, so they stay aligned sample-for-sample.
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (FS, EDGE_CROP_SECONDS, BP_LOW_HZ, BP_HIGH_HZ,   # noqa: E402
                    FILTER_ORDER, HP_CUTOFF_HZ, DC_LP_HZ)


def bandpass_filter(signal: np.ndarray, fs: int = FS, low: float = BP_LOW_HZ,
                    high: float = BP_HIGH_HZ, order: int = FILTER_ORDER) -> np.ndarray:
    signal = pd.Series(signal).dropna().to_numpy()
    nyq = fs / 2
    b, a = butter(order, [low / nyq, high / nyq], btype="bandpass")
    return filtfilt(b, a, signal)


def highpass_filter(signal: np.ndarray, fs: int = FS, cutoff: float = HP_CUTOFF_HZ,
                    order: int = FILTER_ORDER) -> np.ndarray:
    signal = pd.Series(signal).dropna().to_numpy()
    nyq = fs / 2
    b, a = butter(order, cutoff / nyq, btype="highpass")
    return filtfilt(b, a, signal)


def lowpass_filter(signal: np.ndarray, fs: int = FS, cutoff: float = DC_LP_HZ,
                   order: int = FILTER_ORDER) -> np.ndarray:
    signal = pd.Series(signal).dropna().to_numpy()
    nyq = fs / 2
    b, a = butter(order, cutoff / nyq, btype="lowpass")
    return filtfilt(b, a, signal)


def crop_edges(signal: np.ndarray, crop_seconds: float = EDGE_CROP_SECONDS,
               fs: int = FS) -> np.ndarray:
    crop_samples = int(crop_seconds * fs)
    if len(signal) <= 2 * crop_samples:
        return signal
    return signal[crop_samples:-crop_samples]


def ac_component(raw_signal: np.ndarray) -> np.ndarray:
    x = pd.Series(raw_signal).dropna()
    x = x - x.mean()
    x = bandpass_filter(-x.to_numpy())
    return crop_edges(x)


def dc_component(raw_signal: np.ndarray) -> np.ndarray:
    x = pd.Series(raw_signal).dropna().to_numpy()
    return crop_edges(lowpass_filter(x))


def noise_component(raw_signal: np.ndarray) -> np.ndarray:
    x = pd.Series(raw_signal).dropna().to_numpy()
    return crop_edges(highpass_filter(x))
