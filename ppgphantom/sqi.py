"""
sqi.py -- pulse-level signal quality indices (manuscript Table 1).

Morphology consistency (time + amplitude normalised pulse, ac_norm):
    TMCC, Skewness, Kurtosis, Shannon Entropy, Sample Entropy
Signal strength (non-normalised segments):
    AC/DC, SNR
Spectral characteristics (non-normalised AC segment):
    PSD ratio, Spectral flatness
"""

import os
import sys

import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import kurtosis, skew

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import FS   # noqa: E402

SQI_LABELS = ["TMCC", "AC/DC", "SNR", "Skewness", "PSD ratio", "Shannon Entropy",
              "Kurtosis", "Spec. Flatness", "Sample Entropy"]
SQI_UNITS = ["(-)", "(%)", "(dB)", "(-)", "(-)", "(bits)", "(-)", "(-)", "(nats)"]

_PSD_LF_BAND = (1.0, 2.25)    # cardiac fundamental + lower harmonics (Hz)
_PSD_FULL_BAND = (0.0, 8.0)   # PSD ratio denominator band (Hz)


def _sample_entropy(x: np.ndarray, m: int = 2, r_factor: float = 0.2) -> float:
    """SampEn = -ln(A/B), Chebyshev tolerance r = r_factor * std(x)
    (Richman & Moorman, 2000)."""
    n = len(x)
    r = r_factor * float(np.std(x, ddof=0))
    if r == 0 or n < m + 2:
        return np.nan

    tm = np.lib.stride_tricks.sliding_window_view(x, m)
    tm1 = np.lib.stride_tricks.sliding_window_view(x, m + 1)

    # cdist counts self-matches on the diagonal; subtract them.
    B = int(np.sum(cdist(tm, tm, metric="chebyshev") <= r)) - len(tm)
    A = int(np.sum(cdist(tm1, tm1, metric="chebyshev") <= r)) - len(tm1)

    if B == 0:
        return np.nan
    return float(-np.log(A / B))


def compute_sqi(segments, fs: int = FS):
    """Nine SQIs per pulse.

    Returns (sqi_matrix (n_pulses, 9), template, SQI_LABELS, SQI_UNITS);
    NaN where a metric could not be computed. The template used for TMCC is
    the median of all ac_norm pulses of the same recording.
    """
    if not segments:
        return np.array([]), np.array([]), SQI_LABELS, SQI_UNITS

    n = len(segments)
    norm_stack = np.vstack([s.ac_norm for s in segments])
    template = np.nanmedian(norm_stack, axis=0)

    sqi_matrix = np.full((n, 9), np.nan)

    for i, seg in enumerate(segments):
        ac_norm, ac, dc, noise = seg.ac_norm, seg.ac, seg.dc, seg.noise

        # 1. TMCC
        if not np.any(np.isnan(ac_norm)) and np.std(ac_norm) > 0 and np.std(template) > 0:
            sqi_matrix[i, 0] = float(np.corrcoef(ac_norm, template)[0, 1])

        # 2. AC/DC (%)
        ac_amp = np.max(ac) - np.min(ac)
        if dc.size > 0:
            dc_amp = float(np.abs(np.mean(dc)))
            if dc_amp > 0:
                sqi_matrix[i, 1] = 100.0 * ac_amp / dc_amp

        # 3. SNR (dB): AC (0.5-8 Hz) power vs noise (>8 Hz) power in the same window
        if noise.size > 0:
            p_sig = float(np.sum(ac ** 2))
            p_noise = float(np.sum(noise ** 2))
            if p_noise > 0:
                sqi_matrix[i, 2] = 10.0 * np.log10(p_sig / p_noise)

        # 4. Skewness
        if not np.any(np.isnan(ac_norm)):
            sqi_matrix[i, 3] = skew(ac_norm, bias=True)

        # 5. PSD ratio
        N_fft = 1024
        psd = np.abs(np.fft.fft(ac, N_fft)) ** 2
        freqs = np.arange(N_fft) * fs / N_fft
        lf = float(np.sum(psd[(freqs >= _PSD_LF_BAND[0]) & (freqs <= _PSD_LF_BAND[1])]))
        tp = float(np.sum(psd[(freqs >= _PSD_FULL_BAND[0]) & (freqs <= _PSD_FULL_BAND[1])]))
        if tp > 0:
            sqi_matrix[i, 4] = lf / tp

        # 6. Shannon entropy of the amplitude histogram (32 bins over [0, 1])
        if not np.any(np.isnan(ac_norm)):
            hist, _ = np.histogram(ac_norm, bins=32, range=(0.0, 1.0))
            p_hist = hist.astype(float) / hist.sum()
            p_hist = p_hist[p_hist > 0]
            sqi_matrix[i, 5] = float(-np.sum(p_hist * np.log2(p_hist)))

        # 7. Kurtosis (Pearson, not excess)
        if not np.any(np.isnan(ac_norm)):
            sqi_matrix[i, 6] = float(kurtosis(ac_norm, fisher=False, bias=True))

        # 8. Spectral flatness: geometric / arithmetic mean of the one-sided PSD
        psd_pos = psd[:N_fft // 2]
        psd_pos = psd_pos[psd_pos > 0]
        if len(psd_pos) > 0:
            log_mean = float(np.exp(np.mean(np.log(psd_pos))))
            arith_mean = float(np.mean(psd_pos))
            if arith_mean > 0:
                sqi_matrix[i, 7] = log_mean / arith_mean

        # 9. Sample entropy
        sqi_matrix[i, 8] = _sample_entropy(ac_norm)

    return sqi_matrix, template, SQI_LABELS, SQI_UNITS
