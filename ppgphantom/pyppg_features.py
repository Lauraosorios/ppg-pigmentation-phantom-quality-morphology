"""
pyppg_features.py -- fiducial detection, pulse segmentation and morphology
biomarkers with pyPPG (Goda, Charlton & Behar, Physiol. Meas., 2024).

One pyPPG run per recording x channel produces BOTH:
  - the per-pulse morphology biomarkers (pyPPG's own biomarker functions), and
  - onset-to-onset PulseSegments (AC / DC / noise slices) for the SQIs,
so SQIs and morphology come from the identical segmentation, pulse for pulse.
No additional quality gating is applied.

Requires pyPPG==1.0.73 (see requirements.txt).
"""

import os
import sys

import numpy as np
import pandas as pd
from dotmap import DotMap
from scipy.interpolate import CubicSpline

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import FS, PPG_DIR, recording_id   # noqa: E402

from pyPPG import PPG, Fiducials                # noqa: E402
import pyPPG.preproc as PP                      # noqa: E402
import pyPPG.fiducials as FP                    # noqa: E402
import pyPPG.ppg_bm.ppg_sig as PS               # noqa: E402
import pyPPG.ppg_bm.sig_ratios as SR            # noqa: E402
import pyPPG.ppg_bm.ppg_derivs as PD            # noqa: E402
import pyPPG.ppg_bm.derivs_ratios as DR         # noqa: E402

from ppgphantom.filtering import ac_component, dc_component, noise_component   # noqa: E402
from ppgphantom.segmentation import PulseSegment                               # noqa: E402

# pyPPG biomarker -> name used in the manuscript tables.
RENAMED_BIOMARKERS = {
    "Tsp": "CT",            # divided by Tpi below -> fraction of the pulse interval
    "Ab/Aa": "APG b/a",
    "Tsys/Tdia": "Ts/Td",
    "Adp/Asp": "RI",
    "Tpw50/Tpi": "FWHM",
    "IPA": "IPA",
    "AI": "AI",
    "AGI": "AGI",
    "AUCpi": "AUC",         # already normalised by pyPPG to systolic amplitude x pulse duration
}

# Raw / intermediate pyPPG outputs that are not exported: absolute timings in
# seconds, absolute pulse widths and absolute amplitudes (each kept only through
# its normalised or ratio form), plus IPR (saturated at a ceiling).
EXCLUDED_BIOMARKERS = {
    "Ta", "Tb", "Tc", "Td", "Te", "Tf", "Tu", "Tv", "Tw", "Tp1", "Tp2",
    "Tdp", "Tsys", "Tdia", "Tpi", "Tpp", "deltaT",
    "Tpw10", "Tpw25", "Tpw33", "Tpw50", "Tpw66", "Tpw75", "Tpw90",
    "Tdw10", "Tdw25", "Tdw33", "Tdw50", "Tdw66", "Tdw75", "Tdw90",
    "Tsw10", "Tsw25", "Tsw33", "Tsw50", "Tsw66", "Tsw75", "Tsw90",
    "Asp", "Adp", "Adn", "Aoff", "Au", "Av", "Aw",
    "Tb-c", "Tb-d", "Tp1-dp", "Tp2-dp",   # re-added below divided by Tpi
    "IPR",
}

# Fiducial-time differences without a pyPPG-native normalised form: divided by
# the pulse interval so they stay comparable across heart rates.
NORMALISE_BY_TPI = ["Tb-c", "Tb-d", "Tp1-dp", "Tp2-dp"]


def load_recording(skin: str, rep: int, hr: int, flow: int, ppg_dir: str = PPG_DIR) -> pd.DataFrame | None:
    path = os.path.join(ppg_dir, recording_id(skin, rep, hr, flow) + ".csv")
    return pd.read_csv(path) if os.path.exists(path) else None


def run_pyppg(signal_v: np.ndarray, fs: int, name: str,
              fL: float = 0.5, fH: float = 8.0, order: int = 4):
    """pyPPG preprocessing -> fiducials -> biomarkers on one continuous signal.
    Returns (ppg_obj, fiducials_obj, bm_vals dict)."""
    s = DotMap()
    s.start_sig = 0
    s.end_sig = len(signal_v)
    s.v = signal_v
    s.fs = fs
    s.name = name
    s.filtering = True
    s.fL = fL
    s.fH = fH
    s.order = order
    s.sm_wins = {"ppg": 50, "vpg": 10, "apg": 10, "jpg": 10}

    correction = pd.DataFrame()
    correction.loc[0, ["on", "dn", "dp", "v", "w", "f"]] = True
    s.correction = correction

    prep = PP.Preprocess(fL=s.fL, fH=s.fH, order=s.order, sm_wins=s.sm_wins)
    s.ppg, s.vpg, s.apg, s.jpg = prep.get_signals(s=s)

    ppg_obj = PPG(s=s, check_ppg_len=False)
    fp = Fiducials(fp=FP.FpCollection(s=ppg_obj).get_fiducials(s=ppg_obj))

    bm_vals = {}
    for get_fn in (PS.get_ppg_sig, SR.get_sig_ratios, PD.get_ppg_derivs, DR.get_derivs_ratios):
        _, df_bm, _ = get_fn(ppg_obj, fp)
        bm_vals.update(df_bm)

    return ppg_obj, fp, bm_vals


def _biomarker_table(bm_vals: dict) -> pd.DataFrame:
    """Per-pulse biomarker table (pyPPG pulse index), renamed / normalised."""
    out = pd.DataFrame()
    for key, name in RENAMED_BIOMARKERS.items():
        if key in bm_vals and key != "AUCpi":
            out[name] = pd.Series(bm_vals[key])
    if out.empty:
        return out

    # The first detected pulse is dropped: with no signal before it, pyPPG can
    # anchor its onset at sample 0 and return a merged, over-long first pulse.
    out = out.iloc[1:]
    if out.empty:
        return out

    tpi = pd.Series(bm_vals["Tpi"]) if "Tpi" in bm_vals else None
    if "CT" in out.columns and tpi is not None:
        out["CT"] = out["CT"] / tpi
    if "AUCpi" in bm_vals:
        out["AUC"] = pd.Series(bm_vals["AUCpi"])

    for key, val in bm_vals.items():
        if key in RENAMED_BIOMARKERS or key in EXCLUDED_BIOMARKERS:
            continue
        out[key] = pd.Series(val)
    if tpi is not None:
        for key in NORMALISE_BY_TPI:
            if key in bm_vals:
                out[f"{key}/Tpi"] = pd.Series(bm_vals[key]) / tpi
    return out


def extract_recording_channel(raw: np.ndarray, name: str, fs: int = FS):
    """One channel of one recording -> (morphology_df, segments), row-for-row
    aligned, or (None, None) if pyPPG found no usable pulses."""
    ac, dc, noise = ac_component(raw), dc_component(raw), noise_component(raw)
    n = min(len(ac), len(dc), len(noise))
    ac, dc, noise = ac[:n], dc[:n], noise[:n]

    # pyPPG sees the correctly oriented AC recombined with the DC baseline.
    _, fp, bm_vals = run_pyppg(ac + dc, fs, name)
    fp_df = fp.get_fp()

    out = _biomarker_table(bm_vals)
    if out.empty:
        return None, None

    # Biomarker series are usually one pulse shorter than the fiducial table
    # (pyPPG drops the trailing incomplete pulse): keep only matching pulses.
    fp_valid = fp_df.loc[fp_df.index.isin(out.index)]

    ibis = (fp_valid["off"] - fp_valid["on"]).to_numpy(dtype=float) / fs
    ibis = ibis[np.isfinite(ibis) & (ibis > 0)]
    if len(ibis) == 0:
        return None, None
    med_ibi_samples = max(int(round(float(np.median(ibis)) * fs)), 4)

    segments, keep_idx = [], []
    for idx, row in fp_valid.iterrows():
        onset, offset = row["on"], row["off"]
        if not (np.isfinite(onset) and np.isfinite(offset)):
            continue
        onset, offset = int(onset), int(offset)
        if offset <= onset or offset > n or (offset - onset) < 4:
            continue

        ac_seg, dc_seg, noise_seg = ac[onset:offset], dc[onset:offset], noise[onset:offset]
        try:
            resampled = CubicSpline(np.arange(len(ac_seg), dtype=float), ac_seg)(
                np.linspace(0.0, len(ac_seg) - 1.0, med_ibi_samples))
        except Exception:
            continue

        p_min, p_max = resampled.min(), resampled.max()
        ac_norm = ((resampled - p_min) / (p_max - p_min)
                   if p_max > p_min else np.full(med_ibi_samples, np.nan))

        # Pulse arrays are stored at float32 precision, as in the HDF5 pulse
        # cache the published SQI values were computed from.
        f32 = lambda a: np.asarray(a, dtype=np.float32).astype(np.float64)   # noqa: E731
        segments.append(PulseSegment(
            start=onset, end=offset, ac=f32(ac_seg), dc=f32(dc_seg), noise=f32(noise_seg),
            ac_norm=f32(ac_norm), ibi_s=(offset - onset) / fs, ibi_samples=offset - onset,
        ))
        keep_idx.append(idx)

    if not segments:
        return None, None

    return out.loc[keep_idx].reset_index(drop=True), segments
