"""
03_haemodynamics.py -- mean arterial pressure (aorta and finger-phantom inlet)
per recording, summarised per HR x target flow state (the ~58-150 mmHg MAP
range reported in Methods).

Pressure was logged by a separate instrument and is NOT time-aligned with the
PPG; only per-recording means are used.

Input  : <data>/Pressure/*.csv
Output : results/haemodynamics/{recording,state}_haemodynamics.csv

Usage:
    python scripts/03_haemodynamics.py
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import FS, PRESSURE_DIR, HAEMO_DIR, iter_recordings, recording_id   # noqa: E402
from ppgphantom.filtering import lowpass_filter, crop_edges                     # noqa: E402


def mean_pressure(signal: np.ndarray) -> float:
    """MAP: 0.5 Hz low-pass -> crop edges -> mean."""
    return float(np.mean(crop_edges(lowpass_filter(signal, fs=FS), fs=FS)))


def main():
    os.makedirs(HAEMO_DIR, exist_ok=True)
    rows = []
    for skin, rep, hr, flow in iter_recordings():
        rid = recording_id(skin, rep, hr, flow)
        path = os.path.join(PRESSURE_DIR, rid + ".csv")
        if not os.path.exists(path):
            print(f"  [missing] {rid}")
            continue
        p = pd.read_csv(path)
        rows.append({"recording_id": rid, "skin": skin, "replicate": rep, "hr": hr, "co_lmin": flow,
                     "aorta_map_mmHg": mean_pressure(p["Aorta_Pressure"].to_numpy()),
                     "finger_map_mmHg": mean_pressure(p["Finger_Pressure"].to_numpy())})

    if not rows:
        sys.exit("No pressure files found -- check that the data folder is set (see README).")

    rec = pd.DataFrame(rows)
    rec.round(3).to_csv(os.path.join(HAEMO_DIR, "recording_haemodynamics.csv"), index=False)

    value_cols = ["aorta_map_mmHg", "finger_map_mmHg"]
    state = rec.groupby(["hr", "co_lmin"])[value_cols].agg(["median", "min", "max"])
    state.columns = [f"{c}_{s}" for c, s in state.columns]
    state = state.reset_index()
    state.round(2).to_csv(os.path.join(HAEMO_DIR, "state_haemodynamics.csv"), index=False)

    print(f"Pressure available for {len(rec)} recordings.")
    print(state[["hr", "co_lmin", "aorta_map_mmHg_median", "finger_map_mmHg_median"]]
          .round(1).to_string(index=False))
    for site in ("aorta", "finger"):
        med = state[f"{site}_map_mmHg_median"]
        print(f"{site.capitalize()} MAP across the 9 states (median of recordings): "
              f"{med.min():.0f} to {med.max():.0f} mmHg")


if __name__ == "__main__":
    main()
