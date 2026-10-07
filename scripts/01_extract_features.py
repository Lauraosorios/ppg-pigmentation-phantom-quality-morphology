"""
01_extract_features.py -- PPG recordings -> pulse-level SQIs and morphology
features -> recording-level tables (median over pulses).

Input  : <data>/PPG/*.csv              (81 offset-corrected recordings)
Output : results/tables/recording_sqi_table.csv
         results/tables/recording_morphology_table.csv
         results/tables/pulse_counts.csv
         results/tables/pulse_templates.csv   (median + 5th/95th percentile pulse per state)
         results/beat_level/beat_level_{sqi,morphology}_table.csv   (large, not in git)

Runs pyPPG on 81 recordings x 5 channels; expect a long run time.
The tables it produces are already included in results/tables/, so
scripts/02_tables_and_figures.py can be run without repeating this step.

Usage:
    python scripts/01_extract_features.py
    python scripts/01_extract_features.py --recordings DR1_60HR_5CO --tables-dir results/check
"""

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (CHANNELS, FS, TABLES_DIR, BEAT_LEVEL_DIR,          # noqa: E402
                    iter_recordings, recording_id)
from ppgphantom.pyppg_features import extract_recording_channel, load_recording   # noqa: E402
from ppgphantom.sqi import compute_sqi, SQI_LABELS                                 # noqa: E402
from ppgphantom.templates import TemplateCollector                                 # noqa: E402

META = ["skin", "replicate", "hr", "co_lmin", "channel"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--recordings", nargs="*", help="only these recording IDs, e.g. DR1_60HR_5CO")
    ap.add_argument("--tables-dir", default=TABLES_DIR, help="where to write recording-level tables")
    ap.add_argument("--beat-dir", default=BEAT_LEVEL_DIR, help="where to write beat-level tables")
    args = ap.parse_args()

    morph_frames, sqi_frames, rec_sqi_rows = [], [], []
    templates = TemplateCollector()
    t0 = time.time()

    for skin, rep, hr, flow in iter_recordings():
        rid = recording_id(skin, rep, hr, flow)
        if args.recordings and rid not in args.recordings:
            continue
        df = load_recording(skin, rep, hr, flow)
        if df is None:
            print(f"  [missing] {rid}")
            continue

        print(f"\n  {rid}")
        for channel in CHANNELS:
            name = f"{rid}_{channel.replace(' ', '_')}"
            try:
                morph, segs = extract_recording_channel(df[channel].to_numpy(), name, fs=FS)
            except Exception as e:
                print(f"    {channel:20s} [FAIL] {type(e).__name__}: {e}")
                continue
            if morph is None:
                print(f"    {channel:20s} no pulses")
                continue

            meta = {"skin": skin, "replicate": rep, "hr": hr, "co_lmin": flow, "channel": channel}
            for i, col in enumerate(META):
                morph.insert(i, col, meta[col])
            morph.insert(len(META), "beat", range(len(morph)))
            morph_frames.append(morph)

            sqi_mat = compute_sqi(segs, fs=FS)[0]
            sqi_beats = pd.DataFrame(sqi_mat, columns=SQI_LABELS)
            for i, col in enumerate(META):
                sqi_beats.insert(i, col, meta[col])
            sqi_frames.append(sqi_beats)

            rec_sqi_rows.append({**meta, "n_beats": sqi_mat.shape[0],
                                 **dict(zip(SQI_LABELS, np.nanmedian(sqi_mat, axis=0)))})
            templates.add(skin, hr, flow, channel, [s.ac_norm for s in segs])
            print(f"    {channel:20s} {len(segs):4d} pulses")

    if not morph_frames:
        sys.exit("No recordings processed -- check that the data folder is set (see README).")

    os.makedirs(args.tables_dir, exist_ok=True)
    os.makedirs(args.beat_dir, exist_ok=True)

    beat_morph = pd.concat(morph_frames, ignore_index=True)
    beat_sqi = pd.concat(sqi_frames, ignore_index=True)
    beat_morph.to_csv(os.path.join(args.beat_dir, "beat_level_morphology_table.csv"), index=False)
    beat_sqi.to_csv(os.path.join(args.beat_dir, "beat_level_sqi_table.csv"), index=False)

    rec_morph = (beat_morph.groupby(META, sort=False).median(numeric_only=True)
                 .drop(columns=["beat"]).reset_index())
    rec_sqi = pd.DataFrame(rec_sqi_rows)
    rec_morph.to_csv(os.path.join(args.tables_dir, "recording_morphology_table.csv"), index=False)
    rec_sqi.to_csv(os.path.join(args.tables_dir, "recording_sqi_table.csv"), index=False)

    counts = (rec_sqi.groupby(["skin", "channel"], sort=False)["n_beats"].sum()
              .rename("total_pulses").reset_index())
    counts.to_csv(os.path.join(args.tables_dir, "pulse_counts.csv"), index=False)
    templates.to_table().to_csv(os.path.join(args.tables_dir, "pulse_templates.csv"), index=False)

    print(f"\nDone in {(time.time() - t0) / 60:.1f} min: {len(rec_sqi)} recording x channel rows, "
          f"{len(beat_sqi)} pulses.")
    print(f"Retained pulses per skin tone x channel: {counts['total_pulses'].min()}"
          f" to {counts['total_pulses'].max()}")


if __name__ == "__main__":
    main()
