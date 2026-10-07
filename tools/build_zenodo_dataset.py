"""
build_zenodo_dataset.py -- assemble the Zenodo data deposit from the local
working folders. Files are copied byte-for-byte (only the extension changes
from .txt to .csv), so the deposited data are exactly what the analysis read.

Output layout:
    <out>/README.txt
    <out>/PPG/<id>.csv            81 offset-corrected PPG recordings
    <out>/Pressure/<id>.csv       81 pressure recordings (decimated to 200 Hz)
    <out>/metadata/recordings.csv one row per recording
    <out>/metadata/SHA256SUMS.txt checksums of every data file
    <out>/derived/                recording-level SQI / morphology tables (optional)

Usage:
    python tools/build_zenodo_dataset.py --source "<Data Collection Physiology>/Dataset" --out D:/zenodo_ppg
    python tools/build_zenodo_dataset.py ... --zip      # also write one .zip per folder
"""

import argparse
import hashlib
import os
import shutil
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import FS, SKIN_LABELS, TABLES_DIR, iter_recordings, recording_id   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def _count_rows(path: str) -> int:
    with open(path, "rb") as f:
        return sum(1 for _ in f) - 1   # minus header


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True,
                    help="local Dataset/ folder containing PPG/PPG_Decoded and Pressure/Pressure_Decoded")
    ap.add_argument("--out", required=True, help="output folder for the deposit (must not exist or be empty)")
    ap.add_argument("--no-derived", action="store_true", help="do not include the recording-level tables")
    ap.add_argument("--zip", action="store_true", help="also zip each data folder for upload")
    args = ap.parse_args()

    src = {
        "PPG": os.path.join(args.source, "PPG", "PPG_Decoded"),
        "Pressure": os.path.join(args.source, "Pressure", "Pressure_Decoded"),
    }
    for kind, d in src.items():
        if not os.path.isdir(d):
            sys.exit(f"{kind} source folder not found: {d}")
    if os.path.isdir(args.out) and os.listdir(args.out):
        sys.exit(f"Output folder is not empty: {args.out}")

    for sub in ("PPG", "Pressure", "metadata"):
        os.makedirs(os.path.join(args.out, sub), exist_ok=True)

    rows, copied = [], []
    for skin, rep, hr, flow in iter_recordings():
        rid = recording_id(skin, rep, hr, flow)
        row = {"recording_id": rid, "skin_tone": SKIN_LABELS[skin], "replicate": rep,
               "heart_rate_bpm": hr, "target_flow_l_min": flow}
        for kind, d in src.items():
            s_path = os.path.join(d, rid + ".txt")
            present = os.path.exists(s_path)
            row[f"{kind.lower()}_file"] = f"{kind}/{rid}.csv" if present else ""
            if present:
                d_path = os.path.join(args.out, kind, rid + ".csv")
                shutil.copyfile(s_path, d_path)
                copied.append(d_path)
                n = _count_rows(d_path)
                row[f"{kind.lower()}_samples"] = n
                row[f"{kind.lower()}_duration_s"] = round(n / FS, 2)
        rows.append(row)
        print(f"  {rid}  PPG={'y' if row['ppg_file'] else '-'}  "
              f"Pressure={'y' if row['pressure_file'] else '-'}")

    meta = pd.DataFrame(rows)
    meta.to_csv(os.path.join(args.out, "metadata", "recordings.csv"), index=False)

    if not args.no_derived:
        os.makedirs(os.path.join(args.out, "derived"), exist_ok=True)
        for name in ("recording_sqi_table.csv", "recording_morphology_table.csv"):
            s_path = os.path.join(TABLES_DIR, name)
            if os.path.exists(s_path):
                d_path = os.path.join(args.out, "derived", name)
                shutil.copyfile(s_path, d_path)
                copied.append(d_path)
            else:
                print(f"  [warn] {name} not found in {TABLES_DIR} -- skipped")

    shutil.copyfile(os.path.join(HERE, "zenodo_README.txt"), os.path.join(args.out, "README.txt"))

    print("\nComputing checksums ...")
    with open(os.path.join(args.out, "metadata", "SHA256SUMS.txt"), "w", newline="\n") as f:
        for p in copied:
            f.write(f"{_sha256(p)}  {os.path.relpath(p, args.out).replace(os.sep, '/')}\n")

    if args.zip:
        for sub in ("PPG", "Pressure"):
            shutil.make_archive(os.path.join(args.out, sub), "zip", args.out, sub)
            print(f"  zipped {sub}.zip")

    print(f"\nDone: {meta['ppg_file'].astype(bool).sum()} PPG, "
          f"{meta['pressure_file'].astype(bool).sum()} pressure files -> {args.out}")


if __name__ == "__main__":
    main()
