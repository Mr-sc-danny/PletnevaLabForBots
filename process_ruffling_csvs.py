"""
process_ruffling_csvs.py

Purpose:
    Batch-process every .csv file in a folder of MD "ruffling" analysis
    outputs (e.g. E:\\Desktop\\research\\MD\\slit\\ruffling_MD).

    Each input CSV has:
        - Column 1 ("File"): text frame identifier, e.g. "frame_000000_analysis"
        - Columns 2+: numeric heme-distortion metrics (Doop_exp, Doop_sim,
          Doming, Saddling, Ruffling, WavingX, WavingY, Propellering)

Note on this repo's reference code:
    Checked https://github.com/Mr-sc-danny/PletnevaLabForBots for an existing
    script covering this exact workflow (batch-folder CSV frame-range abs/avg/std).
    Nothing in that repo does this specific task, so there is no existing code
    to reference here -- this script is written from scratch.

Inputs:
    INPUT_DIR   - folder containing the per-run .csv files
    FRAME_START, FRAME_END - inclusive frame range used for the avg/std
                  (frame number is parsed out of the "File" text column,
                  e.g. "frame_001234_analysis" -> 1234)

Calculation:
    1. Read each .csv in INPUT_DIR.
    2. Take the absolute value of every numeric column (the "File" text
       column is left untouched).
    3. Parse the frame number out of the "File" column.
    4. Keep only rows with FRAME_START <= frame <= FRAME_END.
    5. For each numeric column, compute mean and std (sample std, ddof=1)
       over that frame window.

Output:
    One summary csv with:
        - Column 1: "File" = the source csv's filename
        - Remaining columns: one column per original numeric column
          (row 1 of the output header = row 1 of the input csv's header),
          with each cell formatted as "avg ± std", e.g. "3.93 ± 0.22".
"""

import glob
import os
import re
import sys

import numpy as np
import pandas as pd

# ---- Inputs -----------------------------------------------------------
INPUT_DIR = r"E:\Desktop\research\MD\slit\ruffling_MD"
FRAME_START = 1000
FRAME_END = 3000
OUTPUT_CSV = os.path.join(INPUT_DIR, "ruffling_avg_std_summary.csv")

FRAME_ID_COL = "File"          # first column, text frame identifier
FRAME_NUM_RE = re.compile(r"(\d+)")  # first run of digits in the identifier
DECIMALS = 2                   # decimal places shown in "avg ± std"

# Some source csvs use different punctuation/order for the Doop columns,
# e.g. "Doop (exp.)" / "Doop (sim.)" instead of "Doop_exp" / "Doop_sim".
# Normalize any variant to a single canonical name so columns line up
# across files regardless of source format.
_DOOP_EXP_RE = re.compile(r"(?i)doop.*exp")
_DOOP_SIM_RE = re.compile(r"(?i)doop.*sim")


def normalize_column_name(col: str) -> str:
    c = col.strip()
    if _DOOP_EXP_RE.search(c):
        return "Doop_exp"
    if _DOOP_SIM_RE.search(c):
        return "Doop_sim"
    return c


def extract_frame_number(file_id: str) -> int:
    """Pull the integer frame number out of a string like
    'frame_001234_analysis' -> 1234."""
    match = FRAME_NUM_RE.search(str(file_id))
    if not match:
        raise ValueError(f"Could not parse a frame number from: {file_id!r}")
    return int(match.group(1))


def format_avg_std(avg: float, std: float, decimals: int = DECIMALS) -> str:
    """Format a mean/std pair as 'avg ± std', e.g. '3.93 ± 0.22'."""
    if pd.isna(avg) or pd.isna(std):
        return "NaN"
    return f"{avg:.{decimals}f} \u00b1 {std:.{decimals}f}"


def process_one_csv(csv_path: str, frame_start: int, frame_end: int) -> dict:
    """Load one csv, abs() the numeric columns, filter to the frame window,
    and return a dict of {colname: 'avg ± std', ...}."""
    df = pd.read_csv(csv_path, sep=None, engine="python")  # auto-detect , or ;
    df.columns = [normalize_column_name(c) for c in df.columns]

    id_col = df.columns[0]           # first column = text identifier (e.g. "File")
    numeric_cols = df.columns[1:]    # remaining columns = numeric metrics

    # Make sure the numeric columns really are numeric, then take abs()
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors="coerce")
    df[numeric_cols] = df[numeric_cols].abs()

    # Parse frame number from the text identifier column
    frame_numbers = df[id_col].apply(extract_frame_number)

    # Keep only the requested frame window
    mask = (frame_numbers >= frame_start) & (frame_numbers <= frame_end)
    windowed = df.loc[mask, numeric_cols]

    if windowed.empty:
        print(f"  WARNING: no rows in frame range {frame_start}-{frame_end} "
              f"for {os.path.basename(csv_path)}; avg/std will be NaN.")

    means = windowed.mean(numeric_only=True)
    stds = windowed.std(numeric_only=True, ddof=1)

    stats = {}
    for col in numeric_cols:
        stats[col] = format_avg_std(means.get(col, np.nan), stds.get(col, np.nan))
    return stats


def main(input_dir: str, frame_start: int, frame_end: int, output_csv: str):
    csv_paths = sorted(glob.glob(os.path.join(input_dir, "*.csv")))
    # Don't re-ingest our own output if this script has already been run
    # once in this folder.
    csv_paths = [p for p in csv_paths if os.path.abspath(p) != os.path.abspath(output_csv)]
    if not csv_paths:
        print(f"No .csv files found in: {input_dir}")
        sys.exit(1)

    rows = []
    for csv_path in csv_paths:
        name = os.path.basename(csv_path)
        print(f"Processing: {name}")
        stats = process_one_csv(csv_path, frame_start, frame_end)
        stats[FRAME_ID_COL] = name  # column 1 = name of the source csv
        rows.append(stats)

    result = pd.DataFrame(rows)

    # Prefer a canonical column order when these known metric names are present;
    # any other/extra columns are appended afterward in their natural order.
    canonical_order = [
        "Doop_exp", "Doop_sim", "Doming", "Saddling",
        "Ruffling", "WavingX", "WavingY", "Propellering",
    ]
    other_cols = [c for c in canonical_order if c in result.columns]
    other_cols += [c for c in result.columns if c not in other_cols and c != FRAME_ID_COL]
    result = result[[FRAME_ID_COL] + other_cols]

    result.to_csv(output_csv, index=False, encoding="utf-8-sig")
    print(f"\nSaved summary for {len(rows)} csv file(s) to: {output_csv}")
    return result


if __name__ == "__main__":
    main(INPUT_DIR, FRAME_START, FRAME_END, OUTPUT_CSV)