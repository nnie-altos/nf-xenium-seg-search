#!/usr/bin/env python3
"""
Filters transcripts.parquet to a spatial bounding box (in microns).

Usage:
    crop_transcripts.py <transcripts.parquet> <x_min_um> <x_max_um> <y_min_um> <y_max_um> <crop_id>

Output:
    <crop_id>_transcripts.parquet
"""
import sys
import pandas as pd
from pathlib import Path


def crop_transcripts(
    parquet_path: str,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    crop_id: str,
) -> None:
    df = pd.read_parquet(parquet_path)

    # Xenium transcripts use x_location and y_location columns (microns)
    x_col = "x_location" if "x_location" in df.columns else "x"
    y_col = "y_location" if "y_location" in df.columns else "y"

    if x_col not in df.columns or y_col not in df.columns:
        sys.exit(f"ERROR: Could not find x/y coordinate columns in {parquet_path}. "
                 f"Available columns: {list(df.columns)}")

    mask = (
        (df[x_col] >= x_min) & (df[x_col] <= x_max) &
        (df[y_col] >= y_min) & (df[y_col] <= y_max)
    )
    cropped = df[mask].reset_index(drop=True)

    out_path = f"{crop_id}_transcripts.parquet"
    cropped.to_parquet(out_path, index=False)

    total = len(df)
    kept = len(cropped)
    print(
        f"Cropped {parquet_path}: {total} → {kept} transcripts "
        f"({100*kept/total:.1f}%) in bbox "
        f"x=[{x_min},{x_max}] y=[{y_min},{y_max}] µm",
        file=sys.stderr,
    )


if __name__ == "__main__":
    if len(sys.argv) != 7:
        sys.exit(
            f"Usage: {sys.argv[0]} <transcripts.parquet> "
            "<x_min_um> <x_max_um> <y_min_um> <y_max_um> <crop_id>"
        )

    crop_transcripts(
        parquet_path=sys.argv[1],
        x_min=float(sys.argv[2]),
        x_max=float(sys.argv[3]),
        y_min=float(sys.argv[4]),
        y_max=float(sys.argv[5]),
        crop_id=sys.argv[6],
    )
