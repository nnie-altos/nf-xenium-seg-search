#!/usr/bin/env python3
"""
Selects N density-stratified spatial crop regions from a Xenium dataset.

Reads nucleus_boundaries.parquet to estimate local cell density across the
tissue, then selects crops spanning low/mid/high density regions.

Usage:
    select_crops.py <nucleus_boundaries.parquet> <n_crops> <crop_size_um> <pixel_size_um>

Output:
    crops.csv  — columns: crop_id, x_min_um, x_max_um, y_min_um, y_max_um
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path


def compute_density_grid(nuclei: pd.DataFrame, tile_um: float = 200.0) -> pd.DataFrame:
    """Divide tissue into tiles and count nuclei per tile."""
    x = nuclei["vertex_x"].values
    y = nuclei["vertex_y"].values

    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()

    x_bins = np.arange(x_min, x_max + tile_um, tile_um)
    y_bins = np.arange(y_min, y_max + tile_um, tile_um)

    counts, xedges, yedges = np.histogram2d(x, y, bins=[x_bins, y_bins])

    tiles = []
    for i in range(len(xedges) - 1):
        for j in range(len(yedges) - 1):
            tiles.append({
                "tile_x_min": xedges[i],
                "tile_x_max": xedges[i + 1],
                "tile_y_min": yedges[j],
                "tile_y_max": yedges[j + 1],
                "count": counts[i, j],
            })

    return pd.DataFrame(tiles)


def select_crops(
    nuclei_parquet: str,
    n_crops: int,
    crop_size_um: float,
    pixel_size_um: float,
) -> pd.DataFrame:
    # Load nucleus boundary coordinates
    df = pd.read_parquet(nuclei_parquet, columns=["vertex_x", "vertex_y"])

    # Compute density grid (200µm tiles)
    tiles = compute_density_grid(df, tile_um=200.0)

    # Remove empty tiles (background / outside tissue)
    tiles = tiles[tiles["count"] > 0].reset_index(drop=True)

    if len(tiles) == 0:
        sys.exit("ERROR: No non-empty density tiles found in nucleus_boundaries.parquet")

    # Stratify into low/mid/high density bins
    tiles["density_bin"] = pd.qcut(tiles["count"], q=3, labels=["low", "mid", "high"])

    # Sample proportionally from each bin
    per_bin = max(1, n_crops // 3)
    remainder = n_crops - per_bin * 3

    sampled = []
    for i, (bin_label, group) in enumerate(tiles.groupby("density_bin", observed=True)):
        n = per_bin + (1 if i < remainder else 0)
        n = min(n, len(group))
        sampled.append(group.sample(n=n, random_state=42))

    selected_tiles = pd.concat(sampled).reset_index(drop=True)

    # Build crop bounding boxes centred on each selected tile
    crops = []
    half = crop_size_um / 2.0
    x_global_min = df["vertex_x"].min()
    x_global_max = df["vertex_x"].max()
    y_global_min = df["vertex_y"].min()
    y_global_max = df["vertex_y"].max()

    for i, row in selected_tiles.iterrows():
        cx = (row["tile_x_min"] + row["tile_x_max"]) / 2.0
        cy = (row["tile_y_min"] + row["tile_y_max"]) / 2.0

        x_min = max(cx - half, x_global_min)
        x_max = min(cx + half, x_global_max)
        y_min = max(cy - half, y_global_min)
        y_max = min(cy + half, y_global_max)

        # Convert micron bbox to pixel bbox
        px_x_min = int(x_min / pixel_size_um)
        px_x_max = int(x_max / pixel_size_um)
        px_y_min = int(y_min / pixel_size_um)
        px_y_max = int(y_max / pixel_size_um)

        crops.append({
            "crop_id": f"crop_{i:03d}",
            "density_bin": row["density_bin"],
            "nuclei_count": int(row["count"]),
            "x_min_um": round(x_min, 2),
            "x_max_um": round(x_max, 2),
            "y_min_um": round(y_min, 2),
            "y_max_um": round(y_max, 2),
            "px_x_min": px_x_min,
            "px_x_max": px_x_max,
            "px_y_min": px_y_min,
            "px_y_max": px_y_max,
        })

    return pd.DataFrame(crops)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(f"Usage: {sys.argv[0]} <nucleus_boundaries.parquet> <n_crops> <crop_size_um> <pixel_size_um>")

    nuclei_parquet = sys.argv[1]
    n_crops = int(sys.argv[2])
    crop_size_um = float(sys.argv[3])
    pixel_size_um = float(sys.argv[4])

    crops = select_crops(nuclei_parquet, n_crops, crop_size_um, pixel_size_um)
    crops.to_csv("crops.csv", index=False)
    print(f"Selected {len(crops)} crops → crops.csv", file=sys.stderr)
    print(crops.to_string(), file=sys.stderr)
