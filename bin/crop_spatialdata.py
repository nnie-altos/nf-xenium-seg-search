#!/usr/bin/env python3
"""
Crop transcripts and morphology image for all crops using SpatialData.

Image crops use SpatialData's Image2DModel + bounding_box_query so the
pixel/micron coordinate transform is handled formally rather than assumed.
Transcript crops use a pandas bounding-box filter to preserve all columns.

Usage:
    crop_spatialdata.py <transcripts.parquet> <morphology.ome.tif> <crops.csv> <pixel_size_um>

Output (per crop in crops.csv):
    <crop_id>_transcripts.parquet
    <crop_id>_morphology.tif
"""
import sys
import json
import itertools
import numpy as np
import pandas as pd
import tifffile
import dask  # noqa: F401 (used via dask.delayed)
import dask.array as da
from math import ceil
from pathlib import Path

from spatialdata import SpatialData, bounding_box_query
from spatialdata.models import Image2DModel
from spatialdata.transformations import Scale


def load_image_as_sdata(morphology_path: str, pixel_size_um: float) -> SpatialData:
    """
    Load an OME-TIFF lazily and wrap it in a SpatialData object with a Scale
    transform so bounding_box_query operates in micron space.

    Uses ZarrTiffStore (zarr 2 MutableMapping) directly rather than zarr.open(),
    which is incompatible with zarr 3.x. Chunks are read on demand via dask.
    """
    tif = tifffile.TiffFile(morphology_path)
    v2_store = tif.aszarr(level=0)  # zarr 2 MutableMapping; zarr 3.x open() rejects it

    meta = json.loads(v2_store['.zarray'])
    shape = tuple(meta['shape'])       # e.g. (C, Y, X)
    chunk_shape = tuple(meta['chunks'])
    dtype = np.dtype(meta['dtype'])
    fill_value = meta.get('fill_value', 0)

    print(f"Morphology image shape: {shape}, dtype: {dtype}", file=sys.stderr)

    ndim = len(shape)
    chunk_grid = tuple(ceil(shape[i] / chunk_shape[i]) for i in range(ndim))

    def get_chunk(*chunk_idx):
        key = '.'.join(str(i) for i in chunk_idx)
        try:
            raw = v2_store[key]
        except KeyError:
            # chunk not stored → fill with fill_value (zarr convention)
            actual = tuple(min(chunk_shape[i], shape[i] - chunk_idx[i] * chunk_shape[i])
                           for i in range(ndim))
            return np.full(actual, fill_value, dtype=dtype)
        actual = tuple(min(chunk_shape[i], shape[i] - chunk_idx[i] * chunk_shape[i])
                       for i in range(ndim))
        return np.frombuffer(raw, dtype=dtype).reshape(actual)

    delayed_chunks = np.empty(chunk_grid, dtype=object)
    for idx in itertools.product(*[range(n) for n in chunk_grid]):
        actual_shape = tuple(min(chunk_shape[i], shape[i] - idx[i] * chunk_shape[i])
                             for i in range(ndim))
        delayed_chunks[idx] = da.from_delayed(
            dask.delayed(get_chunk)(*idx),
            shape=actual_shape,
            dtype=dtype,
        )

    img_dask = da.block(delayed_chunks.tolist())

    shape = img_dask.shape
    print(f"Morphology image shape: {shape}, dtype: {img_dask.dtype}", file=sys.stderr)

    # Normalise to (c, y, x) — SpatialData Image2DModel requires this dim order.
    if img_dask.ndim == 2:                         # (Y, X)
        img_dask = img_dask[None]                  # → (1, Y, X)
    elif img_dask.ndim == 3 and shape[0] > 8:      # (Y, X, C) — first dim is spatial
        img_dask = img_dask.transpose(2, 0, 1)     # → (C, Y, X)
    # else: (C, Y, X) already

    image = Image2DModel.parse(
        img_dask,
        dims=("c", "y", "x"),
        transformations={
            "global": Scale([pixel_size_um, pixel_size_um], axes=("y", "x"))
        },
    )
    return SpatialData(images={"morphology": image})


def crop_one(
    sdata: SpatialData,
    transcripts: pd.DataFrame,
    crop_id: str,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
) -> None:
    # ── Transcript crop (pandas bbox filter; preserves every column) ──────────
    mask = (
        (transcripts["x_location"] >= x_min) & (transcripts["x_location"] <= x_max) &
        (transcripts["y_location"] >= y_min) & (transcripts["y_location"] <= y_max)
    )
    tx_crop = transcripts[mask].reset_index(drop=True)
    tx_crop.to_parquet(f"{crop_id}_transcripts.parquet", index=False)

    # ── Image crop (SpatialData bounding_box_query; coordinate-transform-aware)
    cropped_sdata = bounding_box_query(
        sdata,
        axes=("x", "y"),
        min_coordinate=[x_min, y_min],
        max_coordinate=[x_max, y_max],
        target_coordinate_system="global",
    )
    if cropped_sdata is None or "morphology" not in cropped_sdata.images:
        sys.exit(f"ERROR: bounding_box_query returned empty image for crop {crop_id}")

    img_arr = np.array(cropped_sdata.images["morphology"].values)  # (C, Y, X)
    if img_arr.shape[0] == 1:
        img_arr = img_arr[0]   # (Y, X) for single-channel DAPI

    tifffile.imwrite(f"{crop_id}_morphology.tif", img_arr, compression="deflate")

    print(
        f"Crop {crop_id}: {len(tx_crop)} transcripts, image {img_arr.shape} "
        f"(x=[{x_min:.1f},{x_max:.1f}] y=[{y_min:.1f},{y_max:.1f}] µm)",
        file=sys.stderr,
    )


def main() -> None:
    if len(sys.argv) != 5:
        sys.exit(
            f"Usage: {sys.argv[0]} <transcripts.parquet> <morphology.ome.tif> "
            "<crops.csv> <pixel_size_um>"
        )

    tx_path, img_path, crops_csv_path = sys.argv[1], sys.argv[2], sys.argv[3]
    pixel_size_um = float(sys.argv[4])

    print(f"Loading morphology image: {img_path}", file=sys.stderr)
    sdata = load_image_as_sdata(img_path, pixel_size_um)

    print(f"Loading transcripts: {tx_path}", file=sys.stderr)
    transcripts = pd.read_parquet(tx_path)
    x_col = "x_location" if "x_location" in transcripts.columns else "x"
    y_col = "y_location" if "y_location" in transcripts.columns else "y"
    if x_col != "x_location":
        transcripts = transcripts.rename(columns={x_col: "x_location", y_col: "y_location"})

    crops = pd.read_csv(crops_csv_path)
    for _, row in crops.iterrows():
        crop_one(
            sdata=sdata,
            transcripts=transcripts,
            crop_id=row["crop_id"],
            x_min=float(row["x_min_um"]),
            x_max=float(row["x_max_um"]),
            y_min=float(row["y_min_um"]),
            y_max=float(row["y_max_um"]),
        )

    print(f"Done: {len(crops)} crops written.", file=sys.stderr)


if __name__ == "__main__":
    main()
