#!/usr/bin/env python3
"""
Crops an OME-TIFF morphology image to a pixel bounding box.

Uses lazy (zarr) reading to avoid loading the full slide into memory.
Reads only the Level 0 (highest resolution) data.

Usage:
    crop_image.py <morphology.ome.tif> <px_x_min> <px_x_max> <px_y_min> <px_y_max> <crop_id>

Output:
    <crop_id>_morphology.tif  — cropped TIFF (all channels preserved)
"""
import sys
import numpy as np
import tifffile
from pathlib import Path


def crop_image(
    tif_path: str,
    px_x_min: int,
    px_x_max: int,
    px_y_min: int,
    px_y_max: int,
    crop_id: str,
) -> None:
    out_path = f"{crop_id}_morphology.tif"

    with tifffile.TiffFile(tif_path) as tif:
        # Use zarr store for lazy reading — avoids loading full slide into RAM
        store = tif.aszarr(level=0)

        import zarr
        z = zarr.open(store, mode="r")

        # Handle common dimension orders: (C, Y, X), (Y, X, C), or (Y, X)
        shape = z.shape
        ndim = len(shape)

        print(f"Input shape: {shape}, dtype: {z.dtype}", file=sys.stderr)

        # Clamp to valid bounds
        if ndim == 3:
            # Could be (C, Y, X) or (Y, X, C) — detect by comparing dims
            if shape[0] <= 8:  # C, Y, X — first dim is channels (small)
                c, h, w = shape
                px_y_max = min(px_y_max, h)
                px_x_max = min(px_x_max, w)
                crop = z[:, px_y_min:px_y_max, px_x_min:px_x_max]
            else:  # Y, X, C
                h, w, c = shape
                px_y_max = min(px_y_max, h)
                px_x_max = min(px_x_max, w)
                crop = z[px_y_min:px_y_max, px_x_min:px_x_max, :]
        elif ndim == 2:
            h, w = shape
            px_y_max = min(px_y_max, h)
            px_x_max = min(px_x_max, w)
            crop = z[px_y_min:px_y_max, px_x_min:px_x_max]
        else:
            sys.exit(f"ERROR: Unexpected image dimensions: {shape}")

        crop_arr = np.array(crop)

    tifffile.imwrite(out_path, crop_arr, compression="deflate")
    print(
        f"Cropped {tif_path} px[{px_x_min}:{px_x_max}, {px_y_min}:{px_y_max}] "
        f"→ {out_path} shape {crop_arr.shape}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    if len(sys.argv) != 7:
        sys.exit(
            f"Usage: {sys.argv[0]} <morphology.ome.tif> "
            "<px_x_min> <px_x_max> <px_y_min> <px_y_max> <crop_id>"
        )

    crop_image(
        tif_path=sys.argv[1],
        px_x_min=int(sys.argv[2]),
        px_x_max=int(sys.argv[3]),
        px_y_min=int(sys.argv[4]),
        px_y_max=int(sys.argv[5]),
        crop_id=sys.argv[6],
    )
