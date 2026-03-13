#!/usr/bin/env python3
"""
Renders a crop overlay: DAPI morphology image with cell boundary outlines
coloured by transcript count per cell.

Usage:
    render_crop_overlay.py <morphology.tif> <cells_file> <transcripts.parquet>
                           <method> <param_hash> <crop_id> <sample_id>

Output:
    <sample_id>_<crop_id>_<method>_<param_hash>_overlay.png
"""
import sys
import json
import gzip
import numpy as np
import pandas as pd
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection
from pathlib import Path


def load_cell_boundaries(cells_file: str) -> dict:
    """
    Load cell boundary polygons.
    Returns dict: cell_id → list of (x, y) coordinate tuples.
    Handles Xenium parquet (vertex_x/vertex_y) and GeoJSON/GeoJSON.gz formats.
    """
    p = Path(cells_file)
    if p.suffix == ".parquet":
        df = pd.read_parquet(cells_file)
        if "cell_id" in df.columns and "vertex_x" in df.columns:
            cells = {}
            for cell_id, grp in df.groupby("cell_id"):
                cells[cell_id] = list(zip(grp["vertex_x"], grp["vertex_y"]))
            return cells
        return {}

    if ".geojson" in p.name or p.suffix in (".gz", ".json"):
        opener = gzip.open if p.suffix == ".gz" else open
        with opener(cells_file, "rt") as f:
            gj = json.load(f)
        cells = {}
        for feat in gj.get("features", []):
            cell_id = (feat.get("properties") or {}).get("cell_id") or feat.get("id")
            geom = feat.get("geometry") or {}
            if geom.get("type") == "Polygon":
                coords = geom["coordinates"][0]
                cells[cell_id] = [(c[0], c[1]) for c in coords]
        return cells

    print(f"WARNING: unrecognised cell file format: {cells_file}", file=sys.stderr)
    return {}


def compute_transcripts_per_cell(transcripts_parquet: str) -> dict:
    """Returns dict: cell_id → transcript_count for assigned transcripts."""
    df = pd.read_parquet(transcripts_parquet)
    for col in ["cell_id", "cell", "assignment"]:
        if col in df.columns:
            assigned = df[
                df[col].notna()
                & (df[col] != "UNASSIGNED")
                & (df[col] != "")
                & (df[col] != 0)
            ]
            return assigned.groupby(col).size().to_dict()
    return {}


def render_overlay(
    tif_path: str,
    cells_file: str,
    transcripts_parquet: str,
    out_path: str,
    method: str,
    param_hash: str,
    crop_id: str,
) -> None:
    # ── Load and normalise DAPI image ────────────────────────────────────────
    img = tifffile.imread(tif_path)
    if img.ndim == 3:
        img = img[0]  # first channel = DAPI
    p2, p98 = np.percentile(img, (2, 98))
    img_norm = np.clip((img.astype(float) - p2) / max(p98 - p2, 1e-6), 0, 1)

    # ── Load cell boundaries and transcript counts ────────────────────────────
    cells = load_cell_boundaries(cells_file)
    tx_counts = compute_transcripts_per_cell(transcripts_parquet)

    cell_ids = list(cells.keys())
    counts = np.array([tx_counts.get(cid, 0) for cid in cell_ids], dtype=float)

    # ── Render ────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(5, 5), facecolor="black")
    ax.set_facecolor("black")
    ax.imshow(img_norm, cmap="gray", origin="upper", interpolation="nearest")

    if cells:
        vmin = float(np.percentile(counts, 5)) if len(counts) > 1 else 0.0
        vmax = float(np.percentile(counts, 95)) if len(counts) > 1 else max(counts.max(), 1.0)
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.cm.plasma

        patches, patch_counts = [], []
        for cid, count in zip(cell_ids, counts):
            coords = cells[cid]
            if len(coords) >= 3:
                patches.append(MplPolygon(coords, closed=True))
                patch_counts.append(count)

        if patches:
            pc = PatchCollection(
                patches, cmap=cmap, norm=norm,
                facecolor="none", linewidths=0.4, alpha=0.9,
            )
            pc.set_array(np.array(patch_counts))
            ax.add_collection(pc)

            cbar = fig.colorbar(pc, ax=ax, fraction=0.035, pad=0.02)
            cbar.set_label("Transcripts / cell", color="white", fontsize=7)
            cbar.ax.yaxis.set_tick_params(color="white", labelsize=6)
            plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color="white")

    n_cells = len(patches) if cells else 0
    ax.set_title(
        f"{method}  |  {param_hash}  |  {crop_id}\n{n_cells} cells",
        color="white", fontsize=7, pad=4,
    )
    ax.axis("off")
    fig.tight_layout(pad=0.3)
    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="black")
    plt.close(fig)
    print(f"Written: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 8:
        sys.exit(
            f"Usage: {sys.argv[0]} <morphology.tif> <cells_file> <transcripts.parquet> "
            "<method> <param_hash> <crop_id> <sample_id>"
        )

    tif_path, cells_file, tx_parquet, method, param_hash, crop_id, sample_id = sys.argv[1:]
    out_path = f"{sample_id}_{crop_id}_{method}_{param_hash}_overlay.png"
    render_overlay(tif_path, cells_file, tx_parquet, out_path, method, param_hash, crop_id)
