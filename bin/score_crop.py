#!/usr/bin/env python3
"""
Scores a segmented crop by transcript assignment rate and cell yield.

Stage 1 scoring (no MECR — annotation not available on crops):
    composite = 0.6 × assignment_rate + 0.4 × normalized_cell_yield

Cell yield normalization is done relative to a reference cell density
(expected cells per mm² for the tissue type). Default: 1000 cells/mm².

Usage:
    score_crop.py <transcripts.parquet> <cell_polygons_geojson_or_parquet>
                  <method> <param_hash> <crop_id> <crop_area_mm2>

Output:
    score.json  — with keys: method, param_hash, crop_id, assignment_rate,
                  cell_count, cell_yield_norm, composite_score
"""
import sys
import json
import gzip
import pandas as pd
import numpy as np
from pathlib import Path


def load_cell_count(cell_file: str) -> int:
    """Load cell count from geojson.gz or parquet."""
    p = Path(cell_file)
    if p.suffix == ".gz" or "geojson" in p.name:
        import json as json_lib
        opener = gzip.open if p.suffix == ".gz" else open
        with opener(cell_file, "rt") as f:
            gj = json_lib.load(f)
        return len(gj.get("features", []))
    elif p.suffix == ".parquet":
        df = pd.read_parquet(cell_file)
        return len(df)
    else:
        sys.exit(f"ERROR: Unrecognized cell file format: {cell_file}")


def compute_assignment_rate(transcripts_parquet: str) -> tuple[float, int]:
    """Returns (assignment_rate, assigned_count)."""
    df = pd.read_parquet(transcripts_parquet)

    # Cell ID column may be 'cell_id' or 'overlaps_nucleus'
    if "cell_id" in df.columns:
        assigned = df["cell_id"].notna() & (df["cell_id"] != "UNASSIGNED") & (df["cell_id"] != "")
    elif "overlaps_nucleus" in df.columns:
        # fallback: use overlaps_nucleus as proxy
        assigned = df["overlaps_nucleus"] == 1
    else:
        # Last resort: any non-null, non-zero cell assignment
        for col in ["cell", "assignment", "cell_assignment"]:
            if col in df.columns:
                assigned = df[col].notna() & (df[col] != 0) & (df[col] != "UNASSIGNED")
                break
        else:
            sys.exit(f"ERROR: Cannot find cell assignment column in {transcripts_parquet}. "
                     f"Columns: {list(df.columns)}")

    total = len(df)
    assigned_count = int(assigned.sum())
    if total == 0:
        return 0.0, 0
    return float(assigned_count) / total, assigned_count


def score_crop(
    transcripts_parquet: str,
    cell_file: str,
    method: str,
    param_hash: str,
    crop_id: str,
    crop_area_mm2: float,
    ref_density: float = 1000.0,  # expected cells/mm²
) -> dict:
    assignment_rate, assigned_count = compute_assignment_rate(transcripts_parquet)
    cell_count = load_cell_count(cell_file)

    expected_cells = ref_density * crop_area_mm2
    cell_yield_norm = min(cell_count / expected_cells, 1.0) if expected_cells > 0 else 0.0

    mean_transcripts_per_cell = round(assigned_count / cell_count, 2) if cell_count > 0 else 0.0

    composite = 0.6 * assignment_rate + 0.4 * cell_yield_norm

    result = {
        "method": method,
        "param_hash": param_hash,
        "crop_id": crop_id,
        "crop_area_mm2": round(crop_area_mm2, 4),
        "assignment_rate": round(assignment_rate, 4),
        "cell_count": cell_count,
        "cell_yield_norm": round(cell_yield_norm, 4),
        "mean_transcripts_per_cell": mean_transcripts_per_cell,
        "composite_score": round(composite, 4),
    }
    print(json.dumps(result, indent=2), file=sys.stderr)
    return result


if __name__ == "__main__":
    if len(sys.argv) != 7:
        sys.exit(
            f"Usage: {sys.argv[0]} <transcripts.parquet> <cell_file> "
            "<method> <param_hash> <crop_id> <crop_area_mm2>"
        )

    result = score_crop(
        transcripts_parquet=sys.argv[1],
        cell_file=sys.argv[2],
        method=sys.argv[3],
        param_hash=sys.argv[4],
        crop_id=sys.argv[5],
        crop_area_mm2=float(sys.argv[6]),
    )

    Path("score.json").write_text(json.dumps(result, indent=2))
