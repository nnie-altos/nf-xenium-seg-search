#!/usr/bin/env python3
"""
Pairwise IoU-based Average Precision scoring across segmentation methods.

No ground truth required: each method is compared against every other as
both test and reference, producing a symmetric AP matrix.

Algorithm (from Stringer et al. / Cellpose benchmarking, adapted for
cross-method comparison):
  For each ordered pair (A, B):
    1. geopandas spatial join to find intersecting cell pairs (A intersects B)
    2. Compute IoU = intersection.area / union.area for each candidate pair
    3. For each IoU threshold t in np.linspace(0, 1, 11):
         Sort candidates by IoU descending
         1-to-1 greedy matching: drop duplicates on left cell, then right cell
         TP = number of matched pairs with IoU >= t
         precision = TP / |A|
         recall    = TP / |B|
         score(t)  = precision × recall
    4. AP(A→B) = mean score over all thresholds

Output CSV columns:
    sample_id, method_a, method_b, n_cells_a, n_cells_b,
    ap, ap_at_0.5, ap_at_0.8

Usage:
    score_ap.py <sample_id> <method1>=<cells1.parquet> [<method2>=<cells2.parquet> ...]
"""
import gzip
import json
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Polygon

THRESHOLDS = np.round(np.linspace(0, 1, 11), 2)   # [0.0, 0.1, ..., 1.0]
IDX_50  = list(THRESHOLDS).index(0.5)
IDX_80  = list(THRESHOLDS).index(0.8)


def load_cells_as_gdf(cells_file: str) -> gpd.GeoDataFrame:
    """Load cell boundaries → GeoDataFrame with columns [cell_id, geometry].

    Handles:
      - Parquet with vertex_x / vertex_y columns (Xenium-style vertex list)
      - GeoJSON or GeoJSON.gz (any method that exports GeoJSON)
    """
    p = Path(cells_file)

    if p.suffix == ".parquet":
        df = pd.read_parquet(cells_file)
        if "vertex_x" not in df.columns or "vertex_y" not in df.columns:
            print(f"  WARNING: no vertex columns in {cells_file}, returning empty GDF",
                  file=sys.stderr)
            return gpd.GeoDataFrame(columns=["cell_id", "geometry"],
                                    geometry="geometry")
        polys = {}
        for cell_id, grp in df.groupby("cell_id"):
            coords = list(zip(grp["vertex_x"], grp["vertex_y"]))
            if len(coords) >= 3:
                polys[cell_id] = Polygon(coords)
        records = [{"cell_id": cid, "geometry": poly}
                   for cid, poly in polys.items()]
        if not records:
            return gpd.GeoDataFrame(columns=["cell_id", "geometry"],
                                    geometry="geometry")
        gdf = gpd.GeoDataFrame(records, geometry="geometry")
        gdf["geometry"] = gdf["geometry"].buffer(0)   # fix self-intersections
        return gdf.reset_index(drop=True)

    # GeoJSON / GeoJSON.gz
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(cells_file, "rt") as f:
        gj = json.load(f)
    records = []
    for feat in gj.get("features", []):
        cell_id = (feat.get("properties") or {}).get("cell_id") or feat.get("id")
        geom = feat.get("geometry") or {}
        if geom.get("type") == "Polygon" and geom.get("coordinates"):
            coords = geom["coordinates"][0]
            poly = Polygon([(c[0], c[1]) for c in coords])
            records.append({"cell_id": cell_id, "geometry": poly})
    if not records:
        return gpd.GeoDataFrame(columns=["cell_id", "geometry"],
                                geometry="geometry")
    gdf = gpd.GeoDataFrame(records, geometry="geometry")
    gdf["geometry"] = gdf["geometry"].buffer(0)
    return gdf.reset_index(drop=True)


def pairwise_ap(gdf_a: gpd.GeoDataFrame,
                gdf_b: gpd.GeoDataFrame) -> dict:
    """Compute AP of method A (test) vs method B (reference).

    Returns dict: ap, ap_at_0.5, ap_at_0.8, n_a, n_b.
    """
    n_a, n_b = len(gdf_a), len(gdf_b)
    empty = {"ap": 0.0, "ap_at_0.5": 0.0, "ap_at_0.8": 0.0,
             "n_a": n_a, "n_b": n_b}
    if n_a == 0 or n_b == 0:
        return empty

    # Spatial join to find candidate intersecting pairs
    candidates = gpd.sjoin(gdf_a, gdf_b, how="inner", predicate="intersects")
    if len(candidates) == 0:
        return empty

    # Compute IoU for every candidate pair
    ious = []
    for idx_a, row in candidates.iterrows():
        poly_a = gdf_a.at[idx_a, "geometry"]
        poly_b = gdf_b.at[row["index_right"], "geometry"]
        try:
            inter = poly_a.intersection(poly_b).area
            union = poly_a.union(poly_b).area
            ious.append(inter / union if union > 0 else 0.0)
        except Exception:
            ious.append(0.0)

    candidates = candidates.copy()
    candidates["iou"] = ious

    # Per-threshold scoring
    scores = []
    for t in THRESHOLDS:
        matches = (candidates[candidates["iou"] >= t]
                   .sort_values("iou", ascending=False))
        # 1-to-1 greedy matching
        tp = (matches
              .drop_duplicates("cell_id_left")
              .drop_duplicates("index_right")
              .shape[0])
        precision = tp / n_a
        recall    = tp / n_b
        scores.append(precision * recall)

    ap    = float(np.mean(scores))
    ap_50 = float(scores[IDX_50])
    ap_80 = float(scores[IDX_80])

    return {
        "ap":        round(ap, 4),
        "ap_at_0.5": round(ap_50, 4),
        "ap_at_0.8": round(ap_80, 4),
        "n_a": n_a,
        "n_b": n_b,
    }


def main():
    if len(sys.argv) < 3:
        sys.exit(
            f"Usage: {sys.argv[0]} <sample_id> <method>=<cells_file> ..."
        )

    sample_id = sys.argv[1]
    method_cells: dict[str, str] = {}
    for arg in sys.argv[2:]:
        if "=" not in arg:
            sys.exit(f"ERROR: expected method=cells_file, got: {arg}")
        method, cells_file = arg.split("=", 1)
        method_cells[method] = cells_file

    if len(method_cells) < 2:
        # Nothing to compare — write an empty matrix and exit cleanly
        pd.DataFrame(columns=[
            "sample_id", "method_a", "method_b",
            "n_cells_a", "n_cells_b", "ap", "ap_at_0.5", "ap_at_0.8"
        ]).to_csv(f"{sample_id}_ap_matrix.csv", index=False)
        print(f"Only one method provided; empty AP matrix written.", file=sys.stderr)
        return

    # Load all GDFs
    gdfs: dict[str, gpd.GeoDataFrame] = {}
    for method, cells_file in method_cells.items():
        print(f"Loading {method}: {cells_file}", file=sys.stderr)
        gdfs[method] = load_cells_as_gdf(cells_file)
        print(f"  {len(gdfs[method])} cells", file=sys.stderr)

    methods = sorted(gdfs.keys())
    rows = []
    for method_a in methods:
        for method_b in methods:
            if method_a == method_b:
                continue
            print(f"AP: {method_a} → {method_b}", file=sys.stderr)
            result = pairwise_ap(gdfs[method_a], gdfs[method_b])
            rows.append({
                "sample_id": sample_id,
                "method_a":  method_a,
                "method_b":  method_b,
                "n_cells_a": result["n_a"],
                "n_cells_b": result["n_b"],
                "ap":        result["ap"],
                "ap_at_0.5": result["ap_at_0.5"],
                "ap_at_0.8": result["ap_at_0.8"],
            })

    df = pd.DataFrame(rows)
    out_file = f"{sample_id}_ap_matrix.csv"
    df.to_csv(out_file, index=False)
    print(f"Written: {out_file}", file=sys.stderr)
    print(df.to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
