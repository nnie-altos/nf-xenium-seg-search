#!/usr/bin/env python3
"""
Scores a full-scale segmentation result for Stage 2.

Metrics:
  - MECR (Mutually Exclusive Co-expression Rate): lower = better
  - Transcript recovery rate: assigned / total
  - Cell yield: cell count (normalized by XOA3 baseline)

Composite score = 0.5*(1-MECR_norm) + 0.3*recovery + 0.2*yield_norm

Usage:
    score_full.py --h5ad <spatial_with_annotations.h5ad>
                  --transcripts <transcripts.parquet>
                  --cells <cell_boundaries.parquet or cells.parquet>
                  --markers <markers.yaml>
                  --method <method_name>
                  --sample <sample_id>
                  [--baseline-cell-count <n>]
                  [--mecr-weight 0.5]
                  [--recovery-weight 0.3]
                  [--yield-weight 0.2]

Output:
    <sample_id>_<method>_score.csv
"""
import sys
import argparse
import json
import yaml
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from scipy.sparse import issparse


def read_h5ad_safe(path: str):
    """Read h5ad, handling older files that stored uns values as null encoding.

    Older anndata versions wrote uns/log1p/base = None with encoding_type='null',
    which anndata 0.10.x cannot read back. Strip the offending key and retry.
    """
    try:
        return sc.read_h5ad(path)
    except Exception as e:
        if "null" not in str(e):
            raise
        import tempfile, shutil, h5py, os
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
            tmp_path = tmp.name
        shutil.copy2(path, tmp_path)
        try:
            with h5py.File(tmp_path, "r+") as f:
                if "uns/log1p" in f:
                    del f["uns/log1p"]
            return sc.read_h5ad(tmp_path)
        finally:
            os.unlink(tmp_path)


def normalize(s: str) -> str:
    return str(s).lower().replace("-", "").replace("_", "").replace(" ", "")


def load_markers(markers_yaml: str) -> dict:
    with open(markers_yaml) as f:
        data = yaml.safe_load(f)
    return data.get("markers", {})


def get_mecr(adata, marker_dict: dict) -> float:
    """
    MECR: Mutually Exclusive Co-expression Rate.
    Lower score = better segmentation (less transcript bleeding).
    Source: 11.MECR_lung.ipynb get_mecr() function.
    """
    # Fuzzy-match marker genes against adata.var_names
    actual_genes_map = {normalize(g): g for g in adata.var_names}
    available_genes = []
    for gene in marker_dict.keys():
        match = actual_genes_map.get(normalize(gene))
        if match:
            available_genes.append(match)
        else:
            print(f"  WARNING: marker gene '{gene}' not found in adata.var_names (skipped)", file=sys.stderr)

    if len(available_genes) < 4:
        print(f"  WARNING: Only {len(available_genes)} markers found. MECR requires >= 4. Returning 0.", file=sys.stderr)
        return 0.0

    # Map available genes to their cell types (using original marker_dict keys via fuzzy match)
    available_genes_map = {normalize(g): g for g in adata.var_names}
    gene_to_type = {}
    for orig_gene, ctype in marker_dict.items():
        match = available_genes_map.get(normalize(orig_gene))
        if match and match in available_genes:
            gene_to_type[match] = ctype

    # Extract binary expression matrix
    subset = adata[:, available_genes].X
    if issparse(subset):
        binary_mtx = (subset > 0).toarray()
    else:
        binary_mtx = (subset > 0)

    # Pairwise co-expression rate for cross-lineage gene pairs (Jaccard)
    coexp_rates = []
    for i, g1 in enumerate(available_genes):
        for j, g2 in enumerate(available_genes):
            if i >= j:
                continue
            type1 = gene_to_type.get(g1)
            type2 = gene_to_type.get(g2)
            if type1 is None or type2 is None or type1 == type2:
                continue
            vec1 = binary_mtx[:, i]
            vec2 = binary_mtx[:, j]
            intersection = np.logical_and(vec1, vec2).sum()
            union = np.logical_or(vec1, vec2).sum()
            if union > 0:
                coexp_rates.append(intersection / union)

    if not coexp_rates:
        return 0.0

    mecr = float(np.mean(coexp_rates))
    print(f"  MECR: {mecr:.4f} ({len(available_genes)} genes, {len(coexp_rates)} pairs)", file=sys.stderr)
    return mecr


def compute_assignment_rate(transcripts_parquet: str) -> tuple[float, int]:
    """Returns (assignment_rate, assigned_count)."""
    df = pd.read_parquet(transcripts_parquet)
    for col in ["cell_id", "cell", "assignment"]:
        if col in df.columns:
            assigned = df[col].notna() & (df[col] != "UNASSIGNED") & (df[col] != "") & (df[col] != 0)
            assigned_count = int(assigned.sum())
            return float(assigned_count) / max(len(df), 1), assigned_count
    print(f"  WARNING: No cell assignment column found in {transcripts_parquet}", file=sys.stderr)
    return 0.0, 0


def compute_cell_count(cells_path: str) -> int:
    df = pd.read_parquet(cells_path)
    # cell_boundaries.parquet has one row per vertex per cell — get unique cell IDs
    if "cell_id" in df.columns:
        return df["cell_id"].nunique()
    return len(df)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--transcripts", required=True)
    parser.add_argument("--cells", required=True)
    parser.add_argument("--markers", required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--baseline-cell-count", type=int, default=None)
    parser.add_argument("--mecr-weight", type=float, default=0.5)
    parser.add_argument("--recovery-weight", type=float, default=0.3)
    parser.add_argument("--yield-weight", type=float, default=0.2)
    args = parser.parse_args()

    print(f"Scoring {args.method} / {args.sample} ...", file=sys.stderr)

    # Load data
    adata = read_h5ad_safe(args.h5ad)
    marker_dict = load_markers(args.markers)

    # Compute metrics
    mecr = get_mecr(adata, marker_dict)
    assignment_rate, assigned_count = compute_assignment_rate(args.transcripts)
    cell_count = compute_cell_count(args.cells)
    mean_transcripts_per_cell = round(assigned_count / cell_count, 2) if cell_count > 0 else 0.0

    # Yield normalization: relative to baseline if provided, else 1.0
    if args.baseline_cell_count and args.baseline_cell_count > 0:
        yield_norm = min(cell_count / args.baseline_cell_count, 2.0) / 2.0
    else:
        yield_norm = 1.0  # No normalization possible without baseline

    # Composite score (MECR inverted: lower MECR = higher score)
    w_mecr = args.mecr_weight
    w_rec = args.recovery_weight
    w_yield = args.yield_weight

    # Normalize MECR: typical range 0–0.3; cap at 0.3 for normalization
    mecr_norm = min(mecr / 0.3, 1.0)
    mecr_score = 1.0 - mecr_norm

    composite = w_mecr * mecr_score + w_rec * assignment_rate + w_yield * yield_norm

    result = {
        "method": args.method,
        "sample": args.sample,
        "mecr": round(mecr, 4),
        "mecr_score": round(mecr_score, 4),
        "assignment_rate": round(assignment_rate, 4),
        "cell_count": cell_count,
        "mean_transcripts_per_cell": mean_transcripts_per_cell,
        "yield_norm": round(yield_norm, 4),
        "composite_score": round(composite, 4),
        "mecr_weight": w_mecr,
        "recovery_weight": w_rec,
        "yield_weight": w_yield,
    }

    out_path = f"{args.sample}_{args.method}_score.csv"
    pd.DataFrame([result]).to_csv(out_path, index=False)
    print(f"Score written to {out_path}", file=sys.stderr)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
