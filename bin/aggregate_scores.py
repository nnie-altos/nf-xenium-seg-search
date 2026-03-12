#!/usr/bin/env python3
"""
Aggregates per-crop score JSONs and selects optimal parameters per method.

For each (method, param_hash): computes mean composite score across all crops.
Selects the param_hash with the highest mean score per method.
Outputs optimal_params.json and scores_summary.csv.

Usage:
    aggregate_scores.py <score_json_1> [<score_json_2> ...] --param-manifest <param_manifest.json>

param_manifest.json maps param_hash → {method: str, params: dict}
"""
import sys
import json
import argparse
import pandas as pd
from pathlib import Path
from collections import defaultdict


def load_scores(json_files: list) -> pd.DataFrame:
    records = []
    for f in json_files:
        with open(f) as fh:
            records.append(json.load(fh))
    return pd.DataFrame(records)


def aggregate(scores_df: pd.DataFrame, param_manifest: dict) -> dict:
    # Mean composite score per (method, param_hash) across all crops
    summary = (
        scores_df.groupby(["method", "param_hash"])
        .agg(
            mean_composite=("composite_score", "mean"),
            mean_assignment=("assignment_rate", "mean"),
            mean_cell_yield=("cell_yield_norm", "mean"),
            n_crops=("crop_id", "count"),
        )
        .reset_index()
        .sort_values(["method", "mean_composite"], ascending=[True, False])
    )

    summary.to_csv("scores_summary.csv", index=False)
    print(f"Scores summary written to scores_summary.csv", file=sys.stderr)
    print(summary.to_string(), file=sys.stderr)

    # Select optimal (highest mean composite) per method
    optimal = {}
    for method, group in summary.groupby("method"):
        best_row = group.iloc[0]  # already sorted descending
        best_hash = best_row["param_hash"]
        params = param_manifest.get(best_hash, {}).get("params", {})
        optimal[method] = {
            "param_hash": best_hash,
            "mean_composite_score": round(float(best_row["mean_composite"]), 4),
            "mean_assignment_rate": round(float(best_row["mean_assignment"]), 4),
            "mean_cell_yield_norm": round(float(best_row["mean_cell_yield"]), 4),
            "n_crops": int(best_row["n_crops"]),
            "params": params,
        }
        print(
            f"  {method}: best hash={best_hash}, score={best_row['mean_composite']:.4f}, "
            f"params={params}",
            file=sys.stderr,
        )

    return optimal


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("score_jsons", nargs="+")
    parser.add_argument("--param-manifest", required=True)
    args = parser.parse_args()

    with open(args.param_manifest) as f:
        param_manifest = json.load(f)

    scores_df = load_scores(args.score_jsons)
    optimal = aggregate(scores_df, param_manifest)

    Path("optimal_params.json").write_text(json.dumps(optimal, indent=2))
    print(f"\nOptimal params written to optimal_params.json", file=sys.stderr)
    print(json.dumps(optimal, indent=2))
