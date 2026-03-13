#!/usr/bin/env python3
"""
Export the winning segmentation method and its optimal parameters as a
nf-xenium-processing-compatible params YAML (-params-file).

Selects the winner by mean composite score across all samples from Stage 2,
excluding XOA baseline methods (xoa3, xoa4).

Parameter mapping to nf-xenium-processing nextflow.config:
  Cellpose:
    diameter        → cellpose_cell_diameter  (cell mode)
                    → cellpose_nuclei_diameter (nucleus mode)
    flow_threshold  → cellpose_flow_threshold
    sharpen_tiff    → sharpen_tiff
  Cellpose+Baysor:
    diameter        → cellpose_cell_diameter / cellpose_nuclei_diameter
    flow_threshold  → cellpose_flow_threshold
    prior_segmentation_confidence, min_molecules_per_cell → no direct param
      (use custom baysor_config TOML — values noted as YAML comments)
  ProSeg:
    All params use the xenium preset in nf-xenium-processing — not exposed
      (optimal values noted as YAML comments for reference)
  SEGGER:
    min_transcripts_per_cell → min_transcripts_per_cell
    dist_tx, tile_size → no direct param (noted as YAML comments)
  XeniumRanger:
    expansion_distance, dapi_filter, boundary_stain → no direct param
      (noted as YAML comments)

Usage:
    export_params.py \\
        --scores-dir    <dir with *_score.csv files> \\
        --optimal-params optimal_params.json \\
        --nucleus-segmentation-only true|false
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import yaml


# Methods produced by the XOA baseline — not real "resegmentation" options
BASELINE_METHODS = {"xoa3", "xoa4"}

# Canonical nf-xenium-processing method names
NF_METHOD_MAP = {
    "proseg":         "proseg",
    "cellpose":       "cellpose",
    "cellpose_baysor": "cellpose_baysor",
    "segger":         "segger",
    "xr":             "xr",
}


def select_winner(scores_dir: Path) -> tuple[str, float]:
    """Return (method, mean_composite) for the best non-baseline method."""
    files = list(scores_dir.glob("*_score.csv"))
    if not files:
        sys.exit(f"ERROR: no *_score.csv files in {scores_dir}")
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df = df[~df["method"].isin(BASELINE_METHODS)]
    if df.empty:
        sys.exit("ERROR: no non-baseline methods found in score files")
    agg = df.groupby("method")["composite_score"].mean().sort_values(ascending=False)
    winner = agg.index[0]
    return winner, float(agg.iloc[0])


def build_params(method: str, opt_params: dict, nucleus_segmentation_only: bool) -> tuple[dict, list[str]]:
    """Return (yaml_dict, comment_lines) for the winning method.

    yaml_dict   — params that have direct equivalents in nf-xenium-processing
    comment_lines — human-readable notes for params without direct equivalents
    """
    nf_method = NF_METHOD_MAP.get(method, method)
    p = opt_params.get(method, {}).get("params", {})

    params: dict = {
        "segmentation": nf_method,
        "nucleus_segmentation_only": nucleus_segmentation_only,
    }
    comments: list[str] = []

    if method == "cellpose":
        if nucleus_segmentation_only:
            params["cellpose_nuclei_diameter"] = float(p.get("diameter", 15))
        else:
            params["cellpose_cell_diameter"] = float(p.get("diameter", 30))
        params["cellpose_flow_threshold"] = float(p.get("flow_threshold", 0.4))
        if str(p.get("sharpen_tiff", "false")).lower() == "true":
            params["sharpen_tiff"] = True

    elif method == "cellpose_baysor":
        if nucleus_segmentation_only:
            params["cellpose_nuclei_diameter"] = float(p.get("diameter", 15))
        else:
            params["cellpose_cell_diameter"] = float(p.get("diameter", 30))
        params["cellpose_flow_threshold"] = float(p.get("flow_threshold", 0.4))
        # Baysor-specific params have no direct nf-xenium-processing equivalents
        psc  = p.get("prior_segmentation_confidence", "N/A")
        mmpc = p.get("min_molecules_per_cell", "N/A")
        comments += [
            "# Baysor params (no direct nf-xenium-processing equivalents —",
            "# set these in your baysor_config TOML):",
            f"#   prior_segmentation_confidence: {psc}",
            f"#   min_molecules_per_cell:         {mmpc}",
        ]

    elif method == "proseg":
        comments += [
            "# ProSeg uses the --xenium preset in nf-xenium-processing.",
            "# Optimal params from grid search (for reference only —",
            "# not currently exposed as nf-xenium-processing params):",
            f"#   cell_compactness:                {p.get('cell_compactness', 'N/A')}",
            f"#   max_transcript_nucleus_distance: {p.get('max_transcript_nucleus_distance', 'N/A')}",
            f"#   voxel_size:                      {p.get('voxel_size', 'N/A')}",
            f"#   diffusion_probability:           {p.get('diffusion_probability', 'N/A')}",
        ]

    elif method == "segger":
        mttc = p.get("min_transcripts_per_cell")
        if mttc is not None:
            params["min_transcripts_per_cell"] = int(mttc)
        comments += [
            "# SEGGER params not in nf-xenium-processing (for reference):",
            f"#   dist_tx:   {p.get('dist_tx', 'N/A')}",
            f"#   tile_size: {p.get('tile_size', 'N/A')}",
        ]

    elif method == "xr":
        comments += [
            "# XeniumRanger resegment params not in nf-xenium-processing (for reference):",
            f"#   expansion_distance: {p.get('expansion_distance', 'N/A')}",
            f"#   dapi_filter:        {p.get('dapi_filter', 'N/A')}",
            f"#   boundary_stain:     {p.get('boundary_stain', 'N/A')}",
        ]

    return params, comments


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scores-dir",     required=True,
                        help="Directory containing *_score.csv files from Stage 2")
    parser.add_argument("--optimal-params", required=True,
                        help="optimal_params.json from Stage 1")
    parser.add_argument("--nucleus-segmentation-only", default="false",
                        help="true|false — mirrors the pipeline flag")
    args = parser.parse_args()

    nucleus_segmentation_only = args.nucleus_segmentation_only.lower() == "true"

    scores_dir = Path(args.scores_dir)
    winner, score = select_winner(scores_dir)
    print(f"Winner: {winner}  (mean composite = {score:.4f})", file=sys.stderr)

    with open(args.optimal_params) as f:
        opt_params = json.load(f)

    params, comments = build_params(winner, opt_params, nucleus_segmentation_only)

    # Build YAML string with header comment + optional method-specific comments
    header = [
        "# nf-xenium-processing params generated by nf-xenium-seg-search",
        f"# Winning segmentation method: {winner}  (mean composite score: {score:.4f})",
        "# Usage: nextflow run nf-xenium-processing/main.nf -params-file nfxp_params.yaml ...",
        "#",
        "# NOTE: set 'input' and 'outdir' before use.",
        "",
    ]

    yaml_str = "\n".join(header)
    yaml_str += yaml.dump(params, default_flow_style=False, sort_keys=False)
    if comments:
        yaml_str += "\n" + "\n".join(comments) + "\n"

    Path("nfxp_params.yaml").write_text(yaml_str)
    print("Written: nfxp_params.yaml", file=sys.stderr)

    # Also write a machine-readable summary
    summary = {
        "winner": winner,
        "mean_composite_score": round(score, 4),
        "nf_xenium_processing_params": params,
        "params_without_direct_equivalent": comments,
    }
    Path("winner_summary.json").write_text(json.dumps(summary, indent=2))
    print("Written: winner_summary.json", file=sys.stderr)


if __name__ == "__main__":
    main()
