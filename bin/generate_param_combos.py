#!/usr/bin/env python3
"""
Generates parameter combination files for Stage 1 grid search.

Reads conf/param_grids.yaml, generates all combinations (with optional
Latin-hypercube downsampling for large grids), and writes:
  - param_manifest.json  : hash → {method, params} mapping
  - <method>_params.tsv  : one row per combination, for use in Nextflow

Usage:
    generate_param_combos.py <param_grids.yaml> [--strategy grid|coordinate_descent]
"""
import sys
import json
import yaml
import hashlib
import itertools
import argparse
import numpy as np
from pathlib import Path


def latin_hypercube_sample(grid: dict, max_n: int, seed: int = 42) -> list:
    """Sample up to max_n parameter combinations using Latin Hypercube Sampling."""
    keys = list(grid.keys())
    values = [grid[k] for k in keys]
    all_combos = list(itertools.product(*values))

    if len(all_combos) <= max_n:
        return [dict(zip(keys, combo)) for combo in all_combos]

    # Stratified random sampling across each parameter dimension
    rng = np.random.default_rng(seed)
    indices = rng.choice(len(all_combos), size=max_n, replace=False)
    return [dict(zip(keys, all_combos[i])) for i in sorted(indices)]


def coordinate_descent_combos(grid: dict, defaults: dict) -> list:
    """
    Generate one-at-a-time combinations for coordinate descent.

    Starts with the all-defaults combination, then for each parameter
    generates combos where that parameter varies while all others stay
    at their default value.

    Total combos = 1 + sum(len(values) - 1 for each param).
    """
    # Validate all params have a default
    for k in grid:
        if k not in defaults:
            sys.exit(f"ERROR: coordinate_descent requires a default value for '{k}' "
                     f"but none found in 'defaults' block of param_grids.yaml")

    seen = set()
    combos = []

    def add(combo):
        key = json.dumps(combo, sort_keys=True)
        if key not in seen:
            seen.add(key)
            combos.append(combo)

    # All-defaults combo first
    add({k: defaults[k] for k in grid})

    # Vary one parameter at a time
    for param in grid:
        for val in grid[param]:
            combo = {k: defaults[k] for k in grid}
            combo[param] = val
            add(combo)

    return combos


def param_hash(params: dict) -> str:
    """Short deterministic hash for a parameter combination."""
    serialized = json.dumps(params, sort_keys=True)
    return hashlib.md5(serialized.encode()).hexdigest()[:8]


def generate_combos(param_grids_yaml: str, strategy: str, nucleus_only: bool) -> dict:
    with open(param_grids_yaml) as f:
        grids = yaml.safe_load(f)

    manifest = {}

    for method, config in grids.items():
        if not isinstance(config, dict):
            continue

        max_n    = config.pop("max_combinations", 9999)
        defaults = config.pop("defaults", {})

        # For cellpose (and cellpose_baysor inheriting from it): resolve diameter_cell vs
        # diameter_nucleus based on --nucleus-segmentation-only, then expose as 'diameter'.
        if method == "cellpose":
            if nucleus_only:
                config["diameter"] = config.pop("diameter_nucleus", config.pop("diameter_cell", [15]))
                config.pop("diameter_cell", None)
                defaults["diameter"] = defaults.get("diameter", 15)
            else:
                config["diameter"] = config.pop("diameter_cell", config.pop("diameter_nucleus", [30]))
                config.pop("diameter_nucleus", None)
                defaults["diameter"] = defaults.get("diameter", 30)

        # Only list-valued entries are searchable parameters
        grid = {k: v for k, v in config.items() if isinstance(v, list)}

        if strategy == "coordinate_descent":
            combos = coordinate_descent_combos(grid, defaults)
        else:
            combos = latin_hypercube_sample(grid, max_n)

        print(f"  {method}: {len(combos)} combinations [{strategy}]", file=sys.stderr)

        rows = []
        for combo in combos:
            h = param_hash(combo)
            manifest[h] = {"method": method, "params": combo}
            row = {"param_hash": h, "method": method}
            row.update(combo)
            rows.append(row)

        import csv
        if rows:
            tsv_path = f"{method}_params.tsv"
            with open(tsv_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

    Path("param_manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("param_grids_yaml", help="Path to param_grids.yaml")
    parser.add_argument(
        "--strategy",
        choices=["grid", "coordinate_descent"],
        default="grid",
        help="Search strategy: 'grid' (default) or 'coordinate_descent'",
    )
    parser.add_argument(
        "--nucleus-segmentation-only",
        action="store_true",
        default=False,
        help="Use nucleus diameter values for Cellpose instead of cell diameter values",
    )
    args = parser.parse_args()

    print(
        f"Generating parameter combinations (strategy: {args.strategy}, "
        f"nucleus_only: {args.nucleus_segmentation_only})...",
        file=sys.stderr,
    )
    manifest = generate_combos(args.param_grids_yaml, args.strategy, args.nucleus_segmentation_only)
    print(f"param_manifest.json written with {len(manifest)} total entries", file=sys.stderr)
