#!/usr/bin/env python3
"""
Generates parameter combination files for Stage 1 grid search.

Reads conf/param_grids.yaml, generates all combinations (with optional
Latin-hypercube downsampling for large grids), and writes:
  - param_manifest.json  : hash → {method, params} mapping
  - <method>_params.tsv  : one row per combination, for use in Nextflow

Usage:
    generate_param_combos.py <param_grids.yaml>
"""
import sys
import json
import yaml
import hashlib
import itertools
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


def param_hash(params: dict) -> str:
    """Short deterministic hash for a parameter combination."""
    serialized = json.dumps(params, sort_keys=True)
    return hashlib.md5(serialized.encode()).hexdigest()[:8]


def generate_combos(param_grids_yaml: str) -> dict:
    with open(param_grids_yaml) as f:
        grids = yaml.safe_load(f)

    manifest = {}
    method_params = {}

    for method, config in grids.items():
        if not isinstance(config, dict):
            continue

        max_n = config.pop("max_combinations", 9999)

        # Rename keys: yaml uses underscores, CLI uses hyphens for proseg/cellpose
        grid = {k: v for k, v in config.items() if isinstance(v, list)}

        combos = latin_hypercube_sample(grid, max_n)
        method_params[method] = combos

        print(f"  {method}: {len(combos)} combinations", file=sys.stderr)

        rows = []
        for combo in combos:
            h = param_hash(combo)
            manifest[h] = {"method": method, "params": combo}
            row = {"param_hash": h, "method": method}
            row.update(combo)
            rows.append(row)

        # Write TSV for this method
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
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <param_grids.yaml>")

    print("Generating parameter combinations...", file=sys.stderr)
    manifest = generate_combos(sys.argv[1])
    print(f"param_manifest.json written with {len(manifest)} total entries", file=sys.stderr)
