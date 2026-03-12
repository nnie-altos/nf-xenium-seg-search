#!/usr/bin/env python3
"""
Reads experiment.xenium and prints the major XOA version number (3 or 4).
Usage: detect_xoa_version.py <experiment_xenium_path>
Output: writes xoa_version.txt with the integer major version
"""
import json
import sys
from pathlib import Path


def detect_version(experiment_xenium_path: str) -> int:
    path = Path(experiment_xenium_path)
    if not path.exists():
        sys.exit(f"ERROR: {path} does not exist")

    with open(path) as f:
        data = json.load(f)

    sw_version = data.get("analysis_sw_version", "")
    # e.g. "xenium-3.3.0.1" or "xenium-4.0.0"
    if not sw_version:
        sys.exit(f"ERROR: 'analysis_sw_version' not found in {path}")

    # Strip "xenium-" prefix and parse major version
    version_str = sw_version.replace("xenium-", "").strip()
    major = int(version_str.split(".")[0])

    print(f"Detected analysis_sw_version: {sw_version} → major version: {major}", file=sys.stderr)
    return major


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <experiment.xenium>")

    major = detect_version(sys.argv[1])
    Path("xoa_version.txt").write_text(str(major))
    print(major)
