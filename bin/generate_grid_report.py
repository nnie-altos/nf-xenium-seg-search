#!/usr/bin/env python3
"""
Generates an HTML grid search report showing top-3 and bottom-3 parameter
combinations per method per crop, with DAPI+boundary overlays and metrics.

Usage:
    generate_grid_report.py <score_jsons...> --overlays <overlay_pngs...>

Score JSON filenames must follow the convention:
    <sample_id>_<crop_id>_<method>_<param_hash>_score.json

Overlay PNG filenames must follow the convention:
    <sample_id>_<crop_id>_<method>_<param_hash>_overlay.png

Output:
    grid_search_report.html
"""
import sys
import json
import base64
import argparse
from pathlib import Path
from collections import defaultdict


N_TOP = 3
N_BOTTOM = 3


def parse_score_filename(fname: str) -> dict | None:
    """Parse <sample_id>_<crop_id>_<method>_<param_hash>_score.json."""
    stem = Path(fname).stem  # strip .json
    if not stem.endswith("_score"):
        return None
    stem = stem[: -len("_score")]
    # Last token is param_hash (8 chars), second-last is method
    parts = stem.split("_")
    if len(parts) < 4:
        return None
    param_hash = parts[-1]
    method = parts[-2]
    crop_id = parts[-3]
    sample_id = "_".join(parts[:-3])
    return {"sample_id": sample_id, "crop_id": crop_id, "method": method, "param_hash": param_hash}


def load_scores(score_paths: list[str]) -> list[dict]:
    rows = []
    for p in score_paths:
        meta = parse_score_filename(p)
        if meta is None:
            print(f"WARNING: could not parse score filename: {p}", file=sys.stderr)
            continue
        with open(p) as f:
            data = json.load(f)
        data.update(meta)
        rows.append(data)
    return rows


def build_overlay_index(overlay_paths: list[str]) -> dict:
    """Returns dict: (sample_id, crop_id, method, param_hash) → png_path."""
    index = {}
    for p in overlay_paths:
        stem = Path(p).stem
        if not stem.endswith("_overlay"):
            continue
        stem = stem[: -len("_overlay")]
        parts = stem.split("_")
        if len(parts) < 4:
            continue
        param_hash = parts[-1]
        method = parts[-2]
        crop_id = parts[-3]
        sample_id = "_".join(parts[:-3])
        index[(sample_id, crop_id, method, param_hash)] = p
    return index


def img_to_b64(path: str) -> str:
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode()


def metric_row(score: dict) -> str:
    return (
        f"Score: <b>{score.get('composite_score', 'N/A')}</b> &nbsp;|&nbsp; "
        f"Assignment: {score.get('assignment_rate', 'N/A')} &nbsp;|&nbsp; "
        f"Cells: {score.get('cell_count', 'N/A')} &nbsp;|&nbsp; "
        f"Tx/cell: {score.get('mean_transcripts_per_cell', 'N/A')}"
    )


def build_html(scores: list[dict], overlay_index: dict) -> str:
    # Group by (sample_id, method, crop_id)
    groups = defaultdict(list)
    for s in scores:
        key = (s["sample_id"], s["method"], s["crop_id"])
        groups[key].append(s)

    sections = []
    for (sample_id, method, crop_id), combos in sorted(groups.items()):
        combos_sorted = sorted(combos, key=lambda x: x.get("composite_score", 0), reverse=True)
        top = combos_sorted[:N_TOP]
        bottom = combos_sorted[-N_BOTTOM:]
        # Avoid duplicates when total combos <= N_TOP + N_BOTTOM
        shown = top + [c for c in bottom if c not in top]

        cards = []
        for rank_label, group in [("TOP", top), ("BOTTOM", bottom)]:
            for i, combo in enumerate(group, 1):
                key = (combo["sample_id"], combo["crop_id"], combo["method"], combo["param_hash"])
                img_tag = ""
                if key in overlay_index:
                    b64 = img_to_b64(overlay_index[key])
                    img_tag = f'<img src="data:image/png;base64,{b64}" style="width:100%;border-radius:4px;">'
                else:
                    img_tag = '<div style="height:200px;background:#222;display:flex;align-items:center;justify-content:center;color:#666;">No overlay</div>'

                border_color = "#2ecc71" if rank_label == "TOP" else "#e74c3c"
                label_color = border_color
                cards.append(f"""
                <div style="border:2px solid {border_color};border-radius:6px;padding:8px;background:#1a1a1a;">
                    <div style="color:{label_color};font-size:11px;font-weight:bold;margin-bottom:4px;">
                        {rank_label} {i} &nbsp; <span style="color:#aaa;font-weight:normal;">{combo['param_hash']}</span>
                    </div>
                    {img_tag}
                    <div style="color:#ccc;font-size:10px;margin-top:5px;line-height:1.6;">
                        {metric_row(combo)}
                    </div>
                </div>""")

        cards_html = "\n".join(cards)
        sections.append(f"""
        <div style="margin-bottom:40px;">
            <h2 style="color:#fff;border-bottom:1px solid #444;padding-bottom:6px;">
                {method} &nbsp;<span style="color:#888;font-size:14px;">crop: {crop_id} &nbsp;|&nbsp; sample: {sample_id}</span>
            </h2>
            <div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(220px,1fr));gap:12px;">
                {cards_html}
            </div>
        </div>""")

    body = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Stage 1 Grid Search Report</title>
  <style>
    body {{ background:#111; color:#ddd; font-family:sans-serif; padding:30px; }}
    h1   {{ color:#fff; }}
    h2   {{ font-size:16px; }}
  </style>
</head>
<body>
  <h1>Stage 1 Grid Search Report</h1>
  <p style="color:#888;">Showing top-{N_TOP} and bottom-{N_BOTTOM} parameter combinations
     per method per crop, ranked by composite score.</p>
  {body}
</body>
</html>"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("score_jsons", nargs="+", help="Score JSON files")
    parser.add_argument("--overlays", nargs="*", default=[], help="Overlay PNG files")
    args = parser.parse_args()

    scores = load_scores(args.score_jsons)
    if not scores:
        sys.exit("ERROR: no valid score JSON files found")

    overlay_index = build_overlay_index(args.overlays)
    print(f"Loaded {len(scores)} scores, {len(overlay_index)} overlays", file=sys.stderr)

    html = build_html(scores, overlay_index)
    Path("grid_search_report.html").write_text(html)
    print("Written: grid_search_report.html", file=sys.stderr)
