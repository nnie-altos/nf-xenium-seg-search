#!/usr/bin/env python3
"""
Generates the final HTML comparison report for nf-xenium-seg-search.

Inputs:
  --scores-dir       : directory containing <sample>_<method>_score.csv files
  --stage1-summary   : scores_summary.csv from Stage 1 aggregate_scores
  --stage1-manifest  : param_manifest.json from Stage 1
  --bundles          : JSON file mapping method → xenium_bundle_reseg outs/ path
  --morphology-dir   : path to original bundle morphology_focus/ directory
  --outdir           : output directory for report

Output:
    seg_search_report.html
"""
import argparse
import base64
import io
import json
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import geopandas as gpd
import tifffile


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

METHOD_COLORS = {
    "xoa3":     "#4e79a7",
    "xoa4":     "#59a14f",
    "cellpose": "#f28e2b",
    "proseg":   "#e15759",
    "segger":   "#76b7b2",
}

IMAGE_GROUPS = {"xoa3", "xoa4", "cellpose"}
TRANSCRIPT_GROUPS = {"proseg", "segger"}


def fig_to_base64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def method_color(m: str) -> str:
    return METHOD_COLORS.get(m.lower(), "#aaaaaa")


def group_label(m: str) -> str:
    m = m.lower()
    if m in IMAGE_GROUPS:
        return "Image-based"
    if m in TRANSCRIPT_GROUPS:
        return "Transcript-based"
    return "Unknown"


# ─────────────────────────────────────────────────────────────────────────────
# Figure generators
# ─────────────────────────────────────────────────────────────────────────────

def plot_composite_scores(scores_df: pd.DataFrame) -> str:
    """Ranked composite score bar chart with group winner annotations."""
    agg = scores_df.groupby("method")["composite_score"].mean().sort_values(ascending=False).reset_index()

    fig, ax = plt.subplots(figsize=(max(6, len(agg) * 1.2), 4))
    bars = ax.bar(agg["method"], agg["composite_score"],
                  color=[method_color(m) for m in agg["method"]])

    # Annotate group winners
    img_best = agg[agg["method"].isin(IMAGE_GROUPS)].iloc[0]["method"] if not agg[agg["method"].isin(IMAGE_GROUPS)].empty else None
    tx_best = agg[agg["method"].isin(TRANSCRIPT_GROUPS)].iloc[0]["method"] if not agg[agg["method"].isin(TRANSCRIPT_GROUPS)].empty else None

    for bar, method in zip(bars, agg["method"]):
        label = ""
        if method == img_best:
            label = "★ Image winner"
        elif method == tx_best:
            label = "★ Transcript winner"
        if label:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                    label, ha="center", va="bottom", fontsize=8, color="black")

    ax.set_ylabel("Mean Composite Score")
    ax.set_title("Segmentation Method Ranking (Stage 2)")
    ax.set_ylim(0, 1.15)
    ax.axhline(agg["composite_score"].max(), linestyle="--", color="gray", alpha=0.5)
    fig.tight_layout()
    return fig_to_base64(fig)


def plot_mecr_heatmap(scores_df: pd.DataFrame) -> str:
    """MECR heatmap: methods × samples. Lower = greener = better."""
    pivot = scores_df.pivot(index="method", columns="sample", values="mecr")
    fig, ax = plt.subplots(figsize=(max(6, len(pivot.columns) * 1.0), max(3, len(pivot) * 0.7)))
    sns.heatmap(pivot, ax=ax, cmap="RdYlGn_r", annot=True, fmt=".3f",
                linewidths=0.5, vmin=0, vmax=0.3,
                cbar_kws={"label": "MECR (lower = better)"})
    ax.set_title("MECR per Method × Sample")
    fig.tight_layout()
    return fig_to_base64(fig)


def plot_transcript_recovery(scores_df: pd.DataFrame) -> str:
    """Grouped bar: transcript assignment rate per method × sample."""
    methods = scores_df["method"].unique()
    samples = sorted(scores_df["sample"].unique())
    x = np.arange(len(methods))
    width = 0.8 / max(len(samples), 1)

    fig, ax = plt.subplots(figsize=(max(6, len(methods) * 1.5), 4))
    for i, sample in enumerate(samples):
        sub = scores_df[scores_df["sample"] == sample].set_index("method")
        vals = [sub.loc[m, "assignment_rate"] if m in sub.index else 0 for m in methods]
        ax.bar(x + i * width, vals, width, label=sample,
               color=[method_color(m) for m in methods], alpha=0.7 + 0.1 * i)

    ax.set_xticks(x + width * (len(samples) - 1) / 2)
    ax.set_xticklabels(methods, rotation=15)
    ax.set_ylabel("Transcript Assignment Rate")
    ax.set_title("Transcript Recovery per Method and Sample")
    ax.set_ylim(0, 1.1)
    ax.legend(title="Sample", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    return fig_to_base64(fig)


def plot_cell_yield(scores_df: pd.DataFrame) -> str:
    """Cell count per method × sample."""
    methods = scores_df["method"].unique()
    samples = sorted(scores_df["sample"].unique())
    x = np.arange(len(methods))
    width = 0.8 / max(len(samples), 1)

    fig, ax = plt.subplots(figsize=(max(6, len(methods) * 1.5), 4))
    for i, sample in enumerate(samples):
        sub = scores_df[scores_df["sample"] == sample].set_index("method")
        vals = [sub.loc[m, "cell_count"] if m in sub.index else 0 for m in methods]
        ax.bar(x + i * width, vals, width, label=sample,
               color=[method_color(m) for m in methods], alpha=0.7 + 0.1 * i)

    ax.set_xticks(x + width * (len(samples) - 1) / 2)
    ax.set_xticklabels(methods, rotation=15)
    ax.set_ylabel("Cell Count")
    ax.set_title("Cell Yield per Method and Sample")
    ax.legend(title="Sample", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    return fig_to_base64(fig)


def plot_metrics_heatmap(scores_df: pd.DataFrame) -> str:
    """Methods × metrics normalized heatmap."""
    agg = scores_df.groupby("method").agg(
        MECR=("mecr", "mean"),
        Assignment_Rate=("assignment_rate", "mean"),
        Cell_Yield_Norm=("yield_norm", "mean"),
        Composite=("composite_score", "mean"),
    )
    # Normalize each column 0-1
    normed = agg.copy()
    for col in normed.columns:
        rng = normed[col].max() - normed[col].min()
        if rng > 0:
            normed[col] = (normed[col] - normed[col].min()) / rng
    # MECR: lower is better → invert
    normed["MECR"] = 1 - normed["MECR"]

    fig, ax = plt.subplots(figsize=(8, max(3, len(agg) * 0.6)))
    sns.heatmap(normed, ax=ax, cmap="YlGn", annot=True, fmt=".2f",
                linewidths=0.5, vmin=0, vmax=1,
                cbar_kws={"label": "Normalized score (higher = better)"})
    ax.set_title("Method × Metric Heatmap (MECR inverted: higher = better)")
    fig.tight_layout()
    return fig_to_base64(fig)


def plot_ap_heatmap(ap_df: pd.DataFrame, sample_id: str = None) -> str:
    """AP matrix heatmap: method_a (test) × method_b (reference), value = AP.

    Diagonal is NaN (no self-comparison). Higher = more similar segmentation.
    """
    pivot = ap_df.pivot(index="method_a", columns="method_b", values="ap")
    methods = sorted(set(pivot.index) | set(pivot.columns))
    pivot = pivot.reindex(index=methods, columns=methods)

    fig, ax = plt.subplots(figsize=(max(4, len(methods) * 1.0),
                                    max(3, len(methods) * 0.8)))
    mask = np.eye(len(methods), dtype=bool)   # mask diagonal (self-comparison)
    sns.heatmap(pivot, ax=ax, cmap="YlOrRd", annot=True, fmt=".3f",
                linewidths=0.5, vmin=0, vmax=1, mask=mask,
                cbar_kws={"label": "AP (precision × recall, mean over IoU thresholds)"})
    title = "Pairwise AP: test method (row) vs reference method (col)"
    if sample_id:
        title += f"\nSample: {sample_id}"
    ax.set_title(title)
    ax.set_xlabel("Reference method")
    ax.set_ylabel("Test method")
    fig.tight_layout()
    return fig_to_base64(fig)


def plot_stage1_param_heatmap(stage1_df: pd.DataFrame, method: str) -> str:
    """Stage 1 parameter search score heatmap for one method."""
    sub = stage1_df[stage1_df["method"] == method].copy()
    if sub.empty:
        return None

    # Pivot on top-2 most variable params
    non_score_cols = {"method", "param_hash", "mean_composite", "mean_assignment", "mean_cell_yield", "n_crops"}
    param_cols = [c for c in sub.columns if c not in non_score_cols]

    if len(param_cols) < 2:
        return None

    # Pick 2 most variable parameter columns
    variances = sub[param_cols].apply(lambda c: c.astype(str).nunique())
    top2 = variances.nlargest(2).index.tolist()
    p1, p2 = top2[0], top2[1]

    pivot = sub.pivot_table(values="mean_composite", index=p1, columns=p2, aggfunc="mean")
    fig, ax = plt.subplots(figsize=(max(4, len(pivot.columns)), max(3, len(pivot))))
    sns.heatmap(pivot, ax=ax, cmap="YlGn", annot=True, fmt=".3f",
                linewidths=0.5, cbar_kws={"label": "Mean Composite Score"})
    ax.set_title(f"{method.upper()} — Stage 1 Parameter Search ({p1} vs {p2})")
    fig.tight_layout()
    return fig_to_base64(fig)


def plot_roi_overlay(morphology_tif: str, geojson_paths: dict, roi: dict) -> str:
    """
    Overlay segmentation masks from multiple methods on a morphology image crop.
    roi: {x_min_um, x_max_um, y_min_um, y_max_um, pixel_size_um}
    geojson_paths: {method: path_to_cell_boundaries.geojson.gz}
    """
    px = roi.get("pixel_size_um", 0.2125)
    px_x_min = int(roi["x_min_um"] / px)
    px_x_max = int(roi["x_max_um"] / px)
    px_y_min = int(roi["y_min_um"] / px)
    px_y_max = int(roi["y_max_um"] / px)

    try:
        with tifffile.TiffFile(morphology_tif) as tif:
            import zarr
            z = zarr.open(tif.aszarr(level=0), mode="r")
            shape = z.shape
            if len(shape) == 3 and shape[0] <= 8:
                img = np.array(z[0, px_y_min:px_y_max, px_x_min:px_x_max])
            else:
                img = np.array(z[px_y_min:px_y_max, px_x_min:px_x_max])
    except Exception as e:
        print(f"  WARNING: Could not load morphology image: {e}", file=sys.stderr)
        return None

    n_methods = len(geojson_paths)
    fig, axes = plt.subplots(1, n_methods + 1, figsize=((n_methods + 1) * 4, 4))

    # Panel 0: raw morphology
    axes[0].imshow(img, cmap="gray")
    axes[0].set_title("Morphology (DAPI)")
    axes[0].axis("off")

    for ax, (method, gj_path) in zip(axes[1:], geojson_paths.items()):
        ax.imshow(img, cmap="gray")
        try:
            gdf = gpd.read_file(gj_path)
            # Translate polygons to crop coordinates
            from shapely.affinity import translate
            offset_x = -roi["x_min_um"]
            offset_y = -roi["y_min_um"]
            gdf_crop = gdf[gdf.intersects(gdf.total_bounds[0:4])]  # rough filter
            for geom in gdf_crop.geometry:
                if geom is None:
                    continue
                coords = np.array(geom.exterior.coords)
                px_coords = coords / px
                px_coords[:, 0] -= px_x_min
                px_coords[:, 1] -= px_y_min
                ax.plot(px_coords[:, 0], px_coords[:, 1],
                        color=method_color(method), linewidth=0.5, alpha=0.8)
        except Exception as e:
            print(f"  WARNING: Could not overlay {method} boundaries: {e}", file=sys.stderr)

        ax.set_title(method.upper())
        ax.axis("off")

    fig.suptitle(f"ROI overlay — x=[{roi['x_min_um']:.0f},{roi['x_max_um']:.0f}] "
                 f"y=[{roi['y_min_um']:.0f},{roi['y_max_um']:.0f}] µm")
    fig.tight_layout()
    return fig_to_base64(fig)


# ─────────────────────────────────────────────────────────────────────────────
# HTML template
# ─────────────────────────────────────────────────────────────────────────────

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Xenium Segmentation Search Report</title>
<style>
  body {{ font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; background: #f9f9f9; }}
  h1 {{ color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 10px; }}
  h2 {{ color: #34495e; margin-top: 40px; border-left: 4px solid #3498db; padding-left: 12px; }}
  h3 {{ color: #555; }}
  .caption {{ background: #eef6fb; border-left: 4px solid #3498db; padding: 10px 16px;
              margin: 8px 0 20px 0; font-size: 0.92em; color: #2c3e50; border-radius: 0 4px 4px 0; }}
  .winner-box {{ display: inline-block; background: #2ecc71; color: white;
                 border-radius: 4px; padding: 4px 10px; font-weight: bold; margin: 4px; }}
  .winner-box.transcript {{ background: #e67e22; }}
  .winner-box.overall {{ background: #2c3e50; }}
  table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
  th {{ background: #2c3e50; color: white; padding: 8px 12px; text-align: left; }}
  td {{ padding: 7px 12px; border-bottom: 1px solid #ddd; }}
  tr:nth-child(even) {{ background: #f2f2f2; }}
  tr.best {{ background: #d5f5e3 !important; font-weight: bold; }}
  img {{ max-width: 100%; border: 1px solid #ddd; border-radius: 4px; margin: 10px 0; }}
  .figure-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
  .fig-panel {{ background: white; border-radius: 6px; padding: 16px; box-shadow: 0 1px 4px rgba(0,0,0,0.1); }}
  @media (max-width: 900px) {{ .figure-grid {{ grid-template-columns: 1fr; }} }}
</style>
</head>
<body>
<h1>🔬 Xenium Segmentation Search Report</h1>
<p><strong>Generated:</strong> {timestamp}</p>
<p><strong>Samples:</strong> {n_samples} &nbsp;|&nbsp; <strong>Methods evaluated:</strong> {methods_list}</p>

<h2>1. Summary Table</h2>
<div class="caption">
  <strong>How to read:</strong> Methods are ranked by composite score (MECR 50% + transcript recovery 30% + cell yield 20%).
  MECR is the Mutually Exclusive Co-expression Rate — measures how often cell-type-specific marker genes from <em>different</em>
  lineages are co-expressed in the same cell. <strong>Lower MECR = better segmentation</strong> (less transcript bleed-through).
  Green rows = group winners.
</div>
{winners_html}
{summary_table}

<h2>2. Composite Score Ranking</h2>
<div class="caption">
  Mean composite score across all samples. Stars (★) mark the winner of each method group.
  The dashed line shows the top score. Higher = better.
</div>
<img src="data:image/png;base64,{composite_chart}">

<h2>3. MECR Heatmap</h2>
<div class="caption">
  Mutually Exclusive Co-expression Rate per method and sample. Color scale: <strong>green = low MECR = good</strong>,
  red = high MECR = transcript bleeding. A score below 0.05 is considered excellent; above 0.15 indicates
  significant segmentation error.
</div>
<img src="data:image/png;base64,{mecr_heatmap}">

<h2>4. Transcript Recovery</h2>
<div class="caption">
  Fraction of all detected transcripts that are assigned to a cell (vs. left as background/noise).
  Higher = more transcripts captured per cell. Values near 1.0 indicate the segmentation covers most
  of the tissue without leaving transcripts unassigned.
</div>
<img src="data:image/png;base64,{recovery_chart}">

<h2>5. Cell Yield</h2>
<div class="caption">
  Total number of cells detected per method and sample. Higher cell counts are generally better,
  but should be interpreted alongside MECR — a method that over-segments will show high cell count
  but poor MECR. Consider this metric in combination with transcript recovery.
</div>
<img src="data:image/png;base64,{yield_chart}">

<h2>6. Method × Metric Heatmap</h2>
<div class="caption">
  Normalized scores (0–1) across all metrics. MECR is inverted so that higher always means better.
  Use this to compare the overall profile of each method at a glance. A method that scores
  consistently across all metrics is more robust than one that excels at only one.
</div>
<img src="data:image/png;base64,{metrics_heatmap}">

{ap_section}

{roi_section}

<h2>9. Stage 1 Parameter Search</h2>
<div class="caption">
  Heatmaps of Stage 1 grid search scores for each method. Axes show the two most variable
  parameters tested; cell color = mean composite score (assignment rate 60% + cell yield 40%)
  across all 10 density-stratified crops. Brighter green = better. The parameter combination
  with the highest mean score was used for Stage 2.
</div>
{stage1_heatmaps}

</body>
</html>
"""


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores-dir", required=True)
    parser.add_argument("--ap-dir", default=None,
                        help="Directory containing *_ap_matrix.csv files from SCORE_AP")
    parser.add_argument("--stage1-summary", required=True)
    parser.add_argument("--stage1-manifest", required=True)
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    from datetime import datetime
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Load Stage 2 scores
    scores_dir = Path(args.scores_dir)
    score_files = list(scores_dir.glob("*_score.csv"))
    if not score_files:
        sys.exit(f"ERROR: No *_score.csv files found in {scores_dir}")

    scores_df = pd.concat([pd.read_csv(f) for f in score_files], ignore_index=True)
    n_samples = scores_df["sample"].nunique()
    methods = scores_df["method"].unique().tolist()

    # Load Stage 1 summary
    stage1_df = pd.read_csv(args.stage1_summary)
    with open(args.stage1_manifest) as f:
        param_manifest = json.load(f)
    # Merge param info into stage1_df
    for h, info in param_manifest.items():
        for k, v in info.get("params", {}).items():
            mask = stage1_df["param_hash"] == h
            stage1_df.loc[mask, k] = v

    # Generate figures
    composite_chart = plot_composite_scores(scores_df)
    mecr_heatmap = plot_mecr_heatmap(scores_df)
    recovery_chart = plot_transcript_recovery(scores_df)
    yield_chart = plot_cell_yield(scores_df)
    metrics_heatmap = plot_metrics_heatmap(scores_df)

    # Summary table HTML
    agg = scores_df.groupby("method").agg(
        MECR=("mecr", "mean"),
        Assignment_Rate=("assignment_rate", "mean"),
        Cell_Count=("cell_count", "mean"),
        Composite=("composite_score", "mean"),
        Group=("method", lambda x: group_label(x.iloc[0])),
    ).sort_values("Composite", ascending=False).reset_index()

    # Identify winners
    img_winner = agg[agg["method"].isin(IMAGE_GROUPS)].iloc[0]["method"] if not agg[agg["method"].isin(IMAGE_GROUPS)].empty else None
    tx_winner = agg[agg["method"].isin(TRANSCRIPT_GROUPS)].iloc[0]["method"] if not agg[agg["method"].isin(TRANSCRIPT_GROUPS)].empty else None
    overall_winner = agg.iloc[0]["method"]

    winners_html = ""
    if img_winner:
        winners_html += f'<span class="winner-box">🖼 Image-based winner: {img_winner.upper()}</span>'
    if tx_winner:
        winners_html += f'<span class="winner-box transcript">🧬 Transcript-based winner: {tx_winner.upper()}</span>'
    winners_html += f'<span class="winner-box overall">🏆 Overall winner: {overall_winner.upper()}</span>'
    winners_html = f"<p>{winners_html}</p>"

    table_rows = ""
    for _, row in agg.iterrows():
        best = " class='best'" if row["method"] == overall_winner else ""
        table_rows += (
            f"<tr{best}><td>{row['method'].upper()}</td>"
            f"<td>{row['Group']}</td>"
            f"<td>{row['MECR']:.4f}</td>"
            f"<td>{row['Assignment_Rate']:.4f}</td>"
            f"<td>{int(row['Cell_Count'])}</td>"
            f"<td><strong>{row['Composite']:.4f}</strong></td></tr>\n"
        )
    summary_table = (
        "<table><tr><th>Method</th><th>Group</th><th>MECR ↓</th>"
        "<th>Assignment Rate ↑</th><th>Cell Count ↑</th><th>Composite Score ↑</th></tr>\n"
        + table_rows + "</table>"
    )

    # AP heatmaps (one per sample)
    ap_section = ""
    if args.ap_dir:
        ap_dir = Path(args.ap_dir)
        ap_files = list(ap_dir.glob("*_ap_matrix.csv"))
        if ap_files:
            ap_section = (
                "<h2>7. Pairwise AP Scoring</h2>\n"
                "<div class='caption'>"
                "Average Precision (AP) between each pair of segmentation methods. "
                "AP(A→B) measures how well method A recovers the same cells as method B: "
                "for each IoU threshold from 0 to 1, a 1-to-1 greedy match is performed "
                "and <em>precision × recall</em> is computed. AP = mean over all thresholds. "
                "<strong>Higher AP = more similar segmentations.</strong> "
                "Note: AP is not symmetric — AP(A→B) can differ from AP(B→A) when the "
                "cell counts differ."
                "</div>\n"
            )
            for ap_file in sorted(ap_files):
                ap_df = pd.read_csv(ap_file)
                if ap_df.empty:
                    continue
                sample_ids = ap_df["sample_id"].unique()
                for sid in sample_ids:
                    sub = ap_df[ap_df["sample_id"] == sid]
                    b64 = plot_ap_heatmap(sub, sample_id=sid)
                    ap_section += (
                        f"<h3>Sample: {sid}</h3>"
                        f'<img src="data:image/png;base64,{b64}">'
                    )

    # Stage 1 heatmaps
    stage1_heatmaps = ""
    for method in methods:
        b64 = plot_stage1_param_heatmap(stage1_df, method)
        if b64:
            stage1_heatmaps += (
                f"<h3>{method.upper()}</h3>"
                f'<img src="data:image/png;base64,{b64}">'
            )

    # ROI section (placeholder — populated if bundles dir provided)
    roi_section = (
        "<h2>8. ROI Visual Checks</h2>"
        "<div class='caption'>ROI overlays require --bundles and --morphology-dir arguments. "
        "Run with those arguments to generate segmentation mask overlays on morphology images.</div>"
    )

    # Render HTML
    html = HTML_TEMPLATE.format(
        timestamp=timestamp,
        n_samples=n_samples,
        methods_list=", ".join(m.upper() for m in methods),
        winners_html=winners_html,
        summary_table=summary_table,
        composite_chart=composite_chart,
        mecr_heatmap=mecr_heatmap,
        recovery_chart=recovery_chart,
        yield_chart=yield_chart,
        metrics_heatmap=metrics_heatmap,
        ap_section=ap_section,
        roi_section=roi_section,
        stage1_heatmaps=stage1_heatmaps,
    )

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "seg_search_report.html"
    out_path.write_text(html)
    print(f"Report written to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
