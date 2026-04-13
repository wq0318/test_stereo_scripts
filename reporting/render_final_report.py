#!/usr/bin/env python3

from __future__ import annotations

import argparse
import base64
import csv
import html
import os
from collections import defaultdict
from pathlib import Path


def read_metrics(paths: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle)
            rows.extend(reader)
    return rows


def parse_path_list(value: str) -> list[str]:
    return [item for item in value.split(",") if item]


def embed_image(path: str) -> str:
    if not path or not os.path.exists(path):
        return ""
    suffix = Path(path).suffix.lower().lstrip(".") or "png"
    with open(path, "rb") as handle:
        encoded = base64.b64encode(handle.read()).decode("ascii")
    return f"data:image/{suffix};base64,{encoded}"


def aggregate_global_metrics(rows: list[dict[str, str]]) -> dict[str, str]:
    numeric_sums: defaultdict[str, float] = defaultdict(float)
    latest_values: dict[str, str] = {}
    for row in rows:
        metric = row["metric"]
        value = row["value"]
        latest_values[metric] = value
        try:
            numeric_sums[metric] += float(value)
        except ValueError:
            continue

    total_reads = numeric_sums.get("total_reads", 0)
    mapped_reads = numeric_sums.get("mapped_reads", 0)
    exact_reads = numeric_sums.get("barcode_exactly_overlap_reads", 0)
    mis_reads = numeric_sums.get("barcode_mis_overlap_reads", 0)
    chromap_total = numeric_sums.get("chromap_total_reads", 0)
    chromap_hq = numeric_sums.get("chromap_high_quality_reads", 0)

    latest_values["barcode_mapping_rate"] = f"{(mapped_reads * 100 / total_reads) if total_reads else 0:.2f}"
    latest_values["barcode_exact_rate"] = f"{(exact_reads * 100 / total_reads) if total_reads else 0:.2f}"
    latest_values["barcode_mis_rate"] = f"{(mis_reads * 100 / total_reads) if total_reads else 0:.2f}"
    latest_values["chromap_high_quality_rate_total"] = f"{(chromap_hq * 100 / chromap_total) if chromap_total else 0:.2f}"
    latest_values["total_reads_sum"] = f"{int(total_reads)}"
    latest_values["mapped_reads_sum"] = f"{int(mapped_reads)}"
    latest_values["chromap_total_reads_sum"] = f"{int(chromap_total)}"
    latest_values["chromap_high_quality_reads_sum"] = f"{int(chromap_hq)}"
    return latest_values


def render_metric_table(title: str, items: list[tuple[str, str]]) -> str:
    rows = "".join(
        f"<tr><td>{html.escape(label)}</td><td>{html.escape(value)}</td></tr>"
        for label, value in items
    )
    return f"""
    <section class="card">
      <h2>{html.escape(title)}</h2>
      <table>{rows}</table>
    </section>
    """


def render_image(title: str, data_uri: str) -> str:
    if not data_uri:
        return f'<section class="card"><h2>{html.escape(title)}</h2><p>Image not available.</p></section>'
    return f'<section class="card"><h2>{html.escape(title)}</h2><img src="{data_uri}" alt="{html.escape(title)}" /></section>'


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_name")
    parser.add_argument("bin_size")
    parser.add_argument("output_html")
    parser.add_argument("barcode_metrics")
    parser.add_argument("mapping_metrics")
    parser.add_argument("sam2frag_metrics")
    parser.add_argument("merge_metrics")
    parser.add_argument("downstream_metrics")
    parser.add_argument("--saturation-plot", required=True)
    parser.add_argument("--qc-violin", required=True)
    parser.add_argument("--tss-scatter", required=True)
    parser.add_argument("--fragment-size", required=True)
    parser.add_argument("--spatial-qc", required=True)
    parser.add_argument("--cluster-plots", default="")
    args = parser.parse_args()

    global_rows = read_metrics(
        [
            *parse_path_list(args.barcode_metrics),
            *parse_path_list(args.mapping_metrics),
            *parse_path_list(args.sam2frag_metrics),
            args.merge_metrics,
            args.downstream_metrics,
        ]
    )
    metrics = aggregate_global_metrics(global_rows)

    left_items = [
        ("Total reads", metrics.get("total_reads_sum", "0")),
        ("Valid CID mapped reads", metrics.get("mapped_reads_sum", "0")),
        ("Barcode mapping rate", f'{metrics.get("barcode_mapping_rate", "0")} %'),
        ("Exact overlap rate", f'{metrics.get("barcode_exact_rate", "0")} %'),
        ("Mismatch overlap rate", f'{metrics.get("barcode_mis_rate", "0")} %'),
        ("Chromap total reads", metrics.get("chromap_total_reads_sum", "0")),
        ("High-quality chromap reads", metrics.get("chromap_high_quality_reads_sum", "0")),
        ("High-quality chromap rate", f'{metrics.get("chromap_high_quality_rate_total", "0")} %'),
        ("Raw reads after sam2frag", metrics.get("raw_reads", "0")),
        ("Mitochondria ratio", f'{metrics.get("mitochondria_ratio", "0")} %'),
        ("Merged total fragments", metrics.get("total_frag_count", "0")),
        ("Merged unique fragments", metrics.get("total_unique_frag_count", "0")),
    ]
    right_items = [
        ("Bin size", args.bin_size),
        ("Saturation", metrics.get("saturation", "0")),
        ("chrM ratio", f'{metrics.get("chrM_ratio", "0")} %'),
        ("Total barcodes", metrics.get("total_barcodes", "0")),
        ("Valid barcodes", metrics.get("valid_barcodes", "0")),
        ("Filtered fragments", metrics.get("filtered_fragments", "0")),
        ("Median nFrags", metrics.get("median_nFrags", "0")),
        ("Median TSS", metrics.get("median_TSS", "0")),
        ("Median FRiP", metrics.get("median_FRiP", "0")),
        ("Median blacklist", metrics.get("median_blacklist", "0")),
        ("Cell count", metrics.get("n_cells", "0")),
        ("Filtered fragment size", metrics.get("filtered_fragment_size", "0")),
    ]

    cluster_cards = "".join(
        render_image(f"Cluster plot {index + 1}", embed_image(path))
        for index, path in enumerate(parse_path_list(args.cluster_plots))
    )

    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>{html.escape(args.sample_name)} bin{html.escape(args.bin_size)} Stereo-ATAC report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 0; background: #f4f6f8; color: #1f2937; }}
    .page {{ max-width: 1200px; margin: 0 auto; padding: 24px; }}
    .hero {{ background: linear-gradient(135deg, #0f4c81, #3b82f6); color: white; padding: 24px; border-radius: 16px; margin-bottom: 24px; }}
    .metrics {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-bottom: 24px; }}
    .card {{ background: white; border-radius: 14px; padding: 20px; box-shadow: 0 4px 14px rgba(15, 23, 42, 0.08); margin-bottom: 20px; }}
    h1, h2 {{ margin-top: 0; }}
    table {{ width: 100%; border-collapse: collapse; }}
    td {{ padding: 8px 0; border-bottom: 1px solid #e5e7eb; }}
    td:last-child {{ text-align: right; font-weight: 600; }}
    img {{ max-width: 100%; border-radius: 10px; }}
    .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
    @media (max-width: 900px) {{ .metrics, .grid {{ grid-template-columns: 1fr; }} }}
  </style>
</head>
<body>
  <div class="page">
    <section class="hero">
      <h1>{html.escape(args.sample_name)} Stereo-ATAC report</h1>
      <p>Bin size: {html.escape(args.bin_size)}</p>
    </section>
    <section class="metrics">
      {render_metric_table("Quality Control Metrics", left_items)}
      {render_metric_table("Downstream Metrics", right_items)}
    </section>
    {render_image("Saturation", embed_image(args.saturation_plot))}
    {render_image("QC violin", embed_image(args.qc_violin))}
    <section class="grid">
      {render_image("TSS scatter", embed_image(args.tss_scatter))}
      {render_image("Fragment size", embed_image(args.fragment_size))}
    </section>
    {render_image("Spatial QC", embed_image(args.spatial_qc))}
    <section class="card">
      <h2>Cluster plots</h2>
      <div class="grid">
        {cluster_cards}
      </div>
    </section>
  </div>
</body>
</html>
"""

    with open(args.output_html, "w", encoding="utf-8") as handle:
        handle.write(html_text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
