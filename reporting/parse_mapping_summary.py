#!/usr/bin/env python3

from __future__ import annotations

import csv
import os
import sys
from datetime import datetime


def read_summary(path: str) -> tuple[list[str], list[str]]:
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        tokens = [token.strip() for token in handle.read().replace("\n", ",").split(",") if token.strip()]
    if not tokens:
        raise ValueError(f"empty mapping summary: {path}")
    if any(not token.replace(".", "", 1).isdigit() for token in tokens[:5]):
        header = tokens[:5]
        values = tokens[5:10]
    else:
        header = ["barcode", "total", "duplicate", "unmapped", "lowmapq"]
        values = tokens[:5]
    if len(values) < len(header):
        values.extend(["0"] * (len(header) - len(values)))
    return header, values


def main() -> int:
    if len(sys.argv) != 5:
        print("Usage: parse_mapping_summary.py <summary_csv> <sample_name> <lane_label> <output_csv>")
        return 1

    summary_path, sample_name, lane_label, output_path = sys.argv[1:5]
    header, values = read_summary(summary_path)
    parsed = dict(zip(header, values))
    total = float(parsed.get("total", "0") or 0)
    unmapped = float(parsed.get("unmapped", "0") or 0)
    lowmapq = float(parsed.get("lowmapq", "0") or 0)
    duplicate = float(parsed.get("duplicate", "0") or 0)
    mapped = max(total - unmapped, 0)
    high_quality = max(total - unmapped - lowmapq, 0)

    rows = [
        ("mapping", "sample_name", sample_name, "-", lane_label),
        ("mapping", "lane_label", lane_label, "-", lane_label),
        ("mapping", "chromap_total_reads", f"{int(total)}", "reads", os.path.basename(summary_path)),
        ("mapping", "chromap_unmapped_reads", f"{int(unmapped)}", "reads", os.path.basename(summary_path)),
        ("mapping", "chromap_duplicate_reads", f"{int(duplicate)}", "reads", os.path.basename(summary_path)),
        ("mapping", "chromap_lowmapq_reads", f"{int(lowmapq)}", "reads", os.path.basename(summary_path)),
        ("mapping", "chromap_mapped_reads", f"{int(mapped)}", "reads", os.path.basename(summary_path)),
        ("mapping", "chromap_high_quality_reads", f"{int(high_quality)}", "reads", os.path.basename(summary_path)),
        ("mapping", "chromap_mapping_rate", f"{(mapped * 100 / total) if total else 0:.2f}", "%", os.path.basename(summary_path)),
        ("mapping", "chromap_high_quality_rate", f"{(high_quality * 100 / total) if total else 0:.2f}", "%", os.path.basename(summary_path)),
    ]

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["section", "metric", "value", "unit", "source", "timestamp"])
        for row in rows:
            writer.writerow([*row, timestamp])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
