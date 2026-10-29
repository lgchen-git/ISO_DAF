#!/usr/bin/env python3
"""Compute the real FASTQ-level QC stage for a configured Qitan sample."""

from __future__ import annotations

import argparse
import gzip
import statistics
from pathlib import Path

from real_stage_common import load_config, run_root, write_marker, write_tsv


def fastq_handle(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.name.endswith((".gz", ".gzip")) else path.open("r", encoding="utf-8")


def n50(values: list[int]) -> int:
    if not values:
        return 0
    ordered = sorted(values, reverse=True)
    target = sum(ordered) / 2
    running = 0
    for value in ordered:
        running += value
        if running >= target:
            return value
    return ordered[-1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--input", required=False, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tmp", required=True, type=Path)
    args = parser.parse_args()
    config = load_config(args.configfile)
    fastq = Path(config["paths"]["raw_fastq"])
    if not fastq.is_file():
        raise FileNotFoundError(f"configured FASTQ does not exist: {fastq}")

    lengths: list[int] = []
    mean_qualities: list[float] = []
    total_bases = 0
    with fastq_handle(fastq) as handle:
        while True:
            name = handle.readline()
            if not name:
                break
            sequence = handle.readline().rstrip("\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\n")
            if not plus or not quality or not name.startswith("@"):
                raise ValueError(f"malformed FASTQ record near {name[:40]!r}")
            if len(sequence) != len(quality):
                raise ValueError("FASTQ sequence and quality lengths differ")
            lengths.append(len(sequence))
            total_bases += len(sequence)
            mean_qualities.append(sum(max(0, ord(char) - 33) for char in quality) / len(quality) if quality else 0.0)

    rows = [
        {"metric": "reads", "value": str(len(lengths)), "unit": "reads", "definition": "FASTQ records"},
        {"metric": "total_bases", "value": str(total_bases), "unit": "bp", "definition": "sum of FASTQ sequence lengths"},
        {"metric": "median_length", "value": f"{statistics.median(lengths):.0f}", "unit": "bp", "definition": "read-level median sequence length"},
        {"metric": "N50", "value": str(n50(lengths)), "unit": "bp", "definition": "yield-weighted read length N50"},
        {"metric": "median_read_quality", "value": f"{statistics.median(mean_qualities):.2f}", "unit": "Q score", "definition": "median per-read mean Phred quality"},
    ]
    stage_root = run_root(config, args.project_root) / "02_raw_qc"
    write_tsv(stage_root / "raw_read_qc.tsv", rows, ["metric", "value", "unit", "definition"])
    write_marker(args.output, config, "raw_qc", Path(args.tmp), input_path=Path(str(args.output)), extra={"fastq": fastq, "read_count": len(lengths)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
