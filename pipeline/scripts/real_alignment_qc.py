#!/usr/bin/env python3
"""Scan the configured BAM and write alignment-level QC summaries."""

from __future__ import annotations

import argparse
import statistics
from collections import Counter
from pathlib import Path

import pysam

from real_stage_common import load_config, run_root, write_marker, write_tsv


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--input", required=False, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tmp", required=True, type=Path)
    args = parser.parse_args()
    cfg = load_config(args.configfile)
    bam_path = Path(cfg["paths"]["input_bam"])
    canonical = set(cfg["chromosomes"]["canonical"])
    chr_m = str(cfg["chromosomes"]["mitochondrial"])
    mapq = int(cfg["thresholds"]["mapq"])
    counters = Counter()
    aligned_fractions: list[float] = []
    softclip_fractions: list[float] = []
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        for rec in bam.fetch(until_eof=True):
            counters["records"] += 1
            if rec.is_secondary:
                counters["secondary_records"] += 1
                continue
            if rec.is_supplementary:
                counters["supplementary_records"] += 1
                continue
            counters["primary_records"] += 1
            if rec.is_unmapped:
                counters["primary_unmapped"] += 1
                continue
            counters["primary_mapped"] += 1
            if rec.reference_name in canonical:
                counters["canonical_nuclear"] += 1
            elif rec.reference_name == chr_m:
                counters["chrM"] += 1
            else:
                counters["noncanonical"] += 1
            if rec.mapping_quality >= mapq:
                counters["mapq_pass"] += 1
            for threshold in (20, 40, 60):
                if rec.mapping_quality >= threshold:
                    counters[f"mapq_{threshold}"] += 1
            query_length = rec.query_length or 0
            aligned = rec.query_alignment_length or 0
            softclip = sum(length for op, length in (rec.cigartuples or []) if op == 4)
            aligned_fraction = aligned / query_length if query_length else 0.0
            softclip_fraction = softclip / query_length if query_length else 0.0
            aligned_fractions.append(aligned_fraction)
            softclip_fractions.append(softclip_fraction)
            for threshold in (0.90, 0.95):
                if aligned_fraction >= threshold:
                    counters[f"aligned_fraction_{threshold:.2f}"] += 1
            for threshold in (0.10, 0.05):
                if softclip_fraction <= threshold:
                    counters[f"softclip_fraction_le_{threshold:.2f}"] += 1

    primary = counters["primary_records"]
    mapped = counters["primary_mapped"]
    rows = []

    def add(metric: str, count: int | float, denominator: int | float, unit: str, definition: str) -> None:
        rate = float(count) / float(denominator) if denominator else 0.0
        rows.append({"metric": metric, "count": str(int(count)) if float(count).is_integer() else f"{float(count):.4f}", "denominator": str(int(denominator)) if float(denominator).is_integer() else f"{float(denominator):.4f}", "percent": f"{100 * rate:.2f}%", "unit": unit, "definition": definition})

    add("primary_mapping_rate", mapped, primary, "percent", "primary mapped records / primary BAM records")
    add("canonical_nuclear_fraction", counters["canonical_nuclear"], mapped, "percent", "canonical nuclear primary mapped records / primary mapped records")
    for threshold in (20, 40, 60):
        add(f"MAPQ_ge_{threshold}", counters[f"mapq_{threshold}"], mapped, "percent", f"primary mapped records with MAPQ >= {threshold}")
    for threshold in (0.90, 0.95):
        add(f"aligned_fraction_ge_{threshold:.2f}", counters[f"aligned_fraction_{threshold:.2f}"], primary, "percent", f"primary records with aligned query fraction >= {threshold:.2f}")
    for threshold in (0.10, 0.05):
        add(f"softclip_fraction_le_{threshold:.2f}", counters[f"softclip_fraction_le_{threshold:.2f}"], primary, "percent", f"primary records with softclip fraction <= {threshold:.2f}")
    add("supplementary_records_per_primary_record", counters["supplementary_records"], primary, "ratio", "supplementary records / primary records")
    add("secondary_records_per_primary_record", counters["secondary_records"], primary, "ratio", "secondary records / primary records")
    if aligned_fractions:
        rows.extend([
            {"metric": "aligned_fraction_median", "count": f"{statistics.median(aligned_fractions):.4f}", "denominator": str(primary), "percent": "NA", "unit": "ratio", "definition": "median aligned query fraction across primary records"},
            {"metric": "softclip_fraction_median", "count": f"{statistics.median(softclip_fractions):.4f}", "denominator": str(primary), "percent": "NA", "unit": "ratio", "definition": "median softclip fraction across primary records"},
        ])
    stage_root = run_root(cfg, args.project_root) / "03_alignment_qc"
    write_tsv(stage_root / "alignment_quality_summary.tsv", rows, ["metric", "count", "denominator", "percent", "unit", "definition"])
    write_tsv(stage_root / "alignment_record_counts.tsv", [{"metric": k, "count": str(v)} for k, v in sorted(counters.items())], ["metric", "count"])
    write_marker(args.output, cfg, "alignment_qc", Path(args.tmp), input_path=Path(str(args.output)), extra={"bam": bam_path, "primary_records": primary, "primary_mapped": mapped})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
