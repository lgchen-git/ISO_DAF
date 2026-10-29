#!/usr/bin/env python3
"""Add the independently executed chrM FiberHMM scope to delivery tables."""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

import pysam


def main(run_root: Path) -> None:
    bam = run_root / "fiberhmm" / "chrM" / "BC01.fiberhmm.chrM.sorted.bam"
    bed_dir = run_root / "fiberhmm" / "chrM" / "features"
    stats = Counter()
    fibers = set()
    spans = []
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        for rec in handle.fetch(until_eof=True):
            stats["records"] += 1
            if rec.is_unmapped:
                continue
            if rec.is_secondary:
                stats["secondary"] += 1
            if rec.is_supplementary:
                stats["supplementary"] += 1
            if rec.is_duplicate:
                stats["duplicate"] += 1
            if not rec.is_secondary and not rec.is_supplementary:
                stats["primary_mapped"] += 1
            if rec.mapping_quality >= 20:
                stats["mapq20"] += 1
            if rec.is_secondary or rec.is_supplementary or rec.mapping_quality < 20:
                continue
            stats["eligible_records"] += 1
            fibers.add(rec.query_name)
            spans.append((rec.reference_start, rec.reference_end))
    stats["fibers"] = len(fibers)
    spans.sort()
    covered = 0
    if spans:
        start, end = spans[0]
        for left, right in spans[1:]:
            if left <= end:
                end = max(end, right)
            else:
                covered += end - start
                start, end = left, right
        covered += end - start
    stats["covered_bp_union"] = covered
    feature = Counter()
    for path in bed_dir.glob("*.bed"):
        name = path.name
        kind = next((x for x in ("msp", "nucleosome", "tf", "deam") if f"_{x}.bed" in name), None)
        if not kind:
            continue
        with path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    continue
                try:
                    n = int(fields[9])
                    sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
                except ValueError:
                    continue
                feature[kind] += min(n, len(sizes))
                if kind == "msp":
                    feature["msp_gt150"] += sum(x > 150 for x in sizes[:n])
    row = {
        "scope": "chrM_fiberhmm_output",
        "records": stats["records"],
        "scope_records": stats["records"],
        "primary_mapped": stats["primary_mapped"],
        "supplementary": stats["supplementary"],
        "secondary": stats["secondary"],
        "duplicate": stats["duplicate"],
        "mapq20": stats["mapq20"],
        "eligible_records": stats["eligible_records"],
        "fibers": stats["fibers"],
        "covered_bp_union": stats["covered_bp_union"],
        "output_msp": feature["msp"],
        "output_msp_gt150": feature["msp_gt150"],
        "output_nucleosome": feature["nucleosome"],
        "output_tf": feature["tf"],
        "output_deam": feature["deam"],
    }
    path = run_root / "tables" / "CANONICAL_FILTER_SUMMARY.tsv"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
        fields = handle.seek(0) or next(csv.reader(handle, delimiter="\t"))
    rows = [r for r in rows if r.get("scope") != "chrM_fiberhmm_output"]
    rows.append({k: str(row.get(k, "NA")) for k in fields})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    out = run_root / "tables" / "CHRM_FIBERHMM_SUMMARY.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
