#!/usr/bin/env python3
"""Validate the real DAF/FiberHMM artifacts before report rendering."""

from __future__ import annotations

import argparse
from pathlib import Path

from real_stage_common import load_config, run_root, write_marker, write_tsv


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tmp", required=True, type=Path)
    args = parser.parse_args()
    cfg = load_config(args.configfile)
    root = run_root(cfg, args.project_root)
    daf_root = root / "05_daf_signal"
    fiber_root = root / "06_fiberhmm"
    required = [
        daf_root / "tables/approximate_daf_metric.tsv",
        daf_root / "tables/read_flow_reconciliation.tsv",
        daf_root / "tables/event_reconciliation.tsv",
        fiber_root / "DAFseq_fiberhmm_real_report.html",
        fiber_root / "FIBERHMM_EXECUTION_LOG.md",
        fiber_root / "FINAL_VALIDATION.md",
        fiber_root / "fiberhmm/canonical/BC01.fiberhmm.sorted.bam",
        fiber_root / "fiberhmm/canonical/BC01.fiberhmm.sorted.bam.bai",
        fiber_root / "fiberhmm/chrM/BC01.fiberhmm.chrM.sorted.bam",
        fiber_root / "fiberhmm/chrM/BC01.fiberhmm.chrM.sorted.bam.bai",
        fiber_root / "qc/canonical/BC01.fiberhmm.sorted.qc.json",
        fiber_root / "tables/FIBERHMM_QC_SUMMARY.tsv",
        fiber_root / "tables/MSP_SUMMARY.tsv",
        fiber_root / "tables/NUCLEOSOME_SUMMARY.tsv",
        fiber_root / "tables/TF_FOOTPRINT_SUMMARY.tsv",
    ]
    rows = [{"artifact": str(path.relative_to(root)), "status": "PASS" if path.is_file() and path.stat().st_size > 0 else "FAIL", "definition": "required real DAF/FiberHMM delivery artifact"} for path in required]
    write_tsv(root / "07_delivery/DELIVERY_ARTIFACT_VALIDATION.tsv", rows, ["artifact", "status", "definition"])
    failed = [row for row in rows if row["status"] != "PASS"]
    if failed:
        raise RuntimeError("delivery artifact validation failed: " + ", ".join(row["artifact"] for row in failed))
    write_marker(args.output, cfg, "delivery", Path(args.tmp), input_path=args.input, extra={"artifact_count": len(rows), "artifact_validation": "PASS"})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
