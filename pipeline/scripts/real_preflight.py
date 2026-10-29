#!/usr/bin/env python3
"""Validate real Qitan inputs and the external FiberHMM toolchain."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pysam

from real_stage_common import load_config, run_root, write_marker, write_tsv


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tmp", required=True, type=Path)
    args = parser.parse_args()
    cfg = load_config(args.configfile)
    if cfg.get("execution", {}).get("mode") != "real":
        raise ValueError("real_preflight requires execution.mode=real")

    paths = cfg["paths"]
    reference = cfg["reference"]
    runtime = cfg["runtime"]
    required = {
        "raw_fastq": Path(paths["raw_fastq"]),
        "source_root": Path(paths["source_root"]),
        "input_bam": Path(paths["input_bam"]),
        "input_bai": Path(str(paths["input_bam"]) + ".bai"),
        "reference_fasta": Path(reference["fasta"]),
        "reference_fai": Path(str(reference["fasta"]) + ".fai"),
        "reference_gtf": Path(reference["gtf"]),
        "fiberhmm_bin": Path(runtime["fiberhmm_bin"]),
        "fiberhmm_qc": Path(runtime["fiberhmm_qc"]),
        "fiberhmm_extract": Path(runtime["fiberhmm_extract"]),
        "samtools": Path(runtime["samtools"]),
    }
    rows = []
    failures = []
    for label, path in required.items():
        ok = path.exists()
        rows.append({"check": "path_exists", "item": label, "status": "PASS" if ok else "FAIL", "value": str(path), "detail": "configured immutable input/tool"})
        if not ok:
            failures.append(f"{label}: {path}")

    bam_path = required["input_bam"]
    if bam_path.exists():
        with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
            contigs = set(bam.references)
            canonical = set(cfg["chromosomes"]["canonical"])
            allowed = canonical | {str(cfg["chromosomes"]["mitochondrial"]) }
            rows.append({"check": "bam_contigs", "item": "canonical_and_chrM_present", "status": "PASS" if bool(contigs & allowed) else "FAIL", "value": ",".join(sorted(contigs & allowed)), "detail": "input BAM contig names"})
            seen = 0
            md = mm = ml = False
            for rec in bam.fetch(until_eof=True):
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                    continue
                seen += 1
                md = md or rec.has_tag("MD")
                mm = mm or rec.has_tag("MM") or rec.has_tag("Mm")
                ml = ml or rec.has_tag("ML") or rec.has_tag("Ml")
                if seen >= 1000:
                    break
            rows.extend([
                {"check": "bam_tags", "item": "MD", "status": "PASS" if md else "WARN", "value": str(md).lower(), "detail": "FiberHMM DAF fallback tag"},
                {"check": "bam_tags", "item": "MM/Mm", "status": "PASS" if mm else "WARN", "value": str(mm).lower(), "detail": "optional modification tag"},
                {"check": "bam_tags", "item": "ML/Ml", "status": "PASS" if ml else "WARN", "value": str(ml).lower(), "detail": "optional modification probability tag"},
            ])

    try:
        probe = subprocess.run([str(runtime["python"]), "-c", "import importlib.metadata; print(importlib.metadata.version('fiberhmm'))"], capture_output=True, text=True, check=False)
        version = probe.stdout.strip() if probe.returncode == 0 else "NOT_AVAILABLE"
    except OSError:
        version = "NOT_AVAILABLE"
    expected = str(cfg["fiberhmm"]["version"])
    rows.append({"check": "tool_version", "item": "fiberhmm", "status": "PASS" if version == expected else "FAIL", "value": version, "detail": f"expected {expected}"})
    if version != expected:
        failures.append(f"fiberhmm version {version} != {expected}")

    if failures:
        write_tsv(run_root(cfg, args.project_root) / "01_preflight" / "PREFLIGHT_CHECKS.tsv", rows, ["check", "item", "status", "value", "detail"])
        raise RuntimeError("real preflight failed: " + "; ".join(failures))
    write_tsv(run_root(cfg, args.project_root) / "01_preflight" / "PREFLIGHT_CHECKS.tsv", rows, ["check", "item", "status", "value", "detail"])
    write_marker(args.output, cfg, "preflight", Path(args.tmp), extra={"fiberhmm_version": version, "input_tags_checked": "MD,MM/ML"})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
