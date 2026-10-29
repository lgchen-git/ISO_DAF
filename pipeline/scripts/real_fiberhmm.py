#!/usr/bin/env python3
"""Run the real FiberHMM DddA workflow and build its derived delivery."""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

from real_stage_common import load_config, run_command, run_root, write_marker


def sorted_bam(samtools: Path, unsorted: Path, sorted_path: Path, log_path: Path, threads: int) -> None:
    run_command([str(samtools), "sort", "-@", str(threads), "-o", str(sorted_path), str(unsorted)], log_path)
    run_command([str(samtools), "index", "-@", str(max(1, min(threads, 4))), str(sorted_path)], log_path.with_name(log_path.stem + "_index.log"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tmp", required=True, type=Path)
    args = parser.parse_args()
    cfg = load_config(args.configfile)
    if cfg.get("execution", {}).get("mode") != "real":
        raise ValueError("real_fiberhmm requires execution.mode=real")

    stage_root = run_root(cfg, args.project_root) / "06_fiberhmm"
    fiber_root = stage_root / "fiberhmm"
    canonical_dir = fiber_root / "canonical"
    chrm_dir = fiber_root / "chrM"
    qc_root = stage_root / "qc"
    for path in (canonical_dir, chrm_dir, qc_root / "canonical", qc_root / "chrM", stage_root / "logs"):
        path.mkdir(parents=True, exist_ok=True)
    runtime = cfg["runtime"]
    input_bam = Path(cfg["paths"]["input_bam"])
    reference = Path(cfg["reference"]["fasta"])
    samtools = Path(runtime["samtools"])
    call_bin = Path(runtime["fiberhmm_bin"])
    extract_bin = Path(runtime["fiberhmm_extract"])
    qc_bin = Path(runtime["fiberhmm_qc"])
    cores = int(runtime.get("cores", 12))
    io_threads = int(runtime.get("io_threads", 2))
    threshold = int(cfg["fiberhmm"]["probability_threshold"])
    mapq = int(cfg["thresholds"]["mapq"])
    canonical = list(cfg["chromosomes"]["canonical"])
    chrm = str(cfg["chromosomes"]["mitochondrial"])
    env = {
        "NUMBA_DISABLE_JIT": "1" if cfg.get("execution", {}).get("numba_disable_jit", True) else "0",
        "MPLCONFIGDIR": str(Path(args.tmp) / "matplotlib"),
        "TMPDIR": str(Path(args.tmp)),
    }
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    command_lines = []

    canonical_unsorted = canonical_dir / "BC01.fiberhmm.bam"
    canonical_sorted = canonical_dir / "BC01.fiberhmm.sorted.bam"
    canonical_command = [str(call_bin), "-i", str(input_bam), "-o", str(canonical_unsorted), "--enzyme", str(cfg["fiberhmm"]["enzyme"]), "--seq", "nanopore", "--reference", str(reference), "--min-mapq", str(mapq), "--prob-threshold", str(threshold), "--primary", "--no-dedup", "--no-qc", "--region-parallel", "--region-size", str(runtime.get("region_size", 20000000)), "--chroms", *canonical, "--cores", str(cores), "--io-threads", str(io_threads)]
    command_lines.append(" ".join(canonical_command))
    run_command(canonical_command, stage_root / "logs/fiberhmm_call_canonical.log", env=env)
    sorted_bam(samtools, canonical_unsorted, canonical_sorted, stage_root / "logs/samtools_sort_canonical.log", cores)

    chrm_unsorted = chrm_dir / "BC01.fiberhmm.chrM.bam"
    chrm_sorted = chrm_dir / "BC01.fiberhmm.chrM.sorted.bam"
    chrm_command = [str(call_bin), "-i", str(input_bam), "-o", str(chrm_unsorted), "--enzyme", str(cfg["fiberhmm"]["enzyme"]), "--seq", "nanopore", "--reference", str(reference), "--min-mapq", str(mapq), "--prob-threshold", str(threshold), "--primary", "--no-dedup", "--no-qc", "--no-daf-call-snps", "--region-parallel", "--chroms", chrm, "--cores", str(max(1, min(cores, 4))), "--io-threads", str(io_threads)]
    command_lines.append(" ".join(chrm_command))
    run_command(chrm_command, stage_root / "logs/fiberhmm_call_chrM.log", env=env)
    sorted_bam(samtools, chrm_unsorted, chrm_sorted, stage_root / "logs/samtools_sort_chrM.log", max(1, min(cores, 4)))

    extract_common = [str(extract_bin), "-i", str(canonical_sorted), "-o", str(canonical_dir / "features"), "--all", "--bed-only", "--keep-bed", "-q", str(mapq), "-p", str(threshold), "--chroms", ",".join(canonical), "--region-size", str(runtime.get("region_size", 20000000)), "--sort-parallel", str(max(1, min(cores, 8))), "-c", str(max(1, min(cores, 8)))]
    command_lines.append(" ".join(extract_common))
    run_command(extract_common, stage_root / "logs/fiberhmm_extract_canonical.log", env=env)
    extract_chrm = [str(extract_bin), "-i", str(chrm_sorted), "-o", str(chrm_dir / "features"), "--all", "--bed-only", "--keep-bed", "-q", str(mapq), "-p", str(threshold), "--chroms", chrm, "-c", str(max(1, min(cores, 4)))]
    command_lines.append(" ".join(extract_chrm))
    run_command(extract_chrm, stage_root / "logs/fiberhmm_extract_chrM.log", env=env)

    qc_common = [str(qc_bin), "-i", str(canonical_sorted), "-o", str(qc_root / "canonical"), "--mode", "daf", "--enzyme", str(cfg["fiberhmm"]["enzyme"]), "--reference-profile", str(cfg["fiberhmm"]["reference_profile"]), "--reference", str(reference), "--sample-reads", str(runtime.get("qc_sample_reads", 2000)), "--seed", str(runtime.get("qc_seed", 20260824)), "--min-mapq", str(mapq), "--prob-threshold", str(threshold)]
    command_lines.append(" ".join(qc_common))
    run_command(qc_common, stage_root / "logs/fiberhmm_qc_canonical.log", env=env)
    qc_chrm = [str(qc_bin), "-i", str(chrm_sorted), "-o", str(qc_root / "chrM"), "--mode", "daf", "--enzyme", str(cfg["fiberhmm"]["enzyme"]), "--reference-profile", str(cfg["fiberhmm"]["reference_profile"]), "--reference", str(reference), "--sample-reads", str(runtime.get("qc_sample_reads", 2000)), "--seed", str(runtime.get("qc_seed", 20260824)), "--min-mapq", str(mapq), "--prob-threshold", str(threshold)]
    command_lines.append(" ".join(qc_chrm))
    run_command(qc_chrm, stage_root / "logs/fiberhmm_qc_chrM.log", env=env)

    qc_json = next((qc_root / "canonical").glob("*.qc.json"), None)
    if qc_json is None:
        raise FileNotFoundError("fiberhmm-qc did not produce canonical *.qc.json")
    builder = args.project_root / "pipeline" / "scripts" / "build_fiberhmm_delivery.py"
    build_command = [str(runtime.get("python", "python3")), str(builder), "--run-root", str(stage_root), "--raw-bam", str(input_bam), "--canonical-bam", str(canonical_sorted), "--gtf", str(cfg["reference"]["gtf"]), "--qc-json", str(qc_json), "--bed-dir", str(canonical_dir / "features"), "--report-name", "DAFseq_fiberhmm_real_report.html"]
    command_lines.append(" ".join(build_command))
    run_command(build_command, stage_root / "logs/build_delivery.log", env=env)

    finalize = args.project_root / "pipeline" / "scripts" / "finalize_chrm_summary.py"
    run_command([str(runtime.get("python", "python3")), str(finalize), str(stage_root)], stage_root / "logs/finalize_chrm_summary.log", env=env)
    execution_log = stage_root / "FIBERHMM_EXECUTION_LOG.md"
    execution_log.write_text("# FiberHMM 实际执行记录\n\n## 工具与参数\n\n- FiberHMM version: " + str(cfg["fiberhmm"]["version"]) + "\n- mode: DddA DAF mode\n- enzyme: " + str(cfg["fiberhmm"]["enzyme"]) + "\n- reference profile: " + str(cfg["fiberhmm"]["reference_profile"]) + "\n- probability threshold: " + str(threshold) + "\n- primary reads: enabled\n- MAPQ threshold: " + str(mapq) + "\n- canonical chromosomes: " + ",".join(canonical) + "\n- chrM: separate run\n- NUMBA_DISABLE_JIT: " + env["NUMBA_DISABLE_JIT"] + " (runtime compatibility setting; FiberHMM caller remains the executed tool)\n\n## 实际命令\n\n" + "\n".join("```text\n" + command + "\n```" for command in command_lines) + "\n", encoding="utf-8")
    validation = args.project_root / "pipeline" / "scripts" / "validate_fiberhmm_delivery.py"
    run_command([str(runtime.get("python", "python3")), str(validation), str(stage_root)], stage_root / "logs/validate_fiberhmm_delivery.log", env=env)
    write_marker(args.output, cfg, "fiberhmm", Path(args.tmp), input_path=args.input, extra={"fiberhmm_version": cfg["fiberhmm"]["version"], "canonical_bam": canonical_sorted, "chrM_bam": chrm_sorted, "official_qc": qc_json, "computation_status": "COMPLETED"})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
