#!/usr/bin/env python3
"""Finalize and validate a completed FiberHMM run without rerunning the caller."""

from __future__ import annotations

import argparse
from pathlib import Path

from real_stage_common import load_config, run_command, run_root, write_marker


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tmp", required=True, type=Path)
    args = parser.parse_args()
    cfg = load_config(args.configfile)
    stage_root = run_root(cfg, args.project_root) / "06_fiberhmm"
    runtime = cfg["runtime"]
    log_dir = stage_root / "logs"
    command_logs = [
        "fiberhmm_call_canonical.log",
        "samtools_sort_canonical.log",
        "fiberhmm_call_chrM.log",
        "samtools_sort_chrM.log",
        "fiberhmm_extract_canonical.log",
        "fiberhmm_extract_chrM.log",
        "fiberhmm_qc_canonical.log",
        "fiberhmm_qc_chrM.log",
        "build_delivery.log",
        "finalize_chrm_summary.log",
    ]
    commands = []
    for name in command_logs:
        path = log_dir / name
        if path.exists():
            first = path.read_text(encoding="utf-8", errors="replace").splitlines()
            if first and first[0].startswith("$ "):
                commands.append(first[0][2:])
    execution_log = stage_root / "FIBERHMM_EXECUTION_LOG.md"
    canonical = list(cfg["chromosomes"]["canonical"])
    execution_log.write_text(
        "# FiberHMM 实际执行记录\n\n"
        "## 工具与参数\n\n"
        f"- FiberHMM version: {cfg['fiberhmm']['version']}\n"
        "- mode: DddA DAF mode\n"
        f"- enzyme: {cfg['fiberhmm']['enzyme']}\n"
        f"- reference profile: {cfg['fiberhmm']['reference_profile']}\n"
        f"- probability threshold: {cfg['fiberhmm']['probability_threshold']}\n"
        "- primary reads: enabled\n"
        f"- MAPQ threshold: {cfg['thresholds']['mapq']}\n"
        f"- canonical chromosomes: {','.join(canonical)}\n"
        f"- chrM: separate run\n"
        f"- NUMBA_DISABLE_JIT: {'1' if cfg.get('execution', {}).get('numba_disable_jit', True) else '0'}\n\n"
        "## 实际命令\n\n"
        + "\n".join(f"```text\n{command}\n```" for command in commands)
        + "\n",
        encoding="utf-8",
    )
    env = {
        "MPLCONFIGDIR": str(args.tmp / "matplotlib"),
        "TMPDIR": str(args.tmp),
    }
    (args.tmp / "matplotlib").mkdir(parents=True, exist_ok=True)
    validation = args.project_root / "pipeline" / "scripts" / "validate_fiberhmm_delivery.py"
    run_command(
        [str(runtime.get("python", "python3")), str(validation), str(stage_root)],
        stage_root / "logs/validate_fiberhmm_delivery.log",
        env=env,
    )
    write_marker(
        args.output,
        cfg,
        "fiberhmm",
        args.tmp,
        input_path=args.input,
        extra={
            "fiberhmm_version": cfg["fiberhmm"]["version"],
            "canonical_bam": stage_root / "fiberhmm/canonical/BC01.fiberhmm.sorted.bam",
            "chrM_bam": stage_root / "fiberhmm/chrM/BC01.fiberhmm.chrM.sorted.bam",
            "official_qc": stage_root / "qc/canonical/BC01.fiberhmm.sorted.qc.json",
            "computation_status": "COMPLETED",
            "finalized_without_rerun": "true",
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
