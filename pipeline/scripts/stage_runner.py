#!/usr/bin/env python3
"""Execute one standard pipeline-contract stage.

This first migration pass intentionally does not re-run external Qitan data.
It records a reproducible stage contract and validates the configured input
and reference locations. Profile-specific scientific rules can replace the
contract body one stage at a time without changing the Snakemake interface.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import yaml


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tmp", required=True, type=Path)
    parser.add_argument("--input", type=Path)
    return parser.parse_args()


def get_nested(data: dict[str, Any], *keys: str) -> Any:
    value: Any = data
    for key in keys:
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    return value


def external_checks(config: dict[str, Any]) -> list[dict[str, str]]:
    if str(config.get("profile", "unknown")) != "qitan_v2_11":
        return []
    rows = []
    for label, key in (
        ("source_v2_10_root", "source_v2_10_root"),
        ("input_bam", "input_bam"),
        ("fiberhmm_bam", "fiberhmm_bam"),
        ("reference_fasta", "reference.fasta"),
        ("reference_gtf", "reference.gtf"),
    ):
        value = config.get(key) if "." not in key else get_nested(config, *key.split("."))
        path = Path(str(value)) if value else Path("")
        rows.append(
            {
                "resource": label,
                "path": str(path),
                "exists": str(path.exists()).lower() if value else "false",
                "immutable": "true",
            }
        )
    return rows


def write_stage(path: Path, config: dict[str, Any], stage: str, tmp: Path, input_path: Path | None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp.mkdir(parents=True, exist_ok=True)
    profile = str(config.get("profile", "unknown"))
    status = "COMPLETED" if profile == "demo" else "CONFIGURED_NOT_EXECUTED"
    rows = [
        {"field": "stage", "value": stage},
        {"field": "run_id", "value": str(config.get("run_id", ""))},
        {"field": "profile", "value": profile},
        {"field": "status", "value": status},
        {"field": "execution_mode", "value": str(get_nested(config, "execution", "mode") or "unknown")},
        {"field": "input_marker", "value": str(input_path or "")},
        {"field": "tmp_dir", "value": str(tmp)},
        {"field": "scientific_analysis", "value": "not run by layout migration"},
    ]
    checks = external_checks(config)
    if checks:
        rows.append({"field": "external_resource_checks", "value": str(len(checks))})
        rows.append({"field": "external_resource_policy", "value": "immutable; checked only when stage executes"})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["field", "value"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if checks:
        check_path = path.with_name("external_resource_checks.tsv")
        with check_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(checks[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(checks)


def main() -> int:
    cli = get_args()
    config = yaml.safe_load(cli.configfile.read_text(encoding="utf-8")) or {}
    write_stage(cli.output, config, cli.stage, cli.tmp, cli.input)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
