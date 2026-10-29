#!/usr/bin/env python3
"""Shared, side-effect-limited helpers for the real Qitan rules.

The helpers in this module do not orchestrate stages.  They only provide the
common config, output-root, marker, and command-execution conventions used by
the individual Snakemake rule scripts.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
from pathlib import Path
from typing import Any, Iterable

import yaml


def load_config(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def run_root(config: dict[str, Any], project_root: Path) -> Path:
    return project_root / "outputs" / str(config.get("run_id", "run001"))


def log_root(config: dict[str, Any], project_root: Path) -> Path:
    return project_root / "log" / str(config.get("run_id", "run001"))


def tmp_root(config: dict[str, Any], project_root: Path) -> Path:
    return project_root / "tmp" / str(config.get("run_id", "run001"))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_marker(path: Path, config: dict[str, Any], stage: str, tmp: Path, *, input_path: Path | None = None, extra: dict[str, Any] | None = None) -> None:
    rows = [
        {"field": "stage", "value": stage},
        {"field": "run_id", "value": str(config.get("run_id", ""))},
        {"field": "profile", "value": str(config.get("profile", "unknown"))},
        {"field": "status", "value": "COMPLETED"},
        {"field": "execution_mode", "value": str(config.get("execution", {}).get("mode", "unknown"))},
        {"field": "input_marker", "value": str(input_path or "")},
        {"field": "tmp_dir", "value": str(tmp)},
        {"field": "scientific_analysis", "value": "real execution"},
    ]
    for key, value in (extra or {}).items():
        rows.append({"field": str(key), "value": str(value)})
    write_tsv(path, rows, ["field", "value"])


def run_command(command: list[str], log_path: Path, *, env: dict[str, str] | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    with log_path.open("w", encoding="utf-8") as handle:
        handle.write("$ " + " ".join(command) + "\n")
        handle.flush()
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, env=merged_env, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"command failed with exit code {completed.returncode}: {' '.join(command)}; see {log_path}")


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
