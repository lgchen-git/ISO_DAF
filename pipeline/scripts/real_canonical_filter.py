#!/usr/bin/env python3
"""Apply the configured canonical/chrM scope to the audited DAF tables."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from real_stage_common import load_config, run_root, write_marker


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--input", required=False, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tmp", required=True, type=Path)
    args = parser.parse_args()
    cfg = load_config(args.configfile)
    legacy_dir = args.project_root / "pipeline" / "scripts" / "legacy"
    sys.path.insert(0, str(legacy_dir))
    import qitan_daf_v2_7 as v27

    source, read_rows, event_rows, site_rows, fiber_rows = v27.read_source_tables(cfg)
    stage_root = run_root(cfg, args.project_root) / "04_canonical_filter"
    stage_root.mkdir(parents=True, exist_ok=True)
    v27.filter_main_tables(cfg, stage_root, read_rows, event_rows, site_rows, fiber_rows)
    write_marker(args.output, cfg, "canonical_filter", Path(args.tmp), extra={"source_root": source, "canonical_scope": ",".join(cfg["chromosomes"]["canonical"]), "mitochondrial_scope": cfg["chromosomes"]["mitochondrial"]})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
