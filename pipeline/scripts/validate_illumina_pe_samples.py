#!/usr/bin/env python3
"""Validate the two-column mate paths required by the Illumina PE profile."""

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--skip-path-checks", action="store_true")
    args = parser.parse_args()
    with args.samples.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "group", "input_1", "input_2"}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit(f"sample sheet must contain: {', '.join(sorted(required))}")
        rows = list(reader)
    if not rows or len({row["sample_id"] for row in rows}) != len(rows):
        raise SystemExit("sample sheet must contain rows with unique sample_id values")
    for row in rows:
        for key in ("input_1", "input_2"):
            if not row[key]:
                raise SystemExit(f"{row['sample_id']}: {key} is empty")
            if not args.skip_path_checks and not Path(row[key]).is_file():
                raise SystemExit(f"{row['sample_id']}: FASTQ not found: {row[key]}")
    print(f"Validated {len(rows)} paired-end sample row(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
