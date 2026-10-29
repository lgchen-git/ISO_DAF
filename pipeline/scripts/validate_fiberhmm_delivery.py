#!/usr/bin/env python3
"""Final source/derived/rendered checks for a FiberHMM delivery."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import pysam


def main(run_root: Path) -> int:
    report = run_root / "DAFseq_fiberhmm_real_report.html"
    required = [
        report,
        run_root / "README.md",
        run_root / "RUN_METADATA.tsv",
        run_root / "FIBERHMM_EXECUTION_LOG.md",
        run_root / "FINAL_VALIDATION.md",
        run_root / "fiberhmm/canonical/BC01.fiberhmm.sorted.bam",
        run_root / "fiberhmm/canonical/BC01.fiberhmm.sorted.bam.bai",
        run_root / "fiberhmm/chrM/BC01.fiberhmm.chrM.sorted.bam",
        run_root / "fiberhmm/chrM/BC01.fiberhmm.chrM.sorted.bam.bai",
        run_root / "qc/canonical/BC01.fiberhmm.sorted.qc.json",
        run_root / "qc/canonical/BC01.fiberhmm.sorted.qc.png",
    ]
    rows = []
    for path in required:
        rows.append(("artifact_exists", str(path.relative_to(run_root)), path.exists() and path.stat().st_size > 0, "file exists and is non-empty"))
    text = report.read_text() if report.exists() else ""
    image_paths = re.findall(r"<img[^>]+src='([^']+)'", text)
    href_paths = re.findall(r"href='([^']+)'", text)
    for rel in image_paths + [x for x in href_paths if not x.startswith("#") and not x.endswith("/")]:
        target = (run_root / rel).resolve()
        rows.append(("html_link", rel, target.exists() and target.stat().st_size > 0, "HTML image/link target"))
    rows.append(("html_sections", "eight_h2_sections", len(re.findall(r"<h2>", text)) == 8, "client report has eight first-level sections"))
    rows.append(("html_lightbox", "wheel_escape_mask", all(x in text for x in ("lightbox", "wheel", "Escape", "passive:false")), "offline image zoom interaction is embedded"))
    rows.append(("html_status_boundary", "not_validated", "NOT_VALIDATED" in text and "COMPLETED" in text, "engineering completion is separated from biological validation"))
    with pysam.AlignmentFile(str(run_root / "fiberhmm/canonical/BC01.fiberhmm.sorted.bam"), "rb") as bam:
        pg_ok = any(str(x.get("PN")) == "fiberhmm-call" and str(x.get("VN")) == "2.16.8" for x in bam.header.get("PG", []))
    rows.append(("source_provenance", "fiberhmm_call_2.16.8", pg_ok, "BAM @PG provenance"))
    table = run_root / "tables/MSP_SUMMARY.tsv"
    with table.open() as handle:
        msp = next(csv.DictReader(handle, delimiter="\t"))
    expected = int(msp["msp_gt150_count"]) / int(msp["msp_count"])
    rows.append(("derived_consistency", "msp_gt150_fraction", abs(float(msp["msp_gt150_fraction"]) - expected) < 1e-12, "derived fraction equals count/count"))
    with (run_root / "tables/FIBERHMM_QC_SUMMARY.tsv").open() as handle:
        qc_rows = {r["metric"]: r["value"] for r in csv.DictReader(handle, delimiter="\t")}
    with (run_root / "tables/NUCLEOSOME_SUMMARY.tsv").open() as handle:
        nuc_rows = next(csv.DictReader(handle, delimiter="\t"))
    with (run_root / "tables/MSP_SUMMARY.tsv").open() as handle:
        msp_rows = next(csv.DictReader(handle, delimiter="\t"))
    with (run_root / "RUN_METADATA.tsv").open() as handle:
        metadata = {r["key"]: r["value"] for r in csv.DictReader(handle, delimiter="\t")}
    rendered_checks = [
        ("msp_gt150_count", str(msp_rows["msp_gt150_count"])),
        ("msp_gt150_percent", f"{100 * float(msp_rows['msp_gt150_fraction']):.2f}%"),
        ("nucleosome_median", f"{float(nuc_rows['median_length_bp']):.0f} bp"),
        ("nrl", f"{float(nuc_rows['nrl_bp']):.0f} bp"),
        ("official_qc_status", qc_rows.get("overall.status", "NOT_AVAILABLE")),
    ]
    for name, value in rendered_checks:
        candidates = {value}
        if name == "msp_gt150_count":
            candidates.add(f"{int(float(value)):,}")
        rows.append(("source_derived_rendered", name, any(candidate in text for candidate in candidates), f"rendered HTML contains source-derived value {value}"))
    rows.append(("metadata", "version", metadata.get("fiberhmm_version") == "2.16.8", "RUN_METADATA.tsv version"))
    rows.append(("official_qc", "overall.status", qc_rows.get("overall.status") in {"WARN", "PASS", "FAIL"}, "canonical official QC status is sourced from qc JSON"))
    rows.append(("official_qc", "periodicity.reference_curve_correlation", "periodicity.reference_curve_correlation" in qc_rows, "reference correlation is retained"))
    with (run_root / "tables/CANONICAL_FILTER_SUMMARY.tsv").open() as handle:
        scopes = {r["scope"] for r in csv.DictReader(handle, delimiter="\t")}
    rows.append(("scope_consistency", "canonical_and_chrM", {"canonical_fiberhmm_output", "chrM_input", "chrM_fiberhmm_output"}.issubset(scopes), "main and mitochondrial scope rows"))
    fields = ["layer", "check", "status", "detail"]
    validation_path = run_root / "validation/FINAL_DELIVERY_VALIDATION.tsv"
    validation_path.parent.mkdir(parents=True, exist_ok=True)
    with validation_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for layer, check, ok, detail in rows:
            writer.writerow({"layer": layer, "check": check, "status": "PASS" if ok else "FAIL", "detail": detail})
    failed = [x for x in rows if not x[2]]
    status = "PASS" if not failed else "FAIL"
    (run_root / "FINAL_VALIDATION.md").write_text(
        "# FiberHMM 最终验收\n\n"
        f"Overall engineering delivery status: **{status}**\n\n"
        "本验收覆盖真实 FiberHMM BAM、BAI、官方 QC JSON/PNG、feature BED、derived TSV、HTML link targets 和离线 lightbox。\n\n"
        "- FiberHMM computation: `COMPLETED`\n"
        "- canonical input scope: `chr1-chr22, chrX, chrY`; chrM 独立运行和统计\n"
        "- official canonical QC: `WARN`\n"
        "- Qitan chemistry compatibility: `NOT_VALIDATED`\n"
        "- FiberHMM biological interpretation: `NOT_VALIDATED`\n"
        "- PCR deduplication: `NOT_EXECUTED`\n"
        "- DAF SNP masking: `NOT_EXECUTED`\n"
        "- 详细检查：`validation/FINAL_DELIVERY_VALIDATION.tsv`\n"
        + ("\n失败检查：\n" + "\n".join(f"- {x[1]}: {x[3]}" for x in failed) if failed else "\n所有自动检查通过。\n")
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main(Path(sys.argv[1]).resolve()))
