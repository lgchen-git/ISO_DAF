#!/usr/bin/env python3
"""Build protein-coding TSS sites and compute active-fragment TSS profiles."""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import subprocess
import sys
from pathlib import Path

import yaml


def attributes(text: str) -> dict[str, str]:
    return {key: value for key, value in re.findall(r'([\w.]+) "([^"]*)"', text)}


def build_tss(gtf: Path, output: Path, canonical: set[str]) -> int:
    genes: dict[str, tuple[str, int, int, str, str]] = {}
    opener = gzip.open if gtf.suffix == ".gz" else open
    with opener(gtf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene" or fields[0] not in canonical:
                continue
            attrs = attributes(fields[8])
            if attrs.get("gene_type", attrs.get("gene_biotype", "")) != "protein_coding":
                continue
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue
            start, end = int(fields[3]) - 1, int(fields[4])
            strand = fields[6]
            if strand not in {"+", "-"}:
                continue
            genes[gene_id] = (fields[0], start, end, strand, attrs.get("gene_name", gene_id))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        for chrom, start, end, strand, name in sorted(genes.values()):
            pos = start if strand == "+" else end - 1
            handle.write(f"{chrom}\t{pos}\t{pos + 1}\t{name}\t0\t{strand}\n")
    return len(genes)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    args = parser.parse_args()
    project = args.project_root.resolve()
    cfg = yaml.safe_load(args.configfile.read_text(encoding="utf-8"))
    root = project / cfg["paths"].get("results_root", "results") / cfg["run_id"]
    sample_rows = list(csv.DictReader((project / cfg["sample_sheet"]).open(encoding="utf-8"), delimiter="\t"))
    if len(sample_rows) != 1:
        raise ValueError(f"expected one sample row, found {len(sample_rows)}")
    sample = sample_rows[0]["sample_id"]
    out = root / "07_delivery" / "tss"
    out.mkdir(parents=True, exist_ok=True)
    cfg_tss = cfg.get("tss", {})
    tss_bed = out / "protein_coding_gene_tss.bed"
    tss_count = build_tss(Path(cfg["reference"]["gtf"]), tss_bed, set(cfg["chromosomes"]["canonical"]))
    bam = root / "04_canonical_filter" / f"{sample}.canonical.mapq{cfg['thresholds']['mapq']}.bam"
    peak = root / "07_delivery" / "peaks" / f"{sample}_peaks.narrowPeak"
    active_regions = out / "peak_windows_2kb.bed"
    padding = int(cfg_tss.get("active_peak_padding_bp", cfg_tss.get("window_bp", 2000)))
    with peak.open(encoding="utf-8") as source, active_regions.open("w", encoding="utf-8") as target:
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                start, end = int(fields[1]), int(fields[2])
                summit = start + int(fields[9]) if len(fields) >= 10 and int(fields[9]) >= 0 else (start + end) // 2
                target.write(f"{fields[0]}\t{max(0, start-padding)}\t{end+padding}\t{fields[3] if len(fields)>3 else '.'}\t0\t{summit}\n")
    legacy = project / "pipeline/scripts/legacy/tss_enrichment.py"
    paths = {
        "active": out / "active_tss.bed",
        "all_fragments": out / "all_fragments.bed.gz",
        "daf_fragments": out / "daf_positive_fragments.bed.gz",
        "all_table": out / "all_fragments_profile.tsv",
        "daf_table": out / "daf_positive_profile.tsv",
        "all_plot": out / "all_fragments_profile.png",
        "daf_plot": out / "daf_positive_profile.png",
        "summary": out / "TSS_SUMMARY.tsv",
    }
    command = [sys.executable, str(legacy), "--bam", str(bam), "--sample-id", sample, "--tss-bed", str(tss_bed),
               "--active-region-bed", str(active_regions), "--active-region-label", f"MACS3_peaks_plusminus_{padding}bp",
               "--active-tss-out", str(paths["active"]), "--all-fragments-out", str(paths["all_fragments"]),
               "--daf-positive-fragments-out", str(paths["daf_fragments"]), "--all-table-out", str(paths["all_table"]),
               "--daf-positive-table-out", str(paths["daf_table"]), "--all-plot-out", str(paths["all_plot"]),
               "--daf-positive-plot-out", str(paths["daf_plot"]), "--summary-out", str(paths["summary"]),
               "--window-bp", str(cfg_tss.get("window_bp", 2000)), "--flank-bp", str(cfg_tss.get("flank_bp", 100)),
               "--fragment-min-len", str(cfg_tss.get("fragment_min_len", 0)), "--fragment-max-len", str(cfg_tss.get("fragment_max_len", 1000)),
               "--profile-mode", str(cfg_tss.get("profile_mode", "coverage")), "--active-tss-selection", str(cfg_tss.get("active_tss_selection", "all")),
               "--min-mapq", str(cfg["thresholds"]["mapq"]), "--min-base-quality", str(cfg["thresholds"]["base_quality"]),
               "--strand-cutoff", str(cfg["thresholds"]["strand_cutoff"]), "--conversion-types", "CT,GA"]
    if cfg_tss.get("require_proper_pair", True):
        command.append("--require-proper-pair")
    subprocess.run(command, check=True)
    with paths["summary"].open(encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    marker = out / "stage_complete.tsv"
    with marker.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["field", "value"])
        writer.writerows([["stage", "tss"], ["sample_id", sample], ["status", "COMPLETED"], ["protein_coding_gene_tss", tss_count], ["active_tss_count", row["active_tss_count"]], ["active_peak_padding_bp", padding], ["summary", paths["summary"]]])
    print(f"TSS complete: protein-coding TSS={tss_count}; active TSS={row['active_tss_count']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
