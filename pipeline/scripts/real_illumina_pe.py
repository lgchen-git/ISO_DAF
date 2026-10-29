#!/usr/bin/env python3
"""Run the exploratory paired-end Illumina DAF workflow for one sample."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import shlex
import subprocess
from pathlib import Path

import pysam
import yaml


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate_pairs(r1: Path, r2: Path) -> dict[str, int]:
    count = 0
    total_bases = [0, 0]
    with pysam.FastxFile(str(r1)) as left, pysam.FastxFile(str(r2)) as right:
        for pair in itertools.zip_longest(left, right):
            rec1, rec2 = pair
            if rec1 is None or rec2 is None:
                raise ValueError(f"R1/R2 read counts differ at pair {count + 1}")
            name1 = rec1.name.removesuffix("/1").removesuffix("/2")
            name2 = rec2.name.removesuffix("/1").removesuffix("/2")
            if name1 != name2:
                raise ValueError(f"mate names differ at pair {count + 1}: {name1} != {name2}")
            count += 1
            total_bases[0] += len(rec1.sequence)
            total_bases[1] += len(rec2.sequence)
            if count % 1_000_000 == 0:
                print(f"paired FASTQ validation: {count:,} read pairs checked", flush=True)
    if count == 0:
        raise ValueError("FASTQ inputs contain no read pairs")
    return {
        "read_pairs": count,
        "r1_bases": total_bases[0],
        "r2_bases": total_bases[1],
        "r1_mean_length": round(total_bases[0] / count, 4),
        "r2_mean_length": round(total_bases[1] / count, 4),
    }


def run(command: list[str], log: Path, *, shell: bool = False) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    shown = command if shell else shlex.join(command)
    with log.open("w", encoding="utf-8") as handle:
        handle.write(f"$ {shown}\n")
        handle.flush()
        result = subprocess.run(" ".join(command) if shell else command, shell=shell, executable="/bin/bash" if shell else None, stdout=handle, stderr=subprocess.STDOUT, check=False)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}); see {log}")


def marker(path: Path, stage: str, sample: str, status: str = "COMPLETED", **extra: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["field", "value"])
        for key, value in {"stage": stage, "sample_id": sample, "status": status, **extra}.items():
            writer.writerow([key, value])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--tmp", required=True, type=Path)
    args = parser.parse_args()
    cfg = yaml.safe_load(args.configfile.read_text(encoding="utf-8"))
    project = args.project_root.resolve()
    run_id = cfg["run_id"]
    results = project / cfg["paths"].get("results_root", "results") / run_id
    logs = project / cfg["paths"].get("logs_root", "logs") / run_id
    tmp = args.tmp.resolve()
    sample_rows = read_samples(project / cfg["sample_sheet"])
    if len(sample_rows) != 1:
        raise ValueError(f"expected exactly one sample row, found {len(sample_rows)}")
    row = sample_rows[0]
    sample = row["sample_id"]
    r1, r2 = Path(row["input_1"]), Path(row["input_2"])
    fasta = Path(cfg["reference"]["fasta"])
    for path in (r1, r2, fasta):
        if not path.is_file():
            raise FileNotFoundError(path)

    for dirname in ("01_preflight", "02_raw_qc", "03_alignment_qc", "04_canonical_filter", "05_daf_signal", "06_fiberhmm", "07_delivery", "08_report"):
        (results / dirname).mkdir(parents=True, exist_ok=True)
    (tmp / "reference_index").mkdir(parents=True, exist_ok=True)
    marker(results / "01_preflight/stage_complete.tsv", "preflight", sample, profile=cfg["profile"], reference=str(fasta), input_1=str(r1), input_2=str(r2))

    qc = validate_pairs(r1, r2)
    qc_path = results / "02_raw_qc/paired_fastq_qc.tsv"
    with qc_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(qc.items())
    marker(results / "02_raw_qc/stage_complete.tsv", "raw_qc", sample, **qc)

    runtime = cfg["runtime"]
    samtools = runtime["samtools"]
    aligner = str(runtime.get("aligner", "bwa"))
    align_bam = results / "03_alignment_qc" / f"{sample}.sorted.bam"
    if aligner == "minimap2":
        index = Path(str(fasta) + ".mmi")
        if not index.is_file():
            raise FileNotFoundError(f"configured minimap2 index not found: {index}")
        mapper = [runtime["minimap2"], "-ax", "sr", "--MD", "-t", str(runtime.get("cores", 1)), str(index), str(r1), str(r2)]
    elif aligner == "bwa":
        bwa = runtime["bwa"]
        prefix = tmp / "reference_index/genome"
        index_suffixes = (".amb", ".ann", ".bwt", ".pac", ".sa")
        if not all(Path(str(prefix) + suffix).is_file() for suffix in index_suffixes):
            run([bwa, "index", "-p", str(prefix), str(fasta)], logs / "bwa_index.log")
        mapper = [bwa, "mem", "-t", str(runtime.get("cores", 1)), str(prefix), str(r1), str(r2)]
    else:
        raise ValueError(f"unsupported aligner: {aligner}")
    align_cmd = f"set -o pipefail; {shlex.join(mapper)} | {shlex.join([samtools, 'sort', '-@', str(runtime.get('cores', 1)), '-o', str(align_bam), '-'])}"
    run([align_cmd], logs / "03_alignment.log", shell=True)
    run([samtools, "index", str(align_bam)], logs / "03_index.log")
    stats = subprocess.run([samtools, "flagstat", str(align_bam)], text=True, capture_output=True, check=True).stdout
    (results / "03_alignment_qc/alignment_flagstat.txt").write_text(stats, encoding="utf-8")
    marker(results / "03_alignment_qc/stage_complete.tsv", "alignment_qc", sample, bam=str(align_bam), aligner=aligner, preset="sr", prebuilt_index_k=15 if aligner == "minimap2" else "", prebuilt_index_w=10 if aligner == "minimap2" else "", flagstat=str(results / "03_alignment_qc/alignment_flagstat.txt"))

    canonical = set(cfg["chromosomes"]["canonical"])
    min_mapq = int(cfg["thresholds"]["mapq"])
    filtered_bam = results / "04_canonical_filter" / f"{sample}.canonical.mapq{min_mapq}.bam"
    with pysam.AlignmentFile(str(align_bam), "rb") as source, pysam.AlignmentFile(str(filtered_bam), "wb", template=source) as dest:
        for record in source.fetch(until_eof=True):
            if record.is_unmapped or record.is_secondary or record.is_supplementary or record.mapping_quality < min_mapq:
                continue
            if source.get_reference_name(record.reference_id) in canonical:
                dest.write(record)
    run([samtools, "index", str(filtered_bam)], logs / "04_index.log")
    marker(results / "04_canonical_filter/stage_complete.tsv", "canonical_filter", sample, bam=str(filtered_bam), min_mapq=min_mapq, reference_scope="canonical nuclear chromosomes", duplicate_removal="false")

    legacy = project / "pipeline/scripts/legacy"
    daf_dir = results / "05_daf_signal"
    daf_cmd = [str(project / "bin/python"), str(legacy / "daf_stats.py"), "--bam", str(filtered_bam), "--reference", str(fasta), "--sample-id", sample,
               "--per-site-out", str(daf_dir / "per_site.tsv"), "--summary-out", str(daf_dir / "summary.tsv"), "--mismatch-out", str(daf_dir / "mismatch_spectrum.tsv"),
               "--attrition-out", str(daf_dir / "read_attrition.tsv"), "--denominator-out", str(daf_dir / "denominator_audit.tsv"), "--min-mapq", str(min_mapq),
               "--min-base-quality", str(cfg["thresholds"]["base_quality"]), "--conversion-types", "CT,GA", "--strand-cutoff", str(cfg["thresholds"]["strand_cutoff"])]
    run(daf_cmd, logs / "05_daf_signal.log")
    marker(daf_dir / "stage_complete.tsv", "daf_signal", sample, evidence_unit="primary reads; mates not collapsed", method="manuscript_approx_md", strand_cutoff=cfg["thresholds"]["strand_cutoff"])

    marker(results / "06_fiberhmm/stage_complete.tsv", "fiberhmm", sample, status="NOT_APPLICABLE", reason="Illumina short reads are outside this workflow's FiberHMM long-read scope")

    peak_dir = results / "07_delivery/peaks"
    peak_dir.mkdir(parents=True, exist_ok=True)
    run([runtime["macs3"], "callpeak", "-t", str(filtered_bam), "-f", "BAMPE", "-g", "hs", "-n", sample, "--outdir", str(peak_dir), "-q", str(cfg["thresholds"]["peak_qvalue"])], logs / "07_macs3.log")
    bigwig = results / "07_delivery" / f"{sample}.CPM.bw"
    run([runtime["bamcoverage"], "--bam", str(filtered_bam), "--outFileName", str(bigwig), "--normalizeUsing", "CPM", "--numberOfProcessors", str(runtime.get("cores", 1))], logs / "07_bamcoverage.log")
    audit = results / "07_delivery/DELIVERY_VALIDATION.tsv"
    required = [filtered_bam, Path(str(filtered_bam) + ".bai"), bigwig, daf_dir / "summary.tsv"]
    with audit.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["artifact", "exists", "bytes"])
        for artifact in required:
            writer.writerow([str(artifact.relative_to(results)), artifact.is_file(), artifact.stat().st_size if artifact.is_file() else 0])
    if not all(path.is_file() for path in required):
        raise RuntimeError("delivery validation failed")
    marker(results / "07_delivery/stage_complete.tsv", "delivery", sample, status="COMPLETED", peak_scope="exploratory; no input control", peak_duplicate_filter="MACS3 default", bam_duplicate_removal="false", bigwig=str(bigwig))

    if cfg.get("execution", {}).get("tss", False):
        run([str(project / "bin/python"), str(project / "pipeline/scripts/real_illumina_pe_tss.py"),
             "--configfile", str(args.configfile.resolve()), "--project-root", str(project)], logs / "07_tss.log")

    report = results / "08_report" / f"{run_id}_DAFseq_analysis_report.html"
    report_python = cfg.get("report", {}).get("python", str(project / "bin/python"))
    run([report_python, str(project / "pipeline/scripts/render_illumina_pe_report.py"),
         "--configfile", str(args.configfile.resolve()), "--project-root", str(project)],
        logs / "08_report.log")
    print(json.dumps({"sample_id": sample, "read_pairs": qc["read_pairs"], "report": str(report)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
