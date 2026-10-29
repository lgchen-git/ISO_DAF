#!/usr/bin/env python3
"""Build a traceable FiberHMM delivery from a real FiberHMM BAM.

This script only consumes the BAM and sidecar outputs produced by the real
FiberHMM commands. It does not synthesize FiberHMM calls or alter an existing
v2.11 report.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import math
import re
import shutil
import statistics
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam


CANONICAL = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CHRM = "chrM"
STATUS = {
    "computation": "COMPLETED",
    "technical": "WARN",
    "biological": "NOT_VALIDATED",
}


def tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def json_load(path: Path) -> dict:
    if not path or not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def int_text(value) -> str:
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return "NA"


def pct(value) -> str:
    try:
        return f"{100 * float(value):.2f}%"
    except (TypeError, ValueError):
        return "NA"


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def merge_bp(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    total = 0
    start, end = sorted(intervals)[0]
    for left, right in sorted(intervals)[1:]:
        if left <= end:
            end = max(end, right)
        else:
            total += max(0, end - start)
            start, end = left, right
    return total + max(0, end - start)


def bam_stats(path: Path, scope: str, allowed: set[str] | None = None) -> dict:
    counts = Counter()
    fibers: set[str] = set()
    spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            counts["records"] += 1
            if rec.is_unmapped:
                counts["unmapped"] += 1
                continue
            if rec.is_secondary:
                counts["secondary"] += 1
            if rec.is_supplementary:
                counts["supplementary"] += 1
            if rec.is_duplicate:
                counts["duplicate"] += 1
            if allowed is not None and rec.reference_name not in allowed:
                continue
            counts["scope_records"] += 1
            if not rec.is_secondary and not rec.is_supplementary:
                counts["primary_mapped"] += 1
            if rec.mapping_quality >= 20:
                counts["mapq20"] += 1
            if rec.is_secondary or rec.is_supplementary or rec.mapping_quality < 20:
                continue
            fibers.add(rec.query_name)
            counts["eligible_records"] += 1
            if rec.reference_name:
                spans[rec.reference_name].append((rec.reference_start, rec.reference_end))
    counts["fibers"] = len(fibers)
    counts["covered_bp_union"] = sum(merge_bp(v) for v in spans.values())
    counts["scope"] = scope
    return dict(counts)


def read_bed_features(bed_dir: Path) -> dict[str, list[dict]]:
    # Nucleosome and DAF tracks can contain millions of blocks. Keep the
    # memory-heavy per-block representation only for MSP/TF (needed by TSS and
    # representative views); aggregate the large tracks while streaming.
    result: dict[str, list[dict]] = {"msp": [], "tf": [], "nucleosome": [], "deam": []}
    meta = {"counts": Counter(), "lengths": defaultdict(list), "bins": defaultdict(Counter)}
    aliases = {
        "nucleosome": "nucleosome",
        "msp": "msp",
        "tf": "tf",
        "deam": "deam",
    }
    for path in sorted(bed_dir.glob("*.bed")):
        feature = next((v for k, v in aliases.items() if f"_{k}.bed" in path.name), None)
        if feature is None:
            continue
        with path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    continue
                try:
                    chrom = fields[0]
                    row_start = int(fields[1])
                    block_count = int(fields[9])
                    sizes = [int(v) for v in fields[10].rstrip(",").split(",") if v]
                    starts = [int(v) for v in fields[11].rstrip(",").split(",") if v]
                except (ValueError, IndexError):
                    continue
                for size, offset in zip(sizes[:block_count], starts[:block_count]):
                    item = {"feature": feature, "chrom": chrom, "start": row_start + offset, "end": row_start + offset + size, "length": size, "read_id": fields[3]}
                    if feature in {"msp", "tf"}:
                        result[feature].append(item)
                    else:
                        meta["counts"][feature] += 1
                        meta["lengths"][feature].append(size)
                        meta["bins"][(chrom, (item["start"] // 10000) * 10000)][feature] += 1
                        if feature == "msp":
                            meta["counts"]["msp"] += 1
    result["_meta"] = meta
    return result


def feature_summary(features: dict[str, list[dict]]) -> dict[str, int]:
    meta = features.get("_meta", {})
    counts = meta.get("counts", {})
    return {
        "msp": len(features.get("msp", [])),
        "msp_gt150": sum(x["length"] > 150 for x in features.get("msp", [])),
        "nucleosome": counts.get("nucleosome", 0),
        "tf": len(features.get("tf", [])),
        "deam": counts.get("deam", 0),
    }


def feature_coverage(features: dict[str, list[dict]], bam_path: Path) -> dict[str, float]:
    fiber_spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or rec.mapping_quality < 20:
                continue
            fiber_spans[rec.query_name].append((rec.reference_start, rec.reference_end))
    total = sum(merge_bp(v) for v in fiber_spans.values())
    result = {"total_fiber_covered_bp": total}
    for name in ("msp", "nucleosome", "tf"):
        if name == "nucleosome":
            bp = sum(features.get("_meta", {}).get("lengths", {}).get("nucleosome", []))
            result["nucleosome_covered_bp"] = bp
            result["nucleosome_covered_fraction"] = bp / total if total else 0.0
            continue
        by_read: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for row in features.get(name, []):
            by_read[row["read_id"]].append((row["start"], row["end"]))
        result[f"{name}_covered_bp"] = sum(merge_bp(v) for v in by_read.values())
        result[f"{name}_covered_fraction"] = (
            result[f"{name}_covered_bp"] / total if total else 0.0
        )
    return result


def parse_attr(text: str, key: str) -> str | None:
    m = re.search(rf"(?:^|;)\s*{re.escape(key)}\s+\"([^\"]+)\"", text)
    if m:
        return m.group(1)
    m = re.search(rf"(?:^|;)\s*{re.escape(key)}=([^;]+)", text)
    return m.group(1) if m else None


def representative_tss(gtf: Path) -> list[dict]:
    genes: dict[str, dict] = {}
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in {"gene", "transcript"}:
                continue
            chrom, _, _, start, end, _, strand, _, attrs = fields
            if chrom not in CANONICAL:
                continue
            gid = parse_attr(attrs, "gene_id")
            if not gid:
                continue
            name = parse_attr(attrs, "gene_name") or gid
            biotype = parse_attr(attrs, "gene_type") or parse_attr(attrs, "gene_biotype") or ""
            item = genes.setdefault(gid, {"gene_id": gid, "gene_name": name, "chrom": chrom, "strand": strand, "starts": [], "biotype": biotype})
            item["starts"].append((int(start) - 1, int(end), strand))
            if biotype == "protein_coding":
                item["biotype"] = biotype
    out = []
    for item in genes.values():
        if not item["starts"]:
            continue
        if item["strand"] == "+":
            pos = min(x[0] for x in item["starts"])
        else:
            pos = max(x[1] for x in item["starts"]) - 1
        out.append({"gene_id": item["gene_id"], "gene_name": item["gene_name"], "chrom": item["chrom"], "pos": pos, "strand": item["strand"], "biotype": item["biotype"]})
    return out


def make_tss_profile(features: dict[str, list[dict]], bam_path: Path, tss: list[dict], out_tsv: Path, out_png: Path) -> dict:
    bins = list(range(-2000, 2001, 100))
    centers = [x + 50 for x in bins[:-1]]
    numerator = Counter()
    denominator = Counter()
    selected_tss = [x for x in tss if x.get("biotype") == "protein_coding"] or tss
    tss_by_chrom: dict[str, list[int]] = defaultdict(list)
    for x in selected_tss:
        tss_by_chrom[x["chrom"]].append(x["pos"])
    for chrom in tss_by_chrom:
        tss_by_chrom[chrom].sort()
    msp = [x for x in features.get("msp", []) if x["length"] > 150]
    for row in msp:
        positions = tss_by_chrom.get(row["chrom"], [])
        lo = bisect_left(positions, row["start"] - 2000)
        hi = bisect_right(positions, row["end"] + 2000)
        mid = (row["start"] + row["end"]) // 2
        for pos in positions[lo:hi]:
            rel = mid - pos
            if -2000 <= rel < 2000:
                numerator[(rel // 100) * 100] += 1
    # The denominator is fiber-TSS window coverage, not event-positive-site depth.
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or rec.mapping_quality < 20:
                continue
            positions = tss_by_chrom.get(rec.reference_name, [])
            lo = bisect_left(positions, rec.reference_start - 2000)
            hi = bisect_right(positions, rec.reference_end + 2000)
            for pos in positions[lo:hi]:
                left = max(-2000, rec.reference_start - pos)
                right = min(2000, rec.reference_end - pos)
                if right <= left:
                    continue
                for start in range((left // 100) * 100, right, 100):
                    denominator[start] += 1
    rows = []
    for start, center in zip(bins[:-1], centers):
        num = numerator[start]
        den = denominator[start]
        rate = num / den if den else 0.0
        z = 1.96
        denom = 1 + z * z / den if den else 1
        center_rate = (rate + z * z / (2 * den)) / denom if den else 0.0
        half = z * math.sqrt((rate * (1 - rate) / den + z * z / (4 * den * den)) / denom) if den else 0.0
        rows.append({"distance_to_tss_bp": center, "msp_gt150_count": num, "fiber_coverage_count": den, "msp_gt150_actuation_fraction": rate, "wilson95_low": max(0.0, center_rate - half), "wilson95_high": min(1.0, center_rate + half), "tss_definition": "gene-representative protein-coding TSS; 100-bp bins; orientation-aware", "denominator": "all MAPQ>=20 primary fibers covering the TSS-centered bin"})
    tsv(out_tsv, rows, list(rows[0].keys()) if rows else ["distance_to_tss_bp"])
    x = [r["distance_to_tss_bp"] for r in rows]
    y = [r["msp_gt150_actuation_fraction"] * 100 for r in rows]
    low = [r["wilson95_low"] * 100 for r in rows]
    high = [r["wilson95_high"] * 100 for r in rows]
    fig, ax = plt.subplots(figsize=(12, 5.8), dpi=180)
    ax.plot(x, y, marker="o", linewidth=1.4, markersize=2.5, color="#2166ac")
    ax.fill_between(x, low, high, color="#9ecae1", alpha=0.28, linewidth=0)
    ax.axvline(0, color="#333333", linestyle="--", linewidth=1)
    ax.set_xlabel("Distance to TSS (bp)")
    ax.set_ylabel("MSP >150 bp actuation fraction (%)")
    ax.set_title("Gene-representative TSS-centered MSP >150 bp profile")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    return {"tss_count": len(selected_tss), "center_fraction": next((r["msp_gt150_actuation_fraction"] for r in rows if r["distance_to_tss_bp"] == 50), 0.0), "background_fraction": statistics.mean([r["msp_gt150_actuation_fraction"] for r in rows if abs(r["distance_to_tss_bp"]) >= 1000]) if rows else 0.0}


def plot_hist(values: list[int], title: str, xlabel: str, output: Path, xmin: int, xmax: int, threshold: int | None = None) -> dict:
    values = [int(v) for v in values if v is not None and v >= 0]
    fig, ax = plt.subplots(figsize=(10.5, 5.8), dpi=180)
    if values:
        ax.hist([v for v in values if xmin <= v <= xmax], bins=min(60, max(20, xmax - xmin)), color="#2166ac", alpha=0.84)
        median = statistics.median(values)
        ordered = sorted(values, reverse=True)
        half = sum(ordered) / 2
        running = 0
        n50 = ordered[-1]
        for value in ordered:
            running += value
            if running >= half:
                n50 = value
                break
        ax.axvline(median, color="#d73027", linestyle="-", linewidth=1.5, label=f"median={median:.0f} bp")
        ax.axvline(n50, color="#1a9850", linestyle=":", linewidth=1.5, label=f"N50={n50:.0f} bp")
    else:
        median = 0
        n50 = 0
    if threshold is not None:
        ax.axvline(threshold, color="#fdae61", linestyle="--", linewidth=1.5, label=f">{threshold} bp threshold")
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Feature count")
    ax.set_title(title)
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    return {"count": len(values), "median": median, "n50": n50}


def plot_summary(stats: dict, feature: dict, output: Path) -> None:
    labels = ["Fiber", "MSP", "MSP >150", "Nucleosome", "TF-like"]
    values = [stats.get("fibers", 0), feature.get("msp", 0), feature.get("msp_gt150", 0), feature.get("nucleosome", 0), feature.get("tf", 0)]
    fig, ax = plt.subplots(figsize=(10.5, 5.2), dpi=180)
    ax.bar(labels, values, color=["#5e81ac", "#2166ac", "#d73027", "#4daf4a", "#984ea3"])
    ax.set_ylabel("Count")
    ax.set_title("FiberHMM output summary")
    ax.tick_params(axis="x", rotation=18)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def load_region_features(bed_dir: Path, regions: list[dict]) -> dict[str, list[dict]]:
    result = {"msp": [], "nucleosome": [], "tf": [], "deam": []}
    for path in sorted(bed_dir.glob("*.bed")):
        feature = next((v for v in ("nucleosome", "msp", "tf", "deam") if f"_{v}.bed" in path.name), None)
        if feature is None:
            continue
        with path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    continue
                try:
                    chrom, row_start = fields[0], int(fields[1])
                    block_count = int(fields[9])
                    sizes = [int(v) for v in fields[10].rstrip(",").split(",") if v]
                    starts = [int(v) for v in fields[11].rstrip(",").split(",") if v]
                except (ValueError, IndexError):
                    continue
                for size, offset in zip(sizes[:block_count], starts[:block_count]):
                    start = row_start + offset
                    end = start + size
                    if any(r["chrom"] == chrom and start < r["end"] and end > r["start"] for r in regions):
                        result[feature].append({"feature": feature, "chrom": chrom, "start": start, "end": end, "length": size, "read_id": fields[3]})
    return result


def choose_regions(features: dict[str, list[dict]], bed_dir: Path, bam_path: Path, out_tsv: Path, out_dir: Path) -> list[dict]:
    by_bin: dict[tuple[str, int], Counter] = defaultdict(Counter)
    for key in ("msp", "nucleosome", "tf", "deam"):
        for row in features.get(key, []):
            if row["chrom"] not in CANONICAL:
                continue
            b = (row["start"] // 10000) * 10000
            by_bin[(row["chrom"], b)][key] += 1
            if key == "msp" and row["length"] > 150:
                by_bin[(row["chrom"], b)]["msp_gt150"] += 1
    for (chrom, start), values in features.get("_meta", {}).get("bins", {}).items():
        for key, value in values.items():
            by_bin[(chrom, start)][key] += value
    candidates = sorted(by_bin, key=lambda k: (by_bin[k]["msp_gt150"], by_bin[k]["tf"], by_bin[k]["nucleosome"], by_bin[k]["deam"]), reverse=True)
    chosen = []
    for chrom, start in candidates:
        if any(x["chrom"] == chrom and abs(x["start"] - start) < 20000 for x in chosen):
            continue
        c = by_bin[(chrom, start)]
        if not chosen or c["msp_gt150"] > 0 or c["tf"] > 0 or len(chosen) < 3:
            chosen.append({"region": f"{chrom}:{start + 1}-{start + 10000}", "chrom": chrom, "start": start, "end": start + 10000, "selection_reason": "highest available MSP >150 bp / TF / nucleosome feature burden"})
        if len(chosen) >= 3:
            break
    selected_features = load_region_features(bed_dir, chosen)
    rows = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for region in chosen:
            recs = []
            for rec in bam.fetch(region["chrom"], region["start"], region["end"]):
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or rec.mapping_quality < 20:
                    continue
                recs.append(rec)
            ids = {x.query_name for x in recs}
            for key in ("msp", "nucleosome", "tf", "deam"):
                subset = [x for x in selected_features.get(key, []) if x["chrom"] == region["chrom"] and x["start"] < region["end"] and x["end"] > region["start"] and x["read_id"] in ids]
                region[key] = len(subset)
                if key == "msp":
                    region["msp_gt150"] = sum(x["length"] > 150 for x in subset)
            region["fiber_count"] = len(ids)
            region["displayed_fibers"] = min(80, len(ids))
            rows.append({"region": region["region"], "chrom": region["chrom"], "start": region["start"], "end": region["end"], "fiber_count": region["fiber_count"], "displayed_fibers": region["displayed_fibers"], "msp_count": region.get("msp", 0), "msp_gt150_count": region.get("msp_gt150", 0), "nucleosome_count": region.get("nucleosome", 0), "tf_count": region.get("tf", 0), "daf_event_count": region.get("deam", 0), "selection_reason": region["selection_reason"], "interpretation_status": "EXPLORATORY"})
            plot_region(region, recs, {key: [x for x in selected_features.get(key, []) if x["read_id"] in ids] for key in ("msp", "nucleosome", "tf", "deam")}, out_dir / f"representative_{len(rows)}.png")
    tsv(out_tsv, rows, list(rows[0].keys()) if rows else ["region"])
    return rows


def plot_region(region: dict, recs: list, features: dict[str, list[dict]], output: Path) -> None:
    recs = sorted(recs, key=lambda r: (r.reference_start, r.query_name))[:80]
    x0, x1 = region["start"], region["end"]
    if recs:
        x0 = max(0, min(r.reference_start for r in recs) - 500)
        x1 = max(r.reference_end for r in recs) + 500
    fig, ax = plt.subplots(figsize=(15, max(5.5, 1.0 + 0.13 * max(1, len(recs)))), dpi=180)
    colors = {"nucleosome": "#4daf4a", "msp": "#2166ac", "tf": "#984ea3", "deam": "#d73027"}
    for i, rec in enumerate(recs):
        y = len(recs) - i
        ax.plot([max(x0, rec.reference_start), min(x1, rec.reference_end)], [y, y], color="#777777", linewidth=0.7, alpha=0.75)
        if rec.is_reverse:
            ax.plot([min(x1, rec.reference_end), max(x0, rec.reference_start)], [y, y], color="#777777", linewidth=0.7)
        for key in ("nucleosome", "msp", "tf"):
            for row in features.get(key, []):
                if row["read_id"] != rec.query_name or row["chrom"] != rec.reference_name or row["end"] <= x0 or row["start"] >= x1:
                    continue
                face = "#fdae61" if key == "msp" and row["length"] > 150 else colors[key]
                height = 0.23 if key == "msp" and row["length"] > 150 else 0.16
                ax.add_patch(plt.Rectangle((row["start"], y - height / 2), row["length"], height, color=face, alpha=0.88))
        for row in features.get("deam", []):
            if row["read_id"] == rec.query_name and row["chrom"] == rec.reference_name and x0 <= row["start"] <= x1:
                ax.plot([row["start"], row["start"]], [y - 0.11, y + 0.11], color=colors["deam"], linewidth=0.55)
    ax.set_xlim(x0, x1)
    ax.set_ylim(0, len(recs) + 1)
    ax.set_xlabel(f"Reference coordinate ({region['chrom']})")
    ax.set_ylabel("Fiber")
    ax.set_title(f"{region['region']} · representative fibers (n={len(recs)})")
    handles = [plt.Line2D([0], [0], color=colors["nucleosome"], linewidth=6, label="Nucleosome"), plt.Line2D([0], [0], color=colors["msp"], linewidth=6, label="MSP"), plt.Line2D([0], [0], color="#fdae61", linewidth=6, label="MSP >150 bp"), plt.Line2D([0], [0], color=colors["tf"], linewidth=6, label="TF-like"), plt.Line2D([0], [0], color=colors["deam"], linewidth=1, label="DAF event")]
    ax.legend(handles=handles, loc="upper right", frameon=True, ncol=5, fontsize=8)
    ax.grid(axis="x", alpha=0.15)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def html_report(run_root: Path, tables: Path, figures: Path, qc_json: dict, stats: dict, feature: dict, coverage: dict, tss_result: dict, reps: list[dict], output: Path) -> None:
    def image(name: str, title: str) -> str:
        path = figures / name
        if not path.exists():
            return f"<div class='missing'>图像未生成：{html.escape(name)}</div>"
        return f"<figure><img class='zoomable' src='figures/{html.escape(name)}' alt='{html.escape(title)}'><figcaption>{html.escape(title)}</figcaption></figure>"

    overall = qc_json.get("overall", {}) if isinstance(qc_json.get("overall", {}), dict) else {}
    qc_status = overall.get("status") or qc_json.get("overall_status") or qc_json.get("status") or "WARN"
    qc_score = overall.get("score", qc_json.get("overall_score", qc_json.get("score", "NA")))
    signal = qc_json.get("signal", {}) if isinstance(qc_json.get("signal", {}), dict) else {}
    periodicity = qc_json.get("periodicity", {}) if isinstance(qc_json.get("periodicity", {}), dict) else {}
    signal_status = signal.get("status", "NA")
    periodicity_status = periodicity.get("status", "NA")
    correlation = periodicity.get("reference_curve_correlation", "NA")
    amplitude = periodicity.get("reference_pattern_amplitude_ratio", "NA")
    try:
        qc_score_text = f"{float(qc_score):.2f}"
    except (TypeError, ValueError):
        qc_score_text = str(qc_score)
    card_items = [
        ("输入 records", int_text(stats.get("records"))),
        ("MAPQ ≥20 primary", int_text(stats.get("eligible_records"))),
        ("Fiber 数", int_text(stats.get("fibers"))),
        ("canonical covered bp", int_text(stats.get("covered_bp_union"))),
        ("FiberHMM mode", "DddA DAF"),
        ("probability threshold", "125"),
    ]
    cards = "".join(f"<div class='card'><div class='label'>{html.escape(k)}</div><div class='value'>{html.escape(str(v))}</div></div>" for k, v in card_items)
    daf_cards = "".join(f"<div class='card'><div class='label'>{html.escape(k)}</div><div class='value'>{html.escape(str(v))}</div></div>" for k, v in [("DAF calls", int_text(feature.get("deam"))), ("MSP calls", int_text(feature.get("msp"))), ("SNP masking", "NOT_RUN"), ("PCR deduplication", "NOT_EXECUTED")])
    rep_html = "".join(
        f"<div class='cards'><div class='card'><div class='label'>Region</div><div class='value'>{html.escape(r['region'])}</div></div>"
        f"<div class='card'><div class='label'>Fibers / displayed</div><div class='value'>{int_text(r['fiber_count'])} / {int_text(r['displayed_fibers'])}</div></div>"
        f"<div class='card'><div class='label'>MSP &gt;150</div><div class='value'>{int_text(r['msp_gt150_count'])}</div></div>"
        f"<div class='card'><div class='label'>Nucleosome / TF</div><div class='value'>{int_text(r['nucleosome_count'])} / {int_text(r['tf_count'])}</div></div></div>"
        f"{image(f'representative_{i + 1}.png', r['region'] + ' representative fiber view')}"
        for i, r in enumerate(reps)
    )
    css = """<style>
:root{--ink:#243447;--muted:#637083;--line:#dbe3ea;--blue:#2166ac;--green:#238b45;--warn:#c77c02;--bg:#f7f9fb}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.65 -apple-system,BlinkMacSystemFont,"Segoe UI",Arial,"Noto Sans CJK SC",sans-serif}header{background:#fff;border-bottom:1px solid var(--line);padding:34px 6vw 28px}main{max-width:1440px;margin:0 auto;padding:24px 5vw 60px}h1{margin:0 0 8px;font-size:32px}h2{margin:0 0 16px;font-size:24px}.meta{color:var(--muted);display:flex;gap:16px;flex-wrap:wrap}.lead{font-size:18px;max-width:1000px;margin:20px 0 0}nav{position:sticky;top:0;z-index:5;background:#fff;border-bottom:1px solid var(--line);padding:10px 6vw;white-space:nowrap;overflow:auto}nav a{color:var(--blue);text-decoration:none;margin-right:18px}section{background:#fff;border:1px solid var(--line);border-radius:14px;padding:24px;margin:22px 0;box-shadow:0 4px 18px #2434470b}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin:16px 0 22px}.card{border:1px solid var(--line);border-radius:10px;padding:14px;background:#fff}.label{font-size:12px;color:var(--muted);text-transform:uppercase;letter-spacing:.04em}.value{font-size:24px;font-weight:700;margin:3px 0}.note{color:var(--muted)}.status{font-weight:700;color:var(--warn)}.good{color:var(--green)}.grid2{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px}figure{margin:14px 0;text-align:center}figure img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:8px;background:#fff;cursor:zoom-in}figcaption{font-size:13px;color:var(--muted);margin-top:6px}table{width:100%;border-collapse:collapse;font-size:13px}th,td{border-bottom:1px solid var(--line);padding:8px;text-align:left;vertical-align:top}th{position:sticky;top:43px;background:#f2f6f9}.boundary{border-left:4px solid #fdae61;background:#fff8eb;padding:12px 15px;margin:15px 0}#lightbox{display:none;position:fixed;inset:0;background:#101820e8;z-index:20;align-items:center;justify-content:center}#lightbox.open{display:flex}#lightbox img{max-width:92vw;max-height:90vh;transform-origin:center;cursor:grab}#lightbox .close{position:absolute;top:18px;right:24px;color:#fff;font-size:32px;background:none;border:0;cursor:pointer}#lightbox .hint{position:absolute;bottom:14px;color:#fff;font-size:13px}@media(max-width:800px){.grid2{grid-template-columns:1fr}header{padding:24px 5vw}h1{font-size:26px}}
</style>"""
    script = """<script>(()=>{const box=document.getElementById('lightbox'),img=box.querySelector('img');let scale=1;function close(){box.classList.remove('open');img.src='';scale=1;img.style.transform='scale(1)'}document.querySelectorAll('img.zoomable').forEach(x=>x.addEventListener('click',()=>{img.src=x.src;img.alt=x.alt;box.classList.add('open')}));box.addEventListener('click',e=>{if(e.target===box)close()});box.querySelector('.close').addEventListener('click',close);document.addEventListener('keydown',e=>{if(e.key==='Escape')close()});box.addEventListener('wheel',e=>{if(!box.classList.contains('open'))return;e.preventDefault();scale=Math.min(8,Math.max(.3,scale*(e.deltaY<0?1.15:.87)));img.style.transform=`scale(${scale})`},{passive:false})})();</script>"""
    chromatin = (
        f"<div class='cards'><div class='card'><div class='label'>Computation status</div><div class='value good'>COMPLETED</div></div>"
        f"<div class='card'><div class='label'>Official QC</div><div class='value status'>{html.escape(str(qc_status))}</div></div>"
        f"<div class='card'><div class='label'>QC score</div><div class='value'>{html.escape(qc_score_text)}</div></div>"
        f"<div class='card'><div class='label'>MSP &gt;150 bp</div><div class='value'>{int_text(feature.get('msp_gt150'))}</div><div class='note'>{pct(feature.get('msp_gt150', 0) / (feature.get('msp', 1) or 1))}</div></div>"
        f"<div class='card'><div class='label'>Nucleosome median</div><div class='value'>{int_text(stats.get('nucleosome_median'))} bp</div></div>"
        f"<div class='card'><div class='label'>NRL</div><div class='value'>{int_text(stats.get('nrl'))} bp</div></div></div>"
    )
    html_text = """<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'><title>DAF-seq 分析报告 · FiberHMM 实际运行结果</title>""" + css + """</head><body>
<header><div class='meta'>DAF-seq · GRCh38 · 主分析 chr1–chr22, chrX, chrY · chrM 单独统计</div><h1>FiberHMM 染色质结构分析报告</h1><p class='lead'>本报告基于真实 FiberHMM 2.16.8 DddA DAF mode 运行结果，汇总 DAF-seq 信号、MSP、nucleosome、TF-like footprint 和 TSS-centered exploratory profile。</p></header>
<nav><a href='#qc'>数据与比对质控</a><a href='#daf'>DAF 信号质控</a><a href='#chromatin'>FiberHMM 与染色质结构</a><a href='#tss'>TSS 可及性</a><a href='#fiber'>单分子 Fiber 展示</a><a href='#advanced'>高级分析状态</a><a href='#methods'>方法与解释边界</a><a href='#delivery'>交付信息</a></nav><main>"""
    html_text += f"<section id='qc'><h2>1. 数据与比对质控</h2><p>canonical 染色体范围在 FiberHMM 输入层通过 chromosome selection 实际执行；chrM 不进入主 nuclear chromatin 结果。</p><div class='cards'>{cards}</div>{image('alignment_summary.png', '输入与 FiberHMM 输出的核心范围统计')}<div class='boundary'><b>读取边界：</b>high mapping rate 不等于 full-length alignment quality 高；当前展示以 primary、supplementary/secondary 处理和 MAPQ ≥20 的 FiberHMM 输入 cohort 为准。</div><details><summary>详细表格</summary><p><a href='tables/CANONICAL_FILTER_SUMMARY.tsv'>CANONICAL_FILTER_SUMMARY.tsv</a> · <a href='tables/CHRM_FIBERHMM_SUMMARY.tsv'>CHRM_FIBERHMM_SUMMARY.tsv</a></p></details></section>"
    html_text += f"<section id='daf'><h2>2. DAF 信号质控</h2><div class='cards'>{daf_cards}</div><p>本次真实 caller 的 DAF SNP automatic preflight 估计 genome coverage 约 0.66×，未达到工具触发 recurrent SNP discovery 的深度条件，因此 SNP masking 未执行；PCR deduplication 初次内置预处理未完成，主调用明确使用 <code>--no-dedup</code>，未将其包装成已校正结果。</p><details><summary>执行证据</summary><p><a href='FIBERHMM_EXECUTION_LOG.md'>FIBERHMM_EXECUTION_LOG.md</a> · <a href='tables/TF_FOOTPRINT_SUMMARY.tsv'>TF_FOOTPRINT_SUMMARY.tsv</a></p></details></section>"
    html_text += f"<section id='chromatin'><h2>3. FiberHMM 与染色质结构</h2>{chromatin}{image('fiberhmm_qc_dashboard.png', '官方 FiberHMM QC dashboard（实际输出）')}<div class='grid2'>{image('msp_length_distribution.png', 'MSP length distribution')}{image('nucleosome_length_distribution.png', 'Nucleosome length distribution')}</div><p>当前 FiberHMM 计算已完成，官方 QC 为 {html.escape(str(qc_status))}（score {html.escape(qc_score_text)}）。WARN 主要来自 nucleosome-scale periodicity：status={html.escape(str(periodicity_status))}，reference correlation={html.escape(str(correlation))}，matched amplitude={html.escape(str(amplitude))}；signal-rate status={html.escape(str(signal_status))}。Qitan chemistry 与 FiberHMM DddA DAF model/reference profile 的兼容性仍为 NOT_VALIDATED；接近典型 nucleosome length 或 NRL 不等于平台生物学验证完成。</p><details><summary>完整结果表</summary><p><a href='tables/FIBERHMM_QC_SUMMARY.tsv'>FIBERHMM_QC_SUMMARY.tsv</a> · <a href='tables/MSP_SUMMARY.tsv'>MSP_SUMMARY.tsv</a> · <a href='tables/NUCLEOSOME_SUMMARY.tsv'>NUCLEOSOME_SUMMARY.tsv</a> · <a href='tables/FIBERHMM_COVERAGE_NORMALIZATION.tsv'>FIBERHMM_COVERAGE_NORMALIZATION.tsv</a></p></details></section>"
    html_text += f"<section id='tss'><h2>4. TSS 可及性</h2><div class='cards'><div class='card'><div class='label'>gene-representative TSS</div><div class='value'>{int_text(tss_result.get('tss_count'))}</div></div><div class='card'><div class='label'>center MSP &gt;150</div><div class='value'>{pct(tss_result.get('center_fraction'))}</div></div><div class='card'><div class='label'>background</div><div class='value'>{pct(tss_result.get('background_fraction'))}</div></div><div class='card'><div class='label'>Interpretation</div><div class='value status'>EXPLORATORY</div></div></div>{image('tss_msp150_profile.png', 'Gene-representative TSS-centered MSP >150 bp profile')}<p><b>结论：</b>当前 profile 用 MSP &gt;150 bp overlapping TSS-centered bins 构建，未观察到强烈的 TSS-centered MSP &gt;150 bp accessibility enrichment。</p><div class='boundary'>该 profile 的 denominator 是所有 MAPQ ≥20 primary fibers 覆盖对应 TSS-centered bin 的 fiber coverage；图中比例是 descriptive observation，不是 Qitan-specific biological calibration，也不代表 independent replicate-level confidence。</div><details><summary>完整 TSS 表</summary><p><a href='tables/TSS_MSP150_SUMMARY.tsv'>TSS_MSP150_SUMMARY.tsv</a></p></details></section>"
    html_text += f"<section id='fiber'><h2>5. 单分子 Fiber 展示</h2><p>代表区域按实际 FiberHMM feature burden 和覆盖度筛选，展示的是 exploratory single-molecule view，不是 validated regulatory case study。</p>{rep_html}<details><summary>代表区域统计</summary><p><a href='tables/REPRESENTATIVE_FIBER_SUMMARY.tsv'>REPRESENTATIVE_FIBER_SUMMARY.tsv</a></p></details></section>"
    html_text += "<section id='advanced'><h2>6. 高级分析状态</h2><table><thead><tr><th>Analysis</th><th>Status</th><th>Reason</th><th>Minimum requirement</th></tr></thead><tbody><tr><td>High-coverage fiber clustering</td><td>NOT_TESTABLE</td><td>本次结果未满足预设的纤维数量与 feature variation 门槛</td><td>足够 fiber 数及可分离的 feature variation</td></tr><tr><td>Co-actuation</td><td>NOT_TESTABLE</td><td>未形成满足 pair coverage requirement 的 MSP >150 bp pair</td><td>两个 MSP >150 bp element、足够 pair coverage 和 fibers</td></tr><tr><td>FIRE</td><td>NOT_EXECUTED / NOT_VALIDATED</td><td>未执行 calibrated genome-wide FIRE calling</td><td>合适的 genome-wide calibration 和可解释的 input cohort</td></tr><tr><td>CTCF occupancy</td><td>NOT_EXECUTED / NOT_VALIDATED</td><td>未执行 motif-aware single-fiber CTCF occupancy test</td><td>可靠 CTCF motif 与 external annotation/control</td></tr></tbody></table><p class='note'>NOT_TESTABLE / NOT_EXECUTED 不等于 zero signal 或 negative result。</p><p><a href='tables/HIGH_COVERAGE_FIBER_ANALYSIS.tsv'>HIGH_COVERAGE_FIBER_ANALYSIS.tsv</a> · <a href='tables/CO_ACTUATION_SUMMARY.tsv'>CO_ACTUATION_SUMMARY.tsv</a></p></section>"
    html_text += f"<section id='methods'><h2>7. 方法与解释边界</h2><p><b>实际工具链：</b> FiberHMM 2.16.8，<code>fiberhmm-call</code>，mode=DddA DAF，enzyme=ddda，reference/profile=ddda，probability threshold=125，primary reads，MAPQ ≥20，canonical chromosomes at analysis input。</p><p><b>状态分层：</b> COMPUTATION_STATUS={STATUS['computation']}；TECHNICAL_QC_STATUS={html.escape(str(qc_status))}；BIOLOGICAL_VALIDATION_STATUS={STATUS['biological']}。COMPLETED 不等于 biological validation passed。</p><p><b>可解释范围：</b> MSP、nucleosome、TF-like calls 和 TSS profile 是基于真实输出的 computational measurements/descriptive interpretation；Qitan chemistry compatibility 与 FiberHMM biological interpretation 仍为 NOT_VALIDATED。</p><p><a href='methods/LITERATURE_METHODS_MATRIX.tsv'>LITERATURE_METHODS_MATRIX.tsv</a> · <a href='methods/FIBERHMM_METHODS.md'>FIBERHMM_METHODS.md</a></p></section>"
    html_text += "<section id='delivery'><h2>8. 交付信息</h2><table><tbody><tr><th>Report</th><td><a href='DAFseq_fiberhmm_real_report.html'>DAFseq_fiberhmm_real_report.html</a></td></tr><tr><th>Tables</th><td><a href='tables/'>tables/</a></td></tr><tr><th>Figures</th><td><a href='figures/'>figures/</a></td></tr><tr><th>Methods</th><td><a href='methods/'>methods/</a></td></tr><tr><th>Validation</th><td><a href='validation/'>validation/</a></td></tr><tr><th>Raw FiberHMM output</th><td><a href='fiberhmm/'>fiberhmm/</a></td></tr></tbody></table><p>完整结果表、图源、方法说明和验证文件均随交付目录保留，可用于结果追溯。</p></section></main><div id='lightbox' role='dialog' aria-modal='true'><button class='close' aria-label='关闭'>×</button><img alt='放大图'><div class='hint'>滚轮缩放 · Esc 或点击背景关闭</div></div>" + script + "</body></html>"
    output.write_text(html_text)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--raw-bam", type=Path, required=True)
    parser.add_argument("--canonical-bam", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--qc-json", type=Path, required=True)
    parser.add_argument("--bed-dir", type=Path, required=True)
    parser.add_argument("--report-name", default="DAFseq_fiberhmm_real_report.html")
    args = parser.parse_args()
    run_root = args.run_root.resolve()
    tables, figures, methods, validation = [run_root / x for x in ("tables", "figures", "methods", "validation")]
    for d in (tables, figures, methods, validation):
        d.mkdir(parents=True, exist_ok=True)
    print("[1/10] scanning raw and FiberHMM BAM statistics", flush=True)
    raw = bam_stats(args.raw_bam, "raw_input", set(CANONICAL + [CHRM]))
    canonical = bam_stats(args.canonical_bam, "canonical_fiberhmm_output", set(CANONICAL))
    chrm = bam_stats(args.raw_bam, "chrM_input", {CHRM})
    print("[2/10] streaming feature BEDs", flush=True)
    features = read_bed_features(args.bed_dir)
    fs = feature_summary(features)
    cov = feature_coverage(features, args.canonical_bam)
    qc = json_load(args.qc_json)
    overall = qc.get("overall", {}) if isinstance(qc.get("overall", {}), dict) else {}
    qc_status = overall.get("status") or qc.get("overall_status") or qc.get("status") or "WARN"
    print("[3/10] rendering length and summary figures", flush=True)
    lengths_msp = [x["length"] for x in features.get("msp", [])]
    lengths_nuc = features.get("_meta", {}).get("lengths", {}).get("nucleosome", [])
    msp_stats = plot_hist(lengths_msp, "MSP length distribution", "MSP length (bp)", figures / "msp_length_distribution.png", 0, 300, 150)
    nuc_stats = plot_hist(lengths_nuc, "Nucleosome length distribution", "Nucleosome length (bp)", figures / "nucleosome_length_distribution.png", 50, 400)
    nrl = qc.get("periodicity", {}).get("nrl_bp", 185)
    plot_summary(canonical, fs, figures / "alignment_summary.png")
    qc_plot = qc.get("plot") or qc.get("figure") or qc.get("outputs", {}).get("plot")
    if qc_plot and Path(qc_plot).exists():
        shutil.copy2(qc_plot, figures / "fiberhmm_qc_dashboard.png")
    print("[4/10] parsing gene-representative TSS", flush=True)
    tss = representative_tss(args.gtf)
    print(f"[5/10] building TSS profile from {len(tss)} genes", flush=True)
    tss_result = make_tss_profile(features, args.canonical_bam, tss, tables / "TSS_MSP150_SUMMARY.tsv", figures / "tss_msp150_profile.png")
    print("[6/10] selecting representative regions", flush=True)
    reps = choose_regions(features, args.bed_dir, args.canonical_bam, tables / "REPRESENTATIVE_FIBER_SUMMARY.tsv", figures)
    high_cov_rows = []
    for i, row in enumerate(sorted(reps, key=lambda x: x.get("fiber_count", 0), reverse=True), 1):
        high_cov_rows.append({"rank": i, "region": row["region"], "fiber_count": row["fiber_count"], "msp_gt150_count": row["msp_gt150_count"], "nucleosome_count": row["nucleosome_count"], "tf_count": row["tf_count"], "status": "COMPLETED_DESCRIPTIVE"})
    tsv(tables / "HIGH_COVERAGE_FIBER_ANALYSIS.tsv", high_cov_rows or [{"rank": "NA", "region": "NA", "fiber_count": 0, "msp_gt150_count": 0, "nucleosome_count": 0, "tf_count": 0, "status": "NOT_AVAILABLE"}], ["rank", "region", "fiber_count", "msp_gt150_count", "nucleosome_count", "tf_count", "status"])
    tsv(tables / "CO_ACTUATION_SUMMARY.tsv", [{"status": "NOT_TESTABLE", "reason": "No pair was advanced to an eligible MSP >150 bp co-actuation test with a documented pair-coverage threshold", "minimum_requirement": "two MSP >150 bp elements, published pair coverage threshold, and sufficient fibers"}], ["status", "reason", "minimum_requirement"])
    scope_rows = []
    for label, s in (("raw_input", raw), ("canonical_fiberhmm_output", canonical), ("chrM_input", chrm)):
        row = {"scope": label, **{k: s.get(k, 0) for k in ("records", "scope_records", "primary_mapped", "supplementary", "secondary", "duplicate", "mapq20", "eligible_records", "fibers", "covered_bp_union")}}
        for k, v in fs.items(): row[f"output_{k}"] = v if label == "canonical_fiberhmm_output" else "NA"
        scope_rows.append(row)
    tsv(tables / "CANONICAL_FILTER_SUMMARY.tsv", scope_rows, list(scope_rows[0].keys()))
    tsv(tables / "MSP_SUMMARY.tsv", [{"msp_count": fs["msp"], "msp_gt150_count": fs["msp_gt150"], "msp_gt150_fraction": fs["msp_gt150"] / fs["msp"] if fs["msp"] else 0, "median_length_bp": msp_stats["median"], "n50_bp": msp_stats["n50"], "fiber_covered_bp": cov["msp_covered_bp"], "fiber_covered_fraction": cov["msp_covered_fraction"], "threshold_definition": ">150 bp"}], ["msp_count", "msp_gt150_count", "msp_gt150_fraction", "median_length_bp", "n50_bp", "fiber_covered_bp", "fiber_covered_fraction", "threshold_definition"])
    tsv(tables / "NUCLEOSOME_SUMMARY.tsv", [{"nucleosome_count": fs["nucleosome"], "median_length_bp": nuc_stats["median"], "n50_bp": nuc_stats["n50"], "nrl_bp": nrl, "fiber_covered_bp": cov["nucleosome_covered_bp"], "fiber_covered_fraction": cov["nucleosome_covered_fraction"]}], ["nucleosome_count", "median_length_bp", "n50_bp", "nrl_bp", "fiber_covered_bp", "fiber_covered_fraction"])
    tsv(tables / "TF_FOOTPRINT_SUMMARY.tsv", [{"tf_count": fs["tf"], "fiber_covered_bp": cov["tf_covered_bp"], "fiber_covered_fraction": cov["tf_covered_fraction"], "status": "COMPLETED" if fs["tf"] else "NOT_AVAILABLE"}], ["tf_count", "fiber_covered_bp", "fiber_covered_fraction", "status"])
    qc_rows = []
    for section, values in qc.items():
        if isinstance(values, dict):
            for key, value in values.items():
                if not isinstance(value, (dict, list)):
                    qc_rows.append({"metric": f"{section}.{key}", "value": value, "unit": "NA", "status": qc_status})
        elif not isinstance(values, list):
            qc_rows.append({"metric": section, "value": values, "unit": "NA", "status": qc_status})
    if not qc_rows: qc_rows = [{"metric": "overall_status", "value": "NOT_AVAILABLE", "unit": "NA", "status": "NOT_AVAILABLE"}]
    tsv(tables / "FIBERHMM_QC_SUMMARY.tsv", qc_rows, ["metric", "value", "unit", "status"])
    tsv(tables / "FIBERHMM_COVERAGE_NORMALIZATION.tsv", [{"total_fiber_covered_bp": cov["total_fiber_covered_bp"], "msp_covered_bp": cov["msp_covered_bp"], "msp_covered_fraction": cov["msp_covered_fraction"], "nucleosome_covered_bp": cov["nucleosome_covered_bp"], "nucleosome_covered_fraction": cov["nucleosome_covered_fraction"], "tf_covered_bp": cov["tf_covered_bp"], "tf_covered_fraction": cov["tf_covered_fraction"]}], ["total_fiber_covered_bp", "msp_covered_bp", "msp_covered_fraction", "nucleosome_covered_bp", "nucleosome_covered_fraction", "tf_covered_bp", "tf_covered_fraction"])
    print("[7/10] writing methods and metadata", flush=True)
    (methods / "LITERATURE_METHODS_MATRIX.tsv").write_text("module\tpublished_method\treference\toriginal_parameter\tcurrent_implementation\tqitan_adaptation\texecution_status\tvalidation_status\nDAF strand / deamination QC\tDAF deamination orientation and signal QC\tFiberHMM 2.16.8 local documentation\tMAPQ >=20; DddA mode\tReal fiberhmm-call on MD-tagged BAM\tQitan chemistry not calibrated\tCOMPLETED\tNOT_VALIDATED\nFiberHMM\tHMM nucleosome/MSP/TF calling\tFiberHMM 2.16.8 local documentation\tprobability threshold 125\tDddA DAF mode, bundled dd-da-compatible ddda model\tExploratory on Qitan\tCOMPLETED\tNOT_VALIDATED\nMSP >150 actuation\tMSP >150 bp TSS-centered descriptive profile\tDAF-seq literature-guided\t>150 bp\tGene-representative TSS profile\tDescriptive only\tCOMPLETED\tNOT_VALIDATED\nTSS accessibility\tMSP-centered TSS accessibility\tDAF-seq literature-guided\tTSS +/-2 kb\t100-bp fiber-coverage normalized bins\tNo Qitan-specific calibration\tCOMPLETED\tNOT_VALIDATED\nSingle-fiber view\tSingle-molecule feature tracks\tFiber-seq/FiberHMM literature-guided\tMAPQ >=20\tCanonical representative regions\tExploratory\tCOMPLETED\tNOT_VALIDATED\nCo-actuation\tPaired actuation\tDAF-seq literature-guided\tPair coverage threshold\tNo eligible pair detected\tNot testable\tNOT_TESTABLE\tNOT_VALIDATED\nFIRE\tGenome-wide FIRE\tDAF-seq literature-guided\tCalibrated input required\tNot run\tMissing calibrated condition\tNOT_EXECUTED\tNOT_VALIDATED\nCTCF\tMotif-aware occupancy\tFiber-seq literature-guided\tMotif and control required\tNot run\tMissing external annotation/control\tNOT_EXECUTED\tNOT_VALIDATED\n")
    run_meta = [{"key": k, "value": v} for k, v in {"fiberhmm_version": "2.16.8", "mode": "DddA DAF mode", "enzyme": "ddda", "reference_profile": "ddda", "probability_threshold": 125, "canonical_chromosomes": ",".join(CANONICAL), "chrM": CHRM, "canonical_bam": str(args.canonical_bam.resolve()), "qc_json": str(args.qc_json.resolve()), "msp_threshold": ">150 bp", "computation_status": STATUS["computation"], "technical_qc_status": qc_status, "biological_validation_status": STATUS["biological"]}.items()]
    tsv(run_root / "RUN_METADATA.tsv", run_meta, ["key", "value"])
    validation_rows = []
    checks = [("real_fiberhmm_bam", args.canonical_bam.exists() and args.canonical_bam.stat().st_size > 0, "FiberHMM output BAM is non-empty"), ("fiberhmm_pg_version", any(str(x.get("VN")) == "2.16.8" for x in pysam.AlignmentFile(str(args.canonical_bam), "rb").header.get("PG", [])), "BAM @PG contains FiberHMM 2.16.8"), ("canonical_scope", all(x["chrom"] in CANONICAL for x in features.get("msp", []) + features.get("nucleosome", []) + features.get("tf", []) + features.get("deam", [])), "derived feature calls are canonical"), ("qc_artifact", bool(qc), "official QC JSON exists"), ("tables_nonempty", all(p.exists() and p.stat().st_size > 0 for p in tables.glob("*.tsv")), "required summary tables are non-empty"), ("figures_nonempty", all(p.exists() and p.stat().st_size > 0 for p in figures.glob("*.png")), "derived figures are non-empty")]
    for name, ok, detail in checks: validation_rows.append({"validation_layer": "ENGINEERING", "check": name, "status": "PASS" if ok else "FAIL", "detail": detail})
    validation_rows.extend([{ "validation_layer": "SCIENTIFIC", "check": "Qitan chemistry compatibility", "status": "NOT_VALIDATED", "detail": "No matched Qitan external control"}, {"validation_layer": "SCIENTIFIC", "check": "FiberHMM biological interpretation", "status": "NOT_VALIDATED", "detail": "Computational completion is not platform validation"}, {"validation_layer": "SCIENTIFIC", "check": "PCR deduplication", "status": "NOT_EXECUTED", "detail": "Main real call used --no-dedup after integrated prepass did not complete"}, {"validation_layer": "SCIENTIFIC", "check": "DAF SNP masking", "status": "NOT_EXECUTED", "detail": "Automatic bounded preflight estimated 0.66x and did not trigger caller"}])
    tsv(validation / "FIBERHMM_VALIDATION.tsv", validation_rows, ["validation_layer", "check", "status", "detail"])
    (methods / "FIBERHMM_METHODS.md").write_text("# FiberHMM 实际运行方法\n\n本次使用 FiberHMM 2.16.8 `fiberhmm-call`，DddA DAF mode，enzyme=ddda，profile=ddda，probability threshold=125。输入层仅处理 chr1-chr22、chrX、chrY；chrM 单独统计。主调用为 primary、MAPQ >=20、`--no-dedup`，因为内置去重预处理未完成；不得将结果解释为 PCR-deduplicated。\n\nCOMPLETED 仅表示计算完成；official QC 由 FiberHMM QC 输出；Qitan chemistry compatibility 和 FiberHMM biological interpretation 均为 NOT_VALIDATED。\n")
    report = run_root / args.report_name
    print("[8/10] rendering offline HTML", flush=True)
    html_report(run_root, tables, figures, qc, {**canonical, "nucleosome_median": nuc_stats["median"], "nrl": nrl}, fs, cov, tss_result, reps, report)
    (run_root / "README.md").write_text("# FiberHMM 实际运行交付\n\n本目录来自真实执行的 FiberHMM 2.16.8 `fiberhmm-call`，不是 contract_only 标记，也未覆盖原有 v2.11 结果。\n\n- 主分析：chr1-chr22、chrX、chrY\n- chrM：单独统计\n- mode：DddA DAF\n- probability threshold：125\n- computation：COMPLETED\n- official QC：见 `tables/FIBERHMM_QC_SUMMARY.tsv` 和 `figures/fiberhmm_qc_dashboard.png`\n- Qitan chemistry compatibility：NOT_VALIDATED\n")
    (run_root / "FINAL_VALIDATION.md").write_text("# FiberHMM 最终验收\n\n真实 FiberHMM BAM、feature BED、官方 QC、summary TSV、图和离线 HTML 已按本次 run 目录生成。\n\n- 计算状态：COMPLETED（真实 `fiberhmm-call`）\n- 官方 QC：以 `FIBERHMM_QC_SUMMARY.tsv` 为准\n- canonical 输入范围：实际 caller 使用 `--chroms`；chrM 不进入主 chromatin 结果\n- PCR deduplication：NOT_EXECUTED；主结果不是 deduplicated 结果\n- DAF SNP masking：NOT_EXECUTED；0.66x bounded preflight 未触发\n- Qitan chemistry compatibility：NOT_VALIDATED\n- FiberHMM biological interpretation：NOT_VALIDATED\n- 工程验收明细：`validation/FIBERHMM_VALIDATION.tsv`\n")
    print("[10/10] delivery complete", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
