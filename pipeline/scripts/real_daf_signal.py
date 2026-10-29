#!/usr/bin/env python3
"""Rebuild the canonical DAF signal summaries as one Snakemake stage.

The immutable v2.5 DAF event/read tables are the configured source of called
CT/GA observations.  This stage re-derives the canonical filter, event
universes, callable denominators, TSS profiles, and context/site summaries;
it does not silently substitute a contract marker for those calculations.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np

from real_stage_common import load_config, run_root, write_marker, write_tsv


def build_read_flow(cfg: dict, source: Path, data: dict, root: Path) -> None:
    import qitan_daf_v2_7 as v27

    raw = {x["metric"]: x["value"] for x in v27.read_rows(source / "tables/BC01.v2_5_raw_read_qc.tsv")}
    audit = {x["metric"]: x["value"] for x in v27.read_rows(source / "audit/V2_5_BAM_AUDIT.tsv")}
    fastq = int(float(raw["reads"]))
    primary = int(float(audit["primary"]))
    mapped = int(float(audit["primary_mapped"]))
    reads = data["reads"]
    mapq = int(cfg["thresholds"]["mapq"])
    counts = {
        "fastq_reads": fastq,
        "primary_records": primary,
        "primary_mapped": mapped,
        "canonical_nuclear": len(reads),
        "chrM_excluded": len(data["mito_reads"]),
        "noncanonical_excluded": mapped - len(reads) - len(data["mito_reads"]),
        "mapq_pass": sum(int(x["mapq"]) >= mapq for x in reads),
    }
    counts["low_mapq"] = counts["canonical_nuclear"] - counts["mapq_pass"]
    counts["callable_reads"] = sum(int(x["mapq"]) >= mapq and int(x["callable_C"]) + int(x["callable_G"]) > 0 for x in reads)
    counts["not_callable"] = counts["mapq_pass"] - counts["callable_reads"]
    counts["event_evaluable"] = sum(int(x["mapq"]) >= mapq and int(x["ct_events"]) + int(x["ga_events"]) > 0 for x in reads)
    counts["none_unclassified"] = counts["callable_reads"] - counts["event_evaluable"]
    counts["ct_only"] = sum(x["classification"] == "CT" and int(x["mapq"]) >= mapq for x in reads)
    counts["ga_only"] = sum(x["classification"] == "GA" and int(x["mapq"]) >= mapq for x in reads)
    counts["mixed"] = sum(x["classification"] == "mixed_ambiguous" and int(x["mapq"]) >= mapq for x in reads)
    counts["none_after_classification"] = sum(x["classification"] not in ("CT", "GA", "mixed_ambiguous") and int(x["mapq"]) >= mapq and int(x["ct_events"]) + int(x["ga_events"]) > 0 for x in reads)
    parent = {
        "fastq_reads": None,
        "primary_records": "fastq_reads",
        "primary_mapped": "primary_records",
        "canonical_nuclear": "primary_mapped",
        "chrM_excluded": "primary_mapped",
        "noncanonical_excluded": "primary_mapped",
        "mapq_pass": "canonical_nuclear",
        "low_mapq": "canonical_nuclear",
        "callable_reads": "mapq_pass",
        "not_callable": "mapq_pass",
        "event_evaluable": "callable_reads",
        "none_unclassified": "callable_reads",
        "ct_only": "event_evaluable",
        "ga_only": "event_evaluable",
        "mixed": "event_evaluable",
        "none_after_classification": "event_evaluable",
    }
    labels = {
        "fastq_reads": "FASTQ reads", "primary_records": "Primary BAM records", "primary_mapped": "Primary mapped reads",
        "canonical_nuclear": "Canonical nuclear mapped", "chrM_excluded": "chrM separate scope", "noncanonical_excluded": "Noncanonical excluded",
        "mapq_pass": f"MAPQ >= {mapq}", "low_mapq": f"MAPQ < {mapq}", "callable_reads": "Callable reads", "not_callable": "No callable C/G",
        "event_evaluable": "CT/GA event-evaluable", "none_unclassified": "No CT/GA event", "ct_only": "CT-only", "ga_only": "GA-only",
        "mixed": "Mixed / ambiguous", "none_after_classification": "None / unclassified after classification",
    }
    rows = []
    for node, count in counts.items():
        p = parent[node]
        parent_count = counts[p] if p else fastq
        rows.append({"flow_view": "canonical_nuclear", "node_id": node, "label": labels[node], "parent_node": p or "", "count": str(count), "percent_previous": v27.pct(100 * count / parent_count if parent_count else 0), "percent_fastq": v27.pct(100 * count / fastq if fastq else 0), "unit": "reads", "definition": "primary records only; secondary and supplementary records excluded from the conserved read chain"})
    write_tsv(root / "tables/read_flow_reconciliation.tsv", rows, list(rows[0]))
    flow_dir = root / "read_flow"
    write_tsv(flow_dir / "nodes.input.tsv", [{"view_id": "canonical_nuclear", "node_id": n, "label": labels[n], "value": counts[n], "unit": "reads", "scope": "canonical nuclear primary-read flow"} for n in counts], ["view_id", "node_id", "label", "value", "unit", "scope"])
    write_tsv(flow_dir / "views.tsv", [{"view_id": "canonical_nuclear", "display_name": "Canonical nuclear read flow", "order": "0", "is_default": "true"}], ["view_id", "display_name", "order", "is_default"])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tmp", required=True, type=Path)
    args = parser.parse_args()
    cfg = load_config(args.configfile)
    legacy_dir = args.project_root / "pipeline" / "scripts" / "legacy"
    sys.path.insert(0, str(legacy_dir))
    import qitan_daf_v2_7 as v27

    source = Path(cfg["paths"]["source_root"])
    canonical_root = run_root(cfg, args.project_root) / "04_canonical_filter"
    root = run_root(cfg, args.project_root) / "05_daf_signal"
    read_rows = v27.read_rows(canonical_root / "tables/canonical_filtered_read_classification.tsv.gz")
    event_rows = v27.read_rows(canonical_root / "tables/canonical_filtered_events.tsv.gz")
    site_rows = v27.read_rows(source / "tables/BC01.v2_5_event_sites.tsv.gz")
    fiber_rows = [x for x in v27.read_rows(source / "tables/BC01_v2_5_fiber_chromatin.tsv.gz") if x.get("chromosome") in v27.canonical_set(cfg)]
    data = {
        "reads": read_rows,
        "events": event_rows,
        "sites": site_rows,
        "fibers": fiber_rows,
        "mito_reads": v27.read_rows(source / "tables/BC01.v2_5_read_classification.tsv.gz"),
        "canonical": v27.canonical_set(cfg),
    }
    data["mito_reads"] = [x for x in data["mito_reads"] if x.get("reference_name") == cfg["chromosomes"]["mitochondrial"]]
    data["mito_events"] = []
    data["mito_fibers"] = []

    event_data = v27.aggregate_events(cfg, root, data, site_rows)
    tss_data = v27.parse_tss(Path(cfg["reference"]["gtf"]), data["canonical"])
    callable_scan = None
    audited_tss = Path(cfg.get("paths", {}).get("source_v2_10_root", "")) / "tables"
    audited_profile = audited_tss / "tss_profile.tsv.gz"
    if cfg.get("execution", {}).get("reuse_audited_callable_tss") and audited_profile.is_file():
        profile_rows = v27.read_rows(audited_profile)
        callable_scan = {"callable_bases": 0, "context_callable": Counter(), "profiles": {}}
        for name in ("gene_representative", "unique_genomic", "transcript"):
            subset = [x for x in profile_rows if x.get("profile_type") == name]
            callable = np.asarray([int(x.get("callable_opportunity_count", 0)) for x in subset], dtype=np.int64)
            callable_scan["profiles"][name] = {"callable": callable, "assignment": callable.copy()}
        denom_path = audited_tss / "tss_denominator_audit.tsv"
        if denom_path.is_file():
            denom_rows = v27.read_rows(denom_path)
            callable_scan["callable_bases"] = int(denom_rows[0].get("callable_bases_from_bam", 0)) if denom_rows else 0
        for filename in ("tss_profile.tsv.gz", "tss_summary.tsv", "tss_denominator_audit.tsv"):
            source_file = audited_tss / filename
            if source_file.is_file():
                destination = root / "tables" / filename
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(source_file.read_bytes())
        write_tsv(root / "tables/TSS_INPUT_PROVENANCE.tsv", [{"source": str(audited_profile), "status": "REUSED_AUDITED_BAM_DERIVED", "definition": "v2.10 corrected callable C/G TSS profile reused without changing denominator"}], ["source", "status", "definition"])
        tss_assignments = v27.assign_events_to_tss(event_data["candidate"], tss_data, int(cfg["thresholds"]["tss_window_bp"]), int(cfg["thresholds"]["tss_bin_bp"]))
    else:
        callable_scan = v27.scan_alignment_callable(cfg, data, tss_data, int(cfg["thresholds"]["tss_window_bp"]), int(cfg["thresholds"]["tss_bin_bp"]))
        tss_assignments = v27.make_tss_outputs(cfg, root, data, event_data, tss_data, callable_scan)
    metric = v27.write_metric_tables(cfg, root, source, data, event_data)
    v27.alignment_summary(cfg, root, source, data)
    build_read_flow(cfg, source, data, root)
    v27.site_depth_outputs(root, data | {"sites": event_data["sites"]})
    v27.event_reconciliation(cfg, root, data, event_data, tss_assignments)
    v27.event_universe_sensitivity(cfg, root, data, event_data, tss_data, tss_assignments, callable_scan)
    v27.scientific_boundaries(root)
    (root / "methods").mkdir(parents=True, exist_ok=True)
    v27.build_methods(cfg, root, source, data, metric)
    write_tsv(root / "tables/DAF_SIGNAL_SUMMARY.tsv", [{"metric": "candidate_event_count", "value": len(event_data["candidate"]), "unit": "events", "definition": "canonical CT/GA candidate universe after MAPQ/BQ filters"}, {"metric": "selected_event_count", "value": len(event_data["selected"]), "unit": "events", "definition": "orientation-consistent selected universe"}, {"metric": "callable_denominator", "value": sum(int(x["callable_C"]) for x in data["reads"] if x["classification"] == "CT") + sum(int(x["callable_G"]) for x in data["reads"] if x["classification"] == "GA"), "unit": "C/G opportunities", "definition": "assigned-class callable C/G denominator"}, {"metric": "source_root", "value": source, "unit": "path", "definition": "immutable audited v2.5 DAF tables used for event observations"}], ["metric", "value", "unit", "definition"])
    write_marker(args.output, cfg, "daf_signal", Path(args.tmp), input_path=args.input, extra={"source_root": source, "candidate_events": len(event_data["candidate"]), "selected_events": len(event_data["selected"])})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
