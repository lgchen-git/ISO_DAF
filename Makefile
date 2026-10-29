SHELL := /usr/bin/env bash

PYTHON ?= bin/python
CONFIG ?= config/profiles/demo.yaml
SNAKEMAKE := $(PYTHON) -m snakemake
CONFIG_PATH := $(abspath $(CONFIG))
SNAKEMAKE_ARGS := --directory $(CURDIR)/tmp/snakemake-work --snakefile $(CURDIR)/pipeline/Snakefile --configfile $(CONFIG_PATH) --cores 1

.PHONY: help lint validate smoke dry-run run report clean tree layout-audit

help:
	@printf "%s\n" \
	"Available targets:" \
	"  make lint       - Python and Snakefile syntax checks" \
	"  make validate   - validate the standard sample sheet" \
	"  make smoke      - run the minimal Snakemake smoke profile" \
	"  make dry-run    - resolve the configured Snakemake DAG" \
	"  make run        - execute the configured Snakemake profile" \
	"  make report     - rebuild the Illumina PE report from existing results" \
	"  make layout-audit - audit the seven-root layout" \
	"  make clean      - remove only demo smoke outputs" \
	"  make tree       - print the repository layout"

lint:
	$(PYTHON) -m py_compile pipeline/scripts/stage_runner.py pipeline/scripts/run_report.py pipeline/scripts/restyle_fiberhmm_report_with_zoom_template.py pipeline/scripts/real_stage_common.py pipeline/scripts/real_preflight.py pipeline/scripts/real_raw_qc.py pipeline/scripts/real_alignment_qc.py pipeline/scripts/real_canonical_filter.py pipeline/scripts/real_daf_signal.py pipeline/scripts/real_fiberhmm.py pipeline/scripts/finalize_existing_fiberhmm.py pipeline/scripts/real_delivery.py
	$(PYTHON) -m py_compile pipeline/scripts/legacy/check_samples.py pipeline/scripts/legacy/get_config.py pipeline/scripts/legacy/daf_metrics_common.py
	$(PYTHON) -m py_compile pipeline/scripts/real_illumina_pe.py pipeline/scripts/real_illumina_pe_tss.py pipeline/scripts/validate_illumina_pe_samples.py
	$(PYTHON) -m py_compile pipeline/scripts/render_illumina_pe_report.py
	bash -n workflow/run.sh

validate:
	$(PYTHON) pipeline/scripts/legacy/check_samples.py --samples config/samples.tsv --skip-path-checks

smoke:
	bash tests/run_smoke.sh

dry-run:
	CONFIG=$(CONFIG) bash workflow/run.sh dry-run

run:
	CONFIG=$(CONFIG) bash workflow/run.sh run

report:
	CONFIG=$(CONFIG) bash workflow/run.sh report

layout-audit:
	$(PYTHON) /opt/project/.codex/skills/bio-pipeline-layout/scripts/audit_analysis_layout.py \
		--config config/analysis_layout.yaml --project-root . --outputs-root outputs --json

clean:
	rm -rf outputs/demo log/demo tmp/demo

tree:
	find . -maxdepth 3 -print | sort
