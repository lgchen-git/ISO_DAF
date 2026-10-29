#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root_dir"

python_bin="${PYTHON_BIN:-$root_dir/bin/python}"
cache_root="$root_dir/tmp/layout-validation"
snakemake_work="$root_dir/tmp/snakemake-work"
mkdir -p "$cache_root"
mkdir -p "$snakemake_work"
export TMPDIR="$cache_root"
export TMP="$cache_root"
export TEMP="$cache_root"
export XDG_CACHE_HOME="$cache_root/cache"

printf 'Running Snakemake smoke test in %s\n' "$root_dir"

"$python_bin" -m py_compile pipeline/scripts/stage_runner.py pipeline/scripts/run_report.py
"$python_bin" pipeline/scripts/legacy/check_samples.py --samples config/samples.tsv --skip-path-checks
"$python_bin" -m snakemake --directory "$snakemake_work" --snakefile "$root_dir/pipeline/Snakefile" --configfile "$root_dir/config/profiles/demo.yaml" --cores 1 --dry-run
"$python_bin" -m snakemake --directory "$snakemake_work" --snakefile "$root_dir/pipeline/Snakefile" --configfile "$root_dir/config/profiles/demo.yaml" --cores 1 --rerun-incomplete

test -f outputs/demo/08_report/demo_pipeline_report.html
test -f outputs/demo/01_preflight/stage_complete.tsv
test -f log/demo/08_report.log
test -d tmp/demo

grep -qi 'profile.*demo' outputs/demo/08_report/demo_pipeline_report.html
grep -q 'COMPLETED' outputs/demo/08_report/demo_pipeline_report.html

printf 'Snakemake smoke test passed.\n'
