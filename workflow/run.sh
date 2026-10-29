#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mode="${1:-run}"
config_file="${CONFIG:-config/profiles/demo.yaml}"
python_bin="${PYTHON_BIN:-$root_dir/bin/python}"
case "$mode" in
    run) snakemake_mode="--rerun-incomplete" ;;
    dry-run) snakemake_mode="--dry-run" ;;
    report) snakemake_mode="" ;;
    *) echo "Usage: $0 [run|dry-run|report]" >&2; exit 2 ;;
esac
mkdir -p "$root_dir/tmp/layout-validation" "$root_dir/tmp/snakemake-work"
if [[ "$config_file" != /* ]]; then
    config_file="$root_dir/$config_file"
fi
export TMPDIR="$root_dir/tmp/layout-validation" TMP="$root_dir/tmp/layout-validation" TEMP="$root_dir/tmp/layout-validation"
export XDG_CACHE_HOME="$root_dir/tmp/layout-validation/cache"
export MPLCONFIGDIR="$root_dir/tmp/layout-validation/matplotlib"
if [[ "$mode" == report ]]; then
    report_python="$("$python_bin" -c 'import sys,yaml; c=yaml.safe_load(open(sys.argv[1])); print(c.get("report", {}).get("python", sys.executable))' "$config_file")"
    exec "$report_python" "$root_dir/pipeline/scripts/render_illumina_pe_report.py" \
        --configfile "$config_file" --project-root "$root_dir"
fi
exec "$python_bin" -m snakemake --directory "$root_dir/tmp/snakemake-work" \
    --snakefile "$root_dir/pipeline/Snakefile" --configfile "$config_file" \
    --cores 1 "$snakemake_mode"
