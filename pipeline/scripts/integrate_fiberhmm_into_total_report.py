#!/usr/bin/env python3
"""Integrate real FiberHMM result sections into the existing v2.11 report.

This is a presentation-layer operation. It does not rerun any assay analysis.
The existing total report is used as the base; the FiberHMM sections are
copied from the real-run report and their artifact paths are namespaced under
``fiberhmm/`` in a new delivery root.
"""

from __future__ import annotations

import argparse
import base64
import mimetypes
import re
import shutil
from pathlib import Path


SECTION_MAP = {
    "s3": "chromatin",
    "s4": "tss",
    "s5": "fiber",
    "s6": "advanced",
    "s7": "methods",
    "s8": "delivery",
}


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def extract_section(document: str, section_id: str) -> str:
    pattern = re.compile(
        r"<section\b(?=[^>]*\bid=[\"']" + re.escape(section_id) + r"[\"'])[^>]*>.*?</section>",
        re.IGNORECASE | re.DOTALL,
    )
    match = pattern.search(document)
    if not match:
        raise ValueError(f"section not found: {section_id}")
    return match.group(0)


def normalize_section(section: str, report_id: str, base_heading: str) -> str:
    section = re.sub(
        r"(\bid=)[\"'][^\"']+[\"']",
        rf'\1"{report_id}"',
        section,
        count=1,
    )
    heading_match = re.search(r"<h2>.*?</h2>", section, re.IGNORECASE | re.DOTALL)
    if not heading_match:
        raise ValueError(f"section has no h2: {report_id}")
    section = section[: heading_match.start()] + base_heading + section[heading_match.end() :]
    # This exact replacement covers both quote styles without touching
    # external URLs or data URIs.
    for quote in ('"', "'"):
        for directory in ("figures", "tables", "methods", "validation", "logs"):
            section = section.replace(
                f"src={quote}{directory}/", f"src={quote}fiberhmm/{directory}/"
            )
            section = section.replace(
                f"href={quote}{directory}/", f"href={quote}fiberhmm/{directory}/"
            )
    return section


def copy_tree(source: Path, target: Path) -> None:
    if source.exists():
        shutil.copytree(source, target, dirs_exist_ok=True)


def embed_fiberhmm_images(document: str, output_root: Path) -> str:
    pattern = re.compile(
        r"src=(['\"])(fiberhmm/figures/[^'\"]+\.(?:png|jpe?g|svg))\1",
        re.IGNORECASE,
    )

    def replace(match: re.Match[str]) -> str:
        quote, relative = match.group(1), match.group(2)
        image_path = output_root / relative
        if not image_path.is_file():
            raise FileNotFoundError(f"FiberHMM image not found: {image_path}")
        mime = mimetypes.guess_type(image_path.name)[0] or "application/octet-stream"
        encoded = base64.b64encode(image_path.read_bytes()).decode("ascii")
        return f"src={quote}data:{mime};base64,{encoded}{quote}"

    return pattern.sub(replace, document)


def build_report(args: argparse.Namespace) -> Path:
    base_report = Path(args.base_report).resolve()
    fiber_report = Path(args.fiberhmm_report).resolve()
    base_root = Path(args.base_root).resolve()
    fiber_root = Path(args.fiberhmm_root).resolve()
    output_root = Path(args.output_root).resolve()

    if output_root.exists():
        raise FileExistsError(f"refusing to overwrite existing output root: {output_root}")
    output_root.mkdir(parents=True)

    # Preserve the current total-report delivery context.
    for directory in ("tables", "figures", "methods", "validation"):
        copy_tree(base_root / directory, output_root / directory)

    # Copy the real FiberHMM-derived artifacts under a separate namespace.
    fiber_delivery = output_root / "fiberhmm"
    fiber_delivery.mkdir()
    for directory in ("tables", "figures", "methods", "validation", "logs"):
        copy_tree(fiber_root / directory, fiber_delivery / directory)
    for filename in (
        "RUN_METADATA.tsv",
        "FIBERHMM_EXECUTION_LOG.md",
        "FINAL_VALIDATION.md",
        "README.md",
    ):
        source = fiber_root / filename
        if source.exists():
            shutil.copy2(source, fiber_delivery / filename)

    base = read_text(base_report)
    fiber = read_text(fiber_report)

    base_sections = {sid: extract_section(base, sid) for sid in SECTION_MAP}
    fiber_sections = {
        sid: normalize_section(extract_section(fiber, fiber_id), sid, re.search(r"<h2>.*?</h2>", base_sections[sid], re.I | re.S).group(0))
        for sid, fiber_id in SECTION_MAP.items()
    }

    integrated = base
    # Replace longest sections first so section matching remains unambiguous.
    for sid in reversed(list(SECTION_MAP)):
        integrated = integrated.replace(base_sections[sid], fiber_sections[sid], 1)

    # Keep displayed scientific ratios readable while preserving source TSV
    # precision in the copied artifacts.
    integrated = re.sub(
        r"reference correlation=([0-9]+\.[0-9]+)",
        lambda match: f"reference correlation={float(match.group(1)):.3f}",
        integrated,
    )
    integrated = re.sub(
        r"matched amplitude=([0-9]+\.[0-9]+)",
        lambda match: f"matched amplitude={float(match.group(1)):.3f}",
        integrated,
    )

    # The existing hero is retained, but its FiberHMM-derived MSP >150 bp
    # metric must match the real run (3,246 / 282,218 = 1.15%).
    hero_start = integrated.find('<section class="hero">')
    first_section = integrated.find("<section", hero_start + 1)
    if hero_start >= 0 and first_section > hero_start:
        hero = integrated[hero_start:first_section]
        hero = hero.replace("0.97%", "1.15%", 1)
        integrated = integrated[:hero_start] + hero + integrated[first_section:]

    # Keep the total report's customer-facing metadata while identifying the
    # real FiberHMM artifact source in the delivery chapter links.
    integrated = integrated.replace(
        "DAFseq_fiberhmm_real_report.html", args.report_name
    )
    integrated = integrated.replace(
        "FiberHMM 染色质结构分析报告", "DAF-seq 分析报告"
    )

    # The base report supports button/data-image zoom. Add image zoom for the
    # real FiberHMM figures and wheel scaling for both image types.
    integrated = integrated.replace(
        ".zoomable{display:block;width:100%;border:0;background:#fff;padding:0;cursor:zoom-in;text-align:left}",
        ".zoomable{display:block;width:100%;border:0;background:#fff;padding:0;cursor:zoom-in;text-align:left}.report-section figure{width:96%;margin:18px auto;text-align:center}.report-section figure img.zoomable{display:block;width:100%;height:auto;max-width:100%;max-height:none;object-fit:contain}#lightbox{overflow:auto}",
        1,
    )
    script = r"""<script>(()=>{const box=document.getElementById('lightbox'),img=box.querySelector('img');let scale=1;const closeImage=()=>{box.classList.remove('open');img.src='';img.style.transform='scale(1)';scale=1};const imageSource=el=>el.dataset.image||el.currentSrc||el.src;document.querySelectorAll('.zoomable').forEach(el=>el.addEventListener('click',e=>{e.preventDefault();img.src=imageSource(el);img.alt=el.alt||'放大图';box.classList.add('open')}));box.addEventListener('click',e=>{if(e.target===box)closeImage()});box.addEventListener('wheel',e=>{if(!box.classList.contains('open'))return;e.preventDefault();scale=Math.min(8,Math.max(.3,scale*(e.deltaY<0?1.15:.87)));img.style.transform=`scale(${scale})`},{passive:false});const help=document.getElementById('help-popover'),title=help.querySelector('[data-help-title]'),body=help.querySelector('[data-help-body]'),closeHelp=()=>{help.classList.remove('open');document.querySelectorAll('.help-button').forEach(b=>b.setAttribute('aria-expanded','false'))};document.querySelectorAll('.help-button').forEach(btn=>btn.addEventListener('click',e=>{e.stopPropagation();title.textContent=btn.dataset.helpTitle;body.textContent=btn.dataset.helpBody;help.classList.add('open');btn.setAttribute('aria-expanded','true')}));help.querySelector('.help-close').onclick=closeHelp;document.addEventListener('click',e=>{if(help.classList.contains('open')&&!help.contains(e.target)&&!e.target.classList.contains('help-button'))closeHelp()});document.onkeydown=e=>{if(e.key==='Escape'){closeImage();closeHelp()}}})()</script>"""
    integrated = re.sub(r"<script>.*?</script>", script, integrated, count=1, flags=re.DOTALL)

    # The original report embeds its own figures. Do the same for the real
    # FiberHMM figures so the integrated HTML remains usable when copied
    # without its surrounding directory.
    integrated = embed_fiberhmm_images(integrated, output_root)

    report_path = output_root / args.report_name
    report_path.write_text(integrated, encoding="utf-8")

    readme = """# DAF-seq total report with real FiberHMM results

This delivery is a new presentation-layer integration of the existing v2.11
total report and the independently executed FiberHMM result set. The original
v2.11 report and the standalone FiberHMM delivery are preserved unchanged.

- Total-report base: existing v2.11 `DAFseq_analysis_report.html`
- Integrated report: `{report_path.name}`
- FiberHMM source artifacts: `fiberhmm/`
- Main analysis scope: chr1-chr22, chrX, chrY
- chrM: separately summarized and excluded from the main nuclear result
- Scientific boundary: FiberHMM computation is completed; Qitan chemistry
  compatibility and biological interpretation remain NOT_VALIDATED.

The HTML is offline and self-contained for all report figures. It keeps the
existing eight-section customer layout and supports image lightbox zoom with
mouse-wheel scaling and Esc/mask close.
"""
    (output_root / "README.md").write_text(readme, encoding="utf-8")

    checks = []
    checks.append(("report_exists", report_path.exists()))
    checks.append(("eight_sections", all(f'id="{sid}"' in integrated for sid in [f"s{i}" for i in range(1, 9)])))
    checks.append(("real_fiberhmm_status", "COMPUTATION_STATUS" in integrated and "NOT_VALIDATED" in integrated))
    checks.append(("real_qc_score", "57.33" in integrated))
    checks.append(("real_msp150", "3,246" in integrated and "1.15%" in integrated))
    checks.append(("wheel_zoom", "passive:false" in integrated and "deltaY" in integrated))
    checks.append(("help_system", "help-popover" in integrated and "help-button" in integrated))
    checks.append(("fiberhmm_images_embedded", "src='fiberhmm/figures/" not in integrated and 'src="fiberhmm/figures/' not in integrated))
    missing = []
    for match in re.finditer(r"(?:src|href)=[\"'](fiberhmm/[^\"'#]+)", integrated):
        target = output_root / match.group(1)
        if not target.exists():
            missing.append(str(target))
    checks.append(("fiberhmm_links", not missing))
    validation = output_root / "validation" / "INTEGRATED_REPORT_VALIDATION.tsv"
    validation.parent.mkdir(parents=True, exist_ok=True)
    validation.write_text(
        "check\tstatus\tdetail\n"
        + "\n".join(
            f"{name}\t{'PASS' if ok else 'FAIL'}\t{'' if ok else 'integration check failed'}"
            for name, ok in checks
        )
        + "\n",
        encoding="utf-8",
    )
    if missing or not all(ok for _, ok in checks):
        raise RuntimeError(f"integration validation failed: missing={missing}")
    return report_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-report", required=True)
    parser.add_argument("--base-root", required=True)
    parser.add_argument("--fiberhmm-report", required=True)
    parser.add_argument("--fiberhmm-root", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument(
        "--report-name",
        default="DAFseq_analysis_report_with_fiberhmm.html",
    )
    args = parser.parse_args()
    print(build_report(args))


if __name__ == "__main__":
    main()
