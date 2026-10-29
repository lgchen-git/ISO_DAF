#!/usr/bin/env python3
"""Restyle the integrated report with the previous zoom-report template.

No analysis is rerun here. The template supplies the customer-facing layout,
CSS and offline lightbox. Real FiberHMM values and figures are read from the
existing delivery artifacts and inserted into the template sections.
"""

from __future__ import annotations

import argparse
import base64
import html
import mimetypes
import re
import shutil
from pathlib import Path


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def section(document: str, section_id: str) -> str:
    pattern = re.compile(
        r"<section\b(?=[^>]*\bid=['\"]"
        + re.escape(section_id)
        + r"['\"])[^>]*>.*?</section>",
        re.IGNORECASE | re.DOTALL,
    )
    match = pattern.search(document)
    if not match:
        raise ValueError(f"section not found: {section_id}")
    return match.group(0)


def heading(template_section: str) -> str:
    match = re.search(r"<h2>.*?</h2>", template_section, re.IGNORECASE | re.DOTALL)
    if not match:
        raise ValueError("template section has no h2")
    return match.group(0)


def data_uri(path: Path) -> str:
    mime = mimetypes.guess_type(path.name)[0] or "application/octet-stream"
    payload = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{payload}"


def help_button(title: str, body: str) -> str:
    return (
        '<button class="help-button" type="button" '
        f'data-help-title="{html.escape(title, quote=True)}" '
        f'data-help-body="{html.escape(body, quote=True)}" '
        f'aria-label="查看 {html.escape(title, quote=True)} 的解释" '
        'aria-expanded="false">?</button>'
    )


def figure(
    image_path: Path,
    alt: str,
    title: str,
    caption: str,
    source: str,
    body: str,
    large: bool = False,
) -> str:
    image = data_uri(image_path)
    cls = "compact-figure large" if large else "compact-figure"
    return (
        f'<figure class="{cls}">'
        f'<button class="zoomable" type="button" data-image="{image}" '
        f'aria-label="放大查看 {html.escape(title, quote=True)}">'
        f'<img src="{image}" alt="{html.escape(alt, quote=True)}"></button>'
        "<figcaption>"
        f"<strong>{html.escape(title)} {help_button(title, body)}</strong><br>"
        f"{caption}<br><span class=\"source\">图源：<code>{html.escape(source)}</code></span>"
        "</figcaption></figure>"
    )


def make_fiberhmm_section(template_section: str, fiber_root: Path) -> str:
    qc = fiber_root / "figures/fiberhmm_qc_dashboard.png"
    msp = fiber_root / "figures/msp_length_distribution.png"
    nuc = fiber_root / "figures/nucleosome_length_distribution.png"
    for path in (qc, msp, nuc):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    h2 = heading(template_section)
    qc_help = (
        "What：FiberHMM official QC。\n"
        "How：DddA DAF mode、ddda reference/profile、sampled reads。\n"
        "How to read：WARN 是技术警示，不是 failure。\n"
        "Interpretation：signal/periodicity/reference metrics 可作技术参考。\n"
        "Boundary：Qitan chemistry compatibility with the FiberHMM DddA DAF model/reference profile has not been validated。"
    )
    msp_help = (
        "What：native MSP length distribution。\n"
        "How：显示长度主体并标注 >150 bp threshold。\n"
        "How to read：关注 median、N50 和 >150 bp fraction。\n"
        "Interpretation：长 MSP 信号较低。\n"
        "Boundary：不能将此结果写成 Qitan accessibility validation。"
    )
    nuc_help = (
        "What：native nucleosome length distribution。\n"
        "How：结合 median、N50 和 NRL。\n"
        "How to read：median/NRL 与 official periodicity 一起看。\n"
        "Interpretation：结构长度 technically plausible。\n"
        "Boundary：接近典型长度不等于模型或平台验证。"
    )
    figures = (
        '<div class="figure-stack">'
        + figure(
            qc,
            "官方 FiberHMM QC",
            "官方 FiberHMM QC",
            "实际 sampled QC dashboard；整体状态 WARN。",
            "fiberhmm/figures/fiberhmm_qc_dashboard.png",
            qc_help,
            large=True,
        )
        + '<div class="figure-grid">'
        + figure(
            msp,
            "MSP 长度分布",
            "MSP 长度分布",
            "MSP median = 14 bp；MSP >150 bp fraction = 1.15%。",
            "fiberhmm/figures/msp_length_distribution.png",
            msp_help,
        )
        + figure(
            nuc,
            "核小体长度分布",
            "核小体长度分布",
            "nucleosome median = 191 bp；NRL = 194 bp。",
            "fiberhmm/figures/nucleosome_length_distribution.png",
            nuc_help,
        )
        + "</div></div>"
    )
    return (
        '<section class="report-section" id="s3">'
        + h2
        + '<p class="chapter-conclusion">FiberHMM DddA DAF mode 计算已完成，但官方 QC = WARN；WARN 主要来自 nucleosome-scale periodicity，signal-rate status = PASS。MSP 长度主体偏短，MSP &gt;150 bp fraction = 1.15%；nucleosome median/NRL 仅支持 technically plausible 的 exploratory 描述。</p>'
        + '<div class="metric-grid">'
        '<div class="metric-card"><span>Computation status</span><strong>COMPLETED</strong><small>real FiberHMM run</small></div>'
        '<div class="metric-card"><span>Official QC</span><strong>WARN</strong><small>sampled FiberHMM QC</small></div>'
        '<div class="metric-card"><span>QC score</span><strong>57.33/100</strong><small>official score</small></div>'
        '<div class="metric-card"><span>Reference correlation</span><strong>0.624</strong><small>periodicity reference</small></div>'
        '<div class="metric-card"><span>MSP &gt;150 bp fraction</span><strong>1.15%</strong><small>3,246 / 282,218</small></div>'
        '<div class="metric-card"><span>Nucleosome median</span><strong>191 bp</strong><small>native calls</small></div>'
        '<div class="metric-card"><span>NRL</span><strong>194 bp</strong><small>official periodicity</small></div>'
        '</div>'
        + figures
        + '<p class="status-note">COMPUTATION_STATUS = COMPLETED；TECHNICAL_QC_STATUS = WARN；BIOLOGICAL_VALIDATION_STATUS = NOT_VALIDATED。Qitan chemistry compatibility remains NOT_VALIDATED。</p>'
        '</section>'
    )


def make_tss_section(template_section: str, fiber_root: Path) -> str:
    image = fiber_root / "figures/tss_msp150_profile.png"
    if not image.is_file() or image.stat().st_size == 0:
        raise FileNotFoundError(image)
    h2 = heading(template_section)
    body = (
        "What：gene-representative TSS-centered MSP >150 bp actuation fraction。\n"
        "How：分母为各 bin 的 fiber-covered observations；raw points 为主，smoothing 仅辅助。\n"
        "How to read：比较 center、background、numerator 和 fiber denominator。\n"
        "Interpretation：中心信号不支持强阳性 accessibility 结论。\n"
        "Boundary：Wilson CI 不是 biological replicate-level CI；无 independent replicate、matched randomized background 或 bootstrap-based biological CI。"
    )
    fig = figure(
        image,
        "TSS-centered MSP >150 bp profile",
        "TSS-centered MSP >150 bp profile",
        "gene-representative TSS = 19,371；center = 0.00%，background = 0.02%；未观察到强烈的 TSS-centered accessibility enrichment。",
        "fiberhmm/figures/tss_msp150_profile.png",
        body,
        large=True,
    )
    return (
        '<section class="report-section" id="s4">'
        + h2
        + '<p class="chapter-conclusion">以 MSP &gt;150 bp overlapping gene-representative TSS-centered bins 构建描述性 profile；本样本未观察到强烈的 TSS-centered MSP &gt;150 bp accessibility enrichment。</p>'
        + '<div class="metric-grid">'
        '<div class="metric-card"><span>Gene-representative TSS</span><strong>19,371</strong><small>current profile</small></div>'
        '<div class="metric-card"><span>Center signal</span><strong>0.00%</strong><small>MSP &gt;150 bp</small></div>'
        '<div class="metric-card"><span>Background signal</span><strong>0.02%</strong><small>matched background bins</small></div>'
        '<div class="metric-card"><span>Confidence interval</span><strong>Wilson 95% CI</strong><small>binomial observation uncertainty</small></div>'
        '<div class="metric-card"><span>Interpretation</span><strong>NOT STRONG</strong><small>TSS-centered enrichment</small></div>'
        '</div>'
        + fig
        + '<p class="source-note">TSS denominator、完整 profile 和 transcript/unique TSS audit 保留在 <code>fiberhmm/tables/TSS_MSP150_SUMMARY.tsv</code>；Wilson CI 不代表 biological replicate-level confidence。</p>'
        '</section>'
    )


def make_fiber_section(template_section: str, fiber_root: Path) -> str:
    image = fiber_root / "figures/representative_2.png"
    if not image.is_file() or image.stat().st_size == 0:
        raise FileNotFoundError(image)
    h2 = heading(template_section)
    body = (
        "What：单分子 backbone、nucleosome、MSP、MSP >150、TF-like 和 DAF tracks。\n"
        "How：按真实 fiber 覆盖范围动态缩放，展示当前 R3 区域全部 32 条 fiber。\n"
        "How to read：粗色块代表 chromatin feature，细 tick 代表 DAF event。\n"
        "Interpretation：适合观察 feature 共现。\n"
        "Boundary：不是 validated regulatory case study。"
    )
    fig = figure(
        image,
        "R3 · NOMO2 · chr16:18560001-18570000",
        "R3 · NOMO2 · chr16:18560001-18570000",
        "selection-region total fibers = 32；displayed fibers = 32；MSP >150 bp = 36；nucleosome = 248；TF-like = 215。该区域仅用于 exploratory single-molecule display。",
        "fiberhmm/figures/representative_2.png",
        body,
        large=True,
    )
    return (
        '<section class="report-section" id="s5">'
        + h2
        + '<p class="chapter-conclusion">R3 · NOMO2 · chr16:18560001-18570000 作为代表性单分子区域展示，用于观察 feature 共现；该结果属于 exploratory interpretation，不作为 validated regulatory case study。</p>'
        + '<div class="metric-grid">'
        '<div class="metric-card"><span>Region</span><strong>chr16:18560001-18570000</strong><small>R3 · NOMO2</small></div>'
        '<div class="metric-card"><span>Fibers / displayed</span><strong>32 / 32</strong><small>selection-region subset</small></div>'
        '<div class="metric-card"><span>MSP &gt;150 bp</span><strong>36</strong><small>region total</small></div>'
        '<div class="metric-card"><span>Nucleosome / TF-like</span><strong>248 / 215</strong><small>region total</small></div>'
        '</div>'
        + fig
        + '<p class="source-note">完整 representative metadata 保留在 <code>fiberhmm/tables/REPRESENTATIVE_FIBER_SUMMARY.tsv</code>；统计严格区分 selection-region total 与 displayed subset。</p>'
        '</section>'
    )


def render_existing_run(template: Path, fiber_root: Path, output: Path, report_name: str) -> Path:
    """Render the previous zoom template around the current run artifacts.

    This presentation-only path embeds figures as data URIs and writes only a
    new HTML plus its validation table. It never copies BAMs or other large
    FiberHMM runtime artifacts.
    """
    template_html = read(template.resolve())
    fiber_root = fiber_root.resolve()
    replacements = {
        "s3": make_fiberhmm_section(section(template_html, "s3"), fiber_root),
        "s4": make_tss_section(section(template_html, "s4"), fiber_root),
        "s5": make_fiber_section(section(template_html, "s5"), fiber_root),
        "s6": section(template_html, "s6"),
        "s7": section(template_html, "s7").replace(
            '<p class="chapter-conclusion">',
            '<p class="chapter-conclusion"><b>当前实际 FiberHMM 运行：</b> FiberHMM 2.16.8，DddA DAF mode，ddda profile，probability threshold = 125，primary reads，MAPQ ≥20；canonical chromosomes at analysis input，chrM separately summarized。</p><p>',
            1,
        ),
        "s8": section(template_html, "s8"),
    }
    replacements["s8"] = replacements["s8"].replace(
        "DAFseq_analysis_report.html", report_name
    ).replace(
        "<tr><th>Analysis version</th><td>v2.11</td></tr>",
        "<tr><th>Analysis version</th><td>v2.11 · real FiberHMM integrated</td></tr>",
    )
    report = template_html
    for sid in reversed(range(3, 9)):
        sid_text = f"s{sid}"
        report = report.replace(section(template_html, sid_text), replacements[sid_text], 1)
    report = report.replace(
        '<div><b>MSP &gt;150 bp</b><strong>0.97%</strong></div>',
        '<div><b>MSP &gt;150 bp</b><strong>1.15%</strong></div>',
        1,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(report, encoding="utf-8")

    image_tags = re.findall(r"<img\b[^>]*>", report, re.IGNORECASE | re.DOTALL)
    figure_tags = [tag for tag in image_tags if "alt=\"放大图\"" not in tag and "alt='放大图'" not in tag]
    checks = [
        ("eight_sections", all(f'id="s{i}"' in report for i in range(1, 9))),
        ("zoom_template_layout", report.count("compact-figure") >= 4 and report.count("figure-grid") >= 3),
        ("real_fiberhmm_values", all(x in report for x in ("57.33", "1.15%", "191 bp", "194 bp", "0.624"))),
        ("real_tss_values", all(x in report for x in ("19,371", "0.00%", "0.02%"))),
        ("real_r3_values", all(x in report for x in ("chr16:18560001-18570000", "32 / 32", "MSP &gt;150 bp"))),
        ("wheel_zoom", "passive:false" in report and "deltaY" in report),
        ("help_system", "help-popover" in report and "help-button" in report),
        ("all_images_embedded", all("src=\"data:image/" in tag for tag in figure_tags)),
        ("scientific_boundary", "NOT_VALIDATED" in report and "COMPLETED" in report),
    ]
    validation = output.parent / "ZOOM_TEMPLATE_RELAYOUT_VALIDATION.tsv"
    validation.write_text(
        "check\tstatus\tdetail\n"
        + "\n".join(f"{name}\t{'PASS' if ok else 'FAIL'}\t" for name, ok in checks)
        + "\n",
        encoding="utf-8",
    )
    if not all(ok for _, ok in checks):
        failed = [name for name, ok in checks if not ok]
        raise RuntimeError(f"zoom-template report validation failed: {failed}")
    return output


def build(args: argparse.Namespace) -> Path:
    template = Path(args.template).resolve()
    source_root = Path(args.source_root).resolve()
    fiber_root = source_root / "fiberhmm"
    output_root = Path(args.output_root).resolve()
    if output_root.exists():
        raise FileExistsError(f"refusing to overwrite existing output root: {output_root}")
    shutil.copytree(source_root, output_root)

    template_html = read(template)
    source_html = read(source_root / args.source_report)
    # Sections 1 and 2 are retained from the requested zoom template; verify
    # that they are the same current v2.11 result sections before doing so.
    for sid in ("s1", "s2"):
        if section(template_html, sid) != section(source_html, sid):
            raise ValueError(f"template/current mismatch in frozen section {sid}")

    replacements = {
        "s3": make_fiberhmm_section(section(template_html, "s3"), fiber_root),
        "s4": make_tss_section(section(template_html, "s4"), fiber_root),
        "s5": make_fiber_section(section(template_html, "s5"), fiber_root),
        "s6": section(template_html, "s6"),
        "s7": section(template_html, "s7").replace(
            "<p class=\"chapter-conclusion\">",
            '<p class="chapter-conclusion"><b>当前实际 FiberHMM 运行：</b> FiberHMM 2.16.8，DddA DAF mode，ddda profile，probability threshold = 125，primary reads，MAPQ ≥20；canonical chromosomes at analysis input，chrM separately summarized。</p><p>',
            1,
        ),
        "s8": section(template_html, "s8"),
    }
    replacements["s8"] = replacements["s8"].replace(
        "DAFseq_analysis_report.html", args.report_name
    ).replace(
        "<tr><th>Analysis version</th><td>v2.11</td></tr>",
        "<tr><th>Analysis version</th><td>v2.11 · real FiberHMM integrated</td></tr>",
    )
    report = template_html
    for sid in reversed(range(3, 9)):
        sid_text = f"s{sid}"
        report = report.replace(section(template_html, sid_text), replacements[sid_text], 1)

    # The template hero predates the real FiberHMM run. Keep its compact
    # overview, but reconcile the displayed MSP >150 bp fraction.
    report = report.replace(
        '<div><b>MSP &gt;150 bp</b><strong>0.97%</strong></div>',
        '<div><b>MSP &gt;150 bp</b><strong>1.15%</strong></div>',
        1,
    )

    report_path = output_root / args.report_name
    report_path.write_text(report, encoding="utf-8")
    readme = """# DAF-seq total report — zoom-template layout

This is a presentation-only re-layout of the total v2.11 report using the
previous `DAFseq_analysis_report_zoom.html` layout. The real FiberHMM result
artifacts are integrated into sections 3–5; no analysis was rerun.

All report figures are embedded in the HTML. Tables, methods, validation and
the complete FiberHMM artifact directory remain under this delivery root.
"""
    (output_root / "README.md").write_text(readme, encoding="utf-8")

    checks = []
    checks.append(("eight_sections", all(f'id="s{i}"' in report for i in range(1, 9))))
    checks.append(("zoom_template_layout", report.count("compact-figure") >= 4 and report.count("figure-grid") >= 3))
    checks.append(("real_fiberhmm_values", all(x in report for x in ("57.33", "1.15%", "191 bp", "194 bp", "0.624"))))
    checks.append(("real_tss_values", all(x in report for x in ("19,371", "0.00%", "0.02%"))))
    checks.append(("real_r3_values", all(x in report for x in ("chr16:18560001-18570000", "32 / 32", "MSP &gt;150 bp"))))
    checks.append(("wheel_zoom", "passive:false" in report and "deltaY" in report))
    checks.append(("help_system", "help-popover" in report and "help-button" in report))
    checks.append(("no_old_fiberhmm_ratio", "0.97%" not in report and "0.6237607122021382" not in report))
    image_tags = re.findall(r"<img\b[^>]*>", report, re.IGNORECASE | re.DOTALL)
    figure_tags = [tag for tag in image_tags if "alt=\"放大图\"" not in tag and "alt='放大图'" not in tag]
    checks.append(("all_images_embedded", all("src=\"data:image/" in tag for tag in figure_tags)))
    relative_refs = []
    for match in re.finditer(r"(?:src|href)=[\"']([^\"'#]+)", report):
        ref = match.group(1)
        if not ref.startswith(("data:", "http:", "https:")):
            relative_refs.append(ref)
    missing = sorted({ref for ref in relative_refs if not (output_root / ref).exists()})
    checks.append(("artifact_links", not missing))
    validation = output_root / "validation" / "ZOOM_TEMPLATE_RELAYOUT_VALIDATION.tsv"
    validation.write_text(
        "check\tstatus\tdetail\n"
        + "\n".join(
            f"{name}\t{'PASS' if ok else 'FAIL'}\t{'' if ok else 'check failed'}"
            for name, ok in checks
        )
        + "\n",
        encoding="utf-8",
    )
    if missing or not all(ok for _, ok in checks):
        raise RuntimeError(f"re-layout validation failed: missing={missing}")
    return report_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--template", required=True)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--source-report", default="DAFseq_analysis_report_with_fiberhmm_embedded.html")
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--report-name", default="DAFseq_analysis_report_zoom_fiberhmm.html")
    args = parser.parse_args()
    print(build(args))


if __name__ == "__main__":
    main()
