#!/usr/bin/env python3
"""Rebuild the offline Illumina PE report from completed analysis artifacts."""

from __future__ import annotations

import argparse
import base64
import csv
import html
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import yaml


# Layout and colors follow the project's existing DAF zoom report.
CSS = """
:root{--ink:#18202b;--muted:#667085;--line:#dbe3ea;--paper:#f5f7fa;--blue:#3c5488;--green:#00a087;--red:#e64b35}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:var(--paper);color:var(--ink);font:15px/1.65 Arial,"Noto Sans SC","Microsoft YaHei",sans-serif}
.sidebar{position:fixed;inset:0 auto 0 0;width:245px;overflow:auto;background:#15324b;color:white;padding:24px 18px;z-index:5}.sidebar h2{font-size:20px;margin:0 0 10px}.sidebar p{font-size:12px;color:#c8d8e6}.sidebar a{display:block;color:#dceaf5;text-decoration:none;padding:10px 8px;border-radius:6px;font-size:13px}.sidebar a:hover,.sidebar a.active{background:#294b78;color:white}
.content{margin-left:245px;padding:26px 32px 55px;max-width:1780px}.hero{background:var(--blue);color:white;border-radius:12px;padding:24px 28px;margin-bottom:20px}.hero h1{font-size:29px;margin:6px 0}.eyebrow{font-size:12px;letter-spacing:.07em;color:#dceaf5}.hero p{margin:8px 0}.hero-meta{font-size:12px;overflow-wrap:anywhere;color:#e5eef8}
.report-section{scroll-margin-top:18px;background:white;padding:22px 25px;margin:18px 0 24px;border-radius:10px;box-shadow:0 1px 6px #d6dee8}.report-section h2{font-size:22px;margin:0 0 14px;padding-bottom:8px;border-bottom:2px solid var(--line);color:#243b64}.report-section h3{color:#243b64;font-size:17px;margin:22px 0 8px}.chapter-conclusion{font-size:15px;color:#344054}.metric-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(175px,1fr));gap:12px;margin:16px 0}.metric-card{border:1px solid var(--line);border-top:3px solid var(--blue);border-radius:8px;padding:12px}.metric-card span,.metric-card small{display:block;color:var(--muted);font-size:12px}.metric-card strong{display:block;font-size:23px;color:#233d68;margin:5px 0}.finding{border-left:3px solid var(--green);padding:8px 13px;background:#f1fbf8;margin:10px 0}.notice{background:#fff7e6;border-left:4px solid #d9822b;padding:12px 15px;margin:15px 0;border-radius:0 7px 7px 0}.figure-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}.figure-card{border:1px solid var(--line);border-radius:8px;padding:12px;margin:15px 0;background:white}.figure-card img{display:block;width:100%;height:auto}.figure-card figcaption{font-size:13px;color:#475467;padding:10px 2px 0}.figure-card figcaption strong{display:block;color:#243b64;font-size:15px;margin-bottom:4px}.source-note{font-size:11px;color:var(--muted);overflow-wrap:anywhere}.zoomable{display:block;width:100%;border:0;background:white;padding:0;cursor:zoom-in}.export-controls{display:flex;gap:12px;font-size:12px;margin-top:8px}.export-controls a{color:var(--blue)}
.table-scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border:1px solid var(--line);padding:8px 10px;text-align:left;vertical-align:top}th{background:#f3f6f8;color:#344054}td{overflow-wrap:anywhere}thead th{white-space:nowrap}.pager{display:flex;align-items:center;flex-wrap:wrap;gap:12px;margin:12px 0}.pager button,.help-button,select{cursor:pointer}button,select,input{font:inherit}button:focus-visible,a:focus-visible,select:focus-visible,input:focus-visible{outline:3px solid #4dbbd5;outline-offset:2px}.pager button{background:white;border:1px solid var(--line);padding:5px 14px;border-radius:5px}.pager button:disabled{opacity:.4;cursor:default}.filters{display:flex;flex-wrap:wrap;gap:12px;margin:12px 0}.filters input,.filters select,.track-picker select{border:1px solid #b9c5d0;border-radius:5px;padding:6px 8px;max-width:100%}.track-panel[hidden],tr[hidden]{display:none!important}.track-picker{margin:15px 0}details{border:1px solid var(--line);padding:12px 15px;margin:12px 0;border-radius:7px}summary{cursor:pointer;color:#243b64;font-weight:bold}details p:last-child{margin-bottom:0}code{font-size:12px;overflow-wrap:anywhere} .footer{border-top:1px solid var(--line);padding-top:16px;color:var(--muted);font-size:12px;margin-top:28px}
.help-button{width:23px;height:23px;border-radius:50%;border:1px solid #7b8da6;background:white;color:var(--blue);font-weight:bold;font-size:13px;vertical-align:middle;margin-left:8px}.help-popover{display:none;position:fixed;z-index:40;right:24px;bottom:24px;width:min(430px,calc(100vw - 32px));max-height:70vh;overflow:auto;background:white;border:1px solid var(--line);border-top:4px solid var(--blue);border-radius:10px;padding:20px;box-shadow:0 10px 36px #18202b33}.help-popover.open{display:block}.help-popover h3{color:var(--blue);margin:0 26px 10px 0}.help-body{white-space:pre-line;font-size:14px}.help-close{position:absolute;right:8px;top:6px;border:0;background:white;font-size:22px}
#lightbox{display:none;position:fixed;z-index:30;inset:0;background:#101828ed;align-items:center;justify-content:center;padding:45px 25px;overflow:hidden}#lightbox.open{display:flex}#lightbox img{max-width:92vw;max-height:86vh;object-fit:contain;background:white;transform-origin:center;transition:transform .08s;cursor:zoom-in;touch-action:none}.lightbox-close{position:absolute;z-index:31;top:12px;right:18px;border:1px solid #ffffff77;background:#18202b;color:white;border-radius:6px;padding:5px 12px;font-size:23px}.lightbox-note{position:absolute;bottom:10px;color:white;font-size:12px;pointer-events:none}
@media(max-width:1000px){.figure-grid{grid-template-columns:1fr}.content{padding:20px}}
@media(max-width:760px){.sidebar{position:static;width:auto;padding:18px}.sidebar a{display:inline-block;padding:6px 8px}.sidebar p{margin:6px 0}.content{margin:0;padding:14px}.hero{padding:22px 18px}.hero h1{font-size:24px}.report-section{padding:18px 14px}.metric-grid{grid-template-columns:repeat(2,minmax(0,1fr))}.metric-card strong{font-size:19px}}
@media print{body{background:white;font-size:10pt}.sidebar,#lightbox,.help-popover,.help-button,.export-controls,.pager,.filters,.track-picker{display:none!important}.content{margin:0;padding:0;max-width:none}.hero{background:white;color:#18202b;border-bottom:2px solid var(--blue);padding:0 0 12px}.hero-meta,.eyebrow{color:#475467}.report-section{box-shadow:none;padding:12px 0;margin:12px 0}.report-section h2{break-after:avoid}.figure-card,.metric-card{break-inside:avoid}.figure-grid{grid-template-columns:1fr}.track-panel[hidden]{display:block!important}tr[hidden]{display:table-row!important}.table-scroll{overflow:visible}table{font-size:8pt}thead{display:table-header-group}tr{break-inside:avoid}details{display:block}details>*{display:block!important}.footer{font-size:9pt}}
"""

JS = """
(()=>{'use strict';
 const box=document.getElementById('lightbox'), image=box.querySelector('img');let scale=1,lastFigure;
 const closeFigure=()=>{box.classList.remove('open');box.setAttribute('aria-hidden','true');image.style.transform='scale(1)';document.body.style.overflow='';if(lastFigure)lastFigure.focus()};
 document.querySelectorAll('.zoomable').forEach(b=>b.addEventListener('click',()=>{lastFigure=b;image.src=b.querySelector('img').src;image.alt=b.querySelector('img').alt;scale=1;image.style.transform='scale(1)';box.classList.add('open');box.setAttribute('aria-hidden','false');document.body.style.overflow='hidden';box.querySelector('button').focus()}));
 box.querySelector('button').addEventListener('click',closeFigure);box.addEventListener('click',e=>{if(e.target===box)closeFigure()});box.addEventListener('wheel',e=>{if(!box.classList.contains('open'))return;e.preventDefault();scale=Math.min(8,Math.max(.5,scale*(e.deltaY<0?1.15:1/1.15)));image.style.transform=`scale(${scale})`},{passive:false});
 const help=document.getElementById('help-popover');let lastHelp;
 const closeHelp=()=>{help.classList.remove('open');help.setAttribute('aria-hidden','true');if(lastHelp){lastHelp.setAttribute('aria-expanded','false');lastHelp.focus()}};
 document.querySelectorAll('.help-button').forEach(b=>b.addEventListener('click',()=>{if(lastHelp)lastHelp.setAttribute('aria-expanded','false');lastHelp=b;help.querySelector('h3').textContent=b.dataset.helpTitle;help.querySelector('.help-body').textContent=b.dataset.helpBody;help.classList.add('open');help.setAttribute('aria-hidden','false');b.setAttribute('aria-expanded','true');help.querySelector('button').focus()}));help.querySelector('button').addEventListener('click',closeHelp);
 document.addEventListener('keydown',e=>{if(e.key==='Escape'){if(help.classList.contains('open'))closeHelp();else if(box.classList.contains('open'))closeFigure()}});
 document.querySelectorAll('table[data-page-size]').forEach(t=>{const rows=[...t.tBodies[0].rows],pager=document.querySelector(`[data-pager-for="${t.id}"]`),search=document.querySelector(`[data-search-for="${t.id}"]`),chrom=document.querySelector(`[data-chrom-for="${t.id}"]`);let page=1;const size=Number(t.dataset.pageSize);const update=()=>{const query=search?search.value.toLowerCase():'';const filtered=rows.filter(r=>(!query||r.textContent.toLowerCase().includes(query))&&(!chrom||!chrom.value||r.cells[0].textContent===chrom.value));const pages=Math.max(1,Math.ceil(filtered.length/size));page=Math.min(page,pages);rows.forEach(r=>r.hidden=true);filtered.slice((page-1)*size,page*size).forEach(r=>r.hidden=false);pager.querySelector('[data-page-status]').textContent=`第 ${page} / ${pages} 页 · ${filtered.length} 条结果`;pager.querySelector('[data-page-prev]').disabled=page<=1;pager.querySelector('[data-page-next]').disabled=page>=pages};pager.querySelector('[data-page-prev]').onclick=()=>{page--;update()};pager.querySelector('[data-page-next]').onclick=()=>{page++;update()};[search,chrom].filter(Boolean).forEach(x=>x.addEventListener(x===search?'input':'change',()=>{page=1;update()}));update()});
 const picker=document.getElementById('track-select');if(picker)picker.addEventListener('change',()=>document.querySelectorAll('.track-panel').forEach(p=>p.hidden=p.id!==picker.value));
 let openDetails=[];window.addEventListener('beforeprint',()=>{openDetails=[...document.querySelectorAll('details')].filter(x=>!x.open);openDetails.forEach(x=>x.open=true)});window.addEventListener('afterprint',()=>openDetails.forEach(x=>x.open=false));
 if('IntersectionObserver' in window){const obs=new IntersectionObserver(entries=>entries.forEach(e=>{if(e.isIntersecting){document.querySelectorAll('.sidebar a').forEach(a=>a.classList.toggle('active',a.getAttribute('href')==='#'+e.target.id))}}),{rootMargin:'-10% 0px -65% 0px'});document.querySelectorAll('.report-section').forEach(s=>obs.observe(s))}
})();
"""


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fields(path: Path, key: str = "field") -> dict[str, str]:
    return {row[key]: row["value"] for row in rows(path)}


def data_uri(path: Path) -> str:
    mime = {".png": "image/png", ".svg": "image/svg+xml", ".tsv": "text/tab-separated-values"}[path.suffix]
    return f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def write_table(path: Path, headers: list[str], body: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(body)


def table(headers: list[str], body: list[list[object]], ident: str = "", page_size: int = 10) -> str:
    attrs = f' id="{esc(ident)}" class="paginated-table" data-page-size="{page_size}"' if ident else ""
    result = '<div class="table-scroll"><table' + attrs + '><thead><tr>'
    result += "".join(f"<th>{esc(x)}</th>" for x in headers) + "</tr></thead><tbody>"
    result += "".join("<tr>" + "".join(f"<td>{esc(x)}</td>" for x in row) + "</tr>" for row in body)
    result += "</tbody></table></div>"
    if ident:
        result += f'<div class="pager" data-pager-for="{esc(ident)}"><button type="button" data-page-prev>上一页</button><span data-page-status aria-live="polite"></span><button type="button" data-page-next>下一页</button></div>'
    return result


def card(label: str, value: str, note: str) -> str:
    return f'<div class="metric-card"><span>{esc(label)}</span><strong>{esc(value)}</strong><small>{esc(note)}</small></div>'


def section(ident: str, title: str, body: str, help_text: str) -> tuple[str, str, str]:
    button = f'<button type="button" class="help-button" aria-label="查看{esc(title)}的解释" aria-expanded="false" data-help-title="{esc(title)}" data-help-body="{esc(help_text)}">?</button>'
    return ident, title, f'<section class="report-section" id="{ident}"><h2>{esc(title)}{button}</h2>{body}</section>'


class Delivery:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.entries: list[dict[str, str]] = []
        self.figure_number = 0

    def add(self, path: Path, category: str, description: str, source: object) -> None:
        self.entries.append({"category": category, "description": description, "source": str(source), "target": str(path.relative_to(self.root)), "status": "up_to_date"})

    def copy_table(self, source: Path, category: str) -> Path:
        target = self.root / "tables" / category / source.name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        self.add(target, category, source.name, source)
        return target

    def figure(self, fig, name: str, category: str, title: str, caption: str, source: object) -> str:
        import matplotlib.pyplot as plt

        self.figure_number += 1
        prefix = self.root / "figures" / category / name
        prefix.parent.mkdir(parents=True, exist_ok=True)
        png, svg = prefix.with_suffix(".png"), prefix.with_suffix(".svg")
        fig.savefig(png, dpi=450, facecolor="white", bbox_inches="tight")
        fig.savefig(svg, facecolor="white", bbox_inches="tight")
        plt.close(fig)
        self.add(png, category, title + " PNG", source)
        self.add(svg, category, title + " SVG", source)
        return (f'<article class="figure-card expandable"><button type="button" class="zoomable figure-expand" aria-label="放大{esc(title)}"><img src="{data_uri(png)}" alt="{esc(title)}" loading="lazy"></button>'
                f'<figcaption><strong>图 {self.figure_number} · {esc(title)}</strong>{esc(caption)}</figcaption>'
                f'<div class="export-controls"><a download="{esc(png.name)}" href="{data_uri(png)}">下载 PNG</a><a download="{esc(svg.name)}" href="{data_uri(svg)}">下载 SVG</a></div>'
                f'<p class="source-note">数据来源：{esc(source)}</p></article>')


def load_evidence(root: Path, cfg: dict, project: Path) -> dict:
    if cfg.get("execution", {}).get("workflow") != "illumina_pe":
        raise ValueError("report mode requires an illumina_pe profile")
    stages = []
    for dirname in ["01_preflight", "02_raw_qc", "03_alignment_qc", "04_canonical_filter", "05_daf_signal", "06_fiberhmm", "07_delivery"]:
        stage = fields(root / dirname / "stage_complete.tsv")
        expected = "NOT_APPLICABLE" if dirname == "06_fiberhmm" else "COMPLETED"
        if stage["status"] != expected:
            raise ValueError(f"{dirname}: expected {expected}, got {stage['status']}")
        stages.append([dirname, stage["status"], str(root / dirname / "stage_complete.tsv")])
    sample_rows = rows(project / cfg["sample_sheet"])
    if len(sample_rows) != 1:
        raise ValueError("this report requires one paired-end sample")
    sample = sample_rows[0]["sample_id"]
    summary = rows(root / "05_daf_signal/summary.tsv")
    if len(summary) != 1 or summary[0]["sample_id"] != sample:
        raise ValueError("DAF summary/sample sheet mismatch")
    m = summary[0]
    audit = {x["metric"]: x["value"] for x in rows(root / "05_daf_signal/denominator_audit.tsv")}
    qc = fields(root / "02_raw_qc/paired_fastq_qc.tsv", "metric")
    flagtext = (root / "03_alignment_qc/alignment_flagstat.txt").read_text()
    def count(label: str) -> int:
        match = re.search(rf"^(\d+) \+ (\d+) {re.escape(label)}(?: \(|$)", flagtext, re.M)
        if not match:
            raise ValueError(f"missing flagstat metric: {label}")
        return int(match[1]) + int(match[2])
    flag = {k: count(k) for k in ["primary", "primary mapped", "supplementary", "properly paired", "singletons", "with itself and mate mapped"]}
    n = lambda k: int(m[k])
    if flag["primary"] != 2 * int(qc["read_pairs"]):
        raise ValueError("FASTQ pairs and primary BAM records do not reconcile")
    if sum(n(k) for k in ["ct_reads", "ga_reads", "undetermined_reads", "none_reads"]) != n("eligible_reads"):
        raise ValueError("DAF read classes do not reconcile")
    if n("ct_reads") + n("ga_reads") != n("classified_reads"):
        raise ValueError("classified reads do not reconcile")
    if n("da_bp") != n("ct_da_bp") + n("ga_da_bp") or n("tot_bp") != n("da_bp") + n("ct_ref_bp") + n("ga_ref_bp"):
        raise ValueError("DAF denominator does not reconcile")
    for key in ["da_bp", "tot_bp", "ct_da_bp", "ga_da_bp"]:
        if audit[key] != m[key]:
            raise ValueError(f"denominator audit differs from summary: {key}")
    if abs(100 * n("da_bp") / n("tot_bp") - float(m["base_conversion_rate_pct"])) > 0.0001:
        raise ValueError("DAF displayed rate differs from exact denominator")
    chrom_events = defaultdict(Counter)
    event_rows = 0
    for row in iter_site_rows(root / "05_daf_signal/per_site.tsv"):
        if row["sample_id"] != sample:
            raise ValueError("site table/sample mismatch")
        chrom_events[row["chrom"]][row["conversion"]] += int(row["count"])
        event_rows += 1
    for conversion, key in [("CT", "ct_da_bp"), ("GA", "ga_da_bp")]:
        if sum(c[conversion] for c in chrom_events.values()) != n(key):
            raise ValueError(f"per-site selected {conversion} events differ from summary")
    peaks_path = root / "07_delivery/peaks" / f"{sample}_peaks.narrowPeak"
    peaks = []
    if peaks_path.is_file() and peaks_path.stat().st_size:
        with peaks_path.open() as handle:
            for line in csv.reader(handle, delimiter="\t"):
                if len(line) != 10:
                    raise ValueError("narrowPeak row must contain ten fields")
                peaks.append({"chrom": line[0], "start": int(line[1]), "end": int(line[2]), "name": line[3], "signal": float(line[6]), "p": float(line[7]), "q": float(line[8]), "summit": int(line[9])})
    canonical = list(cfg["chromosomes"]["canonical"])
    order = {name: i for i, name in enumerate(canonical)}
    peaks.sort(key=lambda p: (-p["q"], -p["signal"], order.get(p["chrom"], 999), p["start"]))
    tss_path = root / "07_delivery/tss/TSS_SUMMARY.tsv"
    tss = rows(tss_path)[0] if tss_path.is_file() else None
    if tss:
        stages.append(["07_delivery/tss", "COMPLETED", str(tss_path)])
    elif cfg.get("execution", {}).get("tss", False):
        stages.append(["07_delivery/tss", "NOT_RUN", str(tss_path.parent)])
    return dict(sample=sample, sample_row=sample_rows[0], qc=qc, flag=flag, m=m, audit=audit, stages=stages, chrom_events=chrom_events, event_rows=event_rows, peaks=peaks, peaks_path=peaks_path, canonical=canonical, tss=tss)


def iter_site_rows(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def build_report(cfg: dict, project: Path, root: Path, output: Path, evidence: dict) -> Delivery:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#667085", "axes.labelcolor": "#344054", "text.color": "#344054", "xtick.color": "#475467", "ytick.color": "#475467", "svg.fonttype": "none"})
    colors = {"CT": "#E64B35", "GA": "#4DBBD5", "Undetermined": "#00A087", "None": "#8491B4"}
    delivery = Delivery(output.parent)
    m, qc, flag = evidence["m"], evidence["qc"], evidence["flag"]
    n = lambda key: int(m[key])
    fmt = lambda value: f"{int(value):,}"
    primary_rate = 100 * flag["primary mapped"] / flag["primary"]
    category_qc, category_daf, category_peak = "02_qc", "03_daf_signal", "04_coverage_peaks"
    sections = []

    overview = '<p class="chapter-conclusion">本报告分析一个 Illumina 双端样本，展示技术质控、read 层面的 CT/GA 转换及探索性覆盖峰。</p><div class="metric-grid">'
    overview += card("原始 read pairs", fmt(qc["read_pairs"]), "R1 与 R2 是同一样本的配对文件")
    overview += card("Primary 比对率", f"{primary_rate:.2f}%", f"{fmt(flag['primary mapped'])} / {fmt(flag['primary'])} reads")
    overview += card("DAF 分类比例", f"{float(m['classified_read_rate_pct']):.4f}%", "CT + GA reads / 全部 eligible reads")
    overview += card("已分类 reads 转换比例", f"{float(m['base_conversion_rate_pct']):.4f}%", "selected events / classified callable bases")
    overview += '</div><div class="finding">比对及后续统计已完成；DAF 指标属于 MD-based read-level approximation。相关分母与筛选过程见第 2、3 节。</div>'
    if evidence["peaks"]:
        overview += f'<div class="finding">得到 {fmt(len(evidence["peaks"]))} 条探索性峰；峰分布与 CPM 覆盖度见第 4 节。</div>'
    overview += '<div class="notice">本批次无独立重复、匹配 input 对照或平台校准。计算完成不等于实验或生物学验证；双端 reads 未合并为独立分子。</div>'
    sections.append(section("s1", "1 项目概览", overview, "样本：一个 Illumina PE 样本。\n证据：已完成的 QC、DAF 表格和峰文件。\n结论用于技术审阅与探索性描述，不支持组间差异、绝对脱氨效率或染色质机制推断。"))

    qc_headers = ["指标", "结果", "分母 / 说明"]
    qc_body = [
        ["样本", evidence["sample"], "一个样本，两个 mate 文件"],
        ["原始 read pairs", fmt(qc["read_pairs"]), "配对 FASTQ 记录"],
        ["R1 / R2 平均长度", f"{float(qc['r1_mean_length']):g} / {float(qc['r2_mean_length']):g} bp", "对应 FASTQ 的 reads"],
        ["总碱基数", fmt(int(qc["r1_bases"]) + int(qc["r2_bases"])), "R1 + R2"],
        ["Q20 / Q30", "未提供", "当前上游 QC 未输出碱基质量分布"],
        ["Primary reads", fmt(flag["primary"]), "排除 secondary 与 supplementary"],
        ["Primary mapped reads", fmt(flag["primary mapped"]), f"{primary_rate:.2f}% / primary reads"],
        ["Properly paired reads", fmt(flag["properly paired"]), f"{100 * flag['properly paired'] / flag['primary']:.2f}% / primary reads；单位为 reads"],
        ["Singleton reads", fmt(flag["singletons"]), f"{100 * flag['singletons'] / flag['primary']:.2f}% / primary reads"],
        ["Supplementary records", fmt(flag["supplementary"]), "额外比对记录，不计入 primary reads 分母"],
        ["Canonical / MAPQ ≥20 reads", fmt(n("eligible_reads")), f"{100 * n('eligible_reads') / flag['primary']:.2f}% / 原始 primary reads"],
        ["重复率 / 分子饱和度", "未评估", "BAM 未标记或去除 duplicates；flagstat 的零值不代表无重复"],
    ]
    qc_path = output.parent / "tables/02_qc/alignment_qc.tsv"
    write_table(qc_path, qc_headers, qc_body)
    delivery.add(qc_path, category_qc, "原始与比对 QC 汇总", root / "03_alignment_qc/alignment_flagstat.txt")
    delivery.copy_table(root / "02_raw_qc/paired_fastq_qc.tsv", category_qc)
    fig, ax = plt.subplots(figsize=(8.0, 4.5), layout="constrained")
    counts = [flag["primary"], flag["primary mapped"], n("eligible_reads"), n("classified_reads")]
    labels = ["Input primary reads", "Primary mapped", "Canonical + MAPQ >=20", "CT/GA classified"]
    bars = ax.barh(labels[::-1], np.array(counts[::-1]) / 1e6, color=["#00A087", "#3C5488", "#4DBBD5", "#8491B4"])
    for b, value in zip(bars, counts[::-1]): ax.text(b.get_width() + .4, b.get_y() + b.get_height()/2, fmt(value), va="center", fontsize=9)
    ax.set(xlabel="Reads (millions)", xlim=(0, max(counts) / 1e6 * 1.28), title="Read support through the analysis")
    qfigure = delivery.figure(fig, "read_support", category_qc, "Reads 支持与筛选流程", "单位为 reads。原始 reads 为 read pairs ×2；canonical、MAPQ 和主比对过滤后的数量汇总展示，不拆分无法从现有表格确定的排除原因。CT/GA 分类属于后续事件选择。", "paired_fastq_qc.tsv；alignment_flagstat.txt；summary.tsv")
    sections.append(section("s2", "2 数据与比对质控", table(qc_headers, qc_body) + qfigure + '<p class="source-note">参考：GRCh38 / refdata-gex-GRCh38-2024-A。保留 chr1–22、chrX、chrY；排除 chrM、非 canonical contigs、未比对、secondary、supplementary 及 MAPQ &lt;20 的记录。</p>', "Primary mapping rate 以原始 primary reads 为分母。\nDAF summary 的 100% mapping 是已过滤 BAM 内的比率，不能替代原始比对率。\n未做 BAM 去重，无法从 duplicate flags 估计真实重复率。"))

    for name in ["summary.tsv", "mismatch_spectrum.tsv", "read_attrition.tsv", "denominator_audit.tsv"]:
        delivery.copy_table(root / "05_daf_signal" / name, category_daf)
    class_labels = ["CT", "GA", "Undetermined", "None"]
    class_values = [n(k) for k in ["ct_reads", "ga_reads", "undetermined_reads", "none_reads"]]
    fig, ax = plt.subplots(figsize=(7.2, 4.6), layout="constrained")
    bars = ax.bar(class_labels, class_values, color=[colors[x] for x in class_labels], width=.65)
    ax.set_yscale("log")
    ax.set(ylim=(1e4, max(class_values)*3), ylabel="Reads (log scale)", title="Eligible-read conversion classification")
    for b, value in zip(bars, class_values): ax.text(b.get_x()+b.get_width()/2, max(1e4,value)*1.12, fmt(value), ha="center", fontsize=9)
    class_fig = delivery.figure(fig, "read_classes", category_daf, "DAF read 分类", f"分母为 {fmt(n('eligible_reads'))} 条 eligible reads。CT/GA 达到 0.90 方向分类阈值；Undetermined 含转换但未达到阈值，None 没有选定转换。纵轴为对数尺度，柱顶标出原始计数。", "05_daf_signal/summary.tsv")
    spectrum = rows(root / "05_daf_signal/mismatch_spectrum.tsv")
    if sum(int(x["count"]) for x in spectrum) != n("total_mismatches"):
        raise ValueError("mismatch spectrum total differs from summary")
    fig, ax = plt.subplots(figsize=(7.2, 4.6), layout="constrained")
    ax.bar([x["conversion"] for x in spectrum], [float(x["rate_per_10k"]) for x in spectrum], color=[colors.get(x["conversion"], "#8491B4") for x in spectrum], width=.7)
    ax.set(ylabel="Mismatches / 10,000 eligible reference bases", title="Directional mismatch spectrum")
    mismatch_fig = delivery.figure(fig, "mismatch_spectrum", category_daf, "12 类方向错配谱", "每类错配率使用对应参考碱基的 eligible observations 作分母，单位为每万碱基。CT/GA 用固定颜色突出；本图与已分类 reads 内的转换比例使用不同分母。", "05_daf_signal/mismatch_spectrum.tsv")
    denom_headers = ["指标", "数值", "解释"]
    denom_body = [
        ["Eligible reads", fmt(n("eligible_reads")), "canonical / MAPQ ≥20，MD 与序列可用"],
        ["含 CT/GA 转换的 reads", fmt(n("reads_with_conversion")), f"{float(m['conversion_positive_rate_pct']):.4f}% / eligible reads"],
        ["CT / GA classified reads", f"{fmt(n('ct_reads'))} / {fmt(n('ga_reads'))}", "read-level 方向分类"],
        ["Undetermined / None", f"{fmt(n('undetermined_reads'))} / {fmt(n('none_reads'))}", "未分类 / 无选定转换"],
        ["全部候选 CT + GA mismatch events", fmt(n("total_conversions")), f"{fmt(n('ct_count'))} CT + {fmt(n('ga_count'))} GA；包括未分配为对应方向的 events"],
        ["Selected CT / GA events", f"{fmt(n('ct_da_bp'))} / {fmt(n('ga_da_bp'))}", "仅对应已分类 reads 内保留的转换"],
        ["Selected event 总数", fmt(n("da_bp")), f"{fmt(n('ct_da_bp'))} + {fmt(n('ga_da_bp'))}"],
        ["CT / GA callable opportunities", f"{fmt(n('ct_da_bp') + n('ct_ref_bp'))} / {fmt(n('ga_da_bp') + n('ga_ref_bp'))}", "转换事件 + 未转换的 callable 对应参考碱基"],
        ["总 callable opportunities", fmt(n("tot_bp")), "仅已分类 reads 的 C/G opportunities"],
        ["Selected conversion proportion", f"{float(m['base_conversion_rate_pct']):.4f}%", f"{fmt(n('da_bp'))} / {fmt(n('tot_bp'))} ×100"],
        ["CT / GA conversion proportion", f"{float(m['ct_conversion_rate_pct']):.4f}% / {float(m['ga_conversion_rate_pct']):.4f}%", "分别使用对应已分类 reads 的 callable 分母"],
        ["事件键数", fmt(evidence["event_rows"]), "按 chrom/pos/ref/alt 合并的记录；非独立分子数，不具有 site coverage 分母"],
    ]
    chromosome_rows = [[c, evidence["chrom_events"][c]["CT"], evidence["chrom_events"][c]["GA"]] for c in evidence["canonical"]]
    cpath = output.parent / "tables/03_daf_signal/chromosome_selected_events.tsv"
    write_table(cpath, ["chrom", "selected_ct_events", "selected_ga_events"], chromosome_rows)
    delivery.add(cpath, category_daf, "按染色体汇总的 selected events", root / "05_daf_signal/per_site.tsv")
    fig, ax = plt.subplots(figsize=(12, 4.0), layout="constrained")
    x = np.arange(len(chromosome_rows))
    ax.bar(x-.19, [r[1] for r in chromosome_rows], .38, color=colors["CT"], label="CT")
    ax.bar(x+.19, [r[2] for r in chromosome_rows], .38, color=colors["GA"], label="GA")
    ax.set_xticks(x, [r[0].removeprefix("chr") for r in chromosome_rows])
    ax.set(xlabel="Chromosome", ylabel="Selected conversion events", title="Selected CT/GA events by chromosome")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v/1000:g}k"))
    ax.legend(frameon=False, loc="upper right")
    chrom_fig = delivery.figure(fig, "chromosome_events", category_daf, "各染色体 selected CT/GA 事件", f"对 per_site.tsv 的 count 求和；全图 CT 与 GA 分别汇总为 {fmt(n('ct_da_bp'))} 与 {fmt(n('ga_da_bp'))}。未按染色体覆盖深度或长度归一化，不能据此比较染色体的转换效率或可及性。", "05_daf_signal/per_site.tsv")
    dbody = f'<p class="chapter-conclusion">CT/GA 分类 reads 占 eligible reads 的 {float(m["classified_read_rate_pct"]):.4f}%；已分类 reads 内 selected conversion proportion 为 {float(m["base_conversion_rate_pct"]):.4f}%。这是不同层面的两个指标。</p>'
    dbody += table(denom_headers, denom_body) + '<div class="figure-grid">' + class_fig + mismatch_fig + '</div>' + chrom_fig
    dbody += '<div class="notice">MD-based 统计未完成平台特异性校准。存在转换方向的错配不直接证明酶效率或染色质开放；没有位点 callable coverage 分母，不能由事件计数构造位点转换率。</div>'
    spectrum_body = [[x["conversion"], fmt(x["count"]), fmt(x["eligible_ref_bases"]), x["rate_per_10k"]] for x in spectrum]
    dbody += '<details><summary>展开错配谱原始指标</summary>' + table(["转换", "事件数", "Eligible 参考碱基", "每万碱基事件数"], spectrum_body) + '</details>'
    sections.append(section("s3", "3 DAF 转换信号", dbody, "read 分类比例：CT + GA reads / eligible reads。\n已分类转换比例：selected CT + GA events / corresponding classified callable opportunities。\n错配谱：全部 eligible reads 中，对应参考碱基上的各类 mismatch。\n所有分母均是 read/base observations，重叠 mates 未合并。"))

    peaks = evidence["peaks"]
    bigwig = root / "07_delivery" / f"{evidence['sample']}.CPM.bw"
    peakbody = '<div class="notice">CPM bigWig 表示过滤后 BAM 的比对覆盖度，未单独筛选转换事件。MACS3 使用 BAMPE、human effective genome size、q=0.01 且无 input 对照；峰属于探索性计算结果。</div>'
    if peaks:
        delivery.copy_table(evidence["peaks_path"], category_peak)
        pheaders = ["染色体", "起点 (0-based)", "终点 (exclusive)", "长度 (bp)", "峰名称", "signalValue", "−log10(q)"]
        pbody = [[p["chrom"], p["start"], p["end"], p["end"]-p["start"], p["name"], p["signal"], p["q"]] for p in peaks]
        ppath = output.parent / "tables/04_coverage_peaks/peaks_ranked.tsv"
        write_table(ppath, pheaders, pbody)
        delivery.add(ppath, category_peak, "全部峰，按 q score 排序", evidence["peaks_path"])
        chr_counts = Counter(p["chrom"] for p in peaks)
        fig, ax = plt.subplots(figsize=(7.2, 4.4), layout="constrained")
        ax.bar([c.removeprefix("chr") for c in evidence["canonical"]], [chr_counts[c] for c in evidence["canonical"]], color="#3C5488")
        ax.set(xlabel="Chromosome", ylabel="Peak count", title=f"Exploratory peak distribution (n={len(peaks)})")
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        fig_counts = delivery.figure(fig, "peak_chromosomes", category_peak, "探索性峰的染色体分布", f"统计 narrowPeak 中全部 {len(peaks)} 条记录。未按染色体长度或可比对区域归一化。", "07_delivery/peaks/*_peaks.narrowPeak")
        widths = [p["end"]-p["start"] for p in peaks]
        fig, ax = plt.subplots(figsize=(7.2, 4.4), layout="constrained")
        ax.hist(widths, bins=24, color="#00A087", edgecolor="white", linewidth=.5)
        ax.set(xlabel="Peak width (bp)", ylabel="Peak count", title=f"Peak widths (median={np.median(widths):g} bp)")
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        fig_width = delivery.figure(fig, "peak_widths", category_peak, "探索性峰长度分布", "长度 = BED end − start。该分布为覆盖峰宽度，不解释为核小体、MSP 或单分子染色质结构长度。", "07_delivery/peaks/*_peaks.narrowPeak")
        peakbody += '<div class="figure-grid">' + fig_counts + fig_width + '</div>'
        if bigwig.is_file() and bigwig.stat().st_size:
            import pyBigWig
            selected = []
            panels = []
            track_rows = []
            with pyBigWig.open(str(bigwig)) as bw:
                for p in peaks:
                    summit = p["start"] + p["summit"] if p["summit"] >= 0 else (p["start"] + p["end"]) // 2
                    start, end = max(0, summit-2000), min(bw.chroms(p["chrom"]), summit+2000)
                    if any(q["chrom"] == p["chrom"] and start < b and end > a for q, a, b in selected):
                        continue
                    selected.append((p, start, end))
                    if len(selected) == 3:
                        break
                for i, (p, start, end) in enumerate(selected, 1):
                    values = bw.stats(p["chrom"], start, end, type="mean", nBins=100, exact=True)
                    points = np.array([(start + (j+.5)*(end-start)/100) + 1 for j in range(100)])
                    y = np.array([v if v is not None else np.nan for v in values])
                    for j, v in enumerate(values): track_rows.append([p["name"], p["chrom"], start+j*(end-start)/100, start+(j+1)*(end-start)/100, "NA" if v is None else v])
                    fig, ax = plt.subplots(figsize=(12, 3.5), layout="constrained")
                    ax.fill_between(points, y, color="#3C5488", alpha=.2)
                    ax.plot(points, y, color="#3C5488", linewidth=1.3)
                    ax.axvspan(p["start"]+1, p["end"], color="#00A087", alpha=.16, label="Exploratory peak")
                    ax.set(xlabel=f"{p['chrom']} position (1-based)", ylabel="CPM coverage", title=f"{p['name']} | {p['chrom']}:{start+1:,}-{end:,}")
                    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
                    ax.legend(frameon=False, loc="upper right")
                    figure_html = delivery.figure(fig, f"coverage_region_{i}", category_peak, f"覆盖度区域 {i} · {p['chrom']}:{start+1:,}–{end:,}", "按 −log10(q) 排序、排除重叠展示窗口后的前三个区域；峰 summit 两侧各 2 kb，染色体边缘截断。每个窗口分成 100 个等宽 bin，绘制 bigWig 的精确均值；缺失值保留 NA。绿色区域标示覆盖峰，展示区域不代表生物学功能位点。", "07_delivery/*.CPM.bw；peaks/*_peaks.narrowPeak")
                    panels.append((f"track-{i}", p["name"], figure_html))
            track_path = output.parent / "tables/04_coverage_peaks/representative_coverage.tsv"
            write_table(track_path, ["peak", "chrom", "bin_start_0based", "bin_end_exclusive", "mean_cpm"], track_rows)
            delivery.add(track_path, category_peak, "展示区域的 CPM bin 均值", bigwig)
            peakbody += '<h3>代表区域覆盖度</h3><div class="track-picker"><label for="track-select">选择区域 </label><select id="track-select">' + ''.join(f'<option value="{ident}">{esc(name)}</option>' for ident, name, _ in panels) + '</select></div>'
            peakbody += ''.join(f'<div class="track-panel" id="{ident}"' + (' hidden' if i else '') + f'>{body}</div>' for i, (ident, _, body) in enumerate(panels))
        else:
            peakbody += '<p class="notice">未发现可用 CPM bigWig，本报告省略区域覆盖轨道。</p>'
        options = ''.join(f'<option value="{esc(c)}">{esc(c)}</option>' for c in evidence["canonical"] if chr_counts[c])
        peakbody += f'<h3>全部 {len(peaks)} 条峰</h3><p>表格坐标采用 BED 的 0-based、half-open 区间。q score 为调用器模型评分，不作为重复样本层面的显著性证据。</p><div class="filters"><label>染色体 <select data-chrom-for="peaks-table"><option value="">全部</option>{options}</select></label><label>搜索 <input type="search" data-search-for="peaks-table" placeholder="峰名称或坐标"></label></div>'
        peakbody += table(pheaders, pbody, "peaks-table", int(cfg.get("report", {}).get("page_size", 10)))
    else:
        peakbody += '<p>当前没有可报告的峰记录；省略峰图、峰表和代表区域选择。</p>'
    sections.append(section("s4", "4 覆盖度与探索性峰", peakbody, "覆盖度来自 filtered BAM，不是 DAF conversion 专属信号。\nMACS3 默认相同片段过滤仅作用于 peak calling；保留 BAM 没有 deduplication。\n无 input 对照与独立重复，因此峰属于探索性结果。\n表格 BED 坐标与轨道 1-based 坐标已分别标注。"))

    alignment_marker = fields(root / "03_alignment_qc/stage_complete.tsv")
    mapper = alignment_marker.get("aligner", cfg["runtime"].get("aligner", "unknown"))
    mapper_note = ""
    if mapper == "minimap2":
        mapper_note = f"使用 minimap2 -ax sr --MD；预建索引 k={alignment_marker.get('prebuilt_index_k', '未记录')}、w={alignment_marker.get('prebuilt_index_w', '未记录')}，索引参数覆盖 preset 的相应参数。"
    method_body = '<div class="notice">当前结果支持技术质控与探索性 CT/GA 描述。没有独立生物学重复、匹配对照、平台校准或 SNP/重复区域掩码证据；不报告差异、绝对酶效率或生物学显著性。</div>'
    method_body += '<details><summary>分析方法与分母</summary><p>' + esc(mapper_note) + '</p><p>过滤 canonical nuclear chromosomes，primary mapped reads，MAPQ ≥20；DAF base quality ≥20，CT/GA strand cutoff 0.90。Selected conversion proportion = (ct_da_bp + ga_da_bp) / (ct_da_bp + ct_ref_bp + ga_da_bp + ga_ref_bp)。</p><p>Reads 及碱基 observations 是统计单位；重叠 mates 未合并。未执行剪切、BAM 去重、FiberHMM、TSS、差异分析。MACS3 的冗余片段处理不改变保留 BAM 或 CPM 轨道。</p></details>'
    versions = []
    for name, command in [("minimap2", [cfg["runtime"]["minimap2"], "--version"]), ("samtools", [cfg["runtime"]["samtools"], "--version"]), ("MACS3", [cfg["runtime"]["macs3"], "--version"]), ("bamCoverage", [cfg["runtime"]["bamcoverage"], "--version"])]:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        text = (result.stdout or result.stderr).strip().splitlines()
        versions.append([name, text[0] if result.returncode == 0 and text else "未获取", "报告生成时查询的安装版本"])
    versions.append(["报告 Python", sys.version.split()[0], "仅报告生成环境；不代替上游运行版本证据"])
    vpath = output.parent / "tables/05_methods/software_versions.tsv"
    write_table(vpath, ["software", "version", "scope"], versions)
    delivery.add(vpath, "05_methods", "报告时环境版本", "configured runtime executables")
    method_body += '<details><summary>参考、软件与运行信息</summary>' + table(["软件", "当前版本", "来源范围"], versions) + '<p>参考 FASTA：<code>' + esc(cfg["reference"]["fasta"]) + '</code></p><p>参考 GTF：<code>' + esc(cfg["reference"]["gtf"]) + '</code>；本次报告没有重新执行基因注释。</p><p>分析日志：<code>' + esc(project / cfg["paths"].get("logs_root", "logs") / cfg["run_id"]) + '</code></p><p>渲染软件：matplotlib ' + esc(matplotlib.__version__) + '；数据读取自已完成产物。图表采用固定 CT/GA 配色，图源可导出。</p></details>'
    if evidence["tss"]:
        tss = evidence["tss"]
        tss_dir = root / "07_delivery/tss"
        category_tss = "05_tss"
        for filename in ["TSS_SUMMARY.tsv", "active_tss.bed", "all_fragments_profile.tsv", "daf_positive_profile.tsv"]:
            source = tss_dir / filename
            if source.is_file():
                if source.suffix == ".tsv":
                    delivery.copy_table(source, category_tss)
                else:
                    target = output.parent / "tables" / category_tss / source.name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source, target)
                    delivery.add(target, category_tss, "峰区活跃 TSS BED", source)
        figures = []
        for key, title, description in [("all_plot_png", "全部配对片段 TSS 覆盖", "活跃 TSS 窗口内的全部合格 paired fragments。"), ("daf_positive_plot_png", "DAF 阳性片段 TSS 覆盖", "至少一个 mate 满足 CT/GA 阳性分类的局部 paired fragments。")]:
            source = root / "07_delivery" / tss[key]
            if source.is_file():
                target = output.parent / "figures" / category_tss / source.name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
                delivery.add(target, category_tss, title, source)
                uri = data_uri(target)
                figures.append(f'<article class="figure-card expandable"><button type="button" class="zoomable figure-expand" aria-label="放大{title}"><img src="{uri}" alt="{title}" loading="lazy"></button><figcaption><strong>{title}</strong>{description}</figcaption><div class="export-controls"><a download="{esc(target.name)}" href="{uri}">下载 PNG</a></div><p class="source-note">数据来源：07_delivery/tss</p></article>')
        tss_table = [["峰区活跃 TSS", f"{int(tss['active_tss_count']):,}", "蛋白编码基因 TSS 距 MACS3 峰区 ≤2 kb"], ["局部全部片段", f"{int(tss['all_fragments']):,}", "与活跃 TSS ±窗口相交"], ["局部 DAF 阳性片段", f"{int(tss['daf_positive_fragments']):,}", "paired fragment 层面的描述性分类"], ["窗口 / profile", f"±{tss['window_bp']} bp / {tss['profile_mode']}", "coverage；按 TSS 链方向对齐"], ["中心覆盖", f"{int(tss['center_count_all']):,} / {int(tss['center_count_daf_positive']):,}", "全部片段 / DAF 阳性片段"]]
        tss_body = '<p class="chapter-conclusion">在 MACS3 峰区两侧 2 kb 内筛选蛋白编码基因 TSS，并计算其周围 paired-fragment coverage。</p>' + table(["指标", "结果", "口径"], tss_table) + '<div class="figure-grid">' + ''.join(figures) + '</div><p class="notice">活跃 TSS 定义为距 MACS3 峰区不超过 2 kb；只统计活跃 TSS ±窗口内片段。曲线分别按各自外侧 flank 均值归一化，不能比较为组间效应。单样本无随机背景、重复或置信区间，结果仅作描述性探索。</p>'
    else:
        tss_body = '<p class="notice">本批次未执行 TSS 分析。</p>'
    sections.append(section("s5", "5 TSS 富集分析", tss_body, "TSS 来自 GRCh38 GTF 蛋白编码 gene feature；活跃位点距 MACS3 峰区不超过 2 kb。计数要求 proper pair、MAPQ≥20、fragment 长度 0–1000 bp，仅统计 ±窗口内片段。没有独立重复或背景模型，结果为描述性探索。"))
    method_body = method_body.replace("FiberHMM、TSS、差异分析", "FiberHMM、差异分析")
    sections.append(section("s6", "6 方法与解释边界", method_body, "Illumina 短读长不在本流程 FiberHMM 长读长适用范围。缺失模块不等于零信号；原始 QC 未提供 Q30，BAM 未去重。软件表为报告生成时查询，不伪装为历史运行版本记录。"))

    status_rows = evidence["stages"] + [["08_report", "GENERATED", str(output)]]
    status_path = output.parent / "workflow_status.tsv"
    write_table(status_path, ["stage", "status", "source"], status_rows)
    delivery.add(status_path, "07_delivery", "工作流执行状态", "stage_complete.tsv files")
    source_headers = ["文件类型", "产物位置", "说明"]
    source_body = [["最终 HTML", output.name, "全部图片内嵌，可离线浏览"], ["图源与小型表格", "figures/；tables/", "按报告章节组织；PNG / SVG / TSV"], ["交付清单与执行状态", "delivery_manifest.tsv；workflow_status.tsv", "最终 HTML 是清单最后一行"], ["过滤 BAM / BAI", str(root / "04_canonical_filter"), "保留在分析目录；未复制到报告目录"], ["覆盖度与峰", str(root / "07_delivery"), "CPM bigWig 与 MACS3 输出"], ["Selected event 大表", str(root / "05_daf_signal/per_site.tsv"), "保留完整事件明细，不将大表嵌入 HTML"]]
    delivery_body = '<p>本报告由项目内生成器重建，沿用既有 zoom 报告的版式与交互。所有可见指标来自本批次已完成的产物。</p>' + table(source_headers, source_body)
    delivery_body += '<details><summary>查看模块执行状态</summary>' + table(["模块", "状态"], [row[:2] for row in evidence["stages"]] + [["报告重建", "GENERATED；交付验证结果见 REPORT_VALIDATION.json"]]) + '</details>'
    sections.append(section("s7", "7 交付信息", delivery_body, "最终 HTML 自包含，图表和 CSS/JavaScript 无外部资源。大体积 BAM、bigWig 与 per-site 原始结果保留在分析目录。"))

    nav = ''.join(f'<a href="#{ident}">{esc(title)}</a>' for ident, title, _ in sections)
    generated = datetime.now().astimezone().isoformat(timespec="seconds")
    doc = '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>DAF-seq 分析报告 · ' + esc(evidence["sample"]) + '</title><style>' + CSS + '</style></head><body>'
    doc += '<nav class="sidebar"><h2>DAF-seq 分析报告</h2><p>Illumina PE · GRCh38<br>Exploratory read-level analysis</p>' + nav + '</nav><main class="content"><header class="hero"><div class="eyebrow">DAF-seq · Illumina paired-end</div><h1>DAF-seq 分析报告</h1><p>' + esc(evidence["sample"]) + ' · 一个样本 · R1/R2 配对 · GRCh38</p><p class="hero-meta">Run: ' + esc(cfg["run_id"]) + '<br>报告重建：' + esc(generated) + '</p></header>'
    doc += ''.join(body for _, _, body in sections)
    doc += '<footer class="footer">' + esc(cfg.get("report", {}).get("footer_text", "DAF-seq Analysis Report · Illumina PE")) + '<br>离线报告 · 生成器：pipeline/scripts/render_illumina_pe_report.py · ' + esc(generated) + '</footer></main>'
    doc += '<div id="lightbox" role="dialog" aria-modal="true" aria-label="图表放大" aria-hidden="true"><button type="button" class="lightbox-close" aria-label="关闭放大图">×</button><img alt="放大图"><span class="lightbox-note">滚轮：放大/缩小 · Esc：关闭</span></div><aside id="help-popover" class="help-popover" role="dialog" aria-label="指标帮助" aria-hidden="true"><button type="button" class="help-close" aria-label="关闭帮助">×</button><h3></h3><div class="help-body"></div></aside><script>' + JS + '</script></body></html>'
    md = output.parent / "DAFseq_analysis_results.md"
    if md.is_file():
        delivery.add(md, "07_delivery", "中文 Markdown 分析结果", "analysis result narrative")
    delivery.add(output, "07_delivery", "最终离线报告", "pipeline/scripts/render_illumina_pe_report.py")
    write_table(output.parent / "delivery_manifest.tsv", ["category", "description", "source", "target", "status"], [[entry[k] for k in ["category", "description", "source", "target", "status"]] for entry in delivery.entries])
    output.write_text(doc, encoding="utf-8")
    return delivery


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    args = parser.parse_args()
    project = args.project_root.resolve()
    cfg = yaml.safe_load(args.configfile.read_text(encoding="utf-8"))
    root = project / cfg["paths"].get("results_root", "results") / cfg["run_id"]
    output = root / "08_report" / f"{cfg['run_id']}_DAFseq_analysis_report.html"
    scratch = project / "tmp" / cfg["run_id"] / "report-rebuild"
    scratch.mkdir(parents=True, exist_ok=True)
    for key, leaf in [("TMPDIR", "temp"), ("TMP", "temp"), ("TEMP", "temp"), ("MPLCONFIGDIR", "matplotlib"), ("XDG_CACHE_HOME", "cache")]:
        path = scratch / leaf
        path.mkdir(exist_ok=True)
        os.environ[key] = str(path)
    evidence = load_evidence(root, cfg, project)
    validator = Path(cfg["report"]["validator"])
    if not validator.is_file():
        raise FileNotFoundError(f"report validator not found: {validator}")
    page_size = int(cfg.get("report", {}).get("page_size", 10))
    if page_size < 1:
        raise ValueError("report.page_size must be positive")
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        backup = scratch / (datetime.now().strftime("%Y%m%d_%H%M%S_%f") + "_previous_report.html")
        shutil.copy2(output, backup)
    # Do not retain an old success marker if rendering/validation fails.
    stage_marker = output.parent / "stage_complete.tsv"
    write_table(stage_marker, ["field", "value"], [["stage", "report"], ["sample_id", evidence["sample"]], ["status", "GENERATING"], ["report", str(output)]])
    build_report(cfg, project, root, output, evidence)
    command = [sys.executable, str(validator), "--report", str(output), "--manifest", str(output.parent / "delivery_manifest.tsv"), "--require-section", "1 项目概览", "--require-expandable-figures", "--require-print-css", "--forbid-text", "single_sample_report", "--json"]
    validation = subprocess.run(command, capture_output=True, text=True, check=False)
    (output.parent / "REPORT_VALIDATION.json").write_text(validation.stdout, encoding="utf-8")
    if validation.returncode:
        raise RuntimeError(f"report delivery validation failed: {validation.stdout}\n{validation.stderr}")
    statuses = evidence["stages"] + [["08_report", "COMPLETED", str(output)]]
    write_table(output.parent / "workflow_status.tsv", ["stage", "status", "source"], statuses)
    write_table(stage_marker, ["field", "value"], [["stage", "report"], ["sample_id", evidence["sample"]], ["status", "COMPLETED"], ["report", str(output)], ["validation", "PASS"], ["generator", "pipeline/scripts/render_illumina_pe_report.py"]])
    print(json.dumps({"report": str(output), "validation": "PASS", "figures": len(list((output.parent / 'figures').rglob('*.png'))), "peaks": len(evidence["peaks"])}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
