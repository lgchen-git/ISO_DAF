#!/usr/bin/env python3
"""Add display-only FiberHMM figure sizing and wheel zoom to an offline report."""

from __future__ import annotations

import argparse
from pathlib import Path


ZOOM_CSS = """
/* Display-only enhancement for FiberHMM/chromatin figures. */
#s3 .figure-grid { grid-template-columns: 1fr; gap: 22px; }
#s3 .compact-figure, #s3 .compact-figure.large { width: 100%; margin: 15px auto; padding: 14px; }
#s3 .compact-figure img, #s3 .compact-figure.large img { max-height: 920px; }
#lightbox { overflow: hidden; }
#lightbox img { max-width: 92vw; max-height: 90vh; transform-origin: center center; transition: transform .08s linear; cursor: zoom-in; user-select: none; touch-action: none; }
#lightbox.zoomed img { cursor: grab; }
#lightbox.zoomed img:active { cursor: grabbing; }
#lightbox::after { content: "滚轮：放大/缩小 · Esc：关闭"; position: absolute; left: 50%; bottom: 14px; transform: translateX(-50%); color: #fff; font-size: 12px; background: rgba(15,23,42,.72); border-radius: 999px; padding: 5px 10px; pointer-events: none; }
"""


ZOOM_JS = r"""
(()=>{
  const box=document.getElementById('lightbox');
  const img=box.querySelector('img');
  let scale=1;
  const MIN_SCALE=.5, MAX_SCALE=8, STEP=1.15;
  const reset=()=>{
    scale=1;
    img.style.transform='scale(1)';
    box.classList.remove('zoomed');
  };
  const close=()=>{
    box.classList.remove('open');
    box.setAttribute('aria-hidden','true');
    reset();
  };
  document.querySelectorAll('.zoomable').forEach(btn=>btn.addEventListener('click',e=>{
    e.stopPropagation();
    const source=btn.querySelector('img');
    img.src=btn.dataset.image;
    img.alt=source ? source.alt : '';
    box.classList.add('open');
    box.setAttribute('aria-hidden','false');
    reset();
  }));
  box.addEventListener('wheel',e=>{
    if(!box.classList.contains('open')) return;
    e.preventDefault();
    const factor=e.deltaY<0 ? STEP : 1/STEP;
    scale=Math.min(MAX_SCALE,Math.max(MIN_SCALE,scale*factor));
    img.style.transform=`scale(${scale})`;
    box.classList.toggle('zoomed',scale!==1);
  },{passive:false});
  img.addEventListener('click',e=>e.stopPropagation());
  box.addEventListener('click',e=>{if(e.target===box) close();});
  document.addEventListener('keydown',e=>{if(e.key==='Escape'&&box.classList.contains('open')) close();});
})();
"""


def patch_report(source: Path, target: Path) -> None:
    text = source.read_text(encoding="utf-8")
    if "#s3 .figure-grid" not in text:
        text = text.replace("</style>", ZOOM_CSS + "</style>", 1)
    old_js_start = "<script>(()=>{const box=document.getElementById('lightbox'),img=box.querySelector('img');"
    start = text.find(old_js_start)
    if start < 0:
        raise RuntimeError("existing lightbox script was not found")
    end = text.find("</script>", start)
    if end < 0:
        raise RuntimeError("existing lightbox script terminator was not found")
    replacement = "<script>" + ZOOM_JS.strip() + "</script>"
    text = text[:start] + replacement + text[end + len("</script>"):]
    text = text.replace('<div id="lightbox"><img alt="放大图"></div>', '<div id="lightbox" role="dialog" aria-modal="true" aria-hidden="true"><img alt="放大图"></div>', 1)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    patch_report(args.input, args.output)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
