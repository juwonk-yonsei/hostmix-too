"""Read printed properties of the twelve figure PDFs for the A17 figure checklist.

Page size, font names, font sizes, stroke widths and the RGB colors set in the content
streams are read from the PDF bytes; nothing is redrawn. Stroke widths of 0 (hairline
clip paths) are ignored.
"""
from __future__ import annotations

import re
import zlib
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "manuscript/bmc/figures"
OUT = ROOT / "manuscript/checks/figure_pdf_props_A17.tsv"
STEMS = [f"Figure{i}" for i in range(1, 7)] + [f"FigureS{i}" for i in range(1, 7)]
METHOD = {"#0072B2": "Baseline", "#D55E00": "HostMix-TOO", "#CC79A7": "SCOPE", "#009E73": "CUP-AI-Dx"}
NUMBER = rb"(-?[\d.]+)"


def streams(data: bytes):
    for match in re.finditer(rb"stream\r?\n", data):
        end = data.find(b"endstream", match.end())
        chunk = data[match.end():end]
        try:
            yield zlib.decompressobj().decompress(chunk)
        except zlib.error:
            yield chunk


def hexcolor(parts) -> str:
    return "#" + "".join(f"{round(float(p) * 255):02X}" for p in parts)


def props(stem: str) -> dict:
    data = (FIG / f"{stem}.pdf").read_bytes()
    box = re.search(rb"/MediaBox\s*\[\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\]", data)
    width = (float(box.group(3)) - float(box.group(1))) / 72 * 25.4
    height = (float(box.group(4)) - float(box.group(2))) / 72 * 25.4
    fonts = sorted({name.decode().split("+")[-1] for name in re.findall(rb"/BaseFont\s*/([A-Za-z0-9+\-]+)", data)})
    sizes, widths, colors = [], [], set()
    for text in streams(data):
        sizes += [float(v) for v in re.findall(rb"/F\d+\s+([\d.]+)\s+Tf", text)]
        widths += [float(v) for v in re.findall(NUMBER + rb"\s+w\b", text)]
        for parts in re.findall(NUMBER + rb"\s+" + NUMBER + rb"\s+" + NUMBER + rb"\s+(?:rg|RG)\b", text):
            colors.add(hexcolor(parts))
    positive = [w for w in widths if w > 0]
    return {
        "figure": stem, "width_mm": round(width, 2), "height_mm": round(height, 2),
        "fonts": ";".join(fonts), "min_font_pt": min(sizes), "max_font_pt": max(sizes),
        "font_sizes_pt": ";".join(f"{s:g}" for s in sorted(set(sizes))),
        "min_line_pt": min(positive) if positive else None,
        "line_widths_pt": ";".join(f"{w:g}" for w in sorted(set(positive))),
        "method_colors": ";".join(name for code, name in METHOD.items() if code in colors),
        "other_colors": ";".join(sorted(c for c in colors if c not in METHOD)),
    }


def main() -> None:
    table = pd.DataFrame([props(stem) for stem in STEMS])
    table.to_csv(OUT, sep="\t", index=False)
    print(table.to_string())


if __name__ == "__main__":
    main()
