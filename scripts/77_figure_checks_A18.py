"""A18 section 2.3 checks of the twelve redrawn figures; nothing is redrawn.

Writes to manuscript/checks/:
  figure_pdf_props_A18.tsv  page size, fonts, font sizes and stroke widths read from the PDF bytes
                            (reader of 66_figure_pdf_props_A17.py)
  figure_fonts_A18.tsv      pdffonts table per figure (type, embedded, subset, unicode map)
  figure_files_A18.tsv      file sizes, the 10 MB limit, and whether the PNG equals the committed
                            one (pixel count that differs, from git show <rev>:path)
Usage: 77_figure_checks_A18.py [rev]   (default HEAD)
"""
from __future__ import annotations

import importlib
import io
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
P66 = importlib.import_module("66_figure_pdf_props_A17")
ROOT = P66.ROOT
FIG = P66.FIG
CHECKS = ROOT / "manuscript/checks"
STEMS = P66.STEMS
LIMIT_MB = 10.0


def fonts(stem: str) -> list[dict]:
    out = subprocess.run(["pdffonts", str(FIG / f"{stem}.pdf")], capture_output=True, text=True, check=True).stdout
    rows = []
    for line in out.splitlines()[2:]:
        parts = line.split()
        if len(parts) < 7:
            continue
        name, emb, sub, uni = parts[0], parts[-5], parts[-4], parts[-3]
        rows.append({"figure": stem, "font": name.split("+")[-1], "type": " ".join(parts[1:-5]), "embedded": emb,
                     "subset": sub, "unicode_map": uni})
    return rows


def committed_png(rev: str, stem: str):
    path = f"manuscript/bmc/figures/{stem}.png"
    got = subprocess.run(["git", "show", f"{rev}:{path}"], cwd=ROOT, capture_output=True, check=True).stdout
    return np.asarray(Image.open(io.BytesIO(got)).convert("RGB"))


def files(rev: str, stem: str) -> dict:
    pdf, png = FIG / f"{stem}.pdf", FIG / f"{stem}.png"
    new = np.asarray(Image.open(png).convert("RGB"))
    old = committed_png(rev, stem)
    same_shape = new.shape == old.shape
    differing = int((new != old).any(axis=2).sum()) if same_shape else None
    return {"figure": stem, "pdf_mb": round(pdf.stat().st_size / 2**20, 3), "png_mb": round(png.stat().st_size / 2**20, 3),
            "within_10mb": pdf.stat().st_size / 2**20 <= LIMIT_MB and png.stat().st_size / 2**20 <= LIMIT_MB,
            "png_px_new": f"{new.shape[1]}x{new.shape[0]}", "png_px_committed": f"{old.shape[1]}x{old.shape[0]}",
            "differing_px": "size differs" if differing is None else differing,
            "png_equal_to_committed": differing == 0}


def main() -> None:
    rev = sys.argv[1] if len(sys.argv) > 1 else "HEAD"
    props = pd.DataFrame([P66.props(stem) for stem in STEMS])
    props.to_csv(CHECKS / "figure_pdf_props_A18.tsv", sep="\t", index=False)
    font_rows = pd.DataFrame([row for stem in STEMS for row in fonts(stem)])
    font_rows.to_csv(CHECKS / "figure_fonts_A18.tsv", sep="\t", index=False)
    file_rows = pd.DataFrame([files(rev, stem) for stem in STEMS])
    file_rows.to_csv(CHECKS / "figure_files_A18.tsv", sep="\t", index=False)
    print(props[["figure", "width_mm", "height_mm", "fonts", "min_font_pt", "min_line_pt", "line_widths_pt"]].to_string(index=False))
    print(font_rows.to_string(index=False))
    print(file_rows.to_string(index=False))
    print({"min_font_pt": float(props["min_font_pt"].min()), "min_line_pt": float(props["min_line_pt"].min()),
           "fonts_not_embedded": int((font_rows["embedded"] != "yes").sum()),
           "over_10mb": int((~file_rows["within_10mb"]).sum())})


if __name__ == "__main__":
    main()
