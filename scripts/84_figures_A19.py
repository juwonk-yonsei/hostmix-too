"""Draw Figures 1-6 and S1-S6 at their printed size from stored result tables (A19).

A19: 76_figures_A18.py with the main figures renumbered in order of first citation (old 6 -> 4, old 4 -> 5,
old 5 -> 6) by 83_renumber_figures_A19.py; drawing is unchanged. The environment variable A19_FIGURE_OUT moves
all outputs (figures, figure_source, checks) under another root, for comparison with the renamed files.

Every figure is laid out in millimetres on a 170 mm wide page and saved without
cropping. Before a figure is written, each drawn text is checked for font size
(at least 6 pt), position inside the page and inside its axes, and overlap with
other text; the saved PDF is then read back for its page size and font sizes.
Values are read from results/ by key. The simulation means, the six-cohort
auxiliary aggregate and the one-patient microarray layer reuse the loaders of
50_figures_A13.py unchanged.

Every printed value text (bar and cell values, n, arrow and interval texts, I2,
zero marks) is registered with an element key and its figure_source file, row
and column; after drawing, its extent is written to
manuscript/checks/figure_labels.tsv in millimetres from the top-left corner of
the page. Axis ticks are not registered.

Usage: 84_figures_A19.py [Figure1 FigureS3 ...]   (no argument: all figures,
the audit table and the label table)
"""
from __future__ import annotations

import importlib
import os
import re
import sys
import zlib
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Polygon, Rectangle
from matplotlib.text import Text
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
from matplotlib.transforms import blended_transform_factory
import numpy as np
import pandas as pd
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
A13 = importlib.import_module("50_figures_A13")

plt.rcParams.update({
    "font.family": "Liberation Sans",
    "font.size": 7,
    "axes.titlesize": 7.5,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlepad": 3,
    "axes.labelpad": 2.5,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "xtick.major.pad": 2,
    "ytick.major.pad": 2,
    "lines.linewidth": 1.0,
    "hatch.linewidth": 0.5,
    "legend.frameon": False,
    "legend.handlelength": 1.4,
    "legend.handleheight": 0.8,
    "legend.handletextpad": 0.5,
    "legend.columnspacing": 1.4,
    "legend.borderaxespad": 0.3,
    "legend.borderpad": 0.2,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
    "savefig.facecolor": "white",
})

ROOT = A13.ROOT
OUT = Path(os.environ.get("A19_FIGURE_OUT", str(ROOT)))
FIG = OUT / "manuscript/bmc/figures"
SRC = OUT / "manuscript/bmc/figure_source"
CHECKS = OUT / "manuscript/checks"
PACK = ROOT / "manuscript/bmc/review_pack"
MM = 1 / 25.4
WIDTH = 170.0
MAX_HEIGHT = 225.0
MIN_PT = 6.0

BLUE, ORANGE, PINK, GREEN = "#0072B2", "#D55E00", "#CC79A7", "#009E73"
LIGHT = "#E9E9E9"
GUIDE = "#E3E3E3"
ZERO_LINE = "#8C8C8C"
METHODS = ["BASE-Z", "SA-Z", "SCOPE", "CUP-AI-Dx"]
COLOR = {"BASE-Z": BLUE, "SA-Z": ORANGE, "SCOPE": PINK, "CUP-AI-Dx": GREEN}
MARKER = {"BASE-Z": "o", "SA-Z": "s"}
COHORTS3 = ["MET500", "POG570", "aux_rnaseq"]
COHORT3_NAME = {"MET500": "MET500", "POG570": "POG570", "aux_rnaseq": "Auxiliary RNA-seq"}
AUX = ["blca_iatlas_imvigor210_2017", "brca_iatlas_anders_2022", "mel_dfci_2019",
       "paad_iatlas_prince_2022", "GSE50760", "prad_su2c_2019"]
MICROARRAY = ["GSE41258", "GSE14018", "GSE71729", "GSE74685", "prad_fhcrc"]
COHORT = A13.COHORT
label = A13.label
AUDIT: list[dict] = []
LABELS: list[dict] = []
SHORT = {**{f"Figure{i}": f"Fig{i}" for i in range(1, 7)}, **{f"FigureS{i}": f"FigS{i}" for i in range(1, 7)}}


def fmt(value: float, digits: int = 3) -> str:
    number = Decimal(repr(float(value))).quantize(Decimal(1).scaleb(-digits), rounding=ROUND_HALF_UP)
    if number.is_zero():
        number = abs(number)
    return f"{number:.{digits}f}".replace("-", "\u2212")


def ticks(digits: int) -> FuncFormatter:
    return FuncFormatter(lambda value, _pos: "0" if abs(value) < 1e-9 else fmt(value, digits))


def read(path: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / path, sep="\t")


def one(frame: pd.DataFrame, **keys) -> pd.Series:
    mask = np.ones(len(frame), dtype=bool)
    for column, value in keys.items():
        mask &= frame[column].eq(value).to_numpy()
    hit = frame.loc[mask]
    if len(hit) != 1:
        raise SystemExit(f"expected one row for {keys}, found {len(hit)}")
    return hit.iloc[0]


def source(frame: pd.DataFrame, name: str) -> None:
    SRC.mkdir(parents=True, exist_ok=True)
    frame.to_csv(SRC / name, sep="\t", index=False, float_format="%.10g")


def pos(frame: pd.DataFrame, **keys) -> int:
    """Position (0-based data row) of the first row with these keys, as written by source()."""
    mask = np.ones(len(frame), dtype=bool)
    for column, value in keys.items():
        mask &= frame[column].eq(value).to_numpy()
    hits = np.flatnonzero(mask)
    if not len(hits):
        raise SystemExit(f"no row for {keys}")
    return int(hits[0])


class Canvas:
    """A figure of fixed printed size whose children are placed in millimetres."""

    def __init__(self, stem: str, height: float):
        if height > MAX_HEIGHT:
            raise SystemExit(f"{stem}: height {height} mm exceeds {MAX_HEIGHT} mm")
        self.stem = stem
        self.w = WIDTH
        self.h = float(height)
        self.fig = plt.figure(figsize=(self.w * MM, self.h * MM))
        self.geometry: dict = {}
        self.boxes: list = []
        self.marks: list = []

    def mark(self, item, panel: str, element: str, file: str, row, column: str) -> None:
        """Register a printed value text; item is a Text or a callable returning it after drawing."""
        key = f"{SHORT[self.stem]}/{panel}/{element}" if panel else f"{SHORT[self.stem]}/{element}"
        if any(m[1] == key for m in self.marks):
            raise SystemExit(f"duplicate label key {key}")
        self.marks.append((item, key, panel, file, row, column))

    def mark_ticks(self, ax, axis: str, panel: str, elements: list[str], file: str, rows: list, column: str) -> None:
        for i, (element, row) in enumerate(zip(elements, rows)):
            getter = (lambda i=i: ax.get_xticklabels()[i]) if axis == "x" else (lambda i=i: ax.get_yticklabels()[i])
            self.mark(getter, panel, element, file, row, column)

    def axes(self, left: float, top: float, width: float, height: float, **kwargs):
        ax = self.fig.add_axes((left / self.w, 1 - (top + height) / self.h, width / self.w, height / self.h), **kwargs)
        self.geometry[ax] = (left, top, width, height)
        return ax

    def text(self, x: float, y: float, text: str, **kwargs):
        kwargs.setdefault("ha", "left")
        kwargs.setdefault("va", "baseline")
        return self.fig.text(x / self.w, 1 - y / self.h, text, **kwargs)

    def letter(self, x: float, y: float, text: str):
        return self.text(x, y, text, fontsize=10, fontweight="bold")

    def row_text(self, ax, x: float, y: float, text: str, **kwargs):
        """Text at page position x (mm) on data row y of ax, outside the axes."""
        left, _top, width, _height = self.geometry[ax]
        kwargs.setdefault("va", "center")
        kwargs.setdefault("ha", "left")
        item = ax.text((x - left) / width, y, text, clip_on=False,
                       transform=blended_transform_factory(ax.transAxes, ax.transData), **kwargs)
        item.set_gid("free")
        return item

    def legend(self, x: float, y: float, handles, labels, ncol: int, **kwargs):
        """Legend whose upper centre sits at page position (x, y) in mm."""
        kwargs.setdefault("fontsize", 7)
        return self.fig.legend(handles, labels, loc="upper center", ncol=ncol,
                               bbox_to_anchor=(x / self.w, 1 - y / self.h), **kwargs)


def drawn_texts(fig) -> tuple[list, object]:
    """Text artists drawn in the figure, single- and multi-line (the renderer hook sees only single lines)."""
    renderer = fig.canvas.get_renderer()
    drawn: list = []
    original = Text.draw

    def capture(self, renderer):
        if self.get_visible() and str(self.get_text()).strip() and not any(self is item for item in drawn):
            drawn.append(self)
        return original(self, renderer)

    Text.draw = capture
    try:
        fig.draw(renderer)
    finally:
        Text.draw = original
    return drawn, renderer


def audit(canvas: Canvas) -> dict:
    fig = canvas.fig
    drawn, renderer = drawn_texts(fig)
    tol = 0.3
    width, height = fig.bbox.width, fig.bbox.height
    extents = [(item, item.get_window_extent(renderer)) for item in drawn]
    small = sorted({(round(item.get_fontsize(), 2), item.get_text()) for item, _ in extents if item.get_fontsize() < MIN_PT})
    outside_page = [item.get_text() for item, box in extents
                    if box.x0 < -tol or box.y0 < -tol or box.x1 > width + tol or box.y1 > height + tol]
    overlaps = []
    for i in range(len(extents)):
        for j in range(i + 1, len(extents)):
            a, b = extents[i][1], extents[j][1]
            if min(a.x1, b.x1) - max(a.x0, b.x0) > 0.5 and min(a.y1, b.y1) - max(a.y0, b.y0) > 0.5:
                overlaps.append(f"{extents[i][0].get_text()!r} x {extents[j][0].get_text()!r}")
    outside_axes = []
    frames = []
    for ax in fig.axes:
        frame = ax.get_window_extent(renderer)
        if ax.axison:
            frames.append((ax, frame))
        for item in ax.texts:
            if item.get_gid() == "free" or not item.get_visible() or not item.get_text().strip():
                continue
            box = item.get_window_extent(renderer)
            if box.x0 < frame.x0 - tol or box.x1 > frame.x1 + tol or box.y0 < frame.y0 - tol or box.y1 > frame.y1 + tol:
                outside_axes.append(item.get_text())
    into_axes = []
    for item, box in extents:
        own = getattr(item, "axes", None)
        if item.get_gid() != "free" and own is not None:
            continue
        for ax, frame in frames:
            if ax is own and item.get_gid() != "free":
                continue
            if min(box.x1, frame.x1) - max(box.x0, frame.x0) > 0.5 and min(box.y1, frame.y1) - max(box.y0, frame.y0) > 0.5:
                into_axes.append(item.get_text())
    tight_boxes = []
    margin = 0.8 * MM * fig.dpi
    for patch, texts in canvas.boxes:
        frame = patch.get_window_extent(renderer)
        for item in texts:
            box = item.get_window_extent(renderer)
            if (box.x0 < frame.x0 + margin or box.x1 > frame.x1 - margin
                    or box.y0 < frame.y0 + margin * 0.6 or box.y1 > frame.y1 - margin * 0.6):
                tight_boxes.append(item.get_text())
    return {
        "texts": len(extents),
        "min_pt_artists": min(item.get_fontsize() for item, _ in extents),
        "small": small,
        "outside_page": outside_page,
        "outside_axes": outside_axes,
        "text_on_axes": into_axes,
        "overlaps": overlaps,
        "tight_boxes": tight_boxes,
        "labels": label_rows(canvas, drawn, renderer),
    }


def label_rows(canvas: Canvas, drawn: list, renderer) -> list[dict]:
    fig = canvas.fig
    scale = 25.4 / fig.dpi
    rows = []
    for item, key, panel, file, row, column in canvas.marks:
        text = item() if callable(item) else item
        if not any(text is d for d in drawn):
            raise SystemExit(f"{key}: registered text was not drawn")
        box = text.get_window_extent(renderer)
        rows.append({"figure": canvas.stem, "panel": panel, "key": key,
                     "text": " ".join(text.get_text().split()),
                     "x0_mm": round(box.x0 * scale, 3), "y0_mm": round((fig.bbox.height - box.y1) * scale, 3),
                     "x1_mm": round(box.x1 * scale, 3), "y1_mm": round((fig.bbox.height - box.y0) * scale, 3),
                     "rotation": round(float(text.get_rotation()), 1), "fontsize_pt": text.get_fontsize(),
                     "source_file": f"manuscript/bmc/figure_source/{file}", "source_row": row,
                     "source_column": column})
    return rows


def pdf_fonts(path: Path) -> tuple[float, float, float, float]:
    data = path.read_bytes()
    box = re.search(rb"/MediaBox\s*\[\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\]", data)
    width = (float(box.group(3)) - float(box.group(1))) / 72 * 25.4
    height = (float(box.group(4)) - float(box.group(2))) / 72 * 25.4
    sizes: list[float] = []
    for match in re.finditer(rb"stream\r?\n", data):
        end = data.find(b"endstream", match.end())
        chunk = data[match.end():end]
        try:
            text = zlib.decompressobj().decompress(chunk)
        except zlib.error:
            text = chunk
        sizes += [float(value) for value in re.findall(rb"/F\d+\s+([\d.]+)\s+Tf", text)]
    return width, height, min(sizes), max(sizes)


def finish(canvas: Canvas, changes: str) -> None:
    result = audit(canvas)
    FIG.mkdir(parents=True, exist_ok=True)
    pdf = FIG / f"{canvas.stem}.pdf"
    png = FIG / f"{canvas.stem}.png"
    canvas.fig.savefig(pdf)
    canvas.fig.savefig(png, dpi=300)
    plt.close(canvas.fig)
    width, height, low, high = pdf_fonts(pdf)
    with Image.open(png) as image:
        pixels = image.size
    LABELS.extend(result["labels"])
    problems = [key for key in ("small", "outside_page", "outside_axes", "text_on_axes", "overlaps", "tight_boxes") if result[key]]
    if abs(width - WIDTH) > 0.05 or abs(height - canvas.h) > 0.05 or low < MIN_PT:
        problems.append("pdf")
    for key in problems:
        if key != "pdf":
            print(f"  {canvas.stem} {key}: {result[key][:8]}")
    AUDIT.append({
        "figure": canvas.stem, "width_mm": round(width, 2), "height_mm": round(height, 2),
        "png_px": f"{pixels[0]}x{pixels[1]}", "min_font_pt_pdf": low, "max_font_pt_pdf": high,
        "min_font_pt_artists": round(result["min_pt_artists"], 2), "texts": result["texts"],
        "overlaps": len(result["overlaps"]), "outside_page": len(result["outside_page"]),
        "outside_axes": len(result["outside_axes"]), "text_on_axes": len(result["text_on_axes"]),
        "tight_boxes": len(result["tight_boxes"]), "value_labels": len(result["labels"]),
        "pdf_mb": round(pdf.stat().st_size / 2**20, 3), "png_mb": round(png.stat().st_size / 2**20, 3),
        "changes": changes,
    })
    print(f"{canvas.stem}: {width:.1f} x {height:.1f} mm, fonts {low}-{high} pt, "
          f"{result['texts']} texts, problems {problems or 'none'}")


def method_handles(methods: list[str], kind: str = "bar") -> tuple[list, list[str]]:
    if kind == "bar":
        handles = [Patch(facecolor=COLOR[m], edgecolor="none") for m in methods]
    else:
        handles = [Line2D([], [], color=COLOR[m], marker=MARKER.get(m, "o"), ls="none", ms=4) for m in methods]
    return handles, [label(m) for m in methods]


def zero_mark(ax, x: float, value: float, text: str = "0"):
    return ax.annotate(text, (x, value), xytext=(0, 1.2), textcoords="offset points",
                       ha="center", va="bottom", fontsize=6.5, color="black")


# Figure 1 -----------------------------------------------------------------

def flow_box(canvas: Canvas, ax, x: float, y: float, width: float, lines: list[str], *,
             face: str = "white", edge: str = "black", lw: float = 0.6, ls="solid",
             bold: tuple = (), size: float = 7.0) -> float:
    leading = 3.15
    height = (len(lines) - 1) * leading + 2.9 + 2 * 1.9
    patch = FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0,rounding_size=1.0",
                           facecolor=face, edgecolor=edge, linewidth=lw, linestyle=ls,
                           capstyle="round", joinstyle="round", zorder=2)
    ax.add_patch(patch)
    middle = y + height / 2
    texts = []
    for i, line in enumerate(lines):
        offset = (i - (len(lines) - 1) / 2) * leading
        texts.append(ax.text(x + width / 2, middle + offset, line, ha="center", va="center",
                             fontsize=size, fontweight="bold" if i in bold else "normal", zorder=3))
    canvas.boxes.append((patch, texts))
    return y + height


def arrow(ax, x0: float, y0: float, x1: float, y1: float) -> None:
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=1,
                arrowprops=dict(arrowstyle="-|>", color="black", lw=0.8, shrinkA=0.8,
                                shrinkB=0.8, mutation_scale=9))


def fig1() -> None:
    c = Canvas("Figure1", 105)
    ax = c.axes(0, 0, c.w, c.h)
    ax.set_xlim(0, c.w)
    ax.set_ylim(c.h, 0)
    ax.axis("off")
    gap = 4.6

    c.letter(0.5, 4.8, "a")
    c.text(5.2, 4.8, "Example biopsy", fontsize=8)
    x, w = 2.0, 46.0
    y = flow_box(c, ax, x, 8.0, w, ["Liver biopsy:", "colorectal cancer cells", "+ hepatocytes"], face=LIGHT)
    arrow(ax, x + w / 2, y, x + w / 2, y + gap)
    y = flow_box(c, ax, x, y + gap, w, ["Classifier trained", "on primary tumors"])
    arrow(ax, x + w / 2, y, x + w / 2, y + gap)
    y = flow_box(c, ax, x, y + gap, w, ["Call: liver", "True primary: colorectal", "Host attraction"],
                 face=LIGHT, bold=(2,))
    ax.text(x, y + 4.0, "At risk: the true organ lies outside\nthe organs native to the biopsy site.",
            va="top", ha="left", fontsize=7, linespacing=1.3)
    ax.text(x, y + 12.5, "Host-attraction rate: the fraction of\nat-risk biopsies called a native organ.",
            va="top", ha="left", fontsize=7, linespacing=1.3)

    c.letter(52.5, 4.8, "b")
    c.text(57.2, 4.8, "HostMix-TOO training", fontsize=8)
    x, w = 54.0, 60.0
    half = (w - 2.0) / 2
    y_in = flow_box(c, ax, x, 8.0, half, ["TCGA-train:", "7486 primary tumors,", "32 labels"])
    flow_box(c, ax, x + half + 2.0, 8.0, half, ["GTEx-ref:", "1103 normal profiles,", "10 tissues"])
    arrow(ax, x + half / 2, y_in, x + half / 2, y_in + gap)
    arrow(ax, x + half * 1.5 + 2.0, y_in, x + half * 1.5 + 2.0, y_in + gap)
    y = flow_box(c, ax, x, y_in + gap, w, [
        "x = \u03c1t + (1 \u2212 \u03c1)h in linear TPM;",
        "t, tumor; h, host tissue",
        "\u03c1 ~ U(0.15, 1); per tumor 1 pure",
        "+ 4 mixtures (37430 profiles)",
    ], edge=ORANGE, lw=1.3)
    arrow(ax, x + w / 2, y, x + w / 2, y + gap)
    y = flow_box(c, ax, x, y + gap, w, ["Rank-normal scores, 5000 genes;", "standardized"])
    arrow(ax, x + w / 2, y, x + w / 2, y + gap)
    y = flow_box(c, ax, x, y + gap, w, ["Multinomial logistic regression (C = 0.03):",
                                        "32 labels \u2192 26 organs"])
    ax.text(x, y + 4.0, "Baseline: pure tumors only,\nunstandardized, C = 0.1",
            va="top", ha="left", fontsize=7, linespacing=1.3)

    c.letter(118.0, 4.8, "c")
    c.text(122.7, 4.8, "Study stages", fontsize=8)
    x, w = 119.5, 50.0
    thick, dashed, dotted = 1.6, (0, (3.5, 2.0)), (0, (0.1, 2.0))
    stages = [
        (["Development:", "simulation, MET500"], "white", 0.6, "solid", "development, thin solid"),
        (["Lock:", "models, code and plans committed"], LIGHT, 0.6, "solid", "lock, thin solid"),
        (["Confirmation 1:", "POG570 (H1\u2013H3)"], "white", thick, "solid", "preregistered, thick solid"),
        (["Confirmation 2:", "six RNA-seq cohorts", "(AH1\u2013AH3, AS1\u2013AS3)"], "white", thick, "solid",
         "preregistered, thick solid"),
        (["Published classifiers:", "SCOPE, CUP-AI-Dx"], "white", 1.0, dashed, "pre-specified descriptive, dashed"),
        (["Post hoc:", "ablation, error analyses"], "white", 1.1, dotted, "post hoc, dotted"),
    ]
    y = 8.0
    step = 3.6
    for i, (lines, face, lw, ls, _role) in enumerate(stages):
        if i:
            arrow(ax, x + w / 2, y, x + w / 2, y + step)
            y += step
        y = flow_box(c, ax, x, y, w, lines, face=face, lw=lw, ls=ls, bold=(0,))
    key_y = y + 4.2
    entries = [
        ((x, key_y), thick, "solid", "preregistered confirmation"),
        ((x, key_y + 3.6), 1.0, dashed, "pre-specified descriptive"),
        ((x, key_y + 7.2), 0.6, "solid", "development"),
        ((x, key_y + 10.8), 1.1, dotted, "post hoc"),
    ]
    for (kx, ky), lw, ls, text in entries:
        ax.add_line(Line2D([kx, kx + 5.5], [ky, ky], color="black", lw=lw, ls=ls,
                           solid_capstyle="butt", dash_capstyle="round"))
        ax.text(kx + 6.6, ky, text, va="center", ha="left", fontsize=6.5)
    source(pd.DataFrame([{"panel": "c", "stage": " ".join(lines), "border": role}
                         for lines, _f, _lw, _ls, role in stages]), "Fig1_stages.tsv")
    finish(c, "A18: border key words 'preregistered confirmation', 'pre-specified descriptive', 'development', "
              "'post hoc', one per row because the longer words do not fit two per row in panel c.")


# Figure 2 -----------------------------------------------------------------

def fig2() -> None:
    sim = read("results/stage5/sim_ext/full.tsv")
    pool = A13.POOL
    outside = sorted(set(sim["tissue"]) - set(pool))
    if len(outside) != 12:
        raise SystemExit(f"tissues outside the pool: {len(outside)}")
    means = []
    for method, tissues, where in (("BASE-Z", pool, "in"), ("SA-Z", pool, "in"), ("BASE-Z", outside, "out"),
                                   ("SA-Z", outside, "out"), ("SA-pool22", outside, "out")):
        part = A13.tissue_mean(sim, method, tissues)
        part["pool"] = where
        means.append(part)
    means = pd.concat(means, ignore_index=True)
    rows = []
    for tissue in pool:
        hits = {
            "baseline": one(sim, method="BASE-Z", tissue=tissue, rho=0.6),
            "leave_one_host_out": one(sim, method=f"SA-LOHO-{tissue}", tissue=tissue, rho=0.6),
            "hostmix": one(sim, method="SA-Z", tissue=tissue, rho=0.6),
        }
        at_risk = {int(hit["n_at_risk"]) for hit in hits.values()}
        if len(at_risk) != 1:
            raise SystemExit(f"{tissue}: at-risk n differs between models {at_risk}")
        rows.append({"tissue": tissue, "display": A13.display_tissue(tissue), "n_at_risk": at_risk.pop(),
                     **{key: float(hit["host_rate"]) for key, hit in hits.items()}})
    loho = pd.DataFrame(rows).sort_values("baseline", ascending=False, kind="mergesort").reset_index(drop=True)
    source(means, "Fig2_means.tsv")
    source(loho, "Fig2_loho.tsv")

    c = Canvas("Figure2", 142)
    series = {
        "in": [("BASE-Z", BLUE, "o", "solid", label("BASE-Z")), ("SA-Z", ORANGE, "s", "solid", label("SA-Z"))],
        "out": [("BASE-Z", BLUE, "o", "solid", label("BASE-Z")), ("SA-Z", ORANGE, "s", "solid", label("SA-Z")),
                ("SA-pool22", "#6E6E6E", "D", (0, (3, 1.6)), label("SA-pool22"))],
    }
    for left, where, title, letter_x, panel in ((13.0, "in", "Ten host-pool tissues", 0.5, "a"),
                                                (98.0, "out", "Twelve tissues outside the host pool", 85.0, "b")):
        ax = c.axes(left, 11.0, 68.0, 47.0)
        c.letter(letter_x, 8.0, panel)
        for method, color, marker, ls, name in series[where]:
            part = means.loc[means["pool"].eq(where) & means["method"].eq(method)].sort_values("rho")
            ax.plot(part["rho"], part["host_pull_rate"], color=color, marker=marker, ms=3.6, lw=1.0, ls=ls,
                    label=name, clip_on=False, zorder=3)
        ax.set_xlim(1.03, 0.17)
        ax.set_ylim(0, 1)
        ax.xaxis.set_major_locator(FixedLocator([1.0, 0.8, 0.6, 0.4, 0.2]))
        ax.xaxis.set_major_formatter(ticks(1))
        ax.yaxis.set_major_locator(FixedLocator(np.linspace(0, 1, 6)))
        ax.yaxis.set_major_formatter(ticks(1))
        ax.set_xlabel("Tumor RNA fraction \u03c1")
        ax.set_ylabel("Host-attraction rate")
        ax.set_title(title)
        ax.legend(loc="upper left", fontsize=7)

    ax = c.axes(57.0, 76.0, 94.0, 55.0)
    c.letter(0.5, 73.0, "c")
    y = np.arange(len(loho))
    for yi in y:
        ax.axhline(yi, color=GUIDE, lw=0.5, zorder=0)
    ax.scatter(loho["baseline"], y - 0.2, s=17, color=BLUE, marker="o", zorder=3, label=label("BASE-Z"), clip_on=False)
    ax.scatter(loho["leave_one_host_out"], y, s=17, facecolors="white", edgecolors="black", linewidths=0.8,
               marker="o", zorder=3, label="Model trained without that tissue", clip_on=False)
    ax.scatter(loho["hostmix"], y + 0.2, s=15, color=ORANGE, marker="s", zorder=3, label=label("SA-Z"), clip_on=False)
    ax.set_yticks(y)
    ax.set_yticklabels(loho["display"])
    ax.tick_params(axis="y", length=0, pad=7)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(len(loho) - 0.5, -0.5)
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_locator(FixedLocator(np.linspace(0, 1, 6)))
    ax.xaxis.set_major_formatter(ticks(1))
    ax.set_xlabel("Host-attraction rate at \u03c1 = 0.6")
    c.row_text(ax, 166.0, -0.95, "At-risk n", ha="right", fontsize=6.5)
    for yi, n, tissue in zip(y, loho["n_at_risk"], loho["tissue"]):
        c.mark(c.row_text(ax, 166.0, yi, str(int(n)), ha="right", fontsize=6.5), "c", f"{tissue}/n_at_risk",
               "Fig2_loho.tsv", pos(loho, tissue=tissue), "n_at_risk")
    ax.legend(loc="lower right", fontsize=7, handletextpad=0.3, borderaxespad=0.6, frameon=True,
              facecolor="white", edgecolor="none", framealpha=1.0)
    finish(c, "Mean lines only in a and b; new c dot plot at rho = 0.6 sorted by the baseline; "
              "at-risk n per tissue; separate margins for the top row and c.")


# Figure 3 -----------------------------------------------------------------

def fig3() -> None:
    pog = read("results/stage4/confirm/primary.tsv").set_index("hypothesis")
    aux = read("results/stage7/confirm/tables/hypothesis_primary.tsv").set_index("hypothesis")
    overall = read("results/stage4/confirm/secondary_overall.tsv").set_index("method")
    external = read("results/stage8/external/tables/metrics.tsv")
    sets = read("results/stage9/set_metrics.tsv")

    def arm(source_name: str, method: str) -> tuple[float, int]:
        if source_name == "pog_host":
            return float(overall.loc[method, "host_rate"]), int(overall.loc[method, "n_at_risk"])
        if source_name == "pog_top1":
            return float(overall.loc[method, "top1"]), int(overall.loc[method, "n"])
        if source_name == "pog_native":
            return float(overall.loc[method, "native_truth_top1"]), int(overall.loc[method, "n_native_truth"])
        if source_name == "aux_host":
            hit = one(sets, analysis_cohort="aux_rnaseq", set="risk", method=method)
            return float(hit["host_rate"]), int(hit["n"])
        if source_name == "aux_top1":
            hit = one(external, subset="common_label", analysis_cohort="aux_rnaseq", method=method)
            return float(hit["top1"]), int(hit["n"])
        hit = one(sets, analysis_cohort="aux_rnaseq", set="native_truth", method=method)
        return float(hit["top1"]), int(hit["n"])

    spec = [
        ("a", "H1", "POG570", "pog_host"), ("a", "AH1", "Auxiliary RNA-seq", "aux_host"),
        ("b", "H2", "POG570", "pog_top1"), ("b", "AH2", "Auxiliary RNA-seq", "aux_top1"),
        ("c", "H3", "POG570", "pog_native"), ("c", "AH3", "Auxiliary RNA-seq", "aux_native"),
    ]
    rows = []
    for panel, hypothesis, cohort, arm_source in spec:
        if hypothesis.startswith("A"):
            hit = aux.loc[hypothesis]
            diff, n = float(hit["diff_sa_minus_comparator"]), int(hit["n"])
        else:
            hit = pog.loc[hypothesis]
            diff, n = float(hit["diff"]), int(hit["n"])
        base, n_base = arm(arm_source, "BASE-Z")
        mix, n_mix = arm(arm_source, "SA-Z")
        if abs((mix - base) - diff) > 1e-5 or n_base != n or n_mix != n:
            raise SystemExit(f"{hypothesis}: arms {base} {mix} n {n_base} {n_mix} do not match {diff} {n}")
        bound = float(hit["onesided_low"]) if panel == "c" else np.nan
        rows.append({"panel": panel, "hypothesis": hypothesis, "cohort": cohort, "n": n, "diff": diff,
                     "ci_low": float(hit["ci_low"]), "ci_high": float(hit["ci_high"]),
                     "onesided_low": bound, "baseline": base, "hostmix": mix})
    frame = pd.DataFrame(rows)
    source(frame, "Fig3_differences.tsv")

    c = Canvas("Figure3", 64)
    ax = c.axes(48.0, 9.0, 70.0, 44.0)
    titles = {"a": "Host-attraction rate", "b": "Top-1 accuracy", "c": "Native-truth top-1 accuracy"}
    positions = {}
    y = 0.0
    for panel in "abc":
        c.row_text(ax, 0.5, y, panel, fontsize=10, fontweight="bold")
        c.row_text(ax, 5.0, y, titles[panel], fontsize=7, fontweight="bold")
        for row in frame.loc[frame["panel"].eq(panel)].itertuples():
            y += 1.0
            positions[row.hypothesis] = y
            c.row_text(ax, 5.0, y, f"{row.cohort} ({row.hypothesis})")
        y += 1.35
    bottom = y - 0.85
    ax.set_ylim(bottom, -0.6)
    ax.set_xlim(-0.32, 0.36)
    ax.axvline(0, color=ZERO_LINE, lw=0.6, zorder=1)
    margin_top = positions["H3"] - 1.45
    ax.plot([-0.10, -0.10], [margin_top, bottom], color="black", lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.text(-0.112, margin_top + 0.45, "margin", ha="right", va="center", fontsize=6.5)
    for row in frame.itertuples():
        yy = positions[row.hypothesis]
        ax.errorbar(row.diff, yy, xerr=[[row.diff - row.ci_low], [row.ci_high - row.diff]], fmt="o",
                    color="black", ms=4, elinewidth=0.9, capsize=1.8, capthick=0.9, zorder=3)
        if not np.isnan(row.onesided_low):
            ax.plot(row.onesided_low, yy, marker="D", ms=4, mfc="white", mec="black", mew=0.8, ls="none", zorder=4)
        at = pos(frame, hypothesis=row.hypothesis)
        c.mark(c.row_text(ax, 123.0, yy, f"{fmt(row.baseline)} \u2192 {fmt(row.hostmix)}"), row.panel,
               f"{row.hypothesis}/arms", "Fig3_differences.tsv", at, "baseline;hostmix")
        c.mark(c.row_text(ax, 167.0, yy, str(row.n), ha="right"), row.panel, f"{row.hypothesis}/n",
               "Fig3_differences.tsv", at, "n")
    c.text(123.0, 6.6, "Baseline \u2192 HostMix-TOO", fontsize=6.5)
    c.text(167.0, 6.6, "n", fontsize=6.5, ha="right")
    ax.spines["left"].set_visible(False)
    ax.yaxis.set_major_locator(NullLocator())
    ax.xaxis.set_major_locator(FixedLocator([-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]))
    ax.xaxis.set_major_formatter(ticks(2))
    ax.set_xlabel("Difference (HostMix-TOO \u2212 baseline)")
    finish(c, "Read by key from the stored hypothesis tables; one shared difference axis; black points and "
              "intervals; -0.10 tick and margin label; arm values and n in right-hand columns.")


# Figure 5 -----------------------------------------------------------------

def fig5() -> None:
    rna = read("results/stage9/ablation/metrics.tsv")
    micro = read("results/stage10/ablation_microarray/metrics.tsv")
    restriction = read("results/stage10/ablation_microarray/restriction.tsv")
    order = ["BASE-Z", "PURE-Zs", "PURE-Zs-w", "MIX-Z0", "SA-Z"]
    style = {
        "BASE-Z": dict(color=BLUE, edgecolor="none", hatch=None),
        "PURE-Zs": dict(color="#7A7A7A", edgecolor="black", hatch="////"),
        "PURE-Zs-w": dict(color="#B4B4B4", edgecolor="black", hatch="...."),
        "MIX-Z0": dict(color="#DCDCDC", edgecolor="black", hatch="xxxx"),
        "SA-Z": dict(color=ORANGE, edgecolor="none", hatch=None),
    }
    cohorts = COHORTS3 + ["microarray"]
    names = {**COHORT3_NAME, "microarray": "Microarray"}
    rows = []
    for method in order:
        for cohort in COHORTS3:
            hit = one(rna, cohort=cohort, method=method)
            rows.append({"panel": "a", "cohort": cohort, "method": method, "metric": "host_rate",
                         "value": float(hit["host_rate"]), "n": int(hit["n_at_risk"]), "source": "stage9 ablation"})
            rows.append({"panel": "b", "cohort": cohort, "method": method, "metric": "top1",
                         "value": float(hit["top1"]), "n": int(hit["n"]), "source": "stage9 ablation"})
        for panel, set_name, metric in (("a", "at_risk", "host_rate"), ("b", "evaluation", "top1")):
            hit = one(micro, model=method, set=set_name, cohort="all", metric=metric)
            rows.append({"panel": panel, "cohort": "microarray", "method": method, "metric": metric,
                         "value": float(hit["value"]), "n": int(hit["n"]), "source": "stage10 ablation"})
        for cohort in ("MET500", "POG570"):
            part = restriction.loc[restriction["model"].eq(method) & restriction["cohort"].eq(cohort)]
            if len(part) != 4:
                raise SystemExit(f"Figure 5c: {method} {cohort} has {len(part)} platforms")
            rows.append({"panel": "c", "cohort": cohort, "method": method, "metric": "random_removal_mean_drop",
                         "value": float(part["random_mean"].mean()), "n": len(part), "source": "stage10 restriction"})
    frame = pd.DataFrame(rows)
    frame["display"] = frame["cohort"].map(names)
    frame["label"] = frame["method"].map(label)
    source(frame, "Fig5_ablation.tsv")

    c = Canvas("Figure5", 75)
    panels = (
        ("a", cohorts, 12.0, 50.0, 12.0, 49.0, 0.4, 5, "Host-attraction rate", "at risk"),
        ("b", cohorts, 73.0, 50.0, 12.0, 49.0, 1.0, 6, "Top-1 accuracy", "evaluation"),
        ("c", ["MET500", "POG570"], 138.0, 30.0, 12.0, 49.0, 0.1, 6, "Mean top-1 loss,\nrandom gene removal", None),
    )
    width = 0.155
    for panel, groups, left, axis_width, top, axis_height, ymax, nticks, ylabel, n_kind in panels:
        ax = c.axes(left, top, axis_width, axis_height)
        c.letter(left - 11.5, top - 3.0, panel)
        part = frame.loc[frame["panel"].eq(panel)]
        for i, method in enumerate(order):
            values = [float(one(part, cohort=g, method=method)["value"]) for g in groups]
            xs = np.arange(len(groups)) + (i - 2) * width
            ax.bar(xs, values, width=width * 0.9, color=style[method]["color"],
                   edgecolor=style[method]["edgecolor"], hatch=style[method]["hatch"], linewidth=0.5)
            if panel == "c":
                for xp, value, g in zip(xs, values, groups):
                    c.mark(ax.annotate(fmt(value), (xp, value), xytext=(0, 1.5), textcoords="offset points",
                                       ha="center", va="bottom", fontsize=6.5, rotation=90),
                           "c", f"{g}/{method}", "Fig5_ablation.tsv", pos(frame, panel="c", cohort=g, method=method),
                           "value")
        tick_names = []
        for g in groups:
            if n_kind is None:
                tick_names.append(names[g])
                continue
            counts = set(part.loc[part["cohort"].eq(g), "n"])
            if len(counts) != 1:
                raise SystemExit(f"Figure 5{panel}: n differs between models in {g}")
            tick_names.append(f"{names[g].replace('Auxiliary RNA-seq', 'Auxiliary' + chr(10) + 'RNA-seq')}\n(n = {counts.pop()})")
        ax.set_xticks(np.arange(len(groups)))
        ax.set_xticklabels(tick_names, fontsize=6.5)
        ax.tick_params(axis="x", length=0)
        if n_kind is not None:
            c.mark_ticks(ax, "x", panel, [f"{g}/n" for g in groups], "Fig5_ablation.tsv",
                         [pos(frame, panel=panel, cohort=g) for g in groups], "n")
        ax.set_xlim(-0.55, len(groups) - 0.45)
        ax.set_ylim(0, ymax)
        ax.yaxis.set_major_locator(FixedLocator(np.linspace(0, ymax, nticks)))
        ax.yaxis.set_major_formatter(ticks(2 if panel == "c" else 1))
        ax.set_ylabel(ylabel)
    handles = [Patch(facecolor=style[m]["color"], edgecolor=style[m]["edgecolor"], hatch=style[m]["hatch"],
                     linewidth=0.5) for m in order]
    c.legend(85.0, 1.5, handles, [label(m) for m in order], ncol=5, fontsize=6.5, columnspacing=1.0,
             handlelength=1.6)
    finish(c, "A18: a, b and c in one row (axes 50, 50 and 30 mm), height 75 mm; one-row legend on top; "
              "bar edges 0.5 pt; group labels 6.5 pt; values and n unchanged.")


# Figure 6 -----------------------------------------------------------------

def fig6() -> None:
    raw = read("results/stage8/external/tables/metrics.tsv")
    part = raw.loc[raw["subset"].eq("common_label") & raw["method"].isin(METHODS)].copy()
    part["display"] = part["method"].map(label)
    source(part, "Fig6_published.tsv")
    c = Canvas("Figure6", 70)
    width = 0.19
    for left, column, n_column, ylabel, top, panel, letter_x in (
        (14.0, "top1", "n", "Top-1 accuracy", 1.0, "a", 0.5),
        (99.0, "host_rate", "n_at_risk", "Host-attraction rate", 0.4, "b", 86.0),
    ):
        ax = c.axes(left, 12.0, 67.0, 46.0)
        c.letter(letter_x, 9.0, panel)
        for i, method in enumerate(METHODS):
            values = [float(one(part, analysis_cohort=k, method=method)[column]) for k in COHORTS3]
            xs = np.arange(len(COHORTS3)) + (i - 1.5) * width
            ax.bar(xs, values, width=width * 0.9, color=COLOR[method], edgecolor="none")
            if column == "host_rate":
                for xp, value, k in zip(xs, values, COHORTS3):
                    c.mark(ax.annotate(fmt(value), (xp, value), xytext=(0, 1.5), textcoords="offset points",
                                       ha="center", va="bottom", fontsize=6.5, rotation=90),
                           panel, f"{k}/{method}", "Fig6_published.tsv",
                           pos(part, analysis_cohort=k, method=method), "host_rate")
        names = []
        for k in COHORTS3:
            counts = {int(one(part, analysis_cohort=k, method=m)[n_column]) for m in METHODS}
            if len(counts) != 1:
                raise SystemExit(f"{k}: {n_column} differs between methods")
            names.append(f"{COHORT3_NAME[k]}\n(n = {counts.pop()})")
        ax.set_xticks(np.arange(len(COHORTS3)))
        ax.set_xticklabels(names)
        c.mark_ticks(ax, "x", panel, [f"{k}/n" for k in COHORTS3], "Fig6_published.tsv",
                     [pos(part, analysis_cohort=k) for k in COHORTS3], n_column)
        ax.tick_params(axis="x", length=0)
        ax.set_xlim(-0.55, len(COHORTS3) - 0.45)
        ax.set_ylim(0, top)
        ax.yaxis.set_major_locator(FixedLocator(np.linspace(0, top, 5 if top == 0.4 else 6)))
        ax.yaxis.set_major_formatter(ticks(1))
        ax.set_ylabel(ylabel)
    handles, labels = method_handles(METHODS)
    c.legend(85.0, 1.5, handles, labels, ncol=4)
    finish(c, "One legend for both panels; evaluation n under a, at-risk n under b; b on 0-0.4 with values "
              "above the bars.")


# Figure 4 -----------------------------------------------------------------

def fig4() -> None:
    rates = A13.cohort_rates()
    group_name = {"Development": "Development", "Confirmation 1": "Confirmation 1",
                  "Confirmation 2": "Confirmation 2", "Descriptive RNA-seq": "Other RNA-seq",
                  "Descriptive microarray": "Microarray"}
    rates["group_display"] = rates["group"].map(group_name)
    source(rates, "Fig4a_cohorts.tsv")
    site = read("results/stage4/confirm/secondary_by_site.tsv")
    tc = read("results/stage4/confirm/secondary_tc.tsv")
    site = site.loc[site["method"].isin(["BASE-Z", "SA-Z"])]
    tc = tc.loc[tc["method"].isin(["BASE-Z", "SA-Z"])]
    source(site, "Fig4b_site.tsv")
    source(tc, "Fig4b_tc.tsv")
    sub = read("results/stage9/subcohort_metrics.tsv")
    source(sub, "Fig4c_aux_top1.tsv")
    cases = read("results/stage8/diagnostics/error_shift_cases.tsv")
    counts = cases.groupby("analysis").agg(n=("SA_correct", "size"), correct=("SA_correct", "sum")).reset_index()
    counts["other"] = counts["n"] - counts["correct"]
    source(counts, "Fig4d_removed.tsv")

    c = Canvas("Figure4", 218)
    ax = c.axes(36.0, 9.0, 130.0, 88.0)
    c.letter(0.5, 6.0, "a")
    y = 0.0
    rows_y = {}
    for group in group_name.values():
        c.row_text(ax, 1.0, y, group, fontweight="bold")
        y += 1.0
        for cohort in dict.fromkeys(rates.loc[rates["group_display"].eq(group), "cohort"]):
            base = one(rates, cohort=cohort, method="BASE-Z")
            rows_y[cohort] = y
            c.mark(c.row_text(ax, 4.0, y, f"{cohort} ({int(base['n_at_risk'])})"), "a", f"{cohort}/n_at_risk",
                   "Fig4a_cohorts.tsv", pos(rates, cohort=cohort, method="BASE-Z"), "n_at_risk")
            y += 1.0
        y += 0.45
    ax.set_ylim(y - 0.95, -0.6)
    for x in (0.2, 0.4, 0.6, 0.8):
        ax.axvline(x, color=GUIDE, lw=0.5, zorder=0)
    for cohort, yy in rows_y.items():
        base = float(one(rates, cohort=cohort, method="BASE-Z")["host_rate"])
        mix = float(one(rates, cohort=cohort, method="SA-Z")["host_rate"])
        ax.plot([base, mix], [yy, yy], color="#8C8C8C", lw=0.8, zorder=1)
        ax.scatter([base], [yy], s=16, color=BLUE, marker="o", zorder=3, clip_on=False)
        ax.scatter([mix], [yy], s=14, color=ORANGE, marker="s", zorder=3, clip_on=False)
    ax.set_xlim(0, 0.9)
    ax.spines["left"].set_position(("outward", 4))
    ax.xaxis.set_major_locator(FixedLocator([0, 0.2, 0.4, 0.6, 0.8]))
    ax.xaxis.set_major_formatter(ticks(1))
    ax.set_xlabel("Host-attraction rate")
    ax.yaxis.set_major_locator(FixedLocator(list(rows_y.values())))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda *_: ""))
    handles, labels = method_handles(["BASE-Z", "SA-Z"], kind="marker")
    ax.legend(handles, labels, loc="lower right", fontsize=7, borderaxespad=0.6)

    c.letter(0.5, 112.0, "b")
    site_order = [("liver", "Liver"), ("lung", "Lung"), ("lymph_node", "Lymph node"),
                  ("soft_tissue", "Soft tissue"), ("remainder", "Other")]
    tc_order = [("T1", "Low"), ("T2", "Mid"), ("T3", "High")]
    for frame, key, order, left, width_mm, title, kind in (
        (site, "site_group", site_order, 14.0, 72.0, "POG570: biopsy site", "site"),
        (tc, "tc_bin", tc_order, 104.0, 62.0, "POG570: tumor-content tertile", "tc"),
    ):
        ax = c.axes(left, 116.0, width_mm, 34.0)
        xs = np.arange(len(order))
        names = []
        for value, name in order:
            n = {int(one(frame, **{key: value}, method=m)["n_at_risk"]) for m in ("BASE-Z", "SA-Z")}
            names.append(f"{name}\n(n = {n.pop()})")
        for method, offset in (("BASE-Z", -0.17), ("SA-Z", 0.17)):
            values = [float(one(frame, **{key: value}, method=method)["host_rate"]) for value, _ in order]
            ax.bar(xs + offset, values, width=0.32, color=COLOR[method], edgecolor="none")
            for xp, value, (stratum, _name) in zip(xs + offset, values, order):
                if value == 0:
                    c.mark(zero_mark(ax, xp, 0.0), "b", f"{kind}/{stratum}/{method}/zero", f"Fig4b_{kind}.tsv",
                           pos(frame, **{key: stratum}, method=method), "host_rate")
        ax.set_xticks(xs)
        ax.set_xticklabels(names)
        c.mark_ticks(ax, "x", "b", [f"{kind}/{value}/n" for value, _ in order], f"Fig4b_{kind}.tsv",
                     [pos(frame, **{key: value}) for value, _ in order], "n_at_risk")
        ax.tick_params(axis="x", length=0)
        ax.set_xlim(-0.55, len(order) - 0.45)
        ax.set_ylim(0, 0.5)
        ax.yaxis.set_major_locator(FixedLocator(np.linspace(0, 0.5, 6)))
        ax.yaxis.set_major_formatter(ticks(1))
        ax.set_ylabel("Host-attraction rate")
        ax.set_title(title)
        if key == "site_group":
            handles, labels = method_handles(["BASE-Z", "SA-Z"])
            ax.legend(handles, labels, loc="upper right", fontsize=7)

    c.letter(0.5, 171.0, "c")
    ax = c.axes(14.0, 175.0, 96.0, 31.0)
    xs = np.arange(len(AUX))
    width = 0.19
    for i, method in enumerate(METHODS):
        values = [float(one(sub, cohort=k, method=method)["top1"]) for k in AUX]
        ax.bar(xs + (i - 1.5) * width, values, width=width * 0.9, color=COLOR[method], edgecolor="none")
    names = []
    for k in AUX:
        n = {int(one(sub, cohort=k, method=m)["n"]) for m in METHODS}
        names.append(f"{COHORT[k]}\n(n = {n.pop()})")
    ax.set_xticks(xs)
    ax.set_xticklabels(names)
    c.mark_ticks(ax, "x", "c", [f"{k}/n" for k in AUX], "Fig4c_aux_top1.tsv", [pos(sub, cohort=k) for k in AUX], "n")
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.55, len(AUX) - 0.45)
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_locator(FixedLocator(np.linspace(0, 1, 6)))
    ax.yaxis.set_major_formatter(ticks(1))
    ax.set_ylabel("Top-1 accuracy")
    handles, labels = method_handles(METHODS)
    c.legend(62.0, 168.0, handles, labels, ncol=4)

    c.letter(117.0, 171.0, "d")
    ax = c.axes(128.0, 175.0, 38.0, 31.0)
    order = [("MET500", "MET500"), ("POG570", "POG570"), ("aux_rnaseq", "Auxiliary\nRNA-seq")]
    correct = [int(one(counts, analysis=k)["correct"]) for k, _ in order]
    other = [int(one(counts, analysis=k)["other"]) for k, _ in order]
    totals = [int(one(counts, analysis=k)["n"]) for k, _ in order]
    xs = np.arange(len(order))
    ax.bar(xs, correct, width=0.62, color="#4D4D4D", edgecolor="#4D4D4D", linewidth=0.5)
    ax.bar(xs, other, bottom=correct, width=0.62, color="#D9D9D9", edgecolor="#4D4D4D", hatch="////",
           linewidth=0.5)
    for xp, low, high, (analysis, _name) in zip(xs, correct, other, order):
        at = pos(counts, analysis=analysis)
        c.mark(ax.text(xp, low / 2, str(low), ha="center", va="center", fontsize=6.5, color="white"),
               "d", f"{analysis}/correct", "Fig4d_removed.tsv", at, "correct")
        c.mark(ax.text(xp, low + high / 2, str(high), ha="center", va="center", fontsize=6.5, color="black",
                       bbox=dict(boxstyle="square,pad=0.12", facecolor="#D9D9D9", edgecolor="none")),
               "d", f"{analysis}/other", "Fig4d_removed.tsv", at, "other")
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{name}\n(n = {n})" for (_, name), n in zip(order, totals)])
    c.mark_ticks(ax, "x", "d", [f"{k}/n" for k, _ in order], "Fig4d_removed.tsv",
                 [pos(counts, analysis=k) for k, _ in order], "n")
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_ylim(0, 80)
    ax.yaxis.set_major_locator(FixedLocator([0, 20, 40, 60, 80]))
    ax.set_ylabel("Biopsies")
    handles = [Patch(facecolor="#4D4D4D", edgecolor="#4D4D4D", linewidth=0.5),
               Patch(facecolor="#D9D9D9", edgecolor="#4D4D4D", hatch="////", linewidth=0.5)]
    c.legend(147.0, 168.0, handles, ["Correct", "Other error"], ncol=2, columnspacing=1.0)
    finish(c, "A18: bar and legend edges in d 0.5 pt; nothing else changed.")


# Figure S1 ----------------------------------------------------------------

S1_METHODS = ["BASE-Z", "SA-Z", "SC-Z", "NC-Z", "LD-Z", "M1-Z", "IF20-Z", "SC+IF20-Z", "SC+IF40-Z",
              "BASE-K", "SA-K", "V0-K", "SA-G", "SA-pool22", "SA-MLP"]
S1_COHORTS = ["MET500", "POG570", "Auxiliary RNA-seq", "Microarray"]
GREYS = LinearSegmentedColormap.from_list("light_to_dark", ["#F7F7F7", "#262626"])


def s1_status(cohort: str, method: str, status) -> str:
    """Row status by the rule of its result file (POG570: preregistered for the two primary arms)."""
    if cohort == "Auxiliary RNA-seq":
        return "post hoc aggregate"
    if cohort == "POG570" and status == "pre-specified descriptive" and method in ("BASE-Z", "SA-Z"):
        return "preregistered"
    return status


def s1_sizes() -> pd.DataFrame:
    met = read("results/stage2/met500_overall.tsv")
    pog = read("results/stage4/confirm/secondary_overall.tsv")
    aux = read("results/stage10/derived/aux_six_rates.tsv")
    micro = read("results/stage10/derived/microarray_rates.tsv")
    rows = []
    for cohort, frame, key, n_col, risk_col, path in (
        ("MET500", met, {}, "n", "n_at_risk", "results/stage2/met500_overall.tsv"),
        ("POG570", pog, {}, "n", "n_at_risk", "results/stage4/confirm/secondary_overall.tsv"),
        ("Auxiliary RNA-seq", aux, {}, "n_evaluation", "n_at_risk", "results/stage10/derived/aux_six_rates.tsv"),
        ("Microarray", micro, {"cohort": "all"}, "n_evaluation", "n_at_risk", "results/stage10/derived/microarray_rates.tsv"),
    ):
        part = frame
        for column, value in key.items():
            part = part.loc[part[column].eq(value)]
        for panel, column in (("a", n_col), ("b", risk_col)):
            values = set(part[column].astype(int))
            if len(values) != 1:
                raise SystemExit(f"S1 {cohort}: {column} differs between rows {values}")
            rows.append({"panel": panel, "cohort": cohort, "set": "evaluation" if panel == "a" else "at risk",
                         "n": values.pop(), "source_file": path, "source_column": column})
    return pd.DataFrame(rows)


def fig_s1() -> None:
    frame = A13.locked_heatmap()
    frame["status"] = [s1_status(r.cohort, r.method, r.status) for r in frame.itertuples()]
    source(frame, "S1_heatmap.tsv")
    sizes = s1_sizes()
    source(sizes, "S1_n.tsv")
    c = Canvas("FigureS1", 110)
    for left, column, title, panel, letter_x in ((50.0, "top1", "Top-1 accuracy", "a", 0.5),
                                                 (112.0, "host_rate", "Host-attraction rate", "b", 107.0)):
        grid = np.full((len(S1_METHODS), len(S1_COHORTS)), np.nan)
        for i, method in enumerate(S1_METHODS):
            for j, cohort in enumerate(S1_COHORTS):
                hit = frame.loc[frame["method"].eq(method) & frame["cohort"].eq(cohort), column]
                if len(hit) > 1:
                    raise SystemExit(f"S1 duplicate {method} {cohort}")
                if len(hit):
                    grid[i, j] = float(hit.iloc[0])
        vmax = 1.0 if column == "top1" else float(np.nanmax(grid))
        norm = Normalize(0, vmax)
        ax = c.axes(left, 10.0, 54.0, 80.0)
        c.letter(letter_x, 7.0, panel)
        cmap = GREYS.copy()
        cmap.set_bad("white")
        image = ax.imshow(np.ma.masked_invalid(grid), cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
        for i in range(len(S1_METHODS)):
            for j in range(len(S1_COHORTS)):
                value = grid[i, j]
                if np.isnan(value):
                    ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor="white", edgecolor="#B0B0B0",
                                           hatch="////", linewidth=0))
                    ax.text(j, i, "n/c", ha="center", va="center", fontsize=6.5, color="black",
                            bbox=dict(boxstyle="square,pad=0.15", facecolor="white", edgecolor="none"))
                    continue
                r, g, b, _ = cmap(norm(value))
                dark = 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.5
                cohort, method = S1_COHORTS[j], S1_METHODS[i]
                c.mark(ax.text(j, i, fmt(value), ha="center", va="center", fontsize=6.5,
                               color="white" if dark else "black"),
                       panel, f"{cohort}/{method}", "S1_heatmap.tsv", pos(frame, cohort=cohort, method=method), column)
        for i in range(len(S1_METHODS) + 1):
            ax.axhline(i - 0.5, color="white", lw=0.8)
        for j in range(len(S1_COHORTS) + 1):
            ax.axvline(j - 0.5, color="white", lw=0.8)
        for side in ("left", "bottom", "top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_xticks(range(len(S1_COHORTS)))
        n = {cohort: int(one(sizes, panel=panel, cohort=cohort)["n"]) for cohort in S1_COHORTS}
        ax.set_xticklabels([f"{name}\n(n = {n[cohort]})" for cohort, name in
                            zip(S1_COHORTS, ["MET500", "POG570", "Auxiliary\nRNA-seq", "Microarray"])])
        c.mark_ticks(ax, "x", panel, [f"{cohort}/n" for cohort in S1_COHORTS], "S1_n.tsv",
                     [pos(sizes, panel=panel, cohort=cohort) for cohort in S1_COHORTS], "n")
        ax.tick_params(length=0)
        ax.set_yticks(range(len(S1_METHODS)))
        ax.set_yticklabels([label(m) for m in S1_METHODS] if panel == "a" else [])
        ax.set_title(title)
        cax = c.axes(left + 7.0, 100.0, 40.0, 2.2)
        bar = c.fig.colorbar(image, cax=cax, orientation="horizontal")
        bar.outline.set_linewidth(0.5)
        bar.ax.tick_params(labelsize=6.5, length=2, width=0.5, pad=1.5)
        bar.set_ticks([0, 0.5, 1.0] if column == "top1" else [0, 0.1, 0.2])
        bar.ax.xaxis.set_major_formatter(ticks(1))
    finish(c, "A18: n under each column name (a evaluation, b at risk), read by key; nothing else changed.")


# Figure S2 ----------------------------------------------------------------

def fig_s2() -> None:
    sets = read("results/stage9/set_metrics.tsv")
    source(sets, "S2_sets.tsv")
    columns = [("native_truth", "Native truth"), ("risk", "At risk"), ("pool_out_risk", "Pool-out at risk")]
    c = Canvas("FigureS2", 110)
    width = 0.19
    lefts = [15.0, 71.0, 127.0]
    for col, ((set_name, title), left) in enumerate(zip(columns, lefts)):
        c.letter(left - 14.5 if col else 0.5, 10.5, "abc"[col])
        for row, (metric, top_mm, upper) in enumerate((("top1", 13.0, 1.0), ("host_rate", 63.0, 0.4))):
            ax = c.axes(left, top_mm, 40.0, 32.0)
            for i, method in enumerate(METHODS):
                values = [float(one(sets, analysis_cohort=k, set=set_name, method=method)[metric]) for k in COHORTS3]
                xs = np.arange(len(COHORTS3)) + (i - 1.5) * width
                ax.bar(xs, values, width=width * 0.9, color=COLOR[method], edgecolor="none")
                for xp, value, k in zip(xs, values, COHORTS3):
                    if value == 0:
                        c.mark(zero_mark(ax, xp, 0.0), "abc"[col], f"{set_name}/{metric}/{k}/{method}/zero",
                               "S2_sets.tsv", pos(sets, analysis_cohort=k, set=set_name, method=method), metric)
            names = []
            for k in COHORTS3:
                n = {int(one(sets, analysis_cohort=k, set=set_name, method=m)["n"]) for m in METHODS}
                display = "Auxiliary\nRNA-seq" if k == "aux_rnaseq" else k
                names.append(f"{display}\n(n = {n.pop()})")
            ax.set_xticks(np.arange(len(COHORTS3)))
            ax.set_xticklabels(names, fontsize=6.5)
            c.mark_ticks(ax, "x", "abc"[col], [f"{set_name}/{metric}/{k}/n" for k in COHORTS3], "S2_sets.tsv",
                         [pos(sets, analysis_cohort=k, set=set_name, method="BASE-Z") for k in COHORTS3], "n")
            ax.tick_params(axis="x", length=0)
            ax.set_xlim(-0.55, len(COHORTS3) - 0.45)
            ax.set_ylim(0, upper)
            ax.yaxis.set_major_locator(FixedLocator(np.linspace(0, upper, 6 if upper == 1.0 else 5)))
            ax.yaxis.set_major_formatter(ticks(1))
            if row == 0:
                ax.set_title(title)
                if col == 0:
                    ax.set_ylabel("Top-1 accuracy")
                else:
                    ax.yaxis.set_major_formatter(FuncFormatter(lambda *_: ""))
            else:
                if col == 0:
                    ax.set_ylabel("Called another\nnative organ")
                elif col == 1:
                    ax.set_ylabel("Host-attraction rate")
                else:
                    ax.yaxis.set_major_formatter(FuncFormatter(lambda *_: ""))
    handles, labels = method_handles(METHODS)
    c.legend(85.0, 1.5, handles, labels, ncol=4)
    finish(c, "Read from the stored set table; one legend; n under every cell; lower row 0-0.4; "
              "native-truth lower cell relabeled to its own metric.")


# Figure S3 ----------------------------------------------------------------

def fig_s3() -> None:
    absorb = read("results/stage8/diagnostics/esophagus_absorption.tsv")
    layer = read("results/stage10/derived/s3a_microarray_layer.tsv")
    if set(layer["n_evaluation"]) != {536} or set(layer["n_truth_esophagus"]) != {0}:
        raise SystemExit("S3a microarray layer set is not 536 biopsies without esophagus truth")
    layer = layer.rename(columns={"esophagus_fraction": "esophagus_absorption"})
    absorb = pd.concat([absorb.loc[absorb["analysis"].ne("microarray")], layer], ignore_index=True)
    colon = read("results/stage8/diagnostics/gse41258_colon_confusion.tsv")
    mask = read("results/stage9/mask_metrics.tsv")
    random = read("results/stage9/mask_random_summary.tsv")
    platform = mask.loc[mask["rep"].eq("platform") & mask["method"].isin(["BASE-Z", "SA-Z"])]
    control = random.loc[random["method"].isin(["BASE-Z", "SA-Z"]) & random["mask"].ne("unmasked")]
    groups = [("MET500", "MET500", "MET500"), ("POG570", "POG570", "POG570"),
              ("aux_rnaseq", "all", "Auxiliary\nRNA-seq"), ("microarray", "all", "Microarray")]
    methods = [("BASE-Z", BLUE, None, label("BASE-Z")), ("SA-Z", ORANGE, None, label("SA-Z")),
               ("LD-Z", "#9A9A9A", "////", label("LD-Z"))]
    drawn = []
    for analysis, cohort, _ in groups:
        for method, *_rest in methods:
            hit = absorb.loc[absorb["analysis"].eq(analysis) & absorb["cohort"].eq(cohort) & absorb["method"].eq(method)]
            drawn.append({"analysis": analysis, "cohort": cohort, "method": method,
                          "n_truth_not_esophagus": int(hit["n_truth_not_esophagus"].iloc[0]) if len(hit) else np.nan,
                          "esophagus_fraction": float(hit["esophagus_absorption"].iloc[0]) if len(hit) else np.nan,
                          "shown": "value" if len(hit) else "n/c",
                          "source": "results/stage10/derived/s3a_microarray_layer.tsv" if analysis == "microarray"
                          else "results/stage8/diagnostics/esophagus_absorption.tsv"})
    drawn = pd.DataFrame(drawn)
    source(drawn, "S3_esophagus_fraction.tsv")
    slice_rows = colon.loc[colon["method"].eq("SA-Z")].sort_values("n", ascending=False, kind="mergesort")
    if set(slice_rows["n_slice"]) != {183} or int(slice_rows["n"].sum()) != 183:
        raise SystemExit("GSE41258 native-truth slice is not 183 tumors")
    source(slice_rows, "S3_gse41258.tsv")
    source(platform, "S3_platform.tsv")
    source(control, "S3_random.tsv")

    c = Canvas("FigureS3", 120)
    c.letter(0.5, 9.0, "a")
    ax = c.axes(15.0, 12.0, 68.0, 40.0)
    width = 0.26
    for i, (method, color, hatch, _name) in enumerate(methods):
        xs = np.arange(len(groups)) + (i - 1) * width
        for xp, (analysis, cohort, _) in zip(xs, groups):
            hit = one(drawn, analysis=analysis, cohort=cohort, method=method)
            if hit["shown"] == "n/c":
                zero_mark(ax, xp, 0.0, "n/c")
                continue
            ax.bar(xp, hit["esophagus_fraction"], width=width * 0.9, color=color, hatch=hatch,
                   edgecolor="black" if hatch else "none", linewidth=0.5)
    names = []
    for analysis, cohort, display in groups:
        n = set(drawn.loc[drawn["analysis"].eq(analysis) & drawn["cohort"].eq(cohort), "n_truth_not_esophagus"].dropna())
        if len(n) != 1:
            raise SystemExit(f"S3 {analysis}: denominators {n}")
        names.append(f"{display}\n(n = {int(n.pop())})")
    ax.set_xticks(np.arange(len(groups)))
    ax.set_xticklabels(names)
    c.mark_ticks(ax, "x", "a", [f"{analysis}/n" for analysis, _c, _d in groups], "S3_esophagus_fraction.tsv",
                 [pos(drawn, analysis=a, method="BASE-Z") for a, _c, _d in groups],
                 "n_truth_not_esophagus")
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.55, len(groups) - 0.45)
    ax.set_ylim(0, 0.4)
    ax.yaxis.set_major_locator(FixedLocator(np.linspace(0, 0.4, 5)))
    ax.yaxis.set_major_formatter(ticks(1))
    ax.set_ylabel("Fraction called esophagus")
    handles = [Patch(facecolor=color, hatch=hatch, edgecolor="black" if hatch else "none", linewidth=0.5)
               for _m, color, hatch, _n in methods]
    c.legend(49.0, 1.5, handles, [name for *_x, name in methods], ncol=3)

    c.letter(85.0, 9.0, "b")
    ax = c.axes(104.0, 12.0, 62.0, 40.0)
    y = np.arange(len(slice_rows))
    ax.barh(y, slice_rows["n"], height=0.62, color=ORANGE, edgecolor="none")
    for yi, n, predicted in zip(y, slice_rows["n"], slice_rows["predicted"]):
        c.mark(ax.annotate(str(int(n)), (float(n), yi), xytext=(1.5, 0), textcoords="offset points",
                           ha="left", va="center", fontsize=6.5),
               "b", f"{predicted}/n", "S3_gse41258.tsv", pos(slice_rows, predicted=predicted), "n")
    ax.set_yticks(y)
    ax.set_yticklabels(slice_rows["predicted"])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(slice_rows) - 0.5, -0.5)
    ax.set_xlim(0, 120)
    ax.xaxis.set_major_locator(FixedLocator([0, 25, 50, 75, 100]))
    ax.set_xlabel("Tumors")
    c.mark(ax.set_title(f"HostMix-TOO, GSE41258 native truth (n = {int(slice_rows['n_slice'].iloc[0])})"), "b",
           "title/n", "S3_gse41258.tsv", 0, "n_slice")

    c.letter(0.5, 69.0, "c")
    ax = c.axes(15.0, 73.0, 100.0, 37.0)
    platforms = [("GPL96", "GPL96"), ("GPL20769", "GPL20769"), ("GPL15659", "GPL15659"),
                 ("prad_fhcrc_agilent", "FHCRC Agilent")]
    series = [("POG570", "BASE-Z", BLUE, False), ("POG570", "SA-Z", ORANGE, False),
              ("MET500", "BASE-Z", BLUE, True), ("MET500", "SA-Z", ORANGE, True)]
    width = 0.19
    n_cohort = {}
    for i, (cohort, method, color, hatched) in enumerate(series):
        xs = np.arange(len(platforms)) + (i - 1.5) * width
        for xp, (key, _name) in zip(xs, platforms):
            hit = one(platform, cohort=cohort, mask=key, method=method)
            n_cohort.setdefault(cohort, set()).add(int(hit["n"]))
            if hatched:
                ax.bar(xp, hit["top1_change"], width=width * 0.88, facecolor="white", edgecolor=color,
                       hatch="//////", linewidth=0.6)
            else:
                ax.bar(xp, hit["top1_change"], width=width * 0.88, color=color, edgecolor="none")
            rnd = one(control, cohort=cohort, mask=key, method=method)
            origin = float(rnd["top1_mean"]) - float(rnd["top1_change_mean"])
            ax.vlines(xp, float(rnd["top1_min"]) - origin, float(rnd["top1_max"]) - origin, color="black",
                      lw=0.6, zorder=4)
            ax.hlines(float(rnd["top1_change_mean"]), xp - width * 0.3, xp + width * 0.3, color="black",
                      lw=1.1, zorder=5)
    ax.axhline(0, color=ZERO_LINE, lw=0.6, zorder=1)
    ax.set_xticks(np.arange(len(platforms)))
    ax.set_xticklabels([name for _k, name in platforms])
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.55, len(platforms) - 0.45)
    ax.set_ylim(-0.15, 0.03)
    ax.yaxis.set_major_locator(FixedLocator([-0.15, -0.10, -0.05, 0.0]))
    ax.yaxis.set_major_formatter(ticks(2))
    ax.set_ylabel("Change in top-1 accuracy")
    sizes = {cohort: values.pop() for cohort, values in n_cohort.items() if len(values) == 1}
    handles = [
        Patch(facecolor=BLUE, edgecolor="none"), Patch(facecolor=ORANGE, edgecolor="none"),
        Patch(facecolor="white", edgecolor=BLUE, hatch="//////", linewidth=0.6),
        Patch(facecolor="white", edgecolor=ORANGE, hatch="//////", linewidth=0.6),
        Line2D([], [], color="black", marker="_", ms=7, mew=1.1, ls="none"),
        Line2D([], [], color="black", marker="|", ms=7, mew=0.6, ls="none"),
    ]
    labels = [f"POG570 (n = {sizes['POG570']}), baseline", f"POG570 (n = {sizes['POG570']}), HostMix-TOO",
              f"MET500 (n = {sizes['MET500']}), baseline", f"MET500 (n = {sizes['MET500']}), HostMix-TOO",
              "Mean of ten random removals", "Minimum to maximum of ten"]
    key = c.fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(119.0 / c.w, 1 - 74.0 / c.h), ncol=1,
                       fontsize=6.5, labelspacing=0.55)
    for i, (cohort, name) in enumerate([("POG570", "baseline"), ("POG570", "HostMix-TOO"),
                                        ("MET500", "baseline"), ("MET500", "HostMix-TOO")]):
        method = "BASE-Z" if name == "baseline" else "SA-Z"
        c.mark(key.get_texts()[i], "c", f"{cohort}/{name}/n", "S3_platform.tsv",
               pos(platform, cohort=cohort, method=method), "n")
    finish(c, "A18: microarray group of a recounted on the layer set (n = 536, post hoc); bar and legend "
              "edges 0.5 pt; nothing else changed.")


# Figure S4 ----------------------------------------------------------------

def fig_s4() -> None:
    raw = read("results/stage9/imvigor_contribution_fix1.tsv")
    pieces = []
    for method, part in raw.groupby("method"):
        ordered = part.sort_values("mean_contribution", ascending=False, kind="mergesort")
        if (ordered.head(10)["mean_contribution"] <= 0).any() or (ordered.tail(5)["mean_contribution"] >= 0).any():
            raise SystemExit(f"{method}: fewer than ten positive or five negative contributions")
        pieces.append(ordered.head(10).assign(side="positive", order=range(1, 11)))
        pieces.append(ordered.tail(5).assign(side="negative", order=range(1, 6)))
    shown = pd.concat(pieces, ignore_index=True)
    shown["display"] = shown["method"].map(label)
    source(shown, "S4_contributions.tsv")
    c = Canvas("FigureS4", 118)
    for left, method, color, panel, letter_x in ((21.0, "BASE-Z", BLUE, "a", 0.5),
                                                 (106.0, "SA-Z", ORANGE, "b", 87.0)):
        part = shown.loc[shown["method"].eq(method)]
        positive = part.loc[part["side"].eq("positive")].sort_values("mean_contribution", ascending=False)
        negative = part.loc[part["side"].eq("negative")].sort_values("mean_contribution", ascending=False)
        genes = list(positive["gene"]) + list(negative["gene"])
        values = list(positive["mean_contribution"]) + list(negative["mean_contribution"])
        y = list(range(10)) + [11 + i for i in range(5)]
        ax = c.axes(left, 11.0, 60.0, 90.0)
        c.letter(letter_x, 8.0, panel)
        ax.barh(y, values, height=0.68, color=color, edgecolor="none")
        ax.axvline(0, color="black", lw=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels(genes, fontstyle="italic")
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.set_ylim(15.5, -0.6)
        low, high = min(values), max(values)
        pad = 0.06 * (high - low)
        ax.set_xlim(low - pad, high + pad)
        ax.set_xlabel("Mean contribution to the\nesophagus \u2212 bladder logit")
        c.mark(ax.set_title(f"{label(method)} (n = {int(part['n_samples'].iloc[0])})"), panel, f"{method}/n",
               "S4_contributions.tsv", pos(shown, method=method), "n_samples")
    finish(c, "Ten largest positive and five most negative contributions in signed order with a gap; "
              "zero line; no tick at the gap; italic gene symbols.")


# Figure S5 ----------------------------------------------------------------

def fig_s5() -> None:
    raw = read("results/stage9/esophagus_histology_summary.tsv")
    groups = [("Auxiliary RNA-seq", AUX), ("Other RNA-seq", ["GSE209998"]), ("Microarray", MICROARRAY)]
    methods = ["BASE-Z", "SA-Z", "LD-Z"]
    rows = []
    for method in methods:
        for _group, cohorts in groups:
            for cohort in cohorts:
                hit = raw.loc[raw["method"].eq(method) & raw["cohort"].eq(cohort)]
                if len(hit) > 1:
                    raise SystemExit(f"S5 duplicate {method} {cohort}")
                row = hit.iloc[0] if len(hit) else None
                rows.append({"cohort": COHORT.get(cohort, cohort), "method": label(method), "internal": method,
                             "n": int(row["n"]) if row is not None else 0,
                             "adenocarcinoma": int(row["n_closer_adenocarcinoma"]) if row is not None else 0,
                             "squamous": int(row["n_closer_squamous"]) if row is not None else 0,
                             "tie": int(row["n_tie"]) if row is not None else 0})
    table = pd.DataFrame(rows)
    if int(table["n"].sum()) != 904 or (table["adenocarcinoma"] + table["squamous"] + table["tie"] != table["n"]).any():
        raise SystemExit("S5 counts do not add up to the 904 scored calls")
    source(raw, "S5_histology.tsv")
    source(table.drop(columns="internal"), "S5_reconciliation.tsv")

    c = Canvas("FigureS5", 98)
    positions = []
    y = 0.0
    headers = []
    for group, cohorts in groups:
        headers.append((group, y))
        y += 1.0
        for cohort in cohorts:
            positions.append((COHORT.get(cohort, cohort), y))
            y += 1.0
        y += 0.3
    bottom = y - 0.8
    lefts = [26.0, 74.0, 122.0]
    for k, (method, left) in enumerate(zip(methods, lefts)):
        ax = c.axes(left, 13.0, 44.0, 72.0)
        c.letter(left - 5.0 if k else 0.5, 10.5, "abc"[k])
        part = table.loc[table["internal"].eq(method)].set_index("cohort")
        for name, yy in positions:
            adeno, squamous, n = (int(part.loc[name, col]) for col in ("adenocarcinoma", "squamous", "n"))
            ax.barh(yy, adeno, height=0.66, color="#595959", edgecolor="none")
            ax.barh(yy, squamous, left=adeno, height=0.66, color="#D4D4D4", edgecolor="#595959", hatch="////",
                    linewidth=0.5)
            c.mark(ax.annotate(f"n = {n}", (adeno + squamous, yy), xytext=(1.5, 0), textcoords="offset points",
                               ha="left", va="center", fontsize=6.5),
                   "abc"[k], f"{name}/{method}/n", "S5_reconciliation.tsv",
                   pos(table, cohort=name, internal=method), "n")
        if k == 0:
            for name, yy in positions:
                c.row_text(ax, 4.0, yy, name)
            for group, yy in headers:
                c.row_text(ax, 1.0, yy, group, fontweight="bold")
        ax.set_ylim(bottom, -0.6)
        ax.yaxis.set_major_locator(FixedLocator([yy for _n, yy in positions]))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda *_: ""))
        ax.set_xlim(0, 300)
        ax.xaxis.set_major_locator(FixedLocator([0, 100, 200]))
        ax.set_xlabel("Esophagus calls")
        ax.set_title(label(method))
    handles = [Patch(facecolor="#595959", edgecolor="none"),
               Patch(facecolor="#D4D4D4", edgecolor="#595959", hatch="////", linewidth=0.5)]
    c.legend(85.0, 1.5, handles, ["Closer to the adenocarcinoma centroid", "Closer to the squamous centroid"], ncol=2)
    finish(c, "A18: bar and legend edges 0.5 pt; nothing else changed.")


# Figure S6 ----------------------------------------------------------------

def fig_s6() -> None:
    raw = read("results/stage7/confirm/tables/meta_analysis.tsv")
    cohorts = raw.loc[raw["status"].eq("cohort")].copy()
    se = np.sqrt(cohorts["v"].astype(float))
    cohorts["wald_low"] = cohorts["d"] - 1.96 * se
    cohorts["wald_high"] = cohorts["d"] + 1.96 * se
    layers = [("RNA-seq", "RNA-seq"), ("\ub9c8\uc774\ud06c\ub85c\uc5b4\ub808\uc774", "Microarray")]
    drawn = []
    c = Canvas("FigureS6", 112)
    ax = c.axes(42.0, 12.0, 74.0, 86.0)
    y = 0.0
    for layer, title in layers:
        c.row_text(ax, 1.0, y, title, fontweight="bold")
        y += 1.0
        for row in cohorts.loc[cohorts["layer"].eq(layer)].itertuples():
            name = COHORT.get(row.cohort, row.cohort)
            ax.errorbar(row.d, y, xerr=[[row.d - row.wald_low], [row.wald_high - row.d]], fmt="o", color="black",
                        ms=3.6, elinewidth=0.9, capsize=1.6, capthick=0.9, zorder=3)
            at = len(drawn)
            c.mark(c.row_text(ax, 4.0, y, f"{name} ({int(row.n)})"), "", f"{name}/n_at_risk", "S6_meta.tsv", at,
                   "n_at_risk")
            c.mark(c.row_text(ax, 120.0, y, f"{fmt(row.d)} [{fmt(row.wald_low)} to {fmt(row.wald_high)}]"), "",
                   f"{name}/interval", "S6_meta.tsv", at, "difference;ci_low;ci_high")
            drawn.append({"layer": title, "row": "cohort", "cohort": name, "n_at_risk": int(row.n), "difference": row.d,
                          "ci_low": row.wald_low, "ci_high": row.wald_high, "interval": "Wald, recorded variance",
                          "I2": np.nan, "I2_percent": np.nan, "k": np.nan})
            y += 1.0
        overall = raw.loc[raw["status"].eq("\uacc4\uc0b0") & raw["layer"].eq(layer)]
        if len(overall) != 1:
            raise SystemExit(f"meta summary rows for {layer}: {len(overall)}")
        hit = overall.iloc[0]
        estimate, low, high = float(hit["estimate"]), float(hit["ci_low"]), float(hit["ci_high"])
        if int(hit["k"]) != int(cohorts["layer"].eq(layer).sum()):
            raise SystemExit(f"{layer}: k does not match the cohort rows")
        ax.add_patch(Polygon([(low, y), (estimate, y - 0.34), (high, y), (estimate, y + 0.34)], closed=True,
                             facecolor="black", edgecolor="black", linewidth=0.5, zorder=3))
        c.row_text(ax, 4.0, y, "DerSimonian\u2013Laird")
        at = len(drawn)
        percent = float(hit["I2"]) * 100
        c.mark(c.row_text(ax, 120.0, y, f"{fmt(estimate)} [{fmt(low)} to {fmt(high)}]", fontweight="bold"), "",
               f"{title}/pooled/interval", "S6_meta.tsv", at, "difference;ci_low;ci_high")
        c.mark(c.row_text(ax, 167.0, y, fmt(percent, 1), ha="right"), "", f"{title}/pooled/I2", "S6_meta.tsv", at,
               "I2_percent")
        drawn.append({"layer": title, "row": "DerSimonian-Laird", "cohort": "", "n_at_risk": np.nan,
                      "difference": estimate, "ci_low": low, "ci_high": high,
                      "interval": "DerSimonian-Laird 95% interval", "I2": float(hit["I2"]), "I2_percent": percent,
                      "k": int(hit["k"])})
        y += 1.5
    bottom = y - 1.0
    ax.set_ylim(bottom, -0.6)
    ax.axvline(0, color=ZERO_LINE, lw=0.6, zorder=1)
    ax.spines["left"].set_visible(False)
    ax.yaxis.set_major_locator(NullLocator())
    ax.set_xlim(-0.85, 0.25)
    ax.xaxis.set_major_locator(FixedLocator([-0.8, -0.6, -0.4, -0.2, 0.0, 0.2]))
    ax.xaxis.set_major_formatter(ticks(1))
    ax.set_xlabel("Difference in host-attraction rate (HostMix-TOO \u2212 baseline)")
    c.text(120.0, 9.6, "Difference [95% CI]", fontsize=6.5)
    c.text(167.0, 9.6, "I\u00b2 (%)", fontsize=6.5, ha="right")
    source(pd.DataFrame(drawn), "S6_meta.tsv")
    finish(c, "A18: I2 printed as a percentage with one decimal under the column title 'I2 (%)'; "
              "nothing else changed.")


FIGURES = {
    "Figure1": fig1, "Figure2": fig2, "Figure3": fig3, "Figure4": fig4, "Figure5": fig5, "Figure6": fig6,
    "FigureS1": fig_s1, "FigureS2": fig_s2, "FigureS3": fig_s3, "FigureS4": fig_s4, "FigureS5": fig_s5,
    "FigureS6": fig_s6,
}


def main() -> None:
    chosen = sys.argv[1:] or list(FIGURES)
    unknown = [name for name in chosen if name not in FIGURES]
    if unknown:
        raise SystemExit(f"unknown figures {unknown}")
    for name in chosen:
        FIGURES[name]()
    table = pd.DataFrame(AUDIT)
    print(table.drop(columns="changes").to_string(index=False))
    bad = table.loc[(table[["overlaps", "outside_page", "outside_axes", "text_on_axes", "tight_boxes"]].sum(axis=1) > 0)
                    | (table["min_font_pt_pdf"] < MIN_PT) | ((table["width_mm"] - WIDTH).abs() > 0.05)]
    if len(chosen) == len(FIGURES):
        CHECKS.mkdir(parents=True, exist_ok=True)
        table.to_csv(CHECKS / "figure_audit.tsv", sep="\t", index=False)
        labels = pd.DataFrame(LABELS)
        labels.to_csv(CHECKS / "figure_labels.tsv", sep="\t", index=False)
        print("value labels", len(labels), labels.groupby("figure", sort=False).size().to_dict())
    elif not len(bad) and (CHECKS / "figure_audit.tsv").exists():
        old = pd.read_csv(CHECKS / "figure_audit.tsv", sep="\t")
        old = old.loc[~old["figure"].isin(table["figure"])]
        merged = pd.concat([old, table], ignore_index=True)
        merged["_order"] = merged["figure"].map({name: i for i, name in enumerate(FIGURES)})
        merged.sort_values("_order").drop(columns="_order").to_csv(CHECKS / "figure_audit.tsv", sep="\t", index=False)
    if len(bad):
        raise SystemExit(f"audit failed: {bad['figure'].tolist()}")


if __name__ == "__main__":
    main()
