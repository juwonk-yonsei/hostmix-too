#!/usr/bin/env python3
"""A18 section 5: submission DOCX files, figure contact sheet and their checks (A19: outputs _A19).

  manuscript/bmc/A_main.docx                 main text without figures; legends and tables at the end, as in
                                             A_main.md; line numbers, page numbers, double spacing
  manuscript/bmc/Additional_file_1.docx      Additional file 1 with Figures S1-S6, each placed above its legend
  manuscript/bmc/cover_letter.docx           cover_letter.md; single spacing, no line numbers
  manuscript/bmc/figures_contact_sheet.pdf   one figure and its legend per page (PNG of the figure, legend text
                                             from the manuscript)
  manuscript/checks/docx_check_A19.tsv       figure objects, '**', HTML entities, '{{', line numbering, page field,
                                             line spacing, size
  manuscript/checks/contact_sheet_A19.tsv    page, figure, legend words, fonts embedded

pandoc is not installed; the DOCX files are written with python-docx (as 43_bmc_docx.py did).
"""
from __future__ import annotations

import html
import re
import subprocess
import textwrap
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from docx import Document  # noqa: E402
from docx.enum.section import WD_ORIENT  # noqa: E402,F401
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING  # noqa: E402
from docx.oxml import OxmlElement  # noqa: E402
from docx.oxml.ns import qn  # noqa: E402
from docx.shared import Cm, Pt, RGBColor  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BMC = ROOT / "manuscript" / "bmc"
FIGS = BMC / "figures"
CHECKS = ROOT / "manuscript" / "checks"
FONT = "Liberation Serif"
AUTHORS = "Juwon Kang; Junjeong Choi"
INLINE = re.compile(r"(\*\*.+?\*\*|(?<![*\w])\*[^*\s][^*]*?\*(?![*\w])|`[^`]+`|\^(?:\u2212?\d+|[A-Za-z]))")
FIG_TITLE = re.compile(r"^\*\*(Figure (S?\d+))\. (.+?)\*\*$")


# ---------------------------------------------------------------- DOCX writing
def set_section(section) -> None:
    section.page_width, section.page_height = Cm(21.0), Cm(29.7)
    for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(section, side, Cm(2.5))
    sect = section._sectPr
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:restart"), "continuous")
    ln.set(qn("w:distance"), "283")
    anchor = sect.find(qn("w:pgNumType"))
    if anchor is None:
        anchor = sect.find(qn("w:cols"))
    if anchor is not None:
        anchor.addprevious(ln)
    else:
        sect.append(ln)
    footer = section.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for kind, text in (("begin", None), ("instr", " PAGE "), ("separate", None), ("text", "1"), ("end", None)):
        run = footer.add_run()
        if kind == "instr":
            el = OxmlElement("w:instrText")
            el.set(qn("xml:space"), "preserve")
            el.text = text
            run._r.append(el)
        elif kind == "text":
            run.text = text
        else:
            el = OxmlElement("w:fldChar")
            el.set(qn("w:fldCharType"), kind)
            run._r.append(el)


def set_styles(document: Document) -> None:
    for name in ("Normal", "Title", "Heading 1", "Heading 2", "Heading 3", "List Bullet"):
        style = document.styles[name]
        style.font.name = FONT
        style.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
        style.font.color.rgb = RGBColor(0, 0, 0)
        fmt = style.paragraph_format
        fmt.line_spacing_rule = WD_LINE_SPACING.DOUBLE
        fmt.space_before = Pt(0)
        fmt.space_after = Pt(0)
    document.styles["Normal"].font.size = Pt(12)
    document.styles["Title"].font.size = Pt(16)
    document.styles["Title"].font.bold = True
    for name, size in (("Heading 1", 14), ("Heading 2", 13), ("Heading 3", 12)):
        document.styles[name].font.size = Pt(size)
        document.styles[name].font.bold = True
        document.styles[name].font.italic = False


def add_runs(paragraph, text: str, size: float | None = None) -> None:
    for part in INLINE.split(text):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**") and len(part) > 4:
            add_runs_styled(paragraph, part[2:-2], size, bold=True)
        elif part.startswith("*") and part.endswith("*") and len(part) > 2:
            add_runs_styled(paragraph, part[1:-1], size, italic=True)
        elif part.startswith("`") and part.endswith("`"):
            run = paragraph.add_run(part[1:-1])
            run.font.name = "Liberation Mono"
            if size:
                run.font.size = Pt(size)
        elif part.startswith("^") and len(part) > 1:
            run = paragraph.add_run(part[1:])
            run.font.superscript = True
            if size:
                run.font.size = Pt(size)
        else:
            run = paragraph.add_run(part)
            if size:
                run.font.size = Pt(size)


def add_runs_styled(paragraph, text: str, size, bold=False, italic=False) -> None:
    for part in INLINE.split(text):
        if not part:
            continue
        sup = part.startswith("^") and len(part) > 1
        ital = part.startswith("*") and part.endswith("*") and len(part) > 2 and not part.startswith("**")
        run = paragraph.add_run(part[1:] if sup else part[1:-1] if ital else part)
        run.bold = bold or None
        run.italic = (italic or ital) or None
        run.font.superscript = sup or None
        if size:
            run.font.size = Pt(size)


def add_lines(paragraph, lines: list[str]) -> None:
    for i, line in enumerate(lines):
        if i:
            paragraph.add_run().add_break()
        add_runs(paragraph, line)


def add_table(document: Document, lines: list[str]) -> None:
    rows = [[c.strip() for c in line.strip().strip("|").split("|")] for line in lines]
    rows = [r for r in rows if not all(re.fullmatch(r":?-{3,}:?", c) for c in r)]
    width = max(len(r) for r in rows)
    table = document.add_table(rows=len(rows), cols=width)
    table.style = "Table Grid"
    for i, row in enumerate(rows):
        for j in range(width):
            cell = table.rows[i].cells[j]
            paragraph = cell.paragraphs[0]
            paragraph.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            text = row[j] if j < len(row) else ""
            add_runs(paragraph, f"**{text}**" if i == 0 and text and not text.startswith("**") else text, size=9)


def add_markdown(document: Document, text: str, figures: dict[str, Path] | None = None) -> list[str]:
    placed = []
    for block in re.split(r"\n\s*\n", text.strip()):
        lines = [line for line in block.splitlines() if line.strip()]
        if not lines:
            continue
        if all(line.lstrip().startswith("|") for line in lines):
            add_table(document, lines)
            continue
        if all(line.startswith("- ") for line in lines):
            for line in lines:
                add_runs(document.add_paragraph(style="List Bullet"), line[2:])
            continue
        head = re.match(r"^(#{1,3}) (.*)$", lines[0])
        if head and len(lines) == 1:
            level = len(head.group(1))
            title = head.group(2).replace("**", "").strip()
            if level == 1 and not document.paragraphs:
                document.add_paragraph(title, style="Title")
            else:
                document.add_heading(title, level=level)
            continue
        match = FIG_TITLE.match(lines[0])
        if figures and match and match.group(1) in figures:
            picture = document.add_paragraph()
            picture.paragraph_format.keep_with_next = True
            picture.add_run().add_picture(str(figures[match.group(1)]), width=Cm(16.0))
            placed.append(match.group(1))
        add_lines(document.add_paragraph(), lines)
    return placed


def set_properties(document: Document, text: str, title: str | None = None) -> None:
    props = document.core_properties
    props.author = AUTHORS
    props.last_modified_by = AUTHORS
    props.title = title or text.splitlines()[0].lstrip("# ").strip()
    props.comments = ""
    props.created = props.modified = datetime.now(timezone.utc).replace(microsecond=0)


def write_docx(src: Path, dest: Path, figures: dict[str, Path] | None = None) -> list[str]:
    document = Document()
    set_styles(document)
    set_section(document.sections[0])
    text = src.read_text()
    set_properties(document, text)
    placed = add_markdown(document, text, figures)
    document.save(dest)
    return placed


def write_letter(src: Path, dest: Path) -> None:
    document = Document()
    normal = document.styles["Normal"]
    normal.font.name = FONT
    normal.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    normal.font.size = Pt(12)
    normal.paragraph_format.line_spacing = 1.15
    normal.paragraph_format.space_before = Pt(0)
    normal.paragraph_format.space_after = Pt(10)
    section = document.sections[0]
    section.page_width, section.page_height = Cm(21.0), Cm(29.7)
    for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(section, side, Cm(2.5))
    text = src.read_text()
    set_properties(document, text, title="Cover letter")
    for block in re.split(r"\n\s*\n", text.strip()):
        paragraph = document.add_paragraph()
        paragraph.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        add_lines(paragraph, [line for line in block.splitlines() if line.strip()])
    document.save(dest)


# ---------------------------------------------------------------- checks
def docx_check(path: Path, placed: list[str]) -> dict:
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode()
        footer = "".join(z.read(n).decode() for n in z.namelist() if n.startswith("word/footer"))
        media = [n for n in z.namelist() if n.startswith("word/media/")]
    texts = "".join(html.unescape(t) for t in re.findall(r"<w:t(?: [^>]*)?>([^<]*)</w:t>", xml))
    document = Document(path)
    body = [p for p in document.paragraphs if p.text.strip() and p.style.name in {"Normal", "List Bullet"}]
    spacing = {str(p.paragraph_format.line_spacing_rule or p.style.paragraph_format.line_spacing_rule) for p in body}
    return {"file": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
            "size_ok_10MB": path.stat().st_size <= 10 * 1024 * 1024,
            "figure_objects": xml.count("<w:drawing>"), "media_files": len(media), "figures_placed": " ".join(placed),
            "double_asterisk": texts.count("**"), "html_entities": len(re.findall(r"&(?:[A-Za-z]+|#\d+);", texts)),
            "double_brace": texts.count("{{"), "line_numbering": "<w:lnNumType" in xml,
            "page_number_field": "PAGE" in footer, "body_line_spacing": "; ".join(sorted(spacing)),
            "body_paragraphs": len(body), "tables": len(document.tables)}


# ---------------------------------------------------------------- contact sheet
def legends() -> list[tuple[str, str, str, Path]]:
    out = []
    for doc in ("A_main.md", "Additional_file_1.md"):
        text = (BMC / doc).read_text()
        for m in re.finditer(r"^\*\*(Figure (S?\d+))\. (.+?)\*\*\n\n(.+?)(?=\n\n|\Z)", text, re.M | re.S):
            out.append((m.group(1), m.group(3), m.group(4).strip(), FIGS / f"Figure{m.group(2)}.png"))
    return out


def plain(text: str) -> str:
    text = re.sub(r"\*\*(.+?)\*\*", r"\1", text)
    text = re.sub(r"(?<![*\w])\*([^*\s][^*]*?)\*(?![*\w])", r"\1", text)
    return text.replace("$", r"\$")


def contact_sheet(path: Path) -> list[dict]:
    plt.rcParams.update({"pdf.fonttype": 42, "font.family": "Liberation Sans", "text.usetex": False})
    rows = []
    with PdfPages(path) as pdf:
        for page, (name, title, legend, png) in enumerate(legends(), start=1):
            fig = plt.figure(figsize=(210 / 25.4, 297 / 25.4))
            image = plt.imread(png)
            h, w = image.shape[:2]
            width_mm = 170.0
            height_mm = min(width_mm * h / w, 175.0)
            width_mm = height_mm * w / h
            ax = fig.add_axes(((210 - width_mm) / 2 / 210, 1 - (20 + height_mm) / 297, width_mm / 210, height_mm / 297))
            ax.imshow(image)
            ax.axis("off")
            top = 1 - (26 + height_mm) / 297
            fig.text(20 / 210, top, plain(f"{name}. {title}"), fontsize=9, fontweight="bold", va="top")
            wrapped = "\n".join(textwrap.wrap(plain(legend), 118))
            fig.text(20 / 210, top - 6 / 297, wrapped, fontsize=8, va="top", linespacing=1.35)
            fig.text(0.5, 10 / 297, f"{page}", fontsize=8, ha="center")
            pdf.savefig(fig)
            plt.close(fig)
            rows.append({"page": page, "figure": name, "title_words": len(title.split()),
                         "legend_words": len(plain(legend).split()), "png": str(png.relative_to(ROOT))})
    return rows


def pdf_fonts(path: Path) -> tuple[int, int]:
    out = subprocess.run(["pdffonts", str(path)], capture_output=True, text=True, check=True).stdout.splitlines()[2:]
    emb = [line.split()[-5] for line in out if line.strip()]
    return len(emb), sum(e == "yes" for e in emb)


def main() -> None:
    placed_main = write_docx(BMC / "A_main.md", BMC / "A_main.docx")
    figures = {f"Figure S{i}": FIGS / f"FigureS{i}.png" for i in range(1, 7)}
    placed_af1 = write_docx(BMC / "Additional_file_1.md", BMC / "Additional_file_1.docx", figures)
    write_letter(BMC / "cover_letter.md", BMC / "cover_letter.docx")
    rows = [docx_check(BMC / "A_main.docx", placed_main), docx_check(BMC / "Additional_file_1.docx", placed_af1)]
    pd.DataFrame(rows).to_csv(CHECKS / "docx_check_A19.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).T.to_string())
    sheet = BMC / "figures_contact_sheet.pdf"
    pages = contact_sheet(sheet)
    fonts, embedded = pdf_fonts(sheet)
    for r in pages:
        r.update(fonts=fonts, fonts_embedded=embedded, sheet_bytes=sheet.stat().st_size)
    pd.DataFrame(pages).to_csv(CHECKS / "contact_sheet_A19.tsv", sep="\t", index=False)
    print(pd.DataFrame(pages).to_string())


if __name__ == "__main__":
    main()
