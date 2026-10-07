#!/usr/bin/env python3
"""A19 section 2: renumber the main figures in order of first citation (old 6 -> 4, old 4 -> 5, old 5 -> 6).

Every reference to a main figure is first replaced by a placeholder that carries the old number and then by the
new number, so that no number is replaced twice. Patterns (S figures are never matched):
  "Fig. 6a", "Figs. 4", "Figure 4", "Figure4", "Fig4_", "Fig6b", "Fig4/" (label keys), "fig4" (function names).

Steps
  1. SHA-256 of figures/Figure4-6.pdf/.png before the move.
  2. git mv, through temporary names: figures/Figure{4,5,6}.{pdf,png}, figure_source/Fig{4,5,6}*.tsv,
     checks/table_specs/AF2_Fig{4,5,6}*.yaml.
  3. Placeholder renumbering of the text files listed in TEXT (in place) and of
     scripts/76_figures_A18.py into scripts/84_figures_A19.py and links_A18.tsv into links_A19.tsv.
  4. Main text: legends of the figures in the new order 1-6.
  5. Additional file 2: sheet names, row 1 of every sheet, README rows; sheets in the new figure order.
  6. Citation map: every reference to a main figure in the main text, its legends and Additional file 1, with the
     old and new label, written to manuscript/checks/figure_renumber_map_A19.tsv; hashes after the move to
     manuscript/checks/figure_renumber_hash_A19.tsv.
"""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
BMC = ROOT / "manuscript" / "bmc"
CHECKS = ROOT / "manuscript" / "checks"
NEW = {"4": "5", "5": "6", "6": "4"}
OLD_OF = {v: k for k, v in NEW.items()}
PATTERNS = [
    re.compile(r"(Figs?\.\s?)([456])(?![0-9])"),
    re.compile(r"(?<![A-Za-z])(Figures?\s?)([456])(?![0-9])"),
    re.compile(r"(?<![A-Za-z])(Fig)([456])(?![0-9])"),
    re.compile(r"(?<![A-Za-z])(fig)([456])(?![0-9])"),
]
TEXT = [BMC / "A_main.md", BMC / "Additional_file_1.md", BMC / "cover_letter.md",
        CHECKS / "figure_labels.tsv", ROOT / "scripts" / "figcheck_A17.py",
        ROOT / "scripts" / "fig_expected_A17.py", ROOT / "scripts" / "af2_expected_A17.py"]
DERIVED = [(ROOT / "scripts" / "76_figures_A18.py", ROOT / "scripts" / "84_figures_A19.py"),
           (CHECKS / "links_A18.tsv", CHECKS / "links_A19.tsv")]
BACKUP = Path("/tmp/a19/before_renumber")


def renumber(text: str) -> str:
    for rx in PATTERNS:
        text = rx.sub(lambda m: f"{m.group(1)}\x00{m.group(2)}\x00", text)
    return re.sub(r"\x00([456])\x00", lambda m: NEW[m.group(1)], text)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_mv_cycle(pairs: list[tuple[Path, Path]]) -> None:
    temps = []
    for old, new in pairs:
        tmp = old.with_name("__a19_tmp_" + new.name)
        subprocess.run(["git", "mv", str(old), str(tmp)], cwd=ROOT, check=True)
        temps.append((tmp, new))
    for tmp, new in temps:
        subprocess.run(["git", "mv", str(tmp), str(new)], cwd=ROOT, check=True)


def renamed(path: Path) -> Path:
    return path.with_name(renumber(path.name))


def citation_map(main_before: str, af1_before: str, titles: dict[str, str]) -> pd.DataFrame:
    rows = []
    for doc, text in (("A_main.md", main_before), ("Additional_file_1.md", af1_before)):
        legend_start = text.find("## Figure titles and legends")
        for m in re.finditer(r"(Figs?\.\s?|\*\*Figure\s)([1-6])([a-d](?:[,–-]\s?[a-d])*)?(?![0-9])", text):
            line = text.count("\n", 0, m.start()) + 1
            start = text.rfind("\n", 0, m.start()) + 1
            old = m.group(0).lstrip("*")
            new_number = NEW.get(m.group(2), m.group(2))
            where = "legend title" if m.group(1).startswith("**") else (
                "legends" if doc == "A_main.md" and 0 <= legend_start < m.start() else
                ("main text" if doc == "A_main.md" else "Additional file 1"))
            rows.append({"doc": doc, "line_before": line, "where": where, "sentence_start": text[start:start + 60],
                         "old": old, "new": old.replace(m.group(2), new_number, 1),
                         "new_title": titles[new_number]})
    return pd.DataFrame(rows)


def reorder_legends(text: str) -> str:
    head, rest = text.split("## Figure titles and legends\n", 1)
    body, tail = rest.split("\n## Tables", 1)
    blocks = re.split(r"\n(?=\*\*Figure \d\. )", body.strip("\n"))
    blocks = sorted(blocks, key=lambda b: int(re.match(r"\*\*Figure (\d)\. ", b).group(1)))
    return head + "## Figure titles and legends\n\n" + "\n".join(b.rstrip("\n") + "\n" for b in blocks) + "\n## Tables" + tail


def additional_file_2() -> list[dict]:
    path = BMC / "Additional_file_2.xlsx"
    book = load_workbook(path)
    log = []
    for ws in book.worksheets:
        old_title = ws.title
        new_title = renumber(old_title)
        if new_title != old_title:
            ws.title = new_title
        if ws.title == "README":
            for row in ws.iter_rows(min_row=2):
                for cell in row[:2]:
                    if isinstance(cell.value, str) and renumber(cell.value) != cell.value:
                        log.append({"sheet": "README", "cell": cell.coordinate, "before": cell.value,
                                    "after": renumber(cell.value)})
                        cell.value = renumber(cell.value)
        else:
            cell = ws.cell(1, 1)
            if isinstance(cell.value, str) and renumber(cell.value) != cell.value:
                log.append({"sheet": ws.title, "cell": "A1", "before": cell.value, "after": renumber(cell.value)})
                cell.value = renumber(cell.value)
        if new_title != old_title:
            log.append({"sheet": new_title, "cell": "title", "before": old_title, "after": new_title})

    def rank(title: str) -> tuple:
        m = re.match(r"Fig(\d)", title)
        if title == "README":
            return (0, 0)
        if m:
            return (1, int(m.group(1)))
        return (2, 0)

    order = sorted(range(len(book.worksheets)), key=lambda i: (rank(book.worksheets[i].title), i))
    book._sheets = [book.worksheets[i] for i in order]
    readme = book["README"]
    rows = [[c.value for c in r] for r in readme.iter_rows(min_row=2)]
    position = {ws.title: i for i, ws in enumerate(book.worksheets)}
    rows.sort(key=lambda r: position.get(r[0], len(position)))
    for i, values in enumerate(rows, start=2):
        for j, value in enumerate(values, start=1):
            readme.cell(i, j).value = value
    book.save(path)
    return log


def main() -> None:
    if (BMC / "figure_source" / "Fig5_ablation.tsv").exists():
        raise SystemExit("already renumbered")
    BACKUP.mkdir(parents=True, exist_ok=True)
    for p in TEXT + [BMC / "Additional_file_2.xlsx"] + [s for s, _ in DERIVED]:
        shutil.copy2(p, BACKUP / p.name)
    figs = BMC / "figures"
    before = {f"Figure{n}.{ext}": sha(figs / f"Figure{n}.{ext}") for n in "456" for ext in ("pdf", "png")}
    main_before = (BMC / "A_main.md").read_text()
    af1_before = (BMC / "Additional_file_1.md").read_text()

    git_mv_cycle([(figs / f"Figure{n}.{ext}", figs / f"Figure{NEW[n]}.{ext}") for n in "456" for ext in ("pdf", "png")])
    sources = sorted(p for p in (BMC / "figure_source").glob("Fig[456]*.tsv"))
    git_mv_cycle([(p, renamed(p)) for p in sources])
    specs = sorted(p for p in (CHECKS / "table_specs").glob("AF2_Fig[456]*.yaml"))
    git_mv_cycle([(p, renamed(p)) for p in specs])
    for p in specs:
        target = renamed(p)
        target.write_text(renumber(target.read_text()))

    for p in TEXT:
        p.write_text(renumber(p.read_text()))
    for src, dest in DERIVED:
        dest.write_text(renumber(src.read_text()))
    main_after = reorder_legends((BMC / "A_main.md").read_text())
    (BMC / "A_main.md").write_text(main_after)

    titles = {m.group(1): m.group(2) for m in re.finditer(r"^\*\*Figure (\d)\. (.+?)\*\*$", main_after, re.M)}
    table = citation_map(main_before, af1_before, titles)
    table.to_csv(CHECKS / "figure_renumber_map_A19.tsv", sep="\t", index=False)
    log = additional_file_2()
    pd.DataFrame(log).to_csv(CHECKS / "figure_renumber_af2_A19.tsv", sep="\t", index=False)

    hashes = []
    for name, digest in before.items():
        new_name = renumber(name)
        after = sha(figs / new_name)
        hashes.append({"old_name": name, "new_name": new_name, "sha256_before": digest, "sha256_after": after,
                       "equal": digest == after})
    pd.DataFrame(hashes).to_csv(CHECKS / "figure_renumber_hash_A19.tsv", sep="\t", index=False)
    print(pd.DataFrame(hashes)[["old_name", "new_name", "equal"]].to_string(index=False))
    print("moved:", [p.name for p in sources], [p.name for p in specs])
    print("citations:", len(table), "AF2 edits:", len(log))


if __name__ == "__main__":
    main()
