"""A18 section 3.3 and 3.4: edit sheets of Additional files 2 and 3 in place.

Additional file 2 (openpyxl, cell values only; other sheets untouched):
  Power                 row 1 text; columns recorded_in, commit and commit_time of the commit that added
                        the preregistration recording the power table (from git)
  S1_heatmap            status column rewritten per row from figure_source/S1_heatmap.tsv; sentence added to row 1
  S3_esophagus_masking  microarray rows from results/stage10/derived/s3a_microarray_layer.tsv
  S6_meta_analysis      column I2_percent (I2 x 100) after I2; the I2 column is kept
  README                row-1 texts of the edited sheets
Additional file 3: rows of the history sheet sorted by commit time (ties by git topological order);
only cell values of the history rows are moved, the Deviations_en sheet is not touched.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
AF2 = ROOT / "manuscript/bmc/Additional_file_2.xlsx"
AF3 = ROOT / "manuscript/bmc/Additional_file_3.xlsx"
POWER_RECORD = "config/prereg_A6.md"
POWER_TITLE = "Design power, computed before the corresponding unlock; supports the preregistrations."
S1_NOTE = "Status is given per row in the column 'status'."


def git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()


def header(ws, row: int = 2) -> dict[str, int]:
    return {c.value: c.column for c in ws[row] if c.value is not None}


def set_readme(book, sheet: str, text: str) -> None:
    for row in book["README"].iter_rows(min_row=2):
        if row[0].value == sheet:
            row[1].value = text
            return
    raise SystemExit(f"README has no row for {sheet}")


def power(book) -> dict:
    ws = book["Power"]
    added = git("log", "--diff-filter=A", "--format=%H\t%ci", "--", POWER_RECORD).splitlines()
    if len(added) != 1:
        raise SystemExit(f"{POWER_RECORD}: {len(added)} adding commits")
    commit, when = added[0].split("\t")
    cols = header(ws)
    for name in ("recorded_in", "commit", "commit_time"):
        if name not in cols:
            ws.cell(2, ws.max_column + 1, name)
            cols = header(ws)
    for r in range(3, ws.max_row + 1):
        if ws.cell(r, 1).value is None:
            continue
        ws.cell(r, cols["recorded_in"], POWER_RECORD)
        ws.cell(r, cols["commit"], commit)
        ws.cell(r, cols["commit_time"], when)
    ws["A1"] = POWER_TITLE
    set_readme(book, "Power", POWER_TITLE)
    return {"sheet": "Power", "commit": commit, "commit_time": when, "rows": ws.max_row - 2}


def s1_heatmap(book) -> dict:
    ws = book["S1_heatmap"]
    src = pd.read_csv(ROOT / "manuscript/bmc/figure_source/S1_heatmap.tsv", sep="\t")
    status = {(r.cohort, r.method): r.status for r in src.itertuples()}
    cols = header(ws)
    changed = 0
    for r in range(3, ws.max_row + 1):
        key = (ws.cell(r, cols["cohort"]).value, ws.cell(r, cols["method"]).value)
        if key not in status:
            raise SystemExit(f"S1_heatmap row {r}: {key} not in figure_source")
        cell = ws.cell(r, cols["status"])
        changed += cell.value != status[key]
        cell.value = status[key]
    title = ws["A1"].value
    if S1_NOTE not in title:
        ws["A1"] = f"{title} {S1_NOTE}"
    set_readme(book, "S1_heatmap", ws["A1"].value)
    return {"sheet": "S1_heatmap", "status_cells_changed": changed,
            "status_counts": pd.Series(list(status.values())).value_counts().to_dict()}


def s3(book) -> dict:
    ws = book["S3_esophagus_masking"]
    layer = pd.read_csv(ROOT / "results/stage10/derived/s3a_microarray_layer.tsv", sep="\t").set_index("method")
    cols = header(ws)
    before, after = {}, {}
    for r in range(3, ws.max_row + 1):
        if ws.cell(r, cols["analysis"]).value != "microarray":
            continue
        method = ws.cell(r, cols["method"]).value
        hit = layer.loc[method]
        before[method] = (ws.cell(r, cols["n_truth_not_esophagus"]).value, ws.cell(r, cols["esophagus_absorption"]).value)
        ws.cell(r, cols["n_truth_not_esophagus"], int(hit["n_truth_not_esophagus"]))
        ws.cell(r, cols["esophagus_absorption"], float(f"{float(hit['esophagus_fraction']):.10g}"))
        after[method] = (int(hit["n_truth_not_esophagus"]), float(f"{float(hit['esophagus_fraction']):.10g}"))
    if set(before) != {"BASE-Z", "SA-Z", "LD-Z"}:
        raise SystemExit(f"S3 microarray rows: {sorted(before)}")
    return {"sheet": "S3_esophagus_masking", "before": before, "after": after}


def s6(book) -> dict:
    ws = book["S6_meta_analysis"]
    cols = header(ws)
    if "I2_percent" not in cols:
        ws.insert_cols(cols["I2"] + 1)
        ws.cell(2, cols["I2"] + 1, "I2_percent")
        cols = header(ws)
    filled = 0
    for r in range(3, ws.max_row + 1):
        value = ws.cell(r, cols["I2"]).value
        ws.cell(r, cols["I2_percent"], None if value is None else float(value) * 100)
        filled += value is not None
    return {"sheet": "S6_meta_analysis", "I2_percent_rows": filled, "columns": list(header(ws))}


def af3() -> dict:
    book = load_workbook(AF3)
    ws = book["history"]
    cols = header(ws)
    rows = [[c.value for c in row] for row in ws.iter_rows(min_row=3, max_row=ws.max_row)]
    rows = [r for r in rows if any(v is not None for v in r)]
    topo = {h: i for i, h in enumerate(git("rev-list", "--topo-order", "--reverse", "HEAD").splitlines())}
    before = [r[cols["commit"] - 1][:8] for r in rows]
    ordered = sorted(rows, key=lambda r: (pd.Timestamp(r[cols["time"] - 1]), topo[r[cols["commit"] - 1]]))
    for i, values in enumerate(ordered):
        for j, value in enumerate(values):
            ws.cell(3 + i, 1 + j).value = value
    sheets = book.sheetnames
    book.save(AF3)
    return {"file": "Additional_file_3", "rows": len(rows), "before": before,
            "after": [r[cols["commit"] - 1][:8] for r in ordered], "sheets": sheets}


def main() -> None:
    book = load_workbook(AF2)
    out = [power(book), s1_heatmap(book), s3(book), s6(book)]
    book.save(AF2)
    out.append(af3())
    for item in out:
        print(item)


if __name__ == "__main__":
    main()
