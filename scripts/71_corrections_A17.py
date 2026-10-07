"""Apply the A17 §4 corrections found by 69_check_numbers_A17.py and list them.

Idempotent: each correction is computed from its source and written only where the file differs.
Every correction is appended to manuscript/checks/corrections_A17.tsv (doc, location, before, after, basis).
"""
from __future__ import annotations

import csv
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
AF1 = ROOT / "manuscript" / "bmc" / "Additional_file_1.md"
AF2 = ROOT / "manuscript" / "bmc" / "Additional_file_2.xlsx"
OUT = ROOT / "manuscript" / "checks" / "corrections_A17.tsv"
POWER = ROOT / "results" / "stage6" / "power.tsv"


def half_up(text: str, places: int) -> str:
    return f"{Decimal(text).quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP):f}"


def load_log() -> list[list[str]]:
    if not OUT.exists():
        return []
    with OUT.open() as f:
        rows = list(csv.reader(f, delimiter="\t"))
    return rows[1:]


def power_table(log: list[list[str]]) -> int:
    lines = AF1.read_text().splitlines(keepends=True)
    start = next(i for i, l in enumerate(lines) if l.startswith("| hypothesis | scenario | n | p_favor"))
    with POWER.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    changed = 0
    for k, row in enumerate(rows):
        i = start + 2 + k
        cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
        if cells[0] != row["hypothesis"] or cells[1] != row["scenario"]:
            raise SystemExit(f"power table row {k} is {cells[:2]}")
        new = list(cells)
        new[3] = half_up(row["p_favor"], 3)
        new[4] = half_up(row["p_against"], 3)
        new[5] = half_up(row["alpha"], 2 if row["hypothesis"].startswith("AH") else 3)
        for j, name in ((3, "p_favor"), (4, "p_against"), (5, "alpha")):
            if new[j] != cells[j]:
                log.append(["af1", f"Additional file 1, section 10, auxiliary design power table, row "
                                   f"{row['hypothesis']} scenario {row['scenario']}, column {name}", cells[j], new[j],
                            f"results/stage6/power.tsv {name} = {row[name]}; display rule: 3 decimals, ROUND_HALF_UP "
                            "(alpha of the AH rows kept at 2 decimals, 0.05)"])
                changed += 1
        lines[i] = "| " + " | ".join(new) + " |\n"
    AF1.write_text("".join(lines))
    return changed


def readme(log: list[list[str]]) -> int:
    book = load_workbook(AF2)
    sheet = book["README"]
    changed = 0
    for row in sheet.iter_rows(min_row=2):
        name, text = row[0].value, row[1].value
        if not name or name not in book.sheetnames:
            continue
        title = book[name].cell(1, 1).value
        if text != title:
            log.append(["af2", f"Additional file 2, README, row for sheet {name}", text, title,
                        f"row 1 of sheet {name}; the README repeats the status line of each sheet"])
            row[1].value = title
            changed += 1
    if changed:
        book.save(AF2)
    return changed


def main() -> None:
    log = load_log()
    n_power = power_table(log)
    n_readme = readme(log)
    with OUT.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["doc", "location", "before", "after", "basis"])
        w.writerows(log)
    print(f"power cells changed {n_power}; README rows changed {n_readme}; corrections listed {len(log)}")


if __name__ == "__main__":
    main()
