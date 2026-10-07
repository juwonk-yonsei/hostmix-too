#!/usr/bin/env python3
"""A18 section 5: review pack in manuscript/bmc/review_pack/ (A19: A19 check and hygiene outputs).

Contents: A_main.md, Additional_file_1.md, cover_letter.md, figures_contact_sheet.pdf, figures_current/ (the 12
figure PNGs), claims_sample.md (40 random checked positions, seed 20261007, and every abstract position, from the
A19 check), check_summary.json (check_summary_A19.json), hygiene_checklist.md (hygiene_checklist_A19.md) and
open_items.md (written by hand). Files of earlier packs that are not in this list are moved to
manuscript/archive/review_pack_A17/. The text files of the pack are scanned for POG570 patient identifiers.
"""
from __future__ import annotations

import importlib
import random
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
BMC = ROOT / "manuscript" / "bmc"
PACK = BMC / "review_pack"
ARCHIVE = ROOT / "manuscript" / "archive" / "review_pack_A17"
CHECKS = ROOT / "manuscript" / "checks"
SEED = 20261007
COPIES = {"A_main.md": BMC / "A_main.md", "Additional_file_1.md": BMC / "Additional_file_1.md",
          "cover_letter.md": BMC / "cover_letter.md", "figures_contact_sheet.pdf": BMC / "figures_contact_sheet.pdf",
          "check_summary.json": CHECKS / "check_summary_A19.json",
          "hygiene_checklist.md": CHECKS / "hygiene_checklist_A19.md"}
KEEP = set(COPIES) | {"figures_current", "claims_sample.md", "open_items.md"}
FIGURES = [f"Figure{i}.png" for i in range(1, 7)] + [f"FigureS{i}.png" for i in range(1, 7)]


def claims_sample() -> pd.DataFrame:
    sample = importlib.import_module("74_sample_A17")
    check = importlib.import_module("79_check_numbers_A18")
    sample.CHECK = check
    checker = check.Checker()
    checker.run()
    claims = sorted(checker.claims, key=lambda c: c["key"])
    picks = random.Random(SEED).sample(claims, 40)
    abstract = [c for c in checker.claims if c["doc"] == "main" and c["section"].split(" > ")[0] == "Abstract"]
    frame = pd.DataFrame([sample.row(checker, c, "random") for c in picks]
                         + [sample.row(checker, c, "abstract") for c in abstract])
    cols = ["key", "displayed", "context", "link_id", "meaning", "source_file", "selector", "column", "raw", "ok"]
    lines = ["# Sample of the number check (A19)", "",
             f"Random positions: 40 of {len(claims)} (seed {SEED}). Abstract positions: {len(abstract)}.", ""]
    for group, title in (("random", "Random 40"), ("abstract", "Abstract")):
        part = frame.loc[frame["group"].eq(group)]
        lines += [f"## {title}", "", "| # | " + " | ".join(cols) + " |", "|" + "---|" * (len(cols) + 1)]
        for i, r in enumerate(part.itertuples(index=False), 1):
            lines.append(f"| {i} | " + " | ".join(sample.cell(getattr(r, c)) for c in cols) + " |")
        lines.append("")
    (PACK / "claims_sample.md").write_text("\n".join(lines))
    frame.to_csv(CHECKS / "claims_sample_A19.tsv", sep="\t", index=False)
    return frame


def main() -> None:
    PACK.mkdir(exist_ok=True)
    moved = []
    for p in sorted(PACK.iterdir()):
        if p.name not in KEEP:
            ARCHIVE.mkdir(parents=True, exist_ok=True)
            shutil.move(str(p), ARCHIVE / p.name)
            moved.append(p.name)
    for name, src in COPIES.items():
        shutil.copy2(src, PACK / name)
    current = PACK / "figures_current"
    current.mkdir(exist_ok=True)
    for p in current.iterdir():
        if p.name not in FIGURES:
            p.unlink()
    for name in FIGURES:
        shutil.copy2(BMC / "figures" / name, current / name)
    frame = claims_sample()
    hygiene = importlib.import_module("81_hygiene_A18")
    scan = hygiene.pog570_ids_in(PACK)
    print("moved to archive:", moved)
    print("pack:", sorted(p.name for p in PACK.iterdir()), "figures_current:", len(list(current.iterdir())))
    print("claims sample:", len(frame), "not ok:", int((~frame["ok"].astype(bool)).sum()))
    print("POG570 patient identifiers per text file of the pack:", scan)


if __name__ == "__main__":
    main()
