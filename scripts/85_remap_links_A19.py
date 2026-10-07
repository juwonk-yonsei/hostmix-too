#!/usr/bin/env python3
"""A19: move the paragraph anchors of manuscript/checks/links_A19.tsv to the renumbered main text.

The legends of the main figures were put in the new order 1-6, so the paragraph numbers in the legend section
changed. Each linked sentence is read in the main text before the renumbering (backup made by
83_renumber_figures_A19.py), renumbered with the same rule, and looked up in the current main text; the anchor
becomes the paragraph and sentence where the identical sentence now stands. Anchors that cannot be found exactly
once stop the script. Rows of other documents are left unchanged.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import numdoc_A17 as ND  # noqa: E402

spec = importlib.util.spec_from_file_location("renumber", ROOT / "scripts" / "83_renumber_figures_A19.py")
RN = importlib.util.module_from_spec(spec)
spec.loader.exec_module(RN)

LINKS = ROOT / "manuscript" / "checks" / "links_A19.tsv"
BEFORE = Path("/tmp/a19/before_renumber/A_main.md")


def main() -> None:
    old = ND.load("main", BEFORE.read_text())
    new = ND.load("main")
    where = {}
    for key, sentence in new.sentences.items():
        where.setdefault(sentence, []).append(key)
    lines = LINKS.read_text().split("\n")
    moved = []
    for i, line in enumerate(lines[1:], start=1):
        fields = line.split("\t")
        if fields[0] != "main":
            continue
        para, sent = (int(x[1:]) for x in fields[1].split("|"))
        sentence = RN.renumber(old.sentences[(para, sent)])
        hits = where.get(sentence, [])
        if len(hits) != 1:
            raise SystemExit(f"{fields[1]}: {len(hits)} matches for {sentence[:60]!r}")
        anchor = f"p{hits[0][0]}|s{hits[0][1]}"
        if anchor != fields[1]:
            moved.append({"before": fields[1], "after": anchor, "sentence": sentence[:60]})
            fields[1] = anchor
            lines[i] = "\t".join(fields)
    LINKS.write_text("\n".join(lines))
    print(pd.DataFrame(moved).to_string(index=False) if moved else "no anchor moved")


if __name__ == "__main__":
    main()
