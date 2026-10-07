"""Sample report of the number check (A17 section 4.13).

40 random positions (seed 20261007) drawn from every checked position of the documents and figures,
plus every position of the abstract. For each: the display, 60 characters on each side, the link id,
the meaning, the source file, the selector, the column and the raw source value.

Outputs: manuscript/checks/claims_sample_A17.tsv, manuscript/checks/claims_sample_A17.md
"""
from __future__ import annotations

import importlib
import random
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
CHECK = importlib.import_module("69_check_numbers_A17")
NC = CHECK.NC
OUT_TSV = NC.CHECKS / "claims_sample_A17.tsv"
OUT_MD = NC.CHECKS / "claims_sample_A17.md"
SEED = 20261007


def context(checker, claim: dict) -> str:
    key = claim["key"]
    for name, doc in checker.docs.items():
        for t in doc.tokens:
            if t.key == key:
                pos = CHECK.locate(doc, t)
                if pos is None:
                    return t.sentence[max(0, t.start - 60):t.end + 60]
                text = doc.text
                return text[max(0, pos - 60):pos + len(t.text) + 60].replace("\n", " ")
    t = checker.fig_tokens.get(key)
    if t is not None:
        return t.sentence[max(0, t.start - 60):t.end + 60]
    return ""


def raw_value(checker, claim: dict) -> str:
    if claim["link_type"] != "derived":
        return claim["raw"]
    parts = []
    for name in dict.fromkeys(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", claim["raw"])):
        if name in checker.book.rows:
            parts.append(f"{name}={checker.book.get(name)['raw']}")
        elif name in checker.constants:
            parts.append(f"{name}={checker.constants[name]['value']}")
    return "; ".join(parts)


def row(checker, claim: dict, group: str) -> dict:
    return {"group": group, "key": claim["key"], "displayed": claim["displayed"], "context": context(checker, claim),
            "link_type": claim["link_type"], "link_id": claim["link_id"], "ok": claim["ok"],
            "meaning": claim["meaning"], "source_file": claim["source_file"], "selector": claim["selector"],
            "column": claim["column"], "raw": raw_value(checker, claim), "status": claim["status"]}


def cell(text) -> str:
    return str(text).replace("|", "\\|").replace("\n", " ")


def main() -> None:
    checker = CHECK.Checker()
    checker.run()
    claims = sorted(checker.claims, key=lambda c: c["key"])
    picks = random.Random(SEED).sample(claims, 40)
    abstract = [c for c in checker.claims if c["doc"] == "main" and c["section"].split(" > ")[0] == "Abstract"]
    rows = [row(checker, c, "random") for c in picks] + [row(checker, c, "abstract") for c in abstract]
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT_TSV, sep="\t", index=False)
    cols = ["key", "displayed", "context", "link_id", "meaning", "source_file", "selector", "column", "raw", "ok"]
    lines = ["# A17 sample report", "", f"Random positions: 40 of {len(claims)} (seed {SEED}). "
             f"Abstract positions: {len(abstract)}.", ""]
    for group, title in (("random", "Random 40"), ("abstract", "Abstract")):
        part = frame.loc[frame["group"].eq(group)]
        lines += [f"## {title}", "", "| # | " + " | ".join(cols) + " |", "|" + "---|" * (len(cols) + 1)]
        for i, r in enumerate(part.itertuples(index=False), 1):
            lines.append(f"| {i} | " + " | ".join(cell(getattr(r, c)) for c in cols) + " |")
        lines.append("")
    OUT_MD.write_text("\n".join(lines))
    print(frame.groupby(["group", "link_type"]).size().to_string())
    print("not ok:", int((~frame["ok"].astype(bool)).sum()))


if __name__ == "__main__":
    main()
