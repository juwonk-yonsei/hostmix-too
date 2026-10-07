"""Working lists for writing manuscript/checks/links_A17.tsv (A17 section 4.6).

Writes, outside the repository,
  <out>/unlinked_<doc>.txt   every sentence that still has unlinked numbers: anchor, section, starts,
                             the live tokens in order and the sentence
  <out>/fact_catalog.tsv     fact id, display, meaning, keys and status of every fact, sorted by id
  <out>/cell_catalog.tsv     table cell ids with display and meaning
The lists only show what exists; the link of each token is chosen by reading the sentence.
Usage: 72_link_helper_A17.py <out dir> [doc ...]
"""
from __future__ import annotations

import importlib
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
CHECK = importlib.import_module("69_check_numbers_A17")


def main() -> None:
    out = Path(sys.argv[1])
    wanted = set(sys.argv[2:])
    out.mkdir(parents=True, exist_ok=True)
    checker = CHECK.Checker()
    checker.run()
    by_sentence = defaultdict(list)
    for c in checker.claims:
        if c["link_type"] == "unlinked" and not c["doc"].startswith("fig:"):
            by_sentence[(c["doc"], c["key"].rsplit("|", 1)[0])].append(c)
    tokens = {t.key: t for d in checker.docs.values() for t in d.tokens}
    for name, doc in checker.docs.items():
        if wanted and name not in wanted:
            continue
        lines = []
        for (para, sent), sentence in sorted(doc.sentences.items()):
            claims = by_sentence.get((name, f"{name}|p{para}|s{sent}"))
            if not claims:
                continue
            live = [t for t in doc.tokens if not t.table and (t.para, t.sent) == (para, sent) and not t.exempt]
            section = doc.sections.get(para, "")
            lines.append(f"## {name}\tp{para}|s{sent}\t{section}\n"
                         f"starts: {sentence[:40]!r}\n"
                         f"tokens: {' '.join(f'[{t.text}]' for t in live)}\n"
                         f"{sentence}\n")
        table_claims = [c for c in checker.claims if c["doc"] == name and c["link_type"] == "unlinked"
                        and tokens.get(c["key"]) is not None and tokens[c["key"]].table]
        for c in table_claims:
            t = tokens[c["key"]]
            lines.append(f"## {name}\ttable {t.table} r{t.row} c{t.col}\t{c['note']}\n{t.sentence}\n")
        (out / f"unlinked_{name}.txt").write_text("\n".join(lines))
    facts = []
    for fid, row in sorted(checker.book.rows.items()):
        res = checker.book.result.get(fid, {})
        facts.append({"fact_id": fid, "display": res.get("display", res.get("error", "")), "meaning": row["meaning"],
                      "cohort": row["cohort"], "method": row["method"], "set": row["set"], "status": row["status"],
                      "transform": row["transform"]})
    pd.DataFrame(facts).to_csv(out / "fact_catalog.tsv", sep="\t", index=False)
    cells = [{"id": k, "display": v.get("display", ""), "meaning": v.get("row", {}).get("meaning", "")}
             for k, v in sorted(checker.cell_facts.items())]
    pd.DataFrame(cells).to_csv(out / "cell_catalog.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
