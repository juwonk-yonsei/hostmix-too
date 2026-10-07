"""Count number displays by document, section and form (A17 section 4.2).

Usage: 67_inventory_A17.py start|end. Writes manuscript/checks/inventory_<label>_A17.tsv with
the number of tokens, the tokens exempted by a fixed rule (by category) and the rest.
At the end the checker's claims map supplies the link of every token.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numdoc_A17 as N  # noqa: E402

CHECKS = N.ROOT / "manuscript/checks"


def main() -> None:
    label = sys.argv[1]
    docs = N.load_all()
    links = {}
    claims = CHECKS / "claims_map.tsv"
    if label == "end" and claims.exists():
        table = pd.read_csv(claims, sep="\t", dtype=str, keep_default_na=False)
        links = dict(zip(table["key"], table["link_type"]))
    rows = Counter()
    for name, doc in docs.items():
        for token in doc.tokens:
            section = token.section.split(" > ")[0]
            if label == "end":
                state = links.get(token.key, "unlinked")
            else:
                state = f"exempt {token.exempt}" if token.exempt else "to link"
            rows[(name, section, token.kind, state)] += 1
    table = pd.DataFrame([{"doc": d, "section": s, "type": k, "state": st, "n": n}
                          for (d, s, k, st), n in sorted(rows.items())])
    out = CHECKS / f"inventory_{label}_A17.tsv"
    table.to_csv(out, sep="\t", index=False)
    print(out, len(table), "rows")
    print(table.groupby(["doc", "state"])["n"].sum().to_string())


if __name__ == "__main__":
    main()
