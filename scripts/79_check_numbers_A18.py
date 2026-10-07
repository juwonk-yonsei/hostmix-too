"""A18 number check: the A17 check (69_check_numbers_A17.py) with the A18 changes.

Every number display is linked to a fact, a table cell, a constant or an exemption, and every
link is recomputed from its source file. Changes from A17:
  - figure value labels are matched by element key and position: the PDF words inside the extent
    recorded in manuscript/checks/figure_labels.tsv must spell the recorded text, and their numbers
    are linked to the expected item with the same key (figcheck_A17.py); ticks keep text matching;
  - a link "F:a&F:b" requires the number to equal each of the facts (sentences that give one
    value for two models);
  - outputs carry the suffix _A19; links are read from links_A19.tsv (A19: main figures renumbered,
    old 6 -> 4, old 4 -> 5, old 5 -> 6; the figure expectations follow the new numbers in figcheck_A17.py,
    fig_expected_A17.py and the table specs; checking logic unchanged).

Inputs
  manuscript/checks/facts.tsv            written by 68_facts_A17.py (rebuilt here in memory)
  manuscript/checks/constants.tsv        design constants with their configuration file and line
  manuscript/checks/table_specs/*.yaml   row and column maps of Tables 1-3, S1-S9, the model card tables
  manuscript/checks/links_A19.tsv        one row per sentence: the link of each token, in token order
  manuscript/checks/figure_labels.tsv    value labels of the figures with their extents (84_figures_A19.py)
Outputs (manuscript/checks/)
  claims_map_A19.tsv, check_issues_A19.tsv, key_flags_A19.tsv, status_flags_A19.tsv,
  rounding_boundaries_A19.tsv, check_summary_A19.json, af2_check_A19.tsv,
  figure_check_A19.tsv (printed numbers of each figure PDF), figure_label_check_A19.tsv (one row per
  value label), figure_source_check_A19.tsv (figure_source files against results)

Modes
  (default)            full check
  --mutate             change the last digit of 15 linked numbers of the documents and of 5 figure value
                       labels in the label table (seed 20261008), in memory, and require exactly those
                       20 to fail
  --override           override the stored value of 5 facts in memory and require every position that
                       reads those cells to fail
  --write-manifest     write source_manifest.tsv (SHA-256 of every file the check reads); the default
                       check compares those hashes with the files on disk
"""
from __future__ import annotations

import hashlib
import importlib
import json
import random
import re
import sys
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numcheck_A17 as NC  # noqa: E402
import numdoc_A17 as N
import af2_expected_A17 as AF2  # noqa: E402

FACTS_MOD = importlib.import_module("68_facts_A17")
SPECS = NC.CHECKS / "table_specs"
LINKS = NC.CHECKS / "links_A19.tsv"
TAG = "A19"
LABEL_TOL_MM = 0.6
CONSTANTS = NC.CHECKS / "constants.tsv"


class Book:
    """Evaluated facts: value, display, raw text and the source cells read."""

    def __init__(self, src: NC.Sources, rows: list[dict], constants: dict | None = None) -> None:
        self.src = src
        self.constants = constants or {}
        self.rows = {r["fact_id"]: r for r in rows}
        self.values: dict[str, Decimal] = {}
        self.result: dict[str, dict] = {}
        self.errors: list[dict] = []
        self.empty: list[str] = []
        for fid in list(self.rows):
            self.get(fid)
        for fid in self.empty:
            del self.rows[fid]
            del self.result[fid]

    def get(self, fid: str) -> dict | None:
        if fid in self.result:
            return self.result[fid]
        row = self.rows.get(fid)
        if row is None:
            return None
        try:
            deps: set = set()
            facts = None
            if row["transform"].startswith("derived:"):
                facts = {}
                names = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", row["transform"][8:]))
                for name in names:
                    if name in self.rows:
                        sub = self.get(name)
                        if sub is None or "value" not in sub:
                            raise NC.SourceError(f"derived input {name} failed")
                        facts[name] = sub["value"]
                        deps |= sub["cells"]
                    elif name in self.constants:
                        const = self.constants[name]
                        facts[name] = Decimal(const["value"])
                        deps.add((const["source_file"], int(const["line"] or 0), "constant"))
            value, raw, index = NC.evaluate(self.src, row["source_file"], row["selector"], row["column"],
                                            row["transform"], facts, FACTS_MOD.LISTS)
            col = row["column"] or row["transform"].split(":", 1)[-1]
            if not row["transform"].startswith("derived:"):
                deps |= {(row["source_file"], int(i), col) for i in index}
                if row["transform"] == "count":
                    deps |= {(row["source_file"], -1, "count:" + row["selector"])}
            places = int(row["decimals"])
            out = {"value": value, "raw": raw, "display": NC.show(value, places, row["format"]),
                   "boundary": NC.boundary(value, places, row["format"]), "cells": deps}
        except NC.SourceError as exc:
            out = {"error": str(exc), "cells": set()}
            if str(exc) == "empty value" and row.get("note") != "specific":
                self.empty.append(fid)
            else:
                self.errors.append({"fact_id": fid, "error": str(exc)})
        self.result[fid] = out
        return out


def load_constants() -> tuple[dict, list]:
    table = pd.read_csv(CONSTANTS, sep="\t", dtype=str, keep_default_na=False)
    issues = []
    out = {}
    for r in table.itertuples(index=False):
        out[r.constant_id] = r._asdict()
        path = NC.ROOT / r.source_file
        if r.check.startswith("literal:"):
            want = r.check[len("literal:"):]
            lines = path.read_text().splitlines()
            line = int(r.line)
            if line < 1 or line > len(lines) or want not in lines[line - 1]:
                issues.append({"kind": "constant source", "where": r.constant_id,
                               "detail": f"line {r.line} of {r.source_file} does not contain {want!r}"})
        elif r.check.startswith("yamllen:"):
            data = yaml.safe_load(path.read_text())
            for part in r.check[len("yamllen:"):].split("."):
                data = data[part]
            if Decimal(len(data)) != Decimal(r.value):
                issues.append({"kind": "constant source", "where": r.constant_id,
                               "detail": f"{r.check} has {len(data)} items, constant says {r.value}"})
        elif r.check.startswith("yamlkeys:"):
            path_part, _, skip = r.check[len("yamlkeys:"):].partition("|")
            data = yaml.safe_load(path.read_text())
            for part in path_part.split("."):
                data = data[part]
            keys = [k for k in data if k not in set(skip.split(",")) - {""}]
            if Decimal(len(keys)) != Decimal(r.value):
                issues.append({"kind": "constant source", "where": r.constant_id,
                               "detail": f"{r.check} has {len(keys)} keys, constant says {r.value}"})
        elif r.check.startswith("linecount:"):
            rx = re.compile(r.check[len("linecount:"):])
            found = sum(bool(rx.search(line)) for line in path.read_text().splitlines())
            if Decimal(found) != Decimal(r.value):
                issues.append({"kind": "constant source", "where": r.constant_id,
                               "detail": f"{r.check} matches {found} lines, constant says {r.value}"})
        elif r.check.startswith("count:"):
            rx = re.compile(r.check[len("count:"):])
            lines = path.read_text().splitlines()
            line = int(r.line)
            found = len(rx.findall(lines[line - 1])) if 1 <= line <= len(lines) else -1
            if Decimal(found) != Decimal(r.value):
                issues.append({"kind": "constant source", "where": r.constant_id,
                               "detail": f"{r.check} matches {found} times on line {r.line}, constant says {r.value}"})
        elif r.check != "none":
            issues.append({"kind": "constant source", "where": r.constant_id, "detail": f"unknown check {r.check}"})
    return out, issues


MANIFEST = NC.CHECKS / "source_manifest.tsv"
REQUIRED = ["results/stage9/imvigor_contribution_fix1.tsv", "results/stage10/**/*", "release/models/*"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def check_manifest(read: dict[str, str]) -> dict:
    """Hashes of the manifest against the files on disk, and files read here that the manifest lacks."""
    if not MANIFEST.exists():
        return {"manifest_files": 0, "manifest_mismatch": 0, "manifest_absent": 0, "manifest_unlisted": len(read),
                "manifest_rows": []}
    table = pd.read_csv(MANIFEST, sep="\t", dtype=str, keep_default_na=False)
    rows = []
    for r in table.itertuples(index=False):
        path = NC.ROOT / r.path
        now = sha256(path) if path.exists() else ""
        rows.append({"path": r.path, "sha256": r.sha256, "now": now, "ok": now == r.sha256})
    listed = set(table["path"])
    return {"manifest_files": len(rows), "manifest_mismatch": sum(not r["ok"] and r["now"] != "" for r in rows),
            "manifest_absent": sum(r["now"] == "" for r in rows),
            "manifest_unlisted": len(set(read) - listed), "manifest_rows": rows}


def token_value(token) -> Decimal | None:
    if token.kind in {"word", "fraction"}:
        return Decimal(str(token.value)) if token.kind == "word" else None
    text = NC.norm_token(token.text)
    if token.kind == "scientific":
        m = re.match(r"^(-?)(\d+(?:\.\d+)?)\s*(?:×\s*10\^|e)(-?\d+)$", text)
        if m:
            return Decimal(m.group(1) + m.group(2)).scaleb(int(m.group(3)))
        m = re.match(r"^10\^(-?\d+)$", text)
        return Decimal(1).scaleb(int(m.group(1))) if m else None
    try:
        return Decimal(text)
    except Exception:  # noqa: BLE001
        return None


def compare(token, display: str, value: Decimal, modifier: str, approx: str) -> tuple[bool, str]:
    if modifier == "abs":
        value = abs(value)
        display = display.lstrip("-")
    if token.kind == "fraction":
        tol = Decimal(approx) if approx else Decimal("0")
        got = Decimal(str(token.value))
        return abs(got - value) <= tol, f"{got:.4f} within {tol} of {value:.4f}"
    if token.kind == "word":
        return Decimal(str(token.value)) == value, str(value)
    shown = NC.norm_token(token.text)
    return shown == NC.norm_token(display), display


class Checker:
    def __init__(self, texts: dict | None = None, overrides: dict | None = None,
                 labels: pd.DataFrame | None = None) -> None:
        self.label_table = labels
        self.label_rows: list[dict] = []
        self.src = NC.Sources()
        if overrides:
            self.src.overrides.update(overrides)
        FACTS_MOD.FACTS.clear()
        FACTS_MOD.SEEN.clear()
        self.fact_rows = FACTS_MOD.build()
        self.constants, self.issues = load_constants()
        self.book = Book(self.src, self.fact_rows, self.constants)
        self.docs = N.load_all(texts)
        self.claims: list[dict] = []
        self.key_flags: list[dict] = []
        self.status_flags: list[dict] = []
        self.cell_facts: dict[str, dict] = {}
        self.fig_cells: dict[str, dict] = {}
        self.fig_rows: list[dict] = []
        self.fig_tokens: dict = {}
        self.figure_source_rows: list[dict] = []
        self.unread_cells: list[dict] = []
        for e in self.book.errors:
            self.issues.append({"kind": "selector", "where": e["fact_id"], "detail": e["error"]})

    # ------------------------------------------------------------ resolution of one link id
    def resolve(self, link: str) -> dict:
        modifier = ""
        if "|" in link:
            link, modifier = link.split("|", 1)
        kind, _, ref = link.partition(":")
        if kind == "F":
            res = self.book.get(ref)
            row = self.book.rows.get(ref)
            if res is None:
                return {"error": f"unknown fact {ref}"}
            if "error" in res:
                return {"error": res["error"]}
            ltype = "derived" if row["transform"].startswith("derived:") else "fact"
            return {"type": ltype, "id": ref, "value": res["value"], "display": res["display"], "raw": res["raw"],
                    "cells": res["cells"], "modifier": modifier, "approx": row.get("approx", ""), "row": row,
                    "boundary": res["boundary"]}
        if kind == "T":
            cell = self.cell_facts.get(ref)
            if cell is None:
                return {"error": f"unknown table cell {ref}"}
            if "error" in cell:
                return {"error": cell["error"]}
            return {"type": "table_cell", "id": ref, "value": cell["value"], "display": cell["display"],
                    "raw": cell["raw"], "cells": cell["cells"], "modifier": modifier, "approx": "",
                    "row": cell["row"], "boundary": cell["boundary"]}
        if kind == "G":
            cell = self.fig_cells.get(ref)
            if cell is None:
                return {"error": f"unknown figure cell {ref}"}
            ltype = "constant" if cell["row"]["column"] == "tick" else "table_cell"
            return {"type": ltype, "id": ref, "value": cell["value"], "display": cell["display"], "raw": cell["raw"],
                    "cells": cell["cells"], "modifier": modifier, "approx": "", "row": cell["row"],
                    "boundary": cell["boundary"]}
        if kind == "C":
            const = self.constants.get(ref)
            if const is None:
                return {"error": f"unknown constant {ref}"}
            value = Decimal(const["value"])
            return {"type": "constant", "id": ref, "value": value, "display": const["value"], "raw": const["value"],
                    "cells": {(const["source_file"], int(const["line"] or 0), "constant")}, "modifier": modifier,
                    "approx": "", "row": {"meaning": const["meaning"], "source_file": const["source_file"],
                                          "selector": f"line {const['line']}", "column": "", "status": "design",
                                          "cohort": "", "method": "", "decimals": ""}, "boundary": False}
        if kind == "X":
            return {"type": "exempt", "id": ref, "cells": set()}
        return {"error": f"bad link {link}"}

    def judge(self, token, link: str, context: str, note: str, where: dict) -> dict:
        if "&" in link:
            parts = [self.judge(token, part, context, note, where) for part in link.split("&")]
            claim = dict(parts[0])
            claim.update(link_id="&".join(p["link_id"] for p in parts), ok=all(p["ok"] for p in parts),
                         expected=" & ".join(str(p["expected"]) for p in parts),
                         meaning=" & ".join(p["meaning"] for p in parts),
                         cells=set().union(*(p["cells"] for p in parts)))
            if len({p["link_type"] for p in parts}) != 1 or len({p["status"] for p in parts}) != 1:
                self.issues.append({"kind": "link", "where": token.key,
                                    "detail": f"{link}: parts differ in type or status"})
            return claim
        res = self.resolve(link)
        claim = {"doc": token.doc, "section": token.section, "sentence": where.get("sentence", ""),
                 "key": token.key, "displayed": token.text, "kind": token.kind, "link_type": res.get("type", "error"),
                 "link_id": res.get("id", link), "expected": "", "ok": False, "note": note, "meaning": "",
                 "source_file": "", "selector": "", "column": "", "raw": "", "status": "", "cells": set()}
        if "error" in res:
            claim["note"] = res["error"]
            self.issues.append({"kind": "link", "where": token.key, "detail": res["error"]})
            return claim
        if res["type"] == "exempt":
            claim.update(ok=True, expected=res["id"])
            return claim
        if res["type"] == "constant" and res["id"] in self.constants and token.kind not in {"word", "fraction"}:
            shown = token_value(token)
            value = abs(res["value"]) if res["modifier"] == "abs" else res["value"]
            ok, expected = shown is not None and shown == value, res["display"]
        else:
            ok, expected = compare(token, res["display"], res["value"], res["modifier"], res["approx"])
        row = res["row"]
        claim.update(ok=ok, expected=expected, meaning=row.get("meaning", ""), source_file=row.get("source_file", ""),
                     selector=row.get("selector", ""), column=row.get("column", "") or row.get("transform", ""),
                     raw=res["raw"], status=row.get("status", ""), cells=res["cells"])
        if not ok:
            self.issues.append({"kind": "mismatch", "where": token.key,
                                "detail": f"shown {token.text} expected {expected} ({res['id']})"})
        if res.get("boundary"):
            self.boundaries.append({"key": token.key, "link_id": res["id"], "raw": res["raw"], "displayed": token.text})
        for kind in ("cohort", "method"):
            key = row.get(kind, "")
            if not key:
                continue
            found = all(NC.key_found(kind, k.strip(), context) for k in key.split(" and "))
            if not found:
                resolved = "keys:" in note
                self.key_flags.append({"key": token.key, "link_id": res["id"], "kind": kind, "fact_key": key,
                                       "displayed": token.text, "resolved": resolved,
                                       "resolution": note.split("keys:", 1)[1].strip() if resolved else "",
                                       "context": context[:300]})
        return claim

    # ------------------------------------------------------------ tables
    def tables(self) -> None:
        self.boundaries: list[dict] = []
        self.table_seen: set[str] = set()
        self.af2_rows: list[dict] = []
        for path in sorted(SPECS.glob("*.yaml")):
            spec = yaml.safe_load(path.read_text())
            if spec.get("doc") == "af2":
                self.af2_sheet(spec)
                continue
            if spec.get("kind", "markdown") != "markdown":
                continue
            self.table(spec, path.name)
        self.af2_readme()

    def af2_sheet(self, spec: dict) -> None:
        name = spec["sheet"]
        try:
            expected = AF2.BUILDERS[spec["builder"]](**spec["params"])
        except (NC.SourceError, KeyError, ValueError, IndexError) as exc:
            self.issues.append({"kind": "af2", "where": name, "detail": f"builder failed: {exc}"})
            return
        title, sheet = AF2.sheet_frame(name, spec["header_row"])
        row = {"sheet": name, "origin": spec.get("origin", ""), "builder": spec["builder"], "title": title,
               "rows": len(sheet), "columns": len(sheet.columns), "cells": 0, "numeric_cells": 0, "mismatch": 0,
               "status_words": ", ".join(w for w, rx in NC.STATUS_WORDS.items() if rx.search(title or "")),
               "first_mismatch": ""}
        if list(sheet.columns) != [str(c) for c in expected.columns]:
            row["mismatch"] += 1
            row["first_mismatch"] = f"header {list(sheet.columns)} vs {list(expected.columns)}"
        elif len(sheet) != len(expected):
            row["mismatch"] += 1
            row["first_mismatch"] = f"rows {len(sheet)} vs {len(expected)}"
        else:
            for i in range(len(sheet)):
                for j, col in enumerate(sheet.columns):
                    a, b = sheet.iat[i, j], expected.iat[i, j]
                    row["cells"] += 1
                    if isinstance(a, (int, float)) and not isinstance(a, bool):
                        row["numeric_cells"] += 1
                    if not AF2.same(a, b):
                        row["mismatch"] += 1
                        if not row["first_mismatch"]:
                            row["first_mismatch"] = f"row {i + spec['header_row'] + 1} col {col}: sheet {a!r} expected {b!r}"
        if row["mismatch"]:
            self.issues.append({"kind": "af2 mismatch", "where": name, "detail": f"{row['mismatch']} cells; {row['first_mismatch']}"})
        self.af2_rows.append(row)

    def af2_readme(self) -> None:
        from openpyxl import load_workbook
        book = load_workbook(AF2.BOOK, read_only=True)
        listed = {r[0]: r[1] for r in book["README"].iter_rows(min_row=2, values_only=True) if r[0]}
        for ws in book.worksheets:
            if ws.title == "README":
                continue
            title = next(ws.iter_rows(max_row=1, values_only=True))[0]
            if listed.get(ws.title) != title:
                self.issues.append({"kind": "af2 readme", "where": ws.title,
                                    "detail": f"README text {listed.get(ws.title)!r} differs from row 1 {title!r}"})
        specced = {r["sheet"] for r in self.af2_rows}
        for ws in book.worksheets:
            if ws.title != "README" and ws.title not in specced:
                self.issues.append({"kind": "af2", "where": ws.title, "detail": "sheet has no spec"})

    def find_table(self, doc: N.Document, caption: str, which: int):
        hits = [tid for tid, t in doc.tables.items() if (t["caption"] or "").startswith(caption)]
        if len(hits) <= which:
            return None
        return hits[which]

    def field_value(self, field: dict, variables: dict) -> dict:
        field = NC.fill(field, variables)
        if "fact" in field:
            res = self.resolve("F:" + field["fact"])
        elif "const" in field:
            res = self.resolve("C:" + field["const"])
        else:
            row = {"fact_id": "", "source_file": field["src"], "selector": field.get("sel", ""),
                   "column": field.get("col", ""), "transform": field.get("tf", "none"),
                   "decimals": field.get("dec", 3), "format": field.get("form", "num")}
            try:
                value, raw, index = NC.evaluate(self.src, row["source_file"], row["selector"], row["column"],
                                                row["transform"], None, {**FACTS_MOD.LISTS, **field.get("lists", {})})
                col = row["column"] or row["transform"].split(":", 1)[-1]
                cells = {(row["source_file"], int(i), col) for i in index}
                if row["transform"] == "count":
                    cells |= {(row["source_file"], -1, "count:" + row["selector"])}
                places = int(row["decimals"])
                row.update(status=NC.status_of(row["source_file"], row["selector"]),
                           meaning=field.get("meaning", ""), cohort=field.get("cohort", ""),
                           method=field.get("method", ""))
                res = {"type": "table_cell", "value": value, "raw": raw, "display": NC.show(value, places, row["format"]),
                       "cells": cells, "row": row, "boundary": NC.boundary(value, places, row["format"]),
                       "modifier": "", "approx": ""}
            except NC.SourceError as exc:
                res = {"error": str(exc)}
        if "error" not in res and field.get("abs"):
            res = dict(res, value=abs(res["value"]), display=res["display"].lstrip("-"))
        if "error" not in res and field.get("meaning") and res["row"] is not None:
            res["row"] = dict(res["row"], meaning=field["meaning"])
        return res

    def table(self, spec: dict, name: str) -> None:
        doc = self.docs[spec["doc"]]
        tid = self.find_table(doc, spec["caption"], spec.get("which", 0))
        if tid is None:
            self.issues.append({"kind": "table spec", "where": name, "detail": f"table {spec['caption']} not found"})
            return
        self.table_seen.add(tid)
        table = doc.tables[tid]
        rows = table["rows"]
        header = rows[0]
        if spec.get("header") and [h.strip() for h in header] != spec["header"]:
            self.issues.append({"kind": "table spec", "where": name, "detail": f"header differs: {header}"})
        tokens = defaultdict(list)
        for t in doc.tokens:
            if t.table == tid:
                tokens[(t.row, t.col)].append(t)
        key_cols = spec.get("key_cols", [0])
        maps = spec.get("maps", {})
        columns = spec.get("columns", {})
        context_base = " ".join(header) + " " + (table["caption"] or "") + " " + spec.get("context", "")
        claimed_before = len(self.claims)
        for c, hspec in (spec.get("header_cells") or {}).items():
            c = int(c)
            self.cell(spec, f"header", {"c0": header[0]}, header, 0, c, hspec, tokens, context_base, name)
        seen_rows = set()
        carry = set(spec.get("carry", []))
        last_labels: dict = {}
        for r in range(2, len(rows)):
            cells = rows[r]
            variables = {f"c{i}": c for i, c in enumerate(cells)}
            ok_row = True
            for kc in key_cols:
                label = cells[kc] if kc < len(cells) else ""
                if kc in carry:
                    label = label or last_labels.get(kc, "")
                    last_labels[kc] = label
                cmap = maps.get(kc, maps.get(str(kc), {}))
                hit = None
                for pattern, values in cmap.items():
                    if pattern == label or (pattern.startswith("re:") and re.search(pattern[3:], label)):
                        hit = values or {}
                        break
                if hit is None:
                    ok_row = False
                    self.issues.append({"kind": "table spec", "where": f"{name} row {r}",
                                        "detail": f"label {label!r} of column {kc} not in the row map"})
                    break
                variables.update(hit)
            if not ok_row:
                continue
            row_id = NC.fill(spec.get("row_id", "{c0}"), variables)
            if row_id in seen_rows:
                self.issues.append({"kind": "table spec", "where": f"{name} row {r}", "detail": f"duplicate row id {row_id}"})
            seen_rows.add(row_id)
            context = context_base + " " + " ".join(cells) + " " + " ".join(last_labels.values())
            for c, cell in enumerate(cells):
                cspec = columns.get(c, columns.get(str(c)))
                self.cell(spec, row_id, variables, cells, r, c, cspec, tokens, context, name)
        keys = {cl["key"] for cl in self.claims[claimed_before:]}
        for t in doc.tokens:
            if t.table == tid and not t.exempt and t.key not in keys:
                self.claims.append(self.unlinked(t, "table token not read by the spec"))
                self.issues.append({"kind": "unlinked", "where": t.key, "detail": f"{t.text} not read"})

    def cell(self, spec, row_id, variables, cells, r, c, cspec, tokens, context, name) -> None:
        cell = cells[c]
        toks = [t for t in tokens[(r, c)] if not t.exempt]
        if True:
            if True:
                if cspec is None:
                    for t in toks:
                        self.issues.append({"kind": "unlinked", "where": t.key, "detail": f"table cell {cell!r} has no spec"})
                        self.claims.append(self.unlinked(t))
                    return
                cid = cspec["id"]
                per_row = (cspec.get("rows") or {}).get(row_id, {})
                if per_row is None or per_row == "skip":
                    for t in toks:
                        self.claims.append(self.unlinked(t))
                    return
                template = per_row.get("t", cspec.get("t"))
                fields = {**(cspec.get("f") or {}), **(per_row.get("f") or {})}
                rx, names = NC.template_regex(NC.fill(template, {k: v for k, v in variables.items() if k not in fields}))
                m = rx.match(cell)
                if not m:
                    self.unread_cells.append({"table": spec["id"], "row": row_id, "col": cid, "cell": cell,
                                              "template": template})
                    for t in toks:
                        self.claims.append(self.unlinked(t, f"cell not readable by template {template!r}"))
                    return
                spans = {nm: (m.start(nm), m.end(nm)) for nm in names}
                for nm in names:
                    if nm not in fields:
                        self.issues.append({"kind": "table spec", "where": f"{spec['id']} {row_id} {cid}",
                                            "detail": f"template field {nm} has no definition"})
                        continue
                    key = f"{spec['id']}[row={row_id};col={cid}]" + (f".{nm}" if len(names) > 1 else "")
                    self.cell_facts[key] = self.field_value(fields[nm], variables)
                covered = set()
                for t in toks:
                    hit = [nm for nm, (a, b) in spans.items() if a <= t.start and t.end <= b]
                    if not hit:
                        self.claims.append(self.unlinked(t, f"number outside template fields of {template!r}"))
                        self.issues.append({"kind": "unlinked", "where": t.key, "detail": f"{t.text} outside fields"})
                        continue
                    nm = hit[0]
                    covered.add(nm)
                    key = f"{spec['id']}[row={row_id};col={cid}]" + (f".{nm}" if len(names) > 1 else "")
                    claim = self.judge(t, "T:" + key, context, per_row.get("note", cspec.get("note", "")),
                                       {"sentence": f"{spec['id']} row {row_id}"})
                    self.claims.append(claim)
                for nm in names:
                    if nm not in covered:
                        self.issues.append({"kind": "table spec", "where": f"{spec['id']} {row_id} {cid}",
                                            "detail": f"field {nm} matched no number token (exempt?)"})

    def unlinked(self, token, note: str = "") -> dict:
        return {"doc": token.doc, "section": token.section, "sentence": "", "key": token.key, "displayed": token.text,
                "kind": token.kind, "link_type": "unlinked", "link_id": "", "expected": "", "ok": False, "note": note,
                "meaning": "", "source_file": "", "selector": "", "column": "", "raw": "", "status": "", "cells": set()}

    # ------------------------------------------------------------ figures (A17 section 4.5, A18 section 2.2)
    def figures(self) -> None:
        import fig_expected_A17 as FE
        import figcheck_A17 as FC
        for name in FE.BUILDERS:
            self.figure_source_rows.append(FE.compare(name))
        table = self.label_table if self.label_table is not None else \
            pd.read_csv(FC.LABELS, sep="\t", dtype={"text": str}, keep_default_na=False)
        for fig in FC.FIGURES:
            e = FC.expected(fig)
            self.fig_cells.update(e.cells)
            words = FC.printed_words(fig)
            lines = [" ".join(w[4] for w in ws) for ws in words]
            doc = f"fig:{fig}"
            number = fig[len("Figure"):]
            source_doc = self.docs["af1" if number.startswith("S") else "main"].text
            legend = re.search(rf"^\*\*Figure {number}\. .*?(?=^\*\*|\Z)", source_doc, re.M | re.S)
            reading = " ".join(lines) + " " + (legend.group(0) if legend else "")
            live = {}
            word_of: dict[str, tuple[int, int]] = {}
            for k, line in enumerate(lines):
                toks = N.tokens_in(doc, fig, 0, k, line, False)
                self.fig_tokens.update({t.key: t for t in toks})
                spans, pos = [], 0
                for w in words[k]:
                    spans.append((pos, pos + len(w[4])))
                    pos += len(w[4]) + 1
                for t in toks:
                    word_of[t.key] = (k, next(i for i, (a, b) in enumerate(spans) if a <= t.start < b))
                    if t.exempt:
                        self.claims.append(self.exempt_claim(t))
                live[k] = [t for t in toks if not t.exempt]
            row = {"figure": fig, "printed_lines": len(lines), "numbers": sum(len(v) for v in live.values()),
                   "expected_items": len(e.items), "value_labels": 0, "label_numbers": 0, "position_mismatch": 0,
                   "labels_without_item": 0, "items_without_label": 0, "keyed_lines": 0, "exact_lines": 0,
                   "tokens": 0, "ticks": 0, "linked": 0, "mismatch": 0, "unlinked": 0, "missing_items": 0,
                   "ambiguous_tick_links": 0}
            done: set[int] = set()
            used: set[str] = set()
            self.figure_labels(fig, e, FC.labels(fig, table), words, live, word_of, reading, row, used)
            partial = {k for k in live if any(t.key in used for t in live[k])}
            for k in live:
                live[k] = [t for t in live[k] if t.key not in used]
            e.items = [it for it in e.items if not it.key]

            def link_line(k: int, item, how: str) -> None:
                done.add(k)
                if len(live[k]) != len(item.links):
                    self.issues.append({"kind": "figure", "where": f"{fig} line {k}",
                                        "detail": f"{len(live[k])} numbers for {len(item.links)} links: {lines[k]}"})
                    for t in live[k]:
                        self.claims.append(self.unlinked(t, "figure line with a different number count"))
                    return
                for t, link in zip(live[k], item.links):
                    claim = self.judge(t, link, lines[k] + " " + reading, how, {"sentence": f"line {k}"})
                    self.claims.append(claim)
                    row["mismatch"] += not claim["ok"]
                    row["linked"] += 1
                    row["ticks"] += ":tick:" in link

            def expected_text(item) -> str:
                return item.template.format(*[self.resolve(link).get("display", "?") for link in item.links])

            count = Counter(it.template for it in e.items if it.whole)
            pending = []
            for item in e.items:
                if not (item.whole and count[item.template] == 1):
                    pending.append(item)
                    continue
                rx = FC.template_rx(item.template)
                hits = [k for k in live if k not in done and k not in partial and live[k] and rx.match(NC.norm_token(lines[k]))]
                if len(hits) != 1:
                    self.issues.append({"kind": "figure missing", "where": fig,
                                        "detail": f"{len(hits)} printed lines fit {item.template!r}"})
                    row["missing_items"] += 1
                    continue
                row["keyed_lines"] += 1
                link_line(hits[0], item, "figure line, keyed by its fixed text")
            standalone = defaultdict(list)
            for item in pending:
                text = NC.norm_token(expected_text(item))
                if not item.whole:
                    standalone[text].append(item)
                    continue
                hits = [k for k in live if k not in done and k not in partial and live[k] and NC.norm_token(lines[k]) == text]
                if not hits:
                    self.issues.append({"kind": "figure missing", "where": fig, "detail": f"no printed line {text!r}"})
                    row["missing_items"] += 1
                    continue
                row["exact_lines"] += 1
                link_line(hits[0], item, "figure line equal to the expected text")
            pools = {text: len({it.links[0] for it in items}) for text, items in standalone.items()}
            for k in live:
                if k in done:
                    continue
                for t in live[k]:
                    text = NC.norm_token(t.text)
                    if not standalone.get(text):
                        self.claims.append(self.unlinked(t, f"figure number with no expected item: {lines[k]}"))
                        row["unlinked"] += 1
                        continue
                    item = standalone[text].pop(0)
                    claim = self.judge(t, item.links[0], lines[k] + " " + reading, "figure tick equal to the expected text",
                                       {"sentence": f"line {k}"})
                    self.claims.append(claim)
                    row["tokens"] += 1
                    row["linked"] += 1
                    row["ticks"] += ":tick:" in item.links[0]
                    row["ambiguous_tick_links"] += pools[text] > 1
            for text, items in standalone.items():
                for item in items:
                    self.issues.append({"kind": "figure missing", "where": fig,
                                        "detail": f"no printed number {text!r} for {item.links[0]}"})
                    row["missing_items"] += 1
            self.fig_rows.append(row)

    def figure_labels(self, fig: str, e, recs: pd.DataFrame, words: list, live: dict, word_of: dict,
                      reading: str, row: dict, used: set) -> None:
        """Value labels: the printed words inside each recorded extent must spell the recorded text;
        their numbers are linked, in printed order, to the expected item with the same element key."""
        def norm(text: str) -> str:
            return " ".join(text.replace("\u2212", "-").split())

        tol = LABEL_TOL_MM
        boxes = recs[["x0_mm", "y0_mm", "x1_mm", "y1_mm"]].astype(float).to_numpy()
        owner: dict[tuple[int, int], int] = {}
        for k, ws in enumerate(words):
            for i, w in enumerate(ws):
                cx, cy = (w[0] + w[2]) / 2, (w[1] + w[3]) / 2
                inside = [r for r, (x0, y0, x1, y1) in enumerate(boxes)
                          if x0 - tol <= cx <= x1 + tol and y0 - tol <= cy <= y1 + tol]
                if inside:
                    owner[(k, i)] = min(inside, key=lambda r: (cx - (boxes[r][0] + boxes[r][2]) / 2) ** 2
                                        + (cy - (boxes[r][1] + boxes[r][3]) / 2) ** 2)
        keyed: dict[str, object] = {}
        for item in e.items:
            if item.key:
                if item.key in keyed:
                    self.issues.append({"kind": "figure", "where": item.key, "detail": "two expected items with one key"})
                keyed[item.key] = item
        seen_keys: set[str] = set()
        for r, rec in recs.iterrows():
            mine = sorted(kw for kw, o in owner.items() if o == r)
            printed = " ".join(words[k][i][4] for k, i in mine)
            toks = [t for k in sorted({k for k, _ in mine}) for t in live[k] if word_of[t.key] in mine]
            used.update(t.key for t in toks)
            row["value_labels"] += 1
            row["label_numbers"] += len(toks)
            out = {"figure": fig, "key": rec["key"], "recorded_text": rec["text"], "printed_text": printed,
                   "x0_mm": rec["x0_mm"], "y0_mm": rec["y0_mm"], "x1_mm": rec["x1_mm"], "y1_mm": rec["y1_mm"],
                   "position_ok": norm(printed) == norm(str(rec["text"])), "numbers": len(toks), "links": "",
                   "linked_ok": False}
            self.label_rows.append(out)
            if rec["key"] in seen_keys:
                self.issues.append({"kind": "figure", "where": rec["key"], "detail": "two label records with one key"})
            seen_keys.add(rec["key"])
            if not out["position_ok"]:
                self.issues.append({"kind": "figure position", "where": rec["key"],
                                    "detail": f"recorded {rec['text']!r}, printed at that place {printed!r}"})
                row["position_mismatch"] += 1
                for t in toks:
                    self.claims.append(self.unlinked(t, f"figure label {rec['key']}: printed text differs from the record"))
                    row["unlinked"] += 1
                continue
            item = keyed.pop(rec["key"], None)
            if item is None:
                self.issues.append({"kind": "figure", "where": rec["key"], "detail": "label record with no expected item"})
                row["labels_without_item"] += 1
                for t in toks:
                    self.claims.append(self.unlinked(t, f"figure label {rec['key']} with no expected item"))
                    row["unlinked"] += 1
                continue
            out["links"] = " ".join(item.links)
            if len(toks) != len(item.links):
                self.issues.append({"kind": "figure", "where": rec["key"],
                                    "detail": f"{len(toks)} numbers for {len(item.links)} links: {printed!r}"})
                for t in toks:
                    self.claims.append(self.unlinked(t, "figure label with a different number count"))
                    row["unlinked"] += 1
                continue
            oks = []
            for t, link in zip(toks, item.links):
                claim = self.judge(t, link, printed + " " + reading, f"figure label {rec['key']}, matched by key and position",
                                   {"sentence": rec["key"]})
                self.claims.append(claim)
                oks.append(claim["ok"])
                row["mismatch"] += not claim["ok"]
                row["linked"] += 1
            out["linked_ok"] = all(oks)
        for key in keyed:
            self.issues.append({"kind": "figure missing", "where": key, "detail": "expected label with no label record"})
            row["items_without_label"] += 1
            row["missing_items"] += 1

    # ------------------------------------------------------------ prose
    def prose(self) -> None:
        links = pd.read_csv(LINKS, sep="\t", dtype=str, keep_default_na=False) if LINKS.exists() else \
            pd.DataFrame(columns=["doc", "anchor", "starts", "ids", "note"])
        by_anchor = {}
        for r in links.itertuples(index=False):
            if (r.doc, r.anchor) in by_anchor:
                self.issues.append({"kind": "links", "where": f"{r.doc} {r.anchor}", "detail": "duplicate anchor"})
            by_anchor[(r.doc, r.anchor)] = r
        used = set()
        for name, doc in self.docs.items():
            groups = defaultdict(list)
            for t in doc.tokens:
                if t.table:
                    if t.table not in self.table_seen and not t.exempt:
                        self.claims.append(self.unlinked(t, "table without spec"))
                    elif t.exempt:
                        self.claims.append(self.exempt_claim(t))
                    continue
                groups[(t.para, t.sent)].append(t)
            order = sorted(doc.sentences)
            for (para, sent), toks in sorted(groups.items()):
                anchor = f"p{para}|s{sent}"
                sentence = doc.sentences[(para, sent)]
                live = [t for t in toks if not t.exempt]
                for t in toks:
                    if t.exempt:
                        self.claims.append(self.exempt_claim(t))
                if not live:
                    continue
                row = by_anchor.get((name, anchor))
                if row is None:
                    for t in live:
                        self.claims.append(self.unlinked(t))
                    continue
                used.add((name, anchor))
                if row.starts and not re.sub(r"\d", "0", sentence).startswith(re.sub(r"\d", "0", row.starts)):
                    self.issues.append({"kind": "links", "where": f"{name}|{anchor}",
                                        "detail": f"sentence no longer starts with {row.starts!r}"})
                    for t in live:
                        self.claims.append(self.unlinked(t, "anchor moved"))
                    continue
                ids = row.ids.split()
                if len(ids) != len(live):
                    self.issues.append({"kind": "links", "where": f"{name}|{anchor}",
                                        "detail": f"{len(ids)} ids for {len(live)} tokens: {[t.text for t in live]}"})
                    for t in live:
                        self.claims.append(self.unlinked(t, "id count differs"))
                    continue
                i = order.index((para, sent))
                prev = doc.sentences[order[i - 1]] if i > 0 else ""
                context = sentence + " " + prev + " " + doc.sections.get(para, "")
                for t, link in zip(live, ids):
                    claim = self.judge(t, link, context, row.note, {"sentence": anchor})
                    self.claims.append(claim)
                self.status_check(name, anchor, sentence, [c for c in self.claims[-len(live):]], row.note)
        for key in by_anchor:
            if key not in used:
                self.issues.append({"kind": "links", "where": f"{key[0]}|{key[1]}", "detail": "link row matches no sentence"})

    def exempt_claim(self, t) -> dict:
        return {"doc": t.doc, "section": t.section, "sentence": "", "key": t.key, "displayed": t.text, "kind": t.kind,
                "link_type": "exempt", "link_id": t.exempt, "expected": t.exempt, "ok": True, "note": t.why,
                "meaning": "", "source_file": "", "selector": "", "column": "", "raw": "", "status": "", "cells": set()}

    def status_check(self, doc: str, anchor: str, sentence: str, claims: list[dict], note: str = "") -> None:
        words = {s for s, rx in NC.STATUS_WORDS.items() if rx.search(sentence)}
        if not words:
            return
        for c in claims:
            if c["link_type"] in {"fact", "derived", "table_cell"} and c["status"] and \
                    c["status"] not in {"design", "data description"} and c["status"] not in words:
                self.status_flags.append({"doc": doc, "anchor": anchor, "key": c["key"], "displayed": c["displayed"],
                                          "link_id": c["link_id"], "fact_status": c["status"],
                                          "sentence_words": "; ".join(sorted(words)),
                                          "resolved": "status:" in note,
                                          "resolution": note[note.index("status:"):] if "status:" in note else "",
                                          "sentence": sentence[:300]})

    def read_files(self) -> dict[str, str]:
        """Every file this check reads, with its role (A17 section 4.12)."""
        import fig_expected_A17 as FE
        import figcheck_A17 as FC
        roles: dict[str, str] = {}
        for path in sorted(self.src.used):
            roles[path] = "fact or table source"
        for const in self.constants.values():
            roles.setdefault(const["source_file"], "constant source")
        for path in sorted(FE.READ):
            roles.setdefault(path, "figure expectation source")
        for name in FE.BUILDERS:
            roles.setdefault(f"manuscript/bmc/figure_source/{name}.tsv", "figure_source file compared with results")
        for fig in FC.FIGURES:
            roles.setdefault(f"manuscript/bmc/figures/{fig}.pdf", "figure PDF read for printed numbers")
        roles.setdefault(str(FC.LABELS.relative_to(NC.ROOT)), "figure value labels with their extents")
        for pattern in REQUIRED:
            for path in sorted(NC.ROOT.glob(pattern)):
                if path.is_file():
                    roles.setdefault(str(path.relative_to(NC.ROOT)), "listed by A17 section 4.12")
        return roles

    def run(self) -> dict:
        self.tables()
        self.figures()
        self.prose()
        self.manifest = check_manifest(self.read_files())
        bad = [c for c in self.claims if c["link_type"] not in {"exempt", "unlinked"} and not c["ok"]]
        unlinked = [c for c in self.claims if c["link_type"] == "unlinked"]
        return {"claims": len(self.claims), "unlinked": len(unlinked), "mismatch": len(bad),
                **{k: v for k, v in self.manifest.items() if k != "manifest_rows"},
                "selector_errors": len(self.book.errors), "issues": len(self.issues),
                "key_flags": len(self.key_flags), "key_flags_unresolved": sum(not f["resolved"] for f in self.key_flags),
                "status_flags": len(self.status_flags),
                "status_flags_unresolved": sum(not f["resolved"] for f in self.status_flags), "unread_cells": len(self.unread_cells),
                "facts": len(self.book.rows), "facts_dropped_empty": len(self.book.empty), "table_cells": len(self.cell_facts),
                "af2_sheets": len(self.af2_rows), "af2_cells": sum(r["cells"] for r in self.af2_rows),
                "af2_mismatch_cells": sum(r["mismatch"] for r in self.af2_rows),
                "figure_source_files": len(self.figure_source_rows),
                "figure_source_mismatch": sum(r["mismatch"] for r in self.figure_source_rows),
                "figure_numbers": sum(r["numbers"] for r in self.fig_rows),
                "figure_unlinked": sum(r["unlinked"] for r in self.fig_rows),
                "figure_missing_items": sum(r["missing_items"] for r in self.fig_rows),
                "figure_mismatch": sum(r["mismatch"] for r in self.fig_rows),
                "figure_value_labels": len(self.label_rows),
                "figure_label_position_mismatch": sum(not r["position_ok"] for r in self.label_rows),
                "figure_labels_without_item": sum(r["labels_without_item"] for r in self.fig_rows),
                "figure_items_without_label": sum(r["items_without_label"] for r in self.fig_rows),
                "figure_label_numbers": sum(r["label_numbers"] for r in self.fig_rows),
                "figure_ambiguous_tick_links": sum(r["ambiguous_tick_links"] for r in self.fig_rows)}

    def write(self, summary: dict) -> None:
        out = NC.CHECKS
        frame = pd.DataFrame([{k: v for k, v in c.items() if k != "cells"} for c in self.claims])
        frame.to_csv(out / f"claims_map_{TAG}.tsv", sep="\t", index=False)
        pd.DataFrame(self.issues, columns=["kind", "where", "detail"]).to_csv(out / f"check_issues_{TAG}.tsv", sep="\t", index=False)
        pd.DataFrame(self.key_flags).to_csv(out / f"key_flags_{TAG}.tsv", sep="\t", index=False)
        pd.DataFrame(self.status_flags).to_csv(out / f"status_flags_{TAG}.tsv", sep="\t", index=False)
        pd.DataFrame(self.boundaries).drop_duplicates().to_csv(out / f"rounding_boundaries_{TAG}.tsv", sep="\t", index=False)
        pd.DataFrame(self.unread_cells).to_csv(out / f"unread_cells_{TAG}.tsv", sep="\t", index=False)
        pd.DataFrame(self.af2_rows).to_csv(out / f"af2_check_{TAG}.tsv", sep="\t", index=False)
        pd.DataFrame(self.fig_rows).to_csv(out / f"figure_check_{TAG}.tsv", sep="\t", index=False)
        pd.DataFrame(self.label_rows).to_csv(out / f"figure_label_check_{TAG}.tsv", sep="\t", index=False)
        pd.DataFrame(self.figure_source_rows).to_csv(out / f"figure_source_check_{TAG}.tsv", sep="\t", index=False)
        pd.DataFrame(self.manifest["manifest_rows"]).to_csv(out / f"manifest_check_{TAG}.tsv", sep="\t", index=False)
        facts = pd.DataFrame([r for r in self.fact_rows if r["fact_id"] in self.book.rows])
        facts.to_csv(out / "facts.tsv", sep="\t", index=False)
        (out / f"check_summary_{TAG}.json").write_text(json.dumps(summary, indent=1))


def remaining_table(claims: list[dict]) -> pd.DataFrame:
    rows = Counter((c["doc"], c["section"].split(" > ")[0]) for c in claims if c["link_type"] == "unlinked")
    return pd.DataFrame([{"doc": d, "section": s, "unlinked": n} for (d, s), n in sorted(rows.items())])


# ---------------------------------------------------------------- mutation tests (A17 section 4.10)

def locate(doc: N.Document, token) -> int | None:
    """Absolute offset of a token in the document text, or None."""
    text = doc.text
    if token.table:
        table = doc.tables[token.table]
        line_no = table["first_line"] + token.row
        lines = text.splitlines(keepends=True)
        start = sum(len(l) for l in lines[:line_no - 1])
        line = lines[line_no - 1]
        cells = N.split_row(line)
        body = line.strip()
        offset = line.index(body) + (1 if body.startswith("|") else 0)
        parts = body[1:].split("|") if body.startswith("|") else body.split("|")
        pos = offset
        for i, part in enumerate(parts):
            if i == token.col:
                lead = len(part) - len(part.lstrip())
                if part.strip() != cells[i]:
                    return None
                return start + pos + lead + token.start
            pos += len(part) + 1
        return None
    sentence = token.sentence
    hits = [m.start() for m in re.finditer(re.escape(sentence), text)]
    if len(hits) != 1:
        return None
    return hits[0] + token.start


def mutate() -> dict:
    base = Checker()
    base.run()
    rng = random.Random(20261008)
    pool = [c for c in base.claims if c["link_type"] in {"fact", "derived", "table_cell", "constant"} and c["ok"]
            and re.search(r"\d", c["displayed"]) and c["kind"] not in {"word", "fraction"}]
    rng.shuffle(pool)
    tokens = {t.key: (name, t) for name, d in base.docs.items() for t in d.tokens}
    pool = [c for c in pool if c["key"] in tokens]
    texts = {name: d.text for name, d in base.docs.items()}
    chosen = []
    for c in pool:
        name, t = tokens[c["key"]]
        pos = locate(base.docs[name], t)
        if pos is None:
            continue
        shown = t.text
        last = max(i for i, ch in enumerate(shown) if ch.isdigit())
        absolute = pos + last
        if texts[name][absolute] != shown[last]:
            continue
        new_digit = str((int(shown[last]) + 1) % 10)
        texts[name] = texts[name][:absolute] + new_digit + texts[name][absolute + 1:]
        chosen.append({"where": "document", "key": c["key"], "before": shown,
                       "after": shown[:last] + new_digit + shown[last + 1:], "link_id": c["link_id"]})
        if len(chosen) == 15:
            break
    import figcheck_A17 as FC
    table = pd.read_csv(FC.LABELS, sep="\t", dtype={"text": str}, keep_default_na=False)
    linked = {r["key"] for r in base.label_rows if r["linked_ok"] and r["numbers"]}
    rows = [i for i in table.index if table.at[i, "key"] in linked]
    for i in rng.sample(rows, 5):
        shown = table.at[i, "text"]
        last = max(j for j, ch in enumerate(shown) if ch.isdigit())
        changed = shown[:last] + str((int(shown[last]) + 1) % 10) + shown[last + 1:]
        table.at[i, "text"] = changed
        chosen.append({"where": "figure label table", "key": table.at[i, "key"], "before": shown, "after": changed,
                       "link_id": ""})
    test = Checker(texts=texts, labels=table)
    test.run()
    failed = {c["key"] for c in test.claims if c["link_type"] not in {"exempt", "unlinked"} and not c["ok"]}
    base_failed = {c["key"] for c in base.claims if c["link_type"] not in {"exempt", "unlinked"} and not c["ok"]}
    moved = {r["key"] for r in test.label_rows if not r["position_ok"]}
    base_moved = {r["key"] for r in base.label_rows if not r["position_ok"]}
    caught = [m["key"] in (failed if m["where"] == "document" else moved) for m in chosen]
    extra = sorted((failed - base_failed - {m["key"] for m in chosen}) | (moved - base_moved - {m["key"] for m in chosen}))
    for m, hit in zip(chosen, caught):
        m["caught"] = hit
        row = m["key"].rsplit("|c", 1)[0] + "|c" if m["where"] == "document" and "|r" in m["key"] else None
        m["same_row_failures"] = sorted(k for k in extra if row and k.startswith(row))
    return {"chosen": chosen, "caught": sum(caught),
            "caught_or_row_failed": sum(m["caught"] or bool(m["same_row_failures"]) for m in chosen),
            "extra_failures": extra,
            "extra_failures_outside_mutated_rows": sorted(set(extra) - {k for m in chosen for k in m["same_row_failures"]})}


def override() -> dict:
    base = Checker()
    base.run()
    rng = random.Random(20261008)
    used = Counter(part for c in base.claims if c["link_type"] == "fact" for part in c["link_id"].split("&"))
    candidates = sorted(fid for fid in used if base.book.rows[fid]["transform"] in {"none", "x100"})
    picks = rng.sample(candidates, 5)
    overrides = {}
    detail = []
    for fid in picks:
        res = base.book.get(fid)
        (path, index, column), = [c for c in res["cells"]]
        raw = Decimal(NC.norm_token(res["raw"]))
        new = raw + Decimal(7) if raw.adjusted() >= 0 or raw == 0 else raw * Decimal("1.5") + Decimal("0.0137")
        overrides[(path, index, column)] = str(new)
        expected = sorted(c["key"] for c in base.claims if (path, index, column) in c["cells"])
        detail.append({"fact_id": fid, "source_file": path, "row": index, "column": column, "raw": res["raw"],
                       "override": str(new), "positions": expected})
    test = Checker(overrides=overrides)
    test.run()
    failed = {c["key"] for c in test.claims if c["link_type"] not in {"exempt", "unlinked"} and not c["ok"]}
    before = {c["key"]: c["expected"] for c in base.claims}
    after = {c["key"]: c["expected"] for c in test.claims}
    total = 0
    reading: set[str] = set()
    for d in detail:
        reading |= set(d["positions"])
        d["changed"] = [k for k in d["positions"] if before.get(k) != after.get(k)]
        d["unchanged"] = [k for k in d["positions"] if k not in d["changed"]]
        d["caught"] = sum(k in failed for k in d["changed"])
        d["unchanged_failed"] = sum(k in failed for k in d["unchanged"])
        total += d["caught"] == len(d["changed"]) and len(d["changed"]) > 0 and d["unchanged_failed"] == 0
    base_failed = {c["key"] for c in base.claims if c["link_type"] not in {"exempt", "unlinked"} and not c["ok"]}
    return {"facts": detail, "facts_fully_caught": total, "extra_failures": sorted(failed - base_failed - reading)}


def main() -> None:
    mode = sys.argv[1] if len(sys.argv) > 1 else "--check"
    if mode == "--write-manifest":
        checker = Checker()
        checker.run()
        rows = [{"path": p, "sha256": sha256(NC.ROOT / p), "bytes": (NC.ROOT / p).stat().st_size, "role": role,
                 "status": NC.status_of(p)} for p, role in sorted(checker.read_files().items())]
        pd.DataFrame(rows).to_csv(MANIFEST, sep="\t", index=False)
        print(len(rows), "files written to", MANIFEST)
        return
    if mode == "--mutate":
        result = mutate()
        (NC.CHECKS / f"mutation_digits_{TAG}.json").write_text(json.dumps(result, indent=1, default=str))
        print("digits caught", result["caught"], "of", len(result["chosen"]), "; caught or failing in their table row",
              result["caught_or_row_failed"], "; extra", len(result["extra_failures"]),
              "; extra outside mutated rows", len(result["extra_failures_outside_mutated_rows"]))
        return
    if mode == "--override":
        result = override()
        (NC.CHECKS / f"mutation_facts_{TAG}.json").write_text(json.dumps(result, indent=1, default=str))
        print("facts fully caught", result["facts_fully_caught"], "of", len(result["facts"]))
        for d in result["facts"]:
            print(d["fact_id"], len(d["positions"]), "positions read the cell,", len(d["changed"]),
                  "change value,", d["caught"], "caught,", d["unchanged_failed"], "unchanged but failed")
        print("extra failures", len(result["extra_failures"]))
        return
    checker = Checker()
    summary = checker.run()
    checker.write(summary)
    table = remaining_table(checker.claims)
    table.to_csv(NC.CHECKS / f"unlinked_by_section_{TAG}.tsv", sep="\t", index=False)
    print(json.dumps(summary, indent=1))
    print(table.groupby("doc")["unlinked"].sum().to_string() if len(table) else "no unlinked")
    kinds = Counter(i["kind"] for i in checker.issues)
    print(dict(kinds))


if __name__ == "__main__":
    main()
