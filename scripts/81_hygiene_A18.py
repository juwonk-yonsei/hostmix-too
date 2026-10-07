#!/usr/bin/env python3
"""A18 section 5: hygiene checklist of the submission files (A19: outputs _A19; Sections 3 and 9 exempt).

Items: abstract length; figure title and legend length; reference order and uncited references;
first-citation order of figures, tables and additional files; abbreviation list against the abbreviations
defined in the text; the standardization terms; forbidden words, internal names, paths and hashes
(Additional file 1 Sections 2, 3, 9, 10 and 11, whose source notes may name files, and Table S7 are exempt);
the figures folder; file sizes and font embedding; placeholders; release dry run.

Outputs (manuscript/checks/): hygiene_A19.tsv (one row per item), hygiene_detail_A19.tsv (every hit),
placeholders_A19.tsv (file, line, text), hygiene_checklist_A19.md, release_dry_run_A19.txt.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BMC = ROOT / "manuscript" / "bmc"
CHECKS = ROOT / "manuscript" / "checks"
MAIN = BMC / "A_main.md"
AF1 = BMC / "Additional_file_1.md"
COVER = BMC / "cover_letter.md"
PLACEHOLDER_FILES = [MAIN, AF1, COVER, ROOT / "release" / "README.md", ROOT / "release" / "models" / "MODEL_CARD.md",
                     ROOT / "release" / "LICENSE", ROOT / "README.md", ROOT / "LICENSE", ROOT / "CITATION.cff",
                     ROOT / ".zenodo.json"]
PLACEHOLDER = re.compile(r"\[[^\[\]]*(?:AUTHOR TO|TO BE ASSIGNED)[^\[\]]*\]")
SCALE = re.compile(r"\b(scaled|unscaled|scaling)\b", re.I)
WORDS = {"stored": re.compile(r"\bstored\b", re.I), "significant": re.compile(r"\bsignifican(?:t|tly|ce)\b", re.I),
         "host-pull": re.compile(r"host-pull", re.I)}
PATH = re.compile(r"\b(?:results|scripts|config|data|release|manuscript|logs)/|\b\w+\.(?:py|tsv|parquet|joblib|npz|yaml|json)\b"
                  r"|\bstage\d+\b")
HASH = re.compile(r"\b[0-9a-f]{40}\b")
EXEMPT_AF1_SECTIONS = {"2", "3", "9", "10", "11", "14"}
TAG = "A19"
FIGURES_EXPECTED = 24
LIMITS = {"abstract": 350, "title": 15, "legend": 300}
TEXT_SUFFIXES = {".tsv", ".csv", ".txt", ".md", ".json", ".yaml"}
POG_MIN = 10
DUMMY_URL, DUMMY_DOI = "https://example.invalid/dry-run", "10.0000/dry-run"


def section(text: str, name: str) -> str:
    m = re.search(rf"^## {re.escape(name)}\n(.*?)(?=^## |\Z)", text, re.M | re.S)
    return m.group(1) if m else ""


def body(text: str) -> str:
    return text[:text.index("## References")]


def plain_words(text: str) -> int:
    return len(re.sub(r"\*\*|\*", "", text).split())


def figure_legends() -> list[dict]:
    rows = []
    for doc in (MAIN, AF1):
        for m in re.finditer(r"^\*\*Figure (S?\d+)\. (.+?)\*\*\n\n(.+?)(?=\n\n|\Z)", doc.read_text(), re.M | re.S):
            rows.append({"file": doc.name, "figure": f"Figure {m.group(1)}", "title_words": plain_words(m.group(2)),
                         "legend_words": plain_words(m.group(3))})
    return rows


def expand(group: str) -> list[int]:
    out = []
    for part in re.split(r"\s*,\s*", group):
        a, _, b = part.partition("–") if "–" in part else part.partition("-")
        out += list(range(int(a), int(b) + 1)) if b else [int(a)]
    return out


def first_order(numbers: list[int]) -> tuple[list[int], int]:
    order = list(dict.fromkeys(numbers))
    return order, sum(n != i + 1 for i, n in enumerate(order))


def references(text: str) -> dict:
    listed = [int(n) for n in re.findall(r"^(\d+)\. ", section(text, "References"), re.M)]
    cites = []
    head = text[:text.index("## References")] + text[text.index("## Figure titles and legends"):]
    for m in re.finditer(r"\[(\d+(?:\s*[,–-]\s*\d+)*)\]", head):
        cites += expand(m.group(1))
    order, deviations = first_order(cites)
    return {"listed": len(listed), "cited": len(set(cites)), "order_deviations": deviations,
            "uncited": sorted(set(listed) - set(cites)), "missing": sorted(set(cites) - set(listed)),
            "listed_sequential": listed == list(range(1, len(listed) + 1))}


def citation_orders(text: str) -> dict:
    main = body(text)
    figs = []
    for m in re.finditer(r"\bFigs?\.\s*((?:\d+[a-z]?(?:[–-][a-z])?(?:,\s*[a-z](?=\W))*(?:\s*(?:,|and|–)\s*)?)+)", main):
        figs += [int(n) for n in re.findall(r"\b(\d+)", m.group(1))]
    tables = [int(n) for n in re.findall(r"\bTables? (\d+)\b", main)]
    files = [int(n) for n in re.findall(r"\bAdditional files? (\d+)\b", main)]
    s_tables = [int(n) for n in re.findall(r"\bTables? S(\d+)", main)]
    s_figs = [int(n) for n in re.findall(r"\b(?:Figs?\.|Figures?) S(\d+)", main)]
    out = {}
    for name, seq in (("figures", figs), ("tables", tables), ("additional_files", files),
                      ("supplementary_tables_in_main", s_tables), ("supplementary_figures_in_main", s_figs)):
        order, dev = first_order(seq)
        out[name] = {"first_citation_order": order, "deviations": dev}
    return out


def abbreviations(text: str) -> dict:
    listed = {}
    for item in section(text, "List of abbreviations").strip().split("; "):
        abbr, _, meaning = item.partition(": ")
        listed[abbr.strip()] = meaning.strip().rstrip(".")
    rest = text.replace(section(text, "List of abbreviations"), "").replace(section(text, "References"), "")
    used = {a: len(re.findall(rf"(?<![\w/-]){re.escape(a)}(?![\w/-])", rest)) for a in listed}
    defined = sorted({m.group(1) for m in re.finditer(r"[a-z][\w\-]* \(([A-Z][A-Za-z0-9/\-]*[A-Z][A-Za-z0-9/\-]*)\)", rest)
                      if sum(c.isupper() for c in m.group(1)) >= 2})
    names = [a for a in defined if re.fullmatch(r"GSE\d+", a) or a == "HostMix-TOO"]
    return {"listed": listed, "used_counts": used, "listed_not_used": [a for a, n in used.items() if n == 0],
            "defined_in_text": defined, "accessions_or_method_name": names,
            "defined_not_listed": [a for a in defined if a not in listed and a not in names]}


def af1_sections(text: str) -> list[tuple[str, str]]:
    out, current = [], "0"
    for line in text.splitlines():
        m = re.match(r"^## (\d+)\.", line)
        if m:
            current = m.group(1)
        out.append((current, line))
    return out


def internal_names(text: str) -> list[str]:
    s7 = text[text.index("### Table S7"):]
    names = []
    for line in s7.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) >= 2 and line.startswith("|") and not set(cells[1]) <= set("-:") and cells[1] != "Internal name":
            names += [n for n in cells[1].split() if re.search(r"[A-Z]-|\d|-[A-Za-z]", n) and not re.fullmatch(r"A?H\d?|AS|AH", n)
                      and n not in {"CUP-AI-Dx"}]
    return sorted(set(names), key=len, reverse=True)


def in_git(sha: str) -> bool:
    return subprocess.run(["git", "cat-file", "-e", f"{sha}^{{commit}}"], cwd=ROOT, capture_output=True).returncode == 0


def forbidden(names: list[str]) -> list[dict]:
    hits = []
    code = re.compile(r"(?<![\w-])(?:" + "|".join(re.escape(n) for n in names) + r")(?![\w-])")
    for path in (MAIN, AF1, COVER):
        lines = af1_sections(path.read_text()) if path == AF1 else [("", l) for l in path.read_text().splitlines()]
        for no, (sec, line) in enumerate(lines, start=1):
            exempt = path == AF1 and sec in EXEMPT_AF1_SECTIONS
            checks = list(WORDS.items()) + [("internal name", code), ("path", PATH), ("hash", HASH)]
            for kind, rx in checks:
                for m in rx.finditer(line):
                    hit = {"file": path.name, "line": no, "kind": kind, "text": m.group(0), "exempt": exempt,
                           "context": line[max(0, m.start() - 50):m.end() + 50]}
                    if kind == "hash":
                        hit["project_commit"] = in_git(m.group(0))
                        hit["exempt"] = exempt or not hit["project_commit"]
                    hits.append(hit)
    return hits


def placeholders() -> list[dict]:
    rows = []
    for path in PLACEHOLDER_FILES:
        heading = ""
        for no, line in enumerate(path.read_text().splitlines(), start=1):
            if line.startswith("#"):
                heading = line.strip("# ")
            for m in PLACEHOLDER.finditer(line):
                rows.append({"file": str(path.relative_to(ROOT)), "line": no, "heading": heading, "text": m.group(0)})
    return rows


def pdf_fonts(path: Path) -> tuple[int, int]:
    out = subprocess.run(["pdffonts", str(path)], capture_output=True, text=True, check=True).stdout.splitlines()[2:]
    emb = [line.split()[-5] for line in out if line.strip()]
    return len(emb), sum(e == "yes" for e in emb)


def files_and_fonts() -> list[dict]:
    rows = []
    paths = sorted((BMC / "figures").iterdir()) + [BMC / "A_main.docx", BMC / "Additional_file_1.docx",
                                                  BMC / "figures_contact_sheet.pdf", BMC / "Additional_file_2.xlsx",
                                                  BMC / "Additional_file_3.xlsx"]
    for p in paths:
        row = {"file": str(p.relative_to(ROOT)), "bytes": p.stat().st_size, "size_ok_10MB": p.stat().st_size <= 10 * 1024 ** 2}
        if p.suffix == ".pdf":
            row["fonts"], row["fonts_embedded"] = pdf_fonts(p)
        rows.append(row)
    return rows


def pog570_ids_in(folder: Path) -> list[tuple[str, int]]:
    """Distinct POG570 patient identifiers (config/pog570_eval_labels.tsv) found as whole tokens in each text file."""
    ids = set(pd.read_csv(ROOT / "config" / "pog570_eval_labels.tsv", sep="\t", dtype=str)["PATIENT_ID"])
    out = []
    for p in sorted(folder.rglob("*")):
        if p.is_file() and "smoke_out" not in p.parts and p.suffix in TEXT_SUFFIXES and p.stat().st_size < 500 * 1024 ** 2:
            out.append((str(p.relative_to(ROOT)), len(ids & set(re.findall(r"\b\d{5}\b", p.read_text(errors="ignore"))))))
    return out


def release_dry_run() -> dict:
    run = subprocess.run(["bash", "release/rebuild_release.sh", "--dry-run", "--url", DUMMY_URL, "--doi", DUMMY_DOI],
                         cwd=ROOT, capture_output=True, text=True)
    text = run.stdout + run.stderr
    scan = pog570_ids_in(ROOT / "release")
    pog_ids = [f"{f}: {n} distinct POG570 patient identifiers" for f, n in scan if n >= POG_MIN]
    named = [str(p.relative_to(ROOT)) for p in (ROOT / "release").rglob("*")
             if re.search("pog", p.name, re.I) and p.name != "SHA256SUMS" and "smoke_out" not in p.parts]
    return {"exit": run.returncode, "output": text.strip(),
            "posthoc_model_included": all(f"would include release/models/MIX_Z0_posthoc.{s}" in text for s in ("joblib", "npz", "json")),
            "pog570_named_files": named, "pog570_sample_ids_in_predictions": pog_ids,
            "pog570_scan": f"{len(scan)} text files scanned; most distinct identifiers in one file: "
                           f"{max((n for _, n in scan), default=0)}"}


def main() -> None:
    text = MAIN.read_text()
    af1 = AF1.read_text()
    items, detail = [], []

    def item(name, value, limit, ok, note=""):
        items.append({"item": name, "value": value, "limit": limit, "pass": bool(ok), "detail": note})

    abstract = section(text, "Abstract")
    with_heads, without = plain_words(abstract), plain_words(re.sub(r"^\*\*\w+\*\*$", "", abstract, flags=re.M))
    item("abstract words (with the three subheadings)", with_heads, LIMITS["abstract"], with_heads <= LIMITS["abstract"],
         f"without subheadings {without}")
    legends = figure_legends()
    for r in legends:
        detail.append({"check": "figure length", **r})
    item("figure titles over 15 words", sum(r["title_words"] > LIMITS["title"] for r in legends), 0,
         all(r["title_words"] <= LIMITS["title"] for r in legends),
         f"{len(legends)} figures; longest title {max(r['title_words'] for r in legends)} words")
    item("figure legends over 300 words", sum(r["legend_words"] > LIMITS["legend"] for r in legends), 0,
         all(r["legend_words"] <= LIMITS["legend"] for r in legends),
         f"longest legend {max(r['legend_words'] for r in legends)} words")
    refs = references(text)
    item("reference first-citation order deviations", refs["order_deviations"], 0, refs["order_deviations"] == 0,
         f"{refs['listed']} listed, {refs['cited']} cited, list numbered 1..n: {refs['listed_sequential']}")
    item("uncited references", len(refs["uncited"]), 0, not refs["uncited"], " ".join(map(str, refs["uncited"])))
    item("cited numbers absent from the list", len(refs["missing"]), 0, not refs["missing"], " ".join(map(str, refs["missing"])))
    orders = citation_orders(text)
    for name in ("figures", "tables", "additional_files"):
        o = orders[name]
        item(f"first-citation order deviations, {name.replace('_', ' ')}", o["deviations"], 0, o["deviations"] == 0,
             "order " + " ".join(map(str, o["first_citation_order"])))
    for name in ("supplementary_tables_in_main", "supplementary_figures_in_main"):
        o = orders[name]
        detail.append({"check": name, "order": " ".join(map(str, o["first_citation_order"])), "deviations": o["deviations"]})
    ab = abbreviations(text)
    item("abbreviations listed but not used", len(ab["listed_not_used"]), 0, not ab["listed_not_used"], " ".join(ab["listed_not_used"]))
    item("abbreviations defined in the text but not listed", len(ab["defined_not_listed"]), 0, not ab["defined_not_listed"],
         " ".join(ab["defined_not_listed"]))
    for a, n in ab["used_counts"].items():
        detail.append({"check": "abbreviation", "text": a, "meaning": ab["listed"][a], "uses_in_body": n,
                       "defined_in_text": a in ab["defined_in_text"]})
    for a in ab["defined_in_text"]:
        if a not in ab["listed"]:
            detail.append({"check": "abbreviation defined, not listed", "text": a,
                           "meaning": "accession or method name, not counted" if a in ab["accessions_or_method_name"] else ""})
    scale = {p.name: len(SCALE.findall(p.read_text())) for p in (MAIN, AF1, COVER)}
    item("'scaled', 'unscaled', 'scaling'", sum(scale.values()), 0, sum(scale.values()) == 0,
         "; ".join(f"{k} {v}" for k, v in scale.items()))
    names = internal_names(af1)
    hits = forbidden(names)
    for h in hits:
        detail.append({"check": "forbidden", **h})
    counted = [h for h in hits if not h["exempt"]]
    for kind in ("stored", "significant", "host-pull", "internal name", "path", "hash"):
        n = sum(h["kind"] == kind for h in counted)
        ex = sum(h["kind"] == kind and h["exempt"] for h in hits)
        item(f"forbidden: {kind}", n, 0, n == 0, f"{ex} in exempt places")
    figs = sorted(p.name for p in (BMC / "figures").iterdir())
    item("files in manuscript/bmc/figures", len(figs), FIGURES_EXPECTED, len(figs) == FIGURES_EXPECTED, " ".join(figs))
    sizes = files_and_fonts()
    for r in sizes:
        detail.append({"check": "file", **r})
    big = [r["file"] for r in sizes if not r["size_ok_10MB"]]
    unembedded = [r["file"] for r in sizes if "fonts" in r and r["fonts"] != r["fonts_embedded"]]
    item("files over 10 MB", len(big), 0, not big, " ".join(big))
    item("PDF files with fonts not embedded", len(unembedded), 0, not unembedded,
         f"{sum('fonts' in r for r in sizes)} PDF files checked")
    ph = placeholders()
    pd.DataFrame(ph).to_csv(CHECKS / f"placeholders_{TAG}.tsv", sep="\t", index=False)
    item("placeholders (listed, not filled)", len(ph), "list", True,
         "; ".join(f"{p['file'].split('/')[-1]}:{p['line']}" for p in ph))
    rel = release_dry_run()
    item("release dry run exit status", rel["exit"], 0, rel["exit"] == 0, "dummy URL and DOI, print only")
    item("release dry run: post hoc model included", rel["posthoc_model_included"], True, rel["posthoc_model_included"])
    n_pog = len(rel["pog570_named_files"]) + len(rel["pog570_sample_ids_in_predictions"])
    item("release: POG570 sample-level files", n_pog, 0, n_pog == 0,
         "; ".join(rel["pog570_named_files"] + rel["pog570_sample_ids_in_predictions"] + [rel["pog570_scan"]]))
    (CHECKS / f"release_dry_run_{TAG}.txt").write_text(rel["output"] + "\n")
    frame = pd.DataFrame(items)
    frame.to_csv(CHECKS / f"hygiene_{TAG}.tsv", sep="\t", index=False)
    pd.DataFrame(detail).to_csv(CHECKS / f"hygiene_detail_{TAG}.tsv", sep="\t", index=False)
    md = [f"# Hygiene checklist ({TAG})", "", "| Item | Value | Limit | Pass | Detail |", "|---|---|---|---|---|"]
    for r in items:
        md.append(f"| {r['item']} | {r['value']} | {r['limit']} | {'yes' if r['pass'] else 'no'} | {str(r['detail']).replace('|', '/')} |")
    md += ["", "## Figure titles and legends", "", "| Figure | Title words | Legend words |", "|---|---|---|"]
    md += [f"| {r['figure']} | {r['title_words']} | {r['legend_words']} |" for r in legends]
    md += ["", "## Placeholders", "", "| File | Line | Section | Text |", "|---|---|---|---|"]
    md += [f"| {p['file']} | {p['line']} | {p['heading']} | {p['text']} |" for p in ph]
    (CHECKS / f"hygiene_checklist_{TAG}.md").write_text("\n".join(md) + "\n")
    print(frame.to_string())


if __name__ == "__main__":
    main()
