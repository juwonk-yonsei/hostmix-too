"""Check manuscript DOIs against Crossref. Does not edit the manuscript."""
from __future__ import annotations

import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DRAFT = ROOT / "manuscript/A_manuscript_draft.md"
OUT = ROOT / "manuscript/checks"
UA = "projA-refcheck/1.0 (mailto:aa@aa.com)"


def get(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(request, timeout=60) as handle:
        return json.loads(handle.read().decode())


def work(doi: str) -> dict:
    return get("https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/()"))["message"]


def year_of(item: dict) -> str:
    parts = (item.get("issued") or {}).get("date-parts") or [[None]]
    return str(parts[0][0]) if parts and parts[0] else ""


def first_author(item: dict) -> str:
    authors = item.get("author") or []
    if not authors:
        return ""
    return str(authors[0].get("family") or "")


def title_of(item: dict) -> str:
    return " ".join(item.get("title") or [])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    text = DRAFT.read_text()
    dois = []
    for match in re.finditer(r"10\.\d{4,9}/[^\s,;]+", text):
        doi = match.group(0).rstrip(".").rstrip(")").lower()
        if doi not in dois:
            dois.append(doi)
    rows = []
    body, _, refs = text.partition("## References")
    mismatches = 0
    for doi in dois:
        item = work(doi)
        author = first_author(item)
        year = year_of(item)
        title = title_of(item)
        journal = "; ".join(item.get("container-title") or [])
        # The manuscript line that contains this DOI should also contain the year and, when present, the author.
        window = ""
        for line in text.splitlines():
            if doi in line.lower() or doi in line:
                window = line
                break
        author_ok = (not author) or (author.lower() in window.lower()) or (author.lower() in body.lower())
        year_ok = (not year) or (year in window) or (year in text)
        title_tokens = [tok for tok in re.findall(r"[A-Za-z]{4,}", title.lower()) if tok not in {"with", "from", "that", "this", "using", "into"}]
        title_hit = sum(tok in text.lower() for tok in title_tokens[:8])
        status = "match" if author_ok and year_ok else "mismatch"
        if status != "match":
            mismatches += 1
        rows.append({
            "doi": doi,
            "crossref_year": year,
            "crossref_first_author": author,
            "crossref_journal": journal,
            "crossref_title": title,
            "manuscript_line": window,
            "author_ok": author_ok,
            "year_ok": year_ok,
            "title_token_hits_in_draft": title_hit,
            "status": status,
        })
        print(status, doi, year, author, flush=True)
    cited = []
    for line in refs.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        surname = re.match(r"([A-Za-z]+)", line.strip())
        name = surname.group(1) if surname else ""
        in_body = name.lower() in body.lower() if name else False
        # "Science" is the journal name on the metastatic-adaptation line, not an author.
        cited.append({"reference_line": line.strip(), "lead_token": name, "lead_token_in_body": in_body})
    verify_left = text.count("[DOI to verify]")
    needed_left = text.count("[CITATION NEEDED]")
    report = {
        "n_doi": len(rows),
        "n_mismatch": mismatches,
        "n_doi_to_verify_left": verify_left,
        "n_citation_needed_left": needed_left,
        "rows": rows,
        "reference_lines": cited,
    }
    (OUT / "refs_check.json").write_text(json.dumps(report, indent=2) + "\n")
    print("dois", len(rows), "mismatch", mismatches, "verify_left", verify_left, flush=True)


if __name__ == "__main__":
    main()
