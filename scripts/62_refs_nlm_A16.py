#!/usr/bin/env python3
"""Journal abbreviations, page ranges and a Crossref recheck for the numbered reference list.

Abbreviation sources, in order: MedAbbr of NLM J_Medline.txt matched by an ISSN of the Crossref
record; the ISO abbreviation of the NLM Catalog (E-utilities) for that ISSN; the Crossref
short-container-title; the original name. Periods are removed from abbreviations. Journals are
never matched by title. Reference numbers, authors, volume(issue) and DOIs are left as they are.

Usage: 62_refs_nlm_A16.py [--write]   (without --write only the check tables are written)
"""
from __future__ import annotations

import html
import json
import re
import sys
import time
import unicodedata
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "manuscript" / "bmc" / "A_main.md"
MEDLINE = ROOT / "data" / "reference" / "nlm" / "J_Medline.txt"
CACHE = ROOT / "data" / "reference" / "crossref_A16"
CHECKS = ROOT / "manuscript" / "checks"
AGENT = "projA-refs/1.0 (mailto:aa@aa.com)"
EXPECTED = ["Br J Cancer", "JAMA Netw Open", "EBioMedicine", "J Mol Diagn", "Cancer Discov", "Sci Rep",
            "Nat Commun", "Nat Med", "Sci Adv", "Cancers (Basel)", "Nat Biotechnol", "Nat Cancer", "Sci Signal",
            "Cancer Res", "Nucleic Acids Res", "J Immunother Cancer", "Mol Oncol", "Proc Natl Acad Sci U S A",
            "Nat Genet", "Clin Exp Metastasis", "Control Clin Trials", "N Engl J Med", "Nat Methods"]
CITATION = re.compile(r"[.?!] (?P<journal>[^.?!]+?)\. (?P<year>\d{4});(?P<volume>[^:(;]+)?(?:\((?P<issue>[^)]+)\))?:"
                      r"(?P<pages>[A-Za-z0-9]+(?:[-\u2013][A-Za-z0-9]+)?)\. ")


def get_json(url: str) -> dict:
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": AGENT, "Accept": "application/json"})
            with urllib.request.urlopen(request, timeout=60) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as error:  # noqa: BLE001
            if attempt == 3:
                raise SystemExit(f"{url}: {error}")
            time.sleep(2 * (attempt + 1))
    raise AssertionError


def crossref(number: int, doi: str) -> dict:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{number:02d}.json"
    if not path.exists():
        path.write_text(json.dumps(get_json("https://api.crossref.org/works/" + urllib.parse.quote(doi))))
        time.sleep(0.2)
    return json.loads(path.read_text())["message"]


def medline() -> dict[str, list[str]]:
    by_issn: dict[str, list[str]] = {}
    for block in MEDLINE.read_text(encoding="utf-8").split("-" * 56):
        fields = dict(line.split(": ", 1) for line in block.strip().splitlines() if ": " in line)
        for key in ("ISSN (Print)", "ISSN (Online)"):
            issn = fields.get(key, "").strip()
            if issn and fields.get("MedAbbr"):
                by_issn.setdefault(issn, [])
                if fields["MedAbbr"] not in by_issn[issn]:
                    by_issn[issn].append(fields["MedAbbr"])
    return by_issn


def catalog_iso(issn: str) -> str:
    found = get_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
                     + urllib.parse.urlencode({"db": "nlmcatalog", "term": f"{issn}[ISSN]", "retmode": "json"}))
    ids = found["esearchresult"]["idlist"]
    time.sleep(0.4)
    if not ids:
        return ""
    summary = get_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
                       + urllib.parse.urlencode({"db": "nlmcatalog", "id": ",".join(ids), "retmode": "json"}))
    time.sleep(0.4)
    values = {summary["result"][i].get("isoabbreviation", "") for i in ids} - {""}
    return values.pop() if len(values) == 1 else ""


def no_periods(text: str) -> str:
    return re.sub(r"\s+", " ", text.replace(".", " ")).strip()


def norm(text: str) -> str:
    text = html.unescape(re.sub(r"<[^>]+>", "", text))
    text = unicodedata.normalize("NFKC", text).casefold()
    text = re.sub(r"[\u2010-\u2015\u2212-]", "-", text)
    text = re.sub(r"[\u2018\u2019\u201c\u201d'\"]", "", text)
    return re.sub(r"\s+", " ", text).strip().rstrip(".")


def references(lines: list[str]) -> tuple[int, int]:
    start = next(i for i, line in enumerate(lines) if line.strip() == "## References")
    end = next(i for i in range(start + 1, len(lines)) if lines[i].startswith("## "))
    return start, end


def main() -> None:
    write = "--write" in sys.argv[1:]
    lines = MAIN.read_text(encoding="utf-8").split("\n")
    start, end = references(lines)
    by_issn = medline()
    mapping, checks = [], []
    for index in range(start + 1, end):
        hit = re.match(r"^(\d+)\. (.*)$", lines[index])
        if not hit:
            continue
        number, text = int(hit.group(1)), hit.group(2)
        doi_hit = re.search(r"doi:(\S+?)\.?$", text)
        record = crossref(number, doi_hit.group(1)) if doi_hit else None
        new = text
        if record and record.get("type") == "proceedings-article":
            names = [name for name in record.get("container-title", []) if re.search(r"\d+(st|nd|rd|th) ", name)]
            editors = record.get("editor") or []
            if len(names) != 1 or editors:
                raise SystemExit(f"{number}: proceedings record needs a decision: {names} {editors}")
            year = record["issued"]["date-parts"][0][0]
            pages = record["page"].replace("-", "\u2013")
            authors_title = text.split(". Proceedings of")[0]
            new = f"{authors_title}. In: {names[0]}. {year}. p. {pages}. doi:{doi_hit.group(1)}."
            mapping.append({"number": number, "original": "Proceedings of the Python in Science Conference",
                            "issn": ";".join(record.get("ISSN", [])), "abbreviation": names[0],
                            "source": "conference format; meeting name from the DOI record, no editors in it"})
        elif record:
            cite = CITATION.search(text)
            if not cite:
                raise SystemExit(f"{number}: citation block not parsed: {text}")
            original = cite.group("journal")
            issns = [item["value"] for item in record.get("issn-type", [])] or record.get("ISSN", [])
            abbreviation, source = "", ""
            medabbr = sorted({abbr for issn in issns for abbr in by_issn.get(issn, [])})
            if len(medabbr) > 1:
                raise SystemExit(f"{number}: ISSNs {issns} give several MedAbbr {medabbr}")
            if medabbr:
                abbreviation, source = medabbr[0], "J_Medline MedAbbr by ISSN"
            if not abbreviation:
                for issn in issns:
                    iso = catalog_iso(issn)
                    if iso:
                        abbreviation, source = no_periods(iso), f"NLM Catalog ISO abbreviation ({issn})"
                        break
            if not abbreviation and record.get("short-container-title"):
                abbreviation, source = no_periods(record["short-container-title"][0]), "Crossref short-container-title"
            if not abbreviation:
                abbreviation, source = original, "original name"
            mapping.append({"number": number, "original": original, "issn": ";".join(issns),
                            "abbreviation": abbreviation, "source": source})
            pages = cite.group("pages").replace("-", "\u2013")
            span = cite.span("journal")
            page_span = cite.span("pages")
            new = text[:span[0]] + abbreviation + text[span[1]:page_span[0]] + pages + text[page_span[1]:]
        else:
            cite = CITATION.search(text)
            if cite:
                page_span = cite.span("pages")
                new = text[:page_span[0]] + cite.group("pages").replace("-", "\u2013") + text[page_span[1]:]
                mapping.append({"number": number, "original": cite.group("journal"), "issn": "",
                                "abbreviation": cite.group("journal"),
                                "source": "no DOI and no Crossref record; name left as in the text"})

        if record:
            cite = CITATION.search(new)
            title = (record.get("title") or [""])[0]
            years = {part.get("date-parts", [[None]])[0][0] for key, part in record.items()
                     if key in ("issued", "published-print", "published-online", "published") and isinstance(part, dict)}
            years.discard(None)
            year_text = cite.group("year") if cite else re.search(r"\. (\d{4})\. p\. ", new).group(1)
            ref_pages = (cite.group("pages") if cite else re.search(r"p\. (\S+)\. doi", new).group(1))
            row = {"number": number, "doi": doi_hit.group(1),
                   "doi_match": norm(record.get("DOI", "")) == norm(doi_hit.group(1)),
                   "title_match": norm(title) in norm(new),
                   "year_text": int(year_text), "year_issued": record["issued"]["date-parts"][0][0],
                   "year_any_crossref_date": int(year_text) in years,
                   "volume_text": (cite.group("volume") if cite else None), "volume_crossref": record.get("volume"),
                   "issue_text": (cite.group("issue") if cite else None), "issue_crossref": record.get("issue"),
                   "pages_text": ref_pages, "pages_crossref": record.get("page") or record.get("article-number"),
                   "crossref_title": title}
            row["volume_match"] = (row["volume_text"] or None) == (record.get("volume") or None)
            row["issue_match"] = (row["issue_text"] or None) == (record.get("issue") or None)
            row["pages_match"] = norm(ref_pages or "") == norm(row["pages_crossref"] or "")
            checks.append(row)
        lines[index] = f"{number}. {new}"

    table = pd.DataFrame(mapping)
    table.to_csv(CHECKS / "refs_nlm_A16.tsv", sep="\t", index=False)
    check = pd.DataFrame(checks)
    check.to_csv(CHECKS / "refs_crossref_A16.tsv", sep="\t", index=False)
    got = set(table["abbreviation"])
    missing = [name for name in EXPECTED if name not in got]
    print(table.to_string(index=False))
    print("expected abbreviations not produced:", missing)
    fields = ["doi_match", "title_match", "year_any_crossref_date", "volume_match", "issue_match", "pages_match"]
    bad = check.loc[~check[fields].all(axis=1)]
    print(len(check), "Crossref records;", len(bad), "with a field that differs")
    if len(bad):
        print(bad.drop(columns=["crossref_title"]).to_string(index=False))
    if write:
        MAIN.write_text("\n".join(lines), encoding="utf-8")
        print("wrote", MAIN)


if __name__ == "__main__":
    main()
