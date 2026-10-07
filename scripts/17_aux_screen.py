#!/usr/bin/env python3
"""Screen cBioPortal and GEO for auxiliary cohorts. No model is applied.

Eligibility was committed before this script runs.
"""
from __future__ import annotations

import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd

from stage5_lib import load_paths, stage5_dir

BASE_GEO = (
    '(metastasis[Title] OR metastases[Title] OR metastatic[Title]) AND "Homo sapiens"[Organism] '
    'AND gse[Entry Type] AND ("expression profiling by high throughput sequencing"[DataSet Type] '
    'OR "expression profiling by array"[DataSet Type])'
)
SITE_WORDS = [
    "liver", '"lymph node"', "lung", "bone", "brain", "peritoneal", "omentum", "adrenal", "skin",
]
SITE_ATTR = re.compile(
    r"biopsy.?site|tissue.?site|tissue.?source|metastatic.?site|site.?of.?met|sample.?site|"
    r"anatomic|anatomical|specimen.?site|collection.?site|tumor.?site|tumour.?site|"
    r"site.?of.?biopsy|distant.?site|organ.?site|biopsy.?location|sample.?location|tissue.?location",
    re.I,
)
EXCLUDE_ID = re.compile(r"tcga|gtex|ccle|met500|pog570|pog_570", re.I)
CELL_TITLE = re.compile(r"cell line|organoid|\bpdx\b|xenograft|single[- ]cell|spatial transcript|pbmc|plasma|serum|whole blood|peripheral blood", re.I)


def get_json(url: str):
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "projA-aux-screen"})
    with urllib.request.urlopen(request, timeout=120) as handle:
        return json.load(handle)


def screen_cbio(rows: list[dict]) -> None:
    profiles = get_json("https://www.cbioportal.org/api/molecular-profiles?projection=SUMMARY&pageSize=20000&pageNumber=0")
    studies = {}
    for profile in profiles:
        if profile.get("molecularAlterationType") != "MRNA_EXPRESSION":
            continue
        if profile.get("datatype") == "Z-SCORE":
            continue
        studies.setdefault(profile["studyId"], []).append(profile["molecularProfileId"])
    catalog = {row["studyId"]: row for row in get_json("https://www.cbioportal.org/api/studies?pageSize=20000&pageNumber=0")}
    for study_id, profile_ids in sorted(studies.items()):
        meta = catalog.get(study_id, {})
        title = meta.get("name", "")
        n = meta.get("allSampleCount")
        memo = "profiles=" + ",".join(profile_ids)
        if EXCLUDE_ID.search(study_id) or EXCLUDE_ID.search(title):
            rows.append(record("cBioPortal", "molecular-profiles MRNA_EXPRESSION datatype!=Z-SCORE", study_id, title, n, ",".join(profile_ids), "제외", "7", memo))
            continue
        if CELL_TITLE.search(title):
            rows.append(record("cBioPortal", "molecular-profiles MRNA_EXPRESSION datatype!=Z-SCORE", study_id, title, n, ",".join(profile_ids), "제외", "2", memo))
            continue
        time.sleep(0.05)
        attrs = get_json(f"https://www.cbioportal.org/api/studies/{urllib.parse.quote(study_id)}/clinical-attributes?projection=SUMMARY")
        matched = []
        for attr in attrs:
            if attr.get("patientAttribute") is True:
                continue
            text = f"{attr.get('clinicalAttributeId','')} {attr.get('displayName','')}"
            if SITE_ATTR.search(text):
                matched.append(attr.get("clinicalAttributeId"))
        if not matched:
            rows.append(record("cBioPortal", "molecular-profiles MRNA_EXPRESSION datatype!=Z-SCORE", study_id, title, n, ",".join(profile_ids), "제외", "3", memo))
            continue
        rows.append(record(
            "cBioPortal", "molecular-profiles MRNA_EXPRESSION datatype!=Z-SCORE", study_id, title, n,
            ",".join(profile_ids), "선별", "", memo + "; site_attributes=" + ",".join(matched),
        ))
        print("cbio", study_id, matched, flush=True)


def screen_geo(rows: list[dict]) -> None:
    seen = set()
    for word in SITE_WORDS:
        term = f"{BASE_GEO} AND {word}[Title]"
        query = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode({
            "db": "gds", "retmode": "json", "retmax": "50", "term": term, "tool": "projA", "email": "aa@aa.com",
        })
        time.sleep(0.4)
        found = get_json(query)
        ids = found["esearchresult"]["idlist"]
        print("geo", word, found["esearchresult"]["count"], len(ids), flush=True)
        for start in range(0, len(ids), 40):
            chunk = ids[start:start + 40]
            if not chunk:
                continue
            time.sleep(0.4)
            summary = get_json(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
                + urllib.parse.urlencode({
                    "db": "gds", "retmode": "json", "id": ",".join(chunk), "tool": "projA", "email": "aa@aa.com",
                })
            )
            result = summary["result"]
            for uid in result.get("uids", []):
                item = result[uid]
                acc = item.get("accession")
                if acc in seen:
                    continue
                seen.add(acc)
                title = item.get("title") or ""
                summary_text = item.get("summary") or ""
                blob = title + " " + summary_text
                n = item.get("n_samples")
                platform = f"GPL{item.get('gpl')} {item.get('ptechtype') or ''}".strip()
                if CELL_TITLE.search(blob):
                    decision, failed, memo = "제외", "2", "title or summary"
                else:
                    decision, failed, memo = "보류", "", "title and summary only; sample characteristics not read"
                rows.append(record("GEO", term, acc, title, n, platform, decision, failed, memo))


def record(path, query, ident, title, n, platform, decision, failed, memo) -> dict:
    return {
        "path": path, "query": query, "id": ident, "title": title, "n": n,
        "platform": platform, "decision": decision, "failed_criterion": failed, "memo": memo,
    }


def main() -> None:
    paths = load_paths()
    root = stage5_dir(paths) / "aux"
    root.mkdir(parents=True, exist_ok=True)
    rows = []
    screen_cbio(rows)
    pd.DataFrame(rows).to_csv(root / "screening_log.tsv", sep="\t", index=False)
    screen_geo(rows)
    pd.DataFrame(rows).to_csv(root / "screening_log.tsv", sep="\t", index=False)
    print("screen rows", len(rows), flush=True)


if __name__ == "__main__":
    main()
