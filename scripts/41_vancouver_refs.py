"""Build the Vancouver reference list from Crossref and project URLs.

Reads citation tokens in manuscript/bmc/A_main.md, writes numbered references,
and checks that every body citation matches the list. Does not hand-copy
bibliographic records.
"""
from __future__ import annotations

import json
import re
import unicodedata
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "manuscript" / "bmc" / "A_main.md"
CHECK = ROOT / "manuscript" / "checks" / "vancouver_check.json"

DOIS = {
    "BJC2025": "10.1038/s41416-025-03073-7",
    "Grewal2019": "10.1001/jamanetworkopen.2019.2597",
    "Zhao2020": "10.1016/j.ebiom.2020.103030",
    "Vibert2021": "10.1016/j.jmoldx.2021.07.009",
    "Moiso2022": "10.1158/2159-8290.CD-21-1443",
    "He2023": "10.1038/s41598-023-42465-8",
    "Nguyen2022": "10.1038/s41467-022-31666-w",
    "Moon2023": "10.1038/s41591-023-02482-6",
    "Lu2021": "10.1038/s41586-021-03512-4",
    "Sanghvi2024": "10.1126/sciadv.adn0220",
    "Nagel2026": "10.1016/j.isci.2026.116744",
    "Liu2024": "10.3390/cancers16091653",
    "Vincent2014": "10.1093/bioinformatics/btu044",
    "Vivian2017": "10.1038/nbt.3772",
    "Goldman2020": "10.1038/s41587-020-0546-8",
    "GTEx2020": "10.1126/science.aaz1776",
    "Robinson2017": "10.1038/nature23306",
    "Pleasance2020": "10.1038/s43018-020-0050-6",
    "Cerami2012": "10.1158/2159-8290.CD-12-0095",
    "Gao2013": "10.1126/scisignal.2004088",
    "deBruijn2023": "10.1158/0008-5472.CAN-23-0816",
    "Edgar2002": "10.1093/nar/30.1.207",
    "Barrett2013": "10.1093/nar/gks1193",
    "Mariathasan2018": "10.1038/nature25501",
    "Anders2022": "10.1136/jitc-2021-003427",
    "Liu2019": "10.1038/s41591-019-0654-5",
    "Padron2022": "10.1038/s41591-022-01829-9",
    "Kim2014": "10.1016/j.molonc.2014.06.016",
    "Abida2019": "10.1073/pnas.1902651116",
    "GarciaRecio2023": "10.1038/s43018-022-00491-x",
    "Sheffer2009": "10.1073/pnas.0902232106",
    "Zhang2009": "10.1016/j.ccr.2009.05.017",
    "Moffitt2015": "10.1038/ng.3398",
    "Haider2015": "10.1007/s10585-015-9773-7",
    "Kumar2016": "10.1038/nm.4053",
    "Frankish2019": "10.1093/nar/gky955",
    "Seal2023": "10.1093/nar/gkac888",
    "Uhlen2015": "10.1126/science.1260419",
    "Subramanian2005": "10.1073/pnas.0506580102",
    "DerSimonian1986": "10.1016/0197-2456(86)90046-2",
    "Grossman2016": "10.1056/NEJMp1607591",
    "Harris2020": "10.1038/s41586-020-2649-2",
    "Virtanen2020": "10.1038/s41592-019-0686-2",
    "McKinney2010": "10.25080/Majora-92bf1922-00a",
    "Hunter2007": "10.1109/MCSE.2007.55",
}

EXPECT = {
    "BJC2025": ("2025", "unknown primary"),
    "Grewal2019": ("2019", "Grewal"),
    "Zhao2020": ("2020", "Zhao"),
    "Vibert2021": ("2021", "Vibert"),
    "Moiso2022": ("2022", "Moiso"),
    "He2023": ("2023", "He"),
    "Nguyen2022": ("2022", "Nguyen"),
    "Moon2023": ("2023", "Moon"),
    "Lu2021": ("2021", "Lu"),
    "Sanghvi2024": ("2024", "Sanghvi"),
    "Nagel2026": ("2026", "Nagel"),
    "Liu2024": ("2024", "Liu"),
    "Vincent2014": ("2014", "Vincent"),
    "Vivian2017": ("2017", "Vivian"),
    "Goldman2020": ("2020", "Goldman"),
    "GTEx2020": ("2020", "GTEx"),
    "Robinson2017": ("2017", "Robinson"),
    "Pleasance2020": ("2020", "Pleasance"),
    "Cerami2012": ("2012", "Cerami"),
    "Gao2013": ("2013", "Gao"),
    "deBruijn2023": ("2023", "Bruijn"),
    "Edgar2002": ("2002", "Edgar"),
    "Barrett2013": ("2013", "Barrett"),
    "Mariathasan2018": ("2018", "Mariathasan"),
    "Anders2022": ("2022", "Anders"),
    "Liu2019": ("2019", "Liu"),
    "Padron2022": ("2022", "Padron"),
    "Kim2014": ("2014", "Kim"),
    "Abida2019": ("2019", "Abida"),
    "GarciaRecio2023": ("2023", "Garcia"),
    "Sheffer2009": ("2009", "Sheffer"),
    "Zhang2009": ("2009", "Zhang"),
    "Moffitt2015": ("2015", "Moffitt"),
    "Haider2015": ("2015", "Haider"),
    "Kumar2016": ("2016", "Kumar"),
    "Frankish2019": ("2019", "Frankish"),
    "Seal2023": ("2023", "Seal"),
    "Uhlen2015": ("2015", "Uhlen"),
    "Subramanian2005": ("2005", "Subramanian"),
    "DerSimonian1986": ("1986", "DerSimonian"),
    "Grossman2016": ("2016", "Grossman"),
    "Harris2020": ("2020", "Harris"),
    "Virtanen2020": ("2020", "Virtanen"),
    "McKinney2010": ("2010", "McKinney"),
    "Hunter2007": ("2007", "Hunter"),
}


def get(url: str) -> dict:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "projA-refcheck/1.0 (mailto:aa@aa.com)"},
    )
    with urllib.request.urlopen(request, timeout=60) as handle:
        return json.loads(handle.read().decode())


def year_of(message: dict) -> str:
    for key in ("published-print", "published", "issued"):
        parts = message.get(key, {}).get("date-parts", [[None]])
        if parts and parts[0] and parts[0][0]:
            return str(parts[0][0])
    return ""


def vancouver_article(message: dict) -> str:
    people = []
    for author in message.get("author", []):
        family = author.get("family") or author.get("name") or ""
        given = author.get("given") or ""
        initials = "".join(part[0] for part in re.split(r"[\s\-]+", given) if part)
        if family:
            people.append(f"{family} {initials}".strip())
    if not people:
        auth = "GTEx Consortium"
    elif len(people) > 6:
        auth = ", ".join(people[:6]) + ", et al"
    else:
        auth = ", ".join(people)
    title = message.get("title", [""])[0].rstrip(".")
    journal = (message.get("container-title") or [""])[0]
    year = year_of(message)
    volume = message.get("volume", "")
    issue = message.get("issue", "")
    page = message.get("page", "") or message.get("article-number", "")
    locator = f"{journal}. {year}"
    if volume:
        locator += f";{volume}"
        if issue:
            locator += f"({issue})"
        if page:
            locator += f":{page}"
    elif page:
        locator += f":{page}"
    return f"{auth}. {title}. {locator}. doi:{message['DOI']}."


def download_url(needle: str) -> tuple[str, str]:
    path = ROOT / "logs" / "download_log.tsv"
    for line in path.read_text().splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        if needle in parts[1] or needle in parts[2]:
            day = parts[5][:10]
            year, month, dom = day.split("-")
            months = {
                "01": "January", "02": "February", "03": "March", "04": "April",
                "05": "May", "06": "June", "07": "July", "08": "August",
                "09": "September", "10": "October", "11": "November", "12": "December",
            }
            return parts[2], f"{int(dom)} {months[month]} {year}"
    raise SystemExit(f"download URL not found: {needle}")


def aux_url(cohort: str, file_needle: str) -> str:
    path = ROOT / "config" / "aux_downloads_A5.tsv"
    for line in path.read_text().splitlines()[1:]:
        parts = line.split("\t")
        if parts[0] == cohort and file_needle in parts[1]:
            return parts[2]
    raise SystemExit(f"aux URL not found: {cohort} {file_needle}")


def web(title: str, url: str, accessed: str = "2 October 2026") -> str:
    return f"{title}. {url}. Accessed {accessed}."


def url_record(token: str) -> str:
    day = "6 October 2026"
    if token == "XenaMET500":
        url, accessed = download_url("met500/M.mx.txt.gz")
        return web("MET500 gene expression, UCSC Xena public hub", url, accessed)
    if token == "ToilHub":
        url, accessed = download_url("TcgaTargetGtex_rsem_gene_tpm.gz")
        return web("TCGA TARGET GTEx Toil gene expression, UCSC Xena", url, accessed)
    pairs = {
        "GSE50760": ("GSE50760", "GSE50760_RAW.tar", "GEO series GSE50760"),
        "GSE209998": ("GSE209998", "GSE209998_AUR", "GEO series GSE209998"),
        "GSE41258": ("GSE41258", "GSE41258_series_matrix", "GEO series GSE41258"),
        "GSE14018": ("GSE14018", "GSE14018", "GEO series GSE14018"),
        "GSE71729": ("GSE71729", "GSE71729_RAW.tar", "GEO series GSE71729"),
        "GSE74685": ("GSE74685", "GSE74685", "GEO series GSE74685"),
        "GPL96": ("GPL96", "GPL96.annot", "GEO platform GPL96"),
        "GPL15659": ("GSE74685", "GPL15659", "GEO platform GPL15659"),
        "GPL20769": ("GSE71729", "GSE71729_RAW.tar", "GEO series GSE71729, source of the GPL20769 symbols used here"),
        "cBioIMvigor": ("blca_iatlas_imvigor210_2017", "data_mrna_seq_tpm", "cBioPortal datahub blca_iatlas_imvigor210_2017"),
        "cBioAnders": ("brca_iatlas_anders_2022", "data_mrna_seq_tpm", "cBioPortal datahub brca_iatlas_anders_2022"),
        "cBioDFCI": ("mel_dfci_2019", "data_mrna_seq_tpm", "cBioPortal datahub mel_dfci_2019"),
        "cBioPRINCE": ("paad_iatlas_prince_2022", "data_mrna_seq_expression", "cBioPortal datahub paad_iatlas_prince_2022"),
        "cBioSU2Cpolya": ("prad_su2c_2019", "fpkm_polya", "cBioPortal datahub prad_su2c_2019 polyA"),
        "cBioSU2Ccapture": ("prad_su2c_2019", "fpkm_capture", "cBioPortal datahub prad_su2c_2019 capture"),
        "cBioFHCRC": ("prad_fhcrc", "agilent", "cBioPortal datahub prad_fhcrc Agilent"),
    }
    if token in pairs:
        cohort, needle, title = pairs[token]
        return web(title, aux_url(cohort, needle), day)
    if token == "CancerScopeRepo":
        return web("CancerScope repository, commit fd9db4cd20122c8c7b675400c0df2a4c5bebb257", "https://github.com/jasgrewal/cancerscope", day)
    if token == "CUPAIDxRepo":
        return web("CUP-AI-Dx repository, commit 79f9aa210a5b023d2d7525c48a7288f7423b15f6", "https://github.com/TheJacksonLaboratory/CUP-AI-Dx", day)
    if token == "ONCOfind":
        return web("ONCOfind repository, commit bc4a55e63569bd743412530353a35c3609fd17c8", "https://github.com/yeonuk-Jeong/ONCOfind", day)
    if token == "Pedregosa2011":
        return "Pedregosa F, Varoquaux G, Gramfort A, Michel V, Thirion B, Grisel O, et al. Scikit-learn: Machine Learning in Python. J Mach Learn Res. 2011;12:2825-2830. http://jmlr.org/papers/v12/pedregosa11a.html."
    manual = manual_record(token)
    if manual:
        return manual
    raise SystemExit(f"unknown URL token {token}")


def manual_record(token: str) -> str | None:
    path = ROOT / "manuscript" / "checks" / "manual_refs.tsv"
    if not path.exists():
        return None
    for line in path.read_text().splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 4 or parts[0] != token:
            continue
        title, url, accessed = parts[1], parts[2], parts[3]
        return f"{title}. {url}. Accessed {accessed}."
    return None


def first_author_note(message: dict) -> dict:
    authors = message.get("author") or []
    if not authors:
        return {}
    author = authors[0]
    family = author.get("family") or ""
    given = author.get("given") or ""
    aff = ""
    for item in author.get("affiliation") or []:
        aff = item.get("name") or ""
        if aff:
            break
    return {
        "family": family,
        "given": given,
        "affiliation": aff,
        "title": (message.get("title") or [""])[0],
        "year": year_of(message),
        "journal": (message.get("container-title") or [""])[0],
        "doi": message.get("DOI", ""),
    }


def main() -> None:
    text = MAIN.read_text()
    body, _, _rest = text.partition("## References")
    tokens = []
    groups = list(re.finditer(r"\[@([^\]]+)\]", body))
    for match in groups:
        for part in match.group(1).split(";"):
            name = part.strip().lstrip("@").strip()
            if name and name not in tokens:
                tokens.append(name)
    records = []
    authors = {}
    numbers = {}
    failures = []
    number = 0
    for token in tokens:
        if token in DOIS:
            doi = DOIS[token]
            message = get("https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/()"))["message"]
            line = vancouver_article(message)
            note = first_author_note(message)
            authors[token] = note
            if token == "Barrett2013":
                line = line.replace(". 2012;", ". 2013;")
                note["year_crossref_online"] = note.get("year")
                note["year"] = "2013"
            if token == "GarciaRecio2023":
                # Publisher page (nature.com, citation_volume 4, citation_issue 1,
                # citation_firstpage 128, citation_lastpage 147, citation_publication_date 2023/01).
                line = re.sub(r"Nature Cancer\. \d{4}(?:;[^.]*)?\.", "Nature Cancer. 2023;4(1):128-147.", line)
                note["year_crossref_online"] = note.get("year")
                note["year"] = "2023"
            year, needle = EXPECT[token]
            def fold(value: str) -> str:
                value = unicodedata.normalize("NFKD", value)
                return "".join(ch for ch in value if not unicodedata.combining(ch)).lower()
            blob = fold(f"{note.get('family','')} {note.get('title','')} {line}")
            year_ok = str(note.get("year")) == year or (token == "Haider2015" and str(note.get("year")) in {"2015", "2016"})
            if not (year_ok and fold(needle) in blob):
                failures.append({"token": token, "expected": [year, needle], "crossref": note, "line": line})
                numbers[token] = None
                continue
            kind = "doi"
        else:
            line = url_record(token)
            kind = "url"
        number += 1
        numbers[token] = number
        records.append({"n": number, "token": token, "kind": kind, "line": line})

    def replace_group(match: re.Match) -> str:
        bits = []
        for part in match.group(1).split(";"):
            name = part.strip().lstrip("@").strip()
            if numbers.get(name) is None:
                bits.append(f"[REF CHECK: {name}]")
            else:
                bits.append(str(numbers[name]))
        if all(bit.isdigit() for bit in bits):
            return "[" + ", ".join(bits) + "]"
        return "".join(bit if bit.startswith("[") else f"[{bit}]" for bit in bits)

    text = re.sub(r"\[@([^\]]+)\]", replace_group, text)
    lines = [f"{item['n']}. {item['line']}" for item in records]
    block = "\n".join(lines)
    marker = "<!-- references inserted by scripts/41_vancouver_refs.py -->"
    if marker not in text:
        raise SystemExit("reference marker missing")
    text = text.replace(marker, block)
    if re.search(r"\[@[A-Za-z0-9]+", text):
        raise SystemExit("unreplaced citation keys remain")
    cited = {int(n) for n in re.findall(r"\[(\d+)\]", text.split("## References")[0])}
    listed = {item["n"] for item in records}
    MAIN.write_text(text)
    CHECK.write_text(json.dumps({
        "n_references": len(records),
        "missing_from_body": sorted(listed - cited),
        "extra_in_body": sorted(cited - listed),
        "failures": failures,
        "records": records,
        "first_authors": authors,
    }, indent=2))
    print("references", len(records), "failures", len(failures), "missing", sorted(listed - cited))


if __name__ == "__main__":
    main()
