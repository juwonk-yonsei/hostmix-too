"""Number extraction for the A17 check: documents, sections, sentences, tokens and exemptions.

Every number display in a document becomes one token with a position (document, section,
paragraph, sentence, token index within the sentence). Tokens are classified by form
(integer, decimal, percent, range, scientific, thousands, word, fraction, alnum) and,
where a fixed rule applies, by exemption category (REF, ID, NAME, VER, DATE, FORMULA,
LEVEL). Exemption rules look only at the form of the token and its neighbouring words,
never at a result value.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = {
    "main": "manuscript/bmc/A_main.md",
    "af1": "manuscript/bmc/Additional_file_1.md",
    "cover": "manuscript/bmc/cover_letter.md",
    "card": "release/models/MODEL_CARD.md",
}

MINUS = "\u2212"
WORDS = ("three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen "
         "seventeen eighteen nineteen twenty thirty forty fifty sixty seventy eighty ninety hundred "
         "thousand million").split()
WORD_VALUE = {w: v for w, v in zip(WORDS, [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                           30, 40, 50, 60, 70, 80, 90, 100, 1000, 1000000])}
SMALL = {"one": 1, "two": 2, **WORD_VALUE}
FRACTIONS = "half halves third thirds quarter quarters tenth tenths".split()
FRACTION_VALUE = {"half": 0.5, "halves": 0.5, "third": 1 / 3, "thirds": 1 / 3, "quarter": 0.25,
                  "quarters": 0.25, "tenth": 0.1, "tenths": 0.1}

SCI = re.compile(r"(?<![\w.])(?P<sign>[−\-])?(?P<m>\d+(?:\.\d+)?)(?:\s*×\s*10\^|e)(?P<e>[−\-]?\d+)(?![\w.])"
                 r"|(?<![\w.×^])(?P<m2>10)\^(?P<e2>[−\-]\d+)(?![\w.])")
NUM = re.compile(r"(?<![\w.^])(?P<sign>[−\-])?(?P<num>\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?)(?![\w])")
ALNUM = re.compile(r"(?<![\w.\-])(?=[A-Za-z0-9_\-]*\d)(?=[A-Za-z0-9_\-]*[A-Za-z])[A-Za-z0-9][A-Za-z0-9_\-]*[A-Za-z0-9]|(?<![\w.\-])[A-Za-z]\d(?![\w])")
WORD_RE = re.compile(r"(?i)\b(?:(?:one|two|three|four|five|six|seven|eight|nine)-)?(?:" + "|".join(FRACTIONS) + r")\b"
                     r"|\b(?:" + "|".join(WORDS) + r")(?:-(?:one|two|three|four|five|six|seven|eight|nine))?\b")

ABBREV = ["et al.", "Fig.", "Figs.", "e.g.", "i.e.", "vs.", "approx.", "No.", "cf.", "Inc.", "Ltd.", "Dr.",
          "U.S.", "St."]


@dataclass
class Token:
    doc: str
    section: str
    para: int
    sent: int
    index: int
    text: str
    kind: str
    start: int
    end: int
    sentence: str
    table: str = ""
    row: int = -1
    col: int = -1
    exempt: str = ""
    why: str = ""
    value: object = None
    meta: dict = field(default_factory=dict)

    @property
    def key(self) -> str:
        if self.table:
            return f"{self.doc}|{self.table}|r{self.row}|c{self.col}|{self.index}"
        return f"{self.doc}|p{self.para}|s{self.sent}|{self.index}"


def norm(text: str) -> str:
    return text.replace(MINUS, "-")


def split_sentences(paragraph: str) -> list[str]:
    protected = paragraph
    for i, word in enumerate(ABBREV):
        protected = protected.replace(word, word.replace(".", f"\x00{i}\x00"))
    protected = re.sub(r"(?<=\b[A-Z])\.(?=\s?[A-Z]\.)", "\x01", protected)
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z\[(*\"“‘0-9−\-])|(?<=[.!?]\*\*)\s+(?=[A-Z*])", protected)
    out = []
    for part in parts:
        for i, word in enumerate(ABBREV):
            part = part.replace(word.replace(".", f"\x00{i}\x00"), word)
        out.append(part.replace("\x01", "."))
    return [p for p in out if p.strip()]


def blocks(text: str):
    """Yield (section, kind, lines, first line number) for paragraphs and tables."""
    section = ""
    sub = ""
    buf: list[str] = []
    table: list[str] = []
    start = 0
    lines = text.splitlines()
    for no, line in enumerate(lines + [""], start=1):
        heading = re.match(r"^(#{1,4})\s+(.*)", line)
        where = f"{section} > {sub}" if sub else section
        if line.startswith("|"):
            if buf:
                yield (where, "para", buf, start)
                buf = []
            if not table:
                start = no
            table.append(line)
            continue
        if table:
            yield (where, "table", table, start)
            table = []
        if heading or not line.strip():
            if buf:
                yield (where, "para", buf, start)
                buf = []
            if heading:
                level = len(heading.group(1))
                title = heading.group(2).strip()
                if level <= 2 and not section_is_quote(section, title, level):
                    section, sub = title, ""
                else:
                    sub = title
                yield (f"{section} > {sub}" if sub else section, "heading", [line], no)
            continue
        if not buf:
            start = no
        buf.append(line)


def section_is_quote(current: str, title: str, level: int) -> bool:
    return current.startswith("9.") and not re.match(r"^\d+\.", title)


def classify_number(sentence: str, start: int, end: int, text: str) -> str:
    after = sentence[end:end + 3]
    before = sentence[max(0, start - 14):start]
    window = sentence[max(0, start - 30):end + 30]
    if "," in text.replace(MINUS, "").lstrip("-") and re.search(r"\d,\d{3}", text):
        return "thousands"
    if after.startswith("%"):
        kind = "percent"
    elif "." in text:
        kind = "decimal"
    else:
        kind = "integer"
    if re.search(r"\d\s*(?:to|–)\s*$", sentence[max(0, start - 12):start]) or \
            re.match(r"^%?\s*(?:to|–)\s*[−\-]?\d", sentence[end:end + 12]):
        return "range"
    return kind


CITE = re.compile(r"\[(?:\d+(?:\s*[–,-]\s*\d+)*)(?:,\s*\d+(?:\s*[–-]\s*\d+)?)*\]")
CODE = re.compile(r"`[^`]*`")
URL = re.compile(r"https?://\S+|doi:\S+|\b10\.\d{4,}/\S+")
DATE = re.compile(r"\b20\d\d-\d\d-\d\d(?:[ T]\d\d:\d\d(?::\d\d)?)?(?:\s*[+−-]\d{4})?|\b\d\d:\d\d:\d\d\b"
                  r"|\b(?:January|February|March|April|May|June|July|August|September|October|November|December)"
                  r"(?:\s+\d{1,2},?)?\s+20\d\d\b|\b\d{1,2}\s+(?:January|February|March|April|May|June|July|August"
                  r"|September|October|November|December)\s+20\d\d\b")
FILE = re.compile(r"\b[\w\-]+(?:\.[\w\-]+)*\.(?:tsv|csv|txt|json|yaml|yml|md|py|parquet|joblib|npz|gz|tar|xlsx|docx|pdf|png|sh|bundle|R)\b")
VERSION = re.compile(r"(?:Python|NumPy|SciPy|pandas|scikit-learn|matplotlib|TensorFlow|Keras|PyTorch|R|openpyxl|PyYAML|statsmodels|joblib|GENCODE|HGNC|Ensembl|release|version|v)\s*v?(\d+(?:\.\d+)*)")
FIGREF = re.compile(r"(?:Fig\.|Figs\.|Figure|Figures|Table|Tables|Additional file|Additional files|Section|Sections|§|Step|step|Part|line|lines|row|rows|column|columns|criterion|criteria|stage|stages|Stage|Confirmation|Supplementary Table|Supplementary Figure)\s+S?\d")


def span_exempt(sentence: str, start: int, end: int, in_refs: bool) -> tuple[str, str]:
    if in_refs:
        return "REF", "reference list"
    for pattern, code, why in ((CODE, "ID", "code span: file, script or path"), (URL, "ID", "URL or DOI"),
                               (DATE, "DATE", "date or time")):
        for match in pattern.finditer(sentence):
            if match.start() <= start and end <= match.end():
                return code, why
    for match in CITE.finditer(sentence):
        if match.start() <= start and end <= match.end():
            return "REF", "citation number"
    for match in FILE.finditer(sentence):
        if match.start() <= start and end <= match.end():
            return "ID", "file name"
    return "", ""


def number_exempt(sentence: str, start: int, end: int, text: str) -> tuple[str, str]:
    before = sentence[max(0, start - 40):start]
    after = sentence[end:end + 20]
    bare = text.lstrip("-−")
    if re.search(r"top-$", before) or re.search(r"top-$", sentence[max(0, start - 4):start]):
        return "LEVEL", "the 1 of top-1 (or top-3)"
    if bare == "95" and re.match(r"^%\s*(?:CI|confidence|percentile|interval|patient|Wald|lower|upper|prediction|bootstrap|two-sided)", after):
        return "LEVEL", "interval level"
    if bare == "95" and re.search(r"(?:one|two)-sided\s+$", before) and re.match(r"^%", after):
        return "LEVEL", "interval level"
    if re.search(r"(?:Fig\.|Figs\.|Figure|Figures|Table|Tables|Additional files?|Sections?|§|Supplementary (?:Table|Figure))\s+(?:S?\d+[a-z]?(?:\s*(?:,|and|–|-|to)\s*)?)*S?$", before):
        return "NAME", "figure, table, file or section number"
    if re.search(r"(?:Fig\.|Figure)\s+$", before) or re.search(r"(?:Fig\.|Figure|Table)\s+S?$", before):
        return "NAME", "figure, table, file or section number"
    if re.search(r"(?:Confirmation|[Cc]onfirmation cohort|Stage|Step)\s+$", before) and "." not in bare:
        return "NAME", "stage name"
    if re.search(r"Ubuntu\s+$", before):
        return "VER", "operating system version"
    if re.search(r"(?:Python|NumPy|SciPy|pandas|scikit-learn|matplotlib|TensorFlow|Keras|openpyxl|PyYAML|statsmodels|joblib|PyMuPDF|R)\s+$", before):
        return "VER", "software version"
    if re.search(r"(?:version|release|v)\s*$", before, re.IGNORECASE) and re.match(r"^(?:\.\d+)*", after):
        return "VER", "version"
    if re.search(r"log$", before) and bare in {"2", "10"}:
        return "FORMULA", "log base"
    if after.startswith("^"):
        return "FORMULA", "base of a power in a formula"
    if re.search(r"\^[−\-]?$", before):
        return "FORMULA", "exponent in a formula"
    if bare == "1" and re.search(r"TPM \+ $", before):
        return "FORMULA", "pseudocount of log2(TPM + 1)"
    if bare == "1" and after.startswith(")") and re.search(r"log2\([^)]*\+ $", sentence[max(0, start - 60):start]):
        return "FORMULA", "pseudocount of a log2 transform"
    if bare == "0" and re.search(r"negative values set to $", before):
        return "FORMULA", "clip of negative values, max(v, 0)"
    if re.match(r"^ [−\-] ρ", after) or re.search(r"\^[vx] [−\-] $", before) or re.search(r"r_g [−\-] $", before):
        return "FORMULA", "symbolic number in a formula"
    if bare == "0" and after.startswith(")") and re.search(r"max\([^)]*, $", before):
        return "FORMULA", "symbolic number in a formula"
    if re.match(r"^-tissue model", after):
        return "NAME", "model name (22-tissue model)"
    if re.fullmatch(r"(?:19|20)\d\d", bare) and after.startswith(")") and re.search(r"\([A-Z][A-Za-z ]+, $", before):
        return "DATE", "publication year"
    if re.fullmatch(r"0\d", bare) and re.search(r"type $", before):
        return "ID", "TCGA sample type code"
    return "", ""


def alnum_exempt(text: str, sentence: str, start: int) -> tuple[str, str]:
    t = text
    before = sentence[max(0, start - 30):start]
    if re.fullmatch(r"(?i)top-\d", t):
        return "LEVEL", "the 1 of top-1 (or top-3)"
    if t == "5th" and sentence[start + 3:start + 15].startswith(" percentile"):
        return "LEVEL", "one-sided bound percentile"
    if re.fullmatch(r"(?i)macro-F1|F1|SHA-256|L[12]|IMvigor210|stage\d+|T[123]|B\d+\w*", t):
        return "NAME", "metric, method, cohort or internal name"
    if re.fullmatch(r"log2-[a-z]+", t):
        return "FORMULA", "log base"
    if re.fullmatch(r"[A-Za-z]+(?:_[A-Za-z0-9]+)+", t):
        return "NAME", "file or internal name"
    if re.fullmatch(r"(?:GSE|GPL|GSM|GDS|phs|PRJNA|PRJEB|SRP|SRR|ERP|EGAS|EGAD|E-MTAB-)\d+[\w.\-]*|PMID\d+|PMC\d+", t):
        return "ID", "accession"
    if re.fullmatch(r"[0-9a-f]{7,40}", t):
        return "ID", "commit hash"
    if re.fullmatch(r"(?:H|AH|AS|AS-|H-)\d+", t):
        return "NAME", "hypothesis number"
    if re.fullmatch(r"S\d+[a-z]?", t) and re.search(r"(?:Table|Figure|Tables|Figures|Fig\.|and|,|–)\s*$", before):
        return "NAME", "supplementary table or figure number"
    if re.fullmatch(r"\d+[a-z]", t) and re.search(r"(?:Fig\.|Figure|Figs\.|and|,|–)\s*$", before):
        return "NAME", "figure panel"
    if re.fullmatch(r"log2|log10|log1p|log2p|ln2", t):
        return "FORMULA", "log base"
    if re.fullmatch(r"v\d+(?:\.\d+)*", t):
        return "VER", "version"
    if re.fullmatch(r"(?:MET500|POG570|SU2C|TCGA\w*|GTEx\w*|U133A?|U133|HG-U133\w*|A549|HER2|ER|PR|CD\d+\w*|HLA-\w+|AH\d|AS\d)", t):
        return "NAME", "cohort, assay or marker name"
    if re.fullmatch(r"[A-Z][A-Z0-9]*\d[A-Z0-9\-]*|[A-Z]+\d+[A-Z]*\d*-\d+", t):
        return "NAME", "gene symbol or other name"
    if re.fullmatch(r"(?:BASE|SA|SC|LD|NC|M1|V0|PURE|MIX|IF|K|G|Z)(?:-[A-Za-z0-9]+)+", t) or re.fullmatch(r"[A-Z]+_[A-Z0-9_]+", t):
        return "NAME", "model or internal name"
    if re.fullmatch(r"\d+(?:st|nd|rd|th)", t):
        return "NAME", "ordinal"
    return "", ""


@dataclass
class Document:
    name: str
    path: str
    text: str
    tokens: list[Token]
    sentences: dict
    sections: dict
    tables: dict


def tokens_in(doc: str, section: str, para: int, sent: int, sentence: str, in_refs: bool,
              table: str = "", row: int = -1, col: int = -1) -> list[Token]:
    found: list[tuple[int, int, str, str, object, dict]] = []
    taken: list[tuple[int, int]] = []

    def free(a: int, b: int) -> bool:
        return all(b <= x or a >= y for x, y in taken)

    for m in SCI.finditer(sentence):
        if m.group("m2"):
            value = 10.0 ** int(norm(m.group("e2")))
            meta = {"mantissa": "", "exp": norm(m.group("e2"))}
        else:
            value = float(norm((m.group("sign") or "") + m.group("m"))) * 10 ** int(norm(m.group("e")))
            meta = {"mantissa": m.group("m"), "exp": norm(m.group("e"))}
        found.append((m.start(), m.end(), m.group(0), "scientific", value, meta))
        taken.append((m.start(), m.end()))
    for m in ALNUM.finditer(sentence):
        if free(m.start(), m.end()) and not re.fullmatch(r"[\d\-]+", m.group(0)):
            if re.fullmatch(r"\d+(?:-\d+)*-[A-Za-z][\w\-]*", m.group(0)):
                continue
            found.append((m.start(), m.end(), m.group(0), "alnum", None, {}))
            taken.append((m.start(), m.end()))
    for m in NUM.finditer(sentence):
        sign = m.group("sign") or ""
        s = m.start()
        if sign == "-" and s > 0 and sentence[s - 1].isalnum():
            sign = ""
            s += 1
        if sign and s > 0 and (sentence[s - 1].isalnum() or sentence[s - 1] in ")]"):
            sign = ""
            s += 1
        if not free(s, m.end()):
            continue
        text = sign + m.group("num")
        value = float(norm(text).replace(",", ""))
        kind = classify_number(sentence, s, m.end(), text)
        found.append((s, m.end(), text, kind, value, {}))
        taken.append((s, m.end()))
    for m in WORD_RE.finditer(sentence):
        if not free(m.start(), m.end()):
            continue
        word = m.group(0).lower()
        if word == "third" and (sentence[m.end():].startswith(" party")
                                or (m.start() == 0 and sentence[m.end():].startswith(","))):
            continue
        if word == "half" and sentence[m.end():].startswith(("-open", "-closed")):
            continue
        if any(word.endswith(f) for f in FRACTIONS):
            parts = word.split("-")
            value = FRACTION_VALUE[parts[-1]] * (SMALL.get(parts[0], 1) if len(parts) > 1 else 1)
            kind = "fraction"
        else:
            parts = word.split("-")
            value = sum(SMALL[p] for p in parts)
            kind = "word"
        found.append((m.start(), m.end(), m.group(0), kind, value, {}))
        taken.append((m.start(), m.end()))
    found.sort()
    out = []
    for i, (a, b, text, kind, value, meta) in enumerate(found):
        code, why = span_exempt(sentence, a, b, in_refs)
        if not code and kind in {"integer", "decimal", "percent", "range", "thousands", "scientific"}:
            code, why = number_exempt(sentence, a, b, text)
        if not code and kind == "alnum":
            code, why = alnum_exempt(text, sentence, a)
        if not code and kind == "word" and text.lower() == "million" and sentence[max(0, a - 4):a] == "per ":
            code, why = "NAME", "unit name (per million)"
        if not code and kind in {"alnum", "integer"} and section.endswith(("Acknowledgements", "Funding")) \
                and (kind == "alnum" or re.search(r"projects? [^)]*$", sentence[:a])):
            code, why = "ID", "grant or award number"
        out.append(Token(doc, section, para, sent, i, text, kind, a, b, sentence, table, row, col,
                         code, why, value, meta))
    return out


def split_row(line: str) -> list[str]:
    body = line.strip()
    if body.startswith("|"):
        body = body[1:]
    if body.endswith("|"):
        body = body[:-1]
    return [cell.strip() for cell in body.split("|")]


def load(name: str, text: str | None = None) -> Document:
    path = DOCS[name]
    text = (ROOT / path).read_text() if text is None else text
    tokens: list[Token] = []
    sentences: dict = {}
    sections: dict = {}
    tables: dict = {}
    para = 0
    last_title = ""
    for section, kind, lines, first in blocks(text):
        in_refs = section.split(" > ")[0] == "References"
        if kind == "heading":
            heading = lines[0]
            title = re.sub(r"^#+\s*", "", heading)
            last_title = title
            para += 1
            sentences[(para, 0)] = heading
            sections[para] = section
            for token in tokens_in(name, section, para, 0, heading, False):
                if not token.exempt and token.kind in {"integer", "decimal"} and \
                        re.match(r"^#+\s*$", heading[:token.start]):
                    token.exempt, token.why = "NAME", "section number"
                tokens.append(token)
            continue
        if kind == "table":
            para += 1
            sections[para] = section
            caption = last_title
            tid = f"{section}::{caption}::{first}"
            rows = [split_row(l) for l in lines]
            tables[tid] = {"section": section, "caption": caption, "first_line": first, "rows": rows, "para": para}
            for r, cells in enumerate(rows):
                if r == 1 and all(re.fullmatch(r":?-{3,}:?", c) for c in cells if c):
                    continue
                for c, cell in enumerate(cells):
                    tokens += tokens_in(name, section, para, 0, cell, in_refs, tid, r, c)
            continue
        para += 1
        sections[para] = section
        paragraph = " ".join(lines)
        bold = re.match(r"^\*\*(Table \d+\.[^*]*|Figure \d+\.[^*]*)\*\*\s*$", paragraph)
        if bold:
            last_title = bold.group(1)
        for s, sentence in enumerate(split_sentences(paragraph), start=1):
            sentences[(para, s)] = sentence
            tokens += tokens_in(name, section, para, s, sentence, in_refs)
    return Document(name, path, text, tokens, sentences, sections, tables)


def load_all(texts: dict | None = None) -> dict[str, Document]:
    texts = texts or {}
    return {name: load(name, texts.get(name)) for name in DOCS}
