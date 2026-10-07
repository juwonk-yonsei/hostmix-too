"""Shared pieces of the A17 number check: sources, values, display, status, keys and table specs.

A value is always read from a named source file by a pandas query and a column. TSV files are
read twice, typed for the query and as strings for the value, so that rounding works on the
stored decimal text (Decimal, ROUND_HALF_UP). Nothing here searches a file for a value.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from decimal import Decimal, ROUND_HALF_UP, ROUND_FLOOR
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
CHECKS = ROOT / "manuscript/checks"
MINUS = "\u2212"

TRANSLATE = {
    "계산": "computed",
    "마이크로어레이": "microarray",
    "사전 지정 서술, 가설 검정 아님": "pre-specified descriptive, not a hypothesis test",
    "사후": "post hoc",
    "확인 검정에 없음": "not in the confirmatory test",
}


class SourceError(Exception):
    pass


class Sources:
    """Cache of source tables. Overrides replace one stored cell in test mode only."""

    def __init__(self) -> None:
        self.typed: dict[str, pd.DataFrame] = {}
        self.text: dict[str, pd.DataFrame] = {}
        self.overrides: dict[tuple[str, int, str], str] = {}
        self.used: set[str] = set()

    def load(self, path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        if path not in self.typed:
            full = ROOT / path
            if not full.exists():
                raise SourceError(f"missing source {path}")
            if path.endswith(".json"):
                data = json.loads(full.read_text())
                typed = pd.json_normalize(data if isinstance(data, list) else [data])
                text = typed.astype(str)
            elif path.endswith(".parquet"):
                import pyarrow.parquet as pq
                n = pq.ParquetFile(full).metadata.num_rows
                typed = pd.DataFrame({"row": range(n)})
                text = typed.astype(str)
            elif path.endswith((".yaml", ".yml")):
                data = yaml.safe_load(full.read_text())
                typed = pd.json_normalize([data])
                text = typed.astype(str)
            else:
                typed = pd.read_csv(full, sep="\t", low_memory=False)
                text = pd.read_csv(full, sep="\t", dtype=str, keep_default_na=False, low_memory=False)
            self.typed[path] = typed
            self.text[path] = text
        self.used.add(path)
        return self.typed[path], self.text[path]

    def rows(self, path: str, selector: str, lists: dict | None = None) -> pd.Index:
        typed, _ = self.load(path)
        if not selector or selector == "all":
            return typed.index
        local = {k: v for k, v in (lists or {}).items()}
        try:
            hit = typed.query(selector, local_dict=local, engine="python")
        except Exception as exc:  # noqa: BLE001
            raise SourceError(f"bad selector on {path}: {selector}: {exc}") from exc
        return hit.index

    def cell(self, path: str, index, column: str) -> str:
        _, text = self.load(path)
        item = re.match(r"^(.*)\[(\d+)\]$", column)
        if item and column not in text.columns:
            whole = self.cell(path, index, item.group(1))
            parts = [p.strip() for p in whole.strip().strip("[]").split(",")]
            return parts[int(item.group(2))]
        if column == "__row__":
            return str(list(text.index).index(index) + 1)
        if column == "__line__":
            return str(list(text.index).index(index) + 2)
        if column not in text.columns:
            raise SourceError(f"{path} has no column {column}")
        key = (path, int(index), column)
        if key in self.overrides:
            return self.overrides[key]
        return text.at[index, column]


def dec(text) -> Decimal:
    s = str(text).strip().replace(MINUS, "-")
    if s in {"", "nan", "NaN", "None"}:
        raise SourceError("empty value")
    if s in {"True", "False"}:
        return Decimal(int(s == "True"))
    return Decimal(s)


def median(*values: Decimal) -> Decimal:
    items = sorted(values)
    n = len(items)
    if n % 2:
        return items[n // 2]
    return (items[n // 2 - 1] + items[n // 2]) / 2


def evaluate(src: Sources, path: str, selector: str, column: str, transform: str,
             facts: dict | None = None, lists: dict | None = None) -> tuple[Decimal, str, list]:
    """Return (value, raw text, row indexes) for one fact definition."""
    if transform.startswith("derived:"):
        expr = transform[len("derived:"):].strip()
        names = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", expr)) - {"min", "max", "abs", "median", "Decimal", "round"}
        scope = {"min": min, "max": max, "abs": abs, "median": median, "Decimal": Decimal,
                 "round": lambda x: Decimal(x).quantize(Decimal(1), rounding=ROUND_HALF_UP)}
        for name in names:
            if facts is None or name not in facts:
                raise SourceError(f"derived expression names unknown fact {name}")
            scope[name] = facts[name]
        expr = re.sub(r"(?<![\w.])(\d+(?:\.\d+)?)(?![\w.])", r"Decimal('\1')", expr)
        value = eval(expr, {"__builtins__": {}}, scope)  # noqa: S307 - expression from our own facts file
        return Decimal(value), expr, []
    index = src.rows(path, selector, lists)
    if transform == "count":
        return Decimal(len(index)), str(len(index)), list(index)
    if transform.startswith("nunique:"):
        col = transform.split(":", 1)[1]
        values = {src.cell(path, i, col) for i in index}
        return Decimal(len(values)), str(len(values)), list(index)
    if transform.startswith(("sum:", "mean:")):
        col = transform.split(":", 1)[1]
        values = [dec(src.cell(path, i, col)) for i in index]
        if not values:
            raise SourceError(f"{transform} over no rows: {path} {selector}")
        total = sum(values, Decimal(0))
        value = total if transform.startswith("sum:") else total / Decimal(len(values))
        return value, ";".join(src.cell(path, i, col) for i in index), list(index)
    if len(index) != 1:
        raise SourceError(f"selector returned {len(index)} rows: {path} :: {selector}")
    raw = src.cell(path, index[0], column)
    value = dec(raw)
    if transform == "x100":
        value = value * 100
    elif transform != "none":
        raise SourceError(f"unknown transform {transform}")
    return value, raw, list(index)


def quant(value: Decimal, places: int) -> Decimal:
    return value.quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP)


def show(value: Decimal, places: int, form: str = "num") -> str:
    """Display text with a hyphen for minus; the comparison normalizes the token the same way."""
    if form == "p":
        if value == 0:
            return "0"
        if value < Decimal("0.001"):
            return show(value, 1, "sci")
        exp = math.floor(value.log10())
        return show(value, max(0, 1 - exp), "num")
    if form == "raw":
        return f"{value:f}"
    if form == "sci":
        if value == 0:
            return "0"
        sign = "-" if value < 0 else ""
        mag = abs(value)
        exp = int((mag.log10()).to_integral_value(rounding=ROUND_FLOOR))
        mant = quant(mag.scaleb(-exp), places)
        if mant >= 10:
            exp += 1
            mant = quant(mag.scaleb(-exp), places)
        return f"{sign}{mant} × 10^{'-' if exp < 0 else ''}{abs(exp)}"
    q = quant(value, places)
    if q == 0:
        q = abs(q)
    return f"{q:f}"


def boundary(value: Decimal, places: int, form: str = "num") -> bool:
    """True when the dropped digits are within 1e-6 of a half unit (a rounding boundary)."""
    if form != "num":
        return False
    scaled = abs(value).scaleb(places)
    frac = scaled - scaled.to_integral_value(rounding=ROUND_FLOOR)
    return abs(frac - Decimal("0.5")) < Decimal("1e-6")


def norm_token(text: str) -> str:
    return text.replace(MINUS, "-").replace(",", "").replace("\u2009", " ").strip()


# ---------------------------------------------------------------- status rule (A14 §2.3 + A17 §4.9)

ARM_METHODS = {"BASE-Z", "SA-Z"}


def status_of(path: str, selector: str = "") -> str:
    p = path
    if p.startswith("results/stage10/") or p.startswith("results/stage8/diagnostics/") or p.startswith("results/stage9/"):
        return "post hoc"
    if p in {"results/stage4/confirm/primary.tsv", "results/stage7/confirm/tables/hypothesis_primary.tsv",
             "results/stage7/confirm/tables/hypothesis_secondary.tsv"}:
        return "preregistered"
    if p == "results/stage4/confirm/secondary_overall.tsv":
        return "preregistered" if re.search(r"method == '(BASE-Z|SA-Z)'", selector) else "pre-specified descriptive"
    if p == "results/stage8/external/tables/metrics.tsv" and re.search(r"analysis_cohort == '(POG570|aux_rnaseq)'", selector) \
            and re.search(r"method == '(BASE-Z|SA-Z)'", selector):
        return "preregistered"
    if p.startswith("results/stage4/confirm/") or p.startswith("results/stage7/confirm/tables/") \
            or p.startswith("results/stage8/external/"):
        return "pre-specified descriptive"
    if p.startswith(("results/stage2/", "results/stage5/")) or re.match(r"results/[^/]+\.(tsv|json)$", p):
        return "exploratory"
    if p.startswith(("results/stage3/", "results/stage6/", "results/stage7/planted/", "results/stage7/synthetic")):
        return "design"
    if p.startswith(("config/", "data/", "results/audit/", "results/stage7/sets_only/")):
        return "data description"
    if p.startswith("results/stage7/posthoc/") or p.startswith("results/stage4/posthoc/"):
        return "post hoc"
    if p.startswith("results/figures/") or p.endswith("predictions.parquet"):
        return "post hoc"
    if p == "manuscript/checks/posthoc_model_check.json":
        return "post hoc"
    return "NA"


STATUS_WORDS = {
    "preregistered": re.compile(r"\bpreregist", re.I),
    "pre-specified descriptive": re.compile(r"\bpre-specified\b", re.I),
    "exploratory": re.compile(r"\bexploratory\b", re.I),
    "post hoc": re.compile(r"\bpost hoc\b", re.I),
}

# ---------------------------------------------------------------- key aliases for the consistency check

COHORT_ALIAS = {
    "MET500": r"MET500|development cohort",
    "POG570": r"POG570|Confirmation 1",
    "aux": r"auxiliary|Auxiliary|six (?:further |auxiliary |confirmatory )?RNA-seq|Confirmation 2",
    "IMvigor210": r"IMvigor210|bladder",
    "Anders": r"Anders|breast cancer cohort|triple-negative|breast tumors",
    "DFCI": r"DFCI|melanoma",
    "PRINCE": r"PRINCE|pancreatic",
    "GSE50760": r"GSE50760|colorectal cancers with matched liver",
    "SU2C": r"SU2C|prostate cancer biopsies",
    "AURORA": r"AURORA|GSE209998",
    "microarray": r"microarray|Microarray",
    "GSE41258": r"GSE41258", "GSE14018": r"GSE14018", "GSE71729": r"GSE71729", "GSE74685": r"GSE74685",
    "FHCRC": r"FHCRC",
    "TCGA-train": r"TCGA-train|TCGA primary tumors|training",
    "TCGA-test": r"TCGA-test|held-out|Held-out",
    "TCGA-met": r"TCGA metastatic|metastatic samples",
    "TCGA": r"TCGA",
    "GTEx": r"GTEx",
    "sim": r"[Ss]imulat|ρ|mixtures",
    "CUP-example": r"CUP-AI-Dx authors|authors' (?:example|metastatic)",
}
METHOD_ALIAS = {
    "baseline": r"[Bb]aseline|BASE-Z",
    "HostMix-TOO": r"HostMix-TOO|SA-Z|SA seed",
    "SCOPE": r"SCOPE",
    "CUP-AI-Dx": r"CUP-AI-Dx",
    "site-specific": r"[Ss]ite-specific",
    "normal classes": r"[Nn]ormal(?:-tissue)? class|normal-class",
    "deconvolution": r"[Dd]econvolution",
    "masking": r"[Mm]asking|masked",
    "gene removal": r"[Gg]ene removal|[Rr]emoval of the 20%",
    "gene sets": r"[Gg]ene-set|[Gg]ene set",
    "gated": r"[Gg]ated",
    "22-tissue": r"22-tissue",
    "3-tissue": r"liver, lung and spleen|3-tissue|three tissues",
    "perceptron": r"perceptron",
    "pure C=0.03": r"standardized (?:pure|models? trained on pure|pure-tumor|model was trained on pure)|pure tumors with standardization|[Ss]tandardized, pure|pure tumors only",
    "pure C=0.15": r"standardized (?:pure|models? trained on pure|pure-tumor)|pure tumors with standardization|[Ss]tandardized, pure",
    "mixture model": r"unstandardized|[Uu]nstandardized",
    "leave-one-host-out": r"left out|leave-one-host-out|[Ll]eave-one-host-out|without (?:liver|brain cortex|one pool tissue|that tissue)",
    "variants": r"variant|mixtures per tumor|lower bound of ρ|C between",
}


def key_found(kind: str, key: str, context: str) -> bool:
    table = COHORT_ALIAS if kind == "cohort" else METHOD_ALIAS
    pattern = table.get(key)
    if pattern is None:
        return False
    return re.search(pattern, context) is not None


# ---------------------------------------------------------------- hashes

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


# ---------------------------------------------------------------- templates

NUMBER_RX = r"[−\-]?\d[\d,]*(?:\.\d+)?(?: × 10\^[−\-]?\d+|e[−\-]?\d+)?"


WORD_RX = r"[A-Za-z]+(?:-[a-z]+)?"


def template_regex(template: str) -> tuple[re.Pattern, list[str]]:
    parts = re.split(r"(\{[^{}]+\})", template)
    out = []
    names = []
    for part in parts:
        if part == "{*}":
            out.append(r".*?")
        elif part.startswith("{") and part.endswith("}"):
            name = part[1:-1]
            names.append(name)
            out.append(f"(?P<{name}>{WORD_RX if name.endswith('_w') else NUMBER_RX})")
        else:
            out.append(re.escape(part))
    return re.compile("^" + "".join(out) + "$"), names


def fill(text, variables: dict):
    if isinstance(text, str):
        def rep(match):
            name = match.group(1)
            if name not in variables:
                return match.group(0)
            return str(variables[name])
        return re.sub(r"\{([A-Za-z_][A-Za-z0-9_]*)\}", rep, text)
    if isinstance(text, dict):
        return {k: fill(v, variables) for k, v in text.items()}
    if isinstance(text, list):
        return [fill(v, variables) for v in text]
    return text
