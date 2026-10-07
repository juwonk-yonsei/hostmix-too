"""Check displayed numbers against an explicit row and column.

A number is paired only when a sentence rule or a table registry names the
fact. The checker never searches a result file for a matching value.
"""
from __future__ import annotations

import random
import re
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import pandas as pd

ROOT = Path("/home/kangjw/WORK/Cisplatin/2026_NewProj/proj_A_host_too")
CHECKS = ROOT / "manuscript/checks"
TOKEN = re.compile(r"(?<![A-Za-z0-9/])([−\-]?\d+(?:\.\d+)?)(?![A-Za-z0-9])")


def r3(value) -> str:
    return str(Decimal(str(value)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))


def r1(value) -> str:
    return str(Decimal(str(float(value) * 100)).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))


def show(value, transform: str) -> str:
    if transform == "none":
        return str(value)
    if transform == "int":
        return str(int(value))
    if transform == "r3":
        return r3(value)
    if transform == "r1pct":
        return r1(value)
    if transform == "p":
        number = float(value)
        if number == 0:
            return "0"
        if number < 0.001:
            text = f"{number:.1e}".replace("e-0", "e-").replace("e-", " × 10^−")
            return text
        return f"{number:.2g}"
    raise ValueError(transform)


class FactBook:
    def __init__(self) -> None:
        self.rows = []
        self.cache = {}

    def frame(self, path: str) -> pd.DataFrame:
        if path not in self.cache:
            full = ROOT / path
            if path.endswith(".json"):
                self.cache[path] = pd.read_json(full)
            else:
                self.cache[path] = pd.read_csv(full, sep="\t")
        return self.cache[path]

    def add(self, fact_id, path, query, column, transform, status, note="") -> None:
        frame = self.frame(path)
        hit = frame.query(query)
        if len(hit) != 1:
            raise SystemExit(f"{fact_id} selector returned {len(hit)} rows: {query}")
        raw = hit.iloc[0][column]
        self.rows.append({
            "fact_id": fact_id, "source_file": path, "selector": query, "column": column,
            "transform": transform, "status": status, "note": note,
            "raw": raw, "displayed": show(raw, transform),
        })

    def get(self, fact_id: str) -> dict:
        for row in self.rows:
            if row["fact_id"] == fact_id:
                return row
        raise KeyError(fact_id)


def build_facts() -> FactBook:
    book = FactBook()
    primary = "results/stage4/confirm/primary.tsv"
    for hypothesis, status in (("H1", "preregistered"), ("H2", "preregistered"), ("H3", "preregistered")):
        query = f"hypothesis == '{hypothesis}'"
        for column, transform in (
            ("n", "int"), ("diff", "r3"), ("ci_low", "r3"), ("ci_high", "r3"),
        ):
            book.add(f"{hypothesis}_{column}", primary, query, column, transform, status)
        if hypothesis != "H3":
            book.add(f"{hypothesis}_n_favor", primary, query, "n_favor", "int", status)
            book.add(f"{hypothesis}_n_against", primary, query, "n_against", "int", status)
            book.add(f"{hypothesis}_p", primary, query, "p_one_sided", "p", status)
        else:
            book.add(f"{hypothesis}_bound", primary, query, "onesided_low", "r3", status)
    aux = "results/stage7/confirm/tables/hypothesis_primary.tsv"
    for hypothesis in ("AH1", "AH2", "AH3"):
        query = f"hypothesis == '{hypothesis}'"
        book.add(f"{hypothesis}_n", aux, query, "n", "int", "preregistered")
        book.add(f"{hypothesis}_diff", aux, query, "diff_sa_minus_comparator", "r3", "preregistered")
        book.add(f"{hypothesis}_low", aux, query, "ci_low", "r3", "preregistered")
        book.add(f"{hypothesis}_high", aux, query, "ci_high", "r3", "preregistered")
        if hypothesis == "AH3":
            book.add(f"{hypothesis}_bound", aux, query, "onesided_low", "r3", "preregistered")
        else:
            book.add(f"{hypothesis}_p", aux, query, "p", "p", "preregistered")
    sec = "results/stage7/confirm/tables/hypothesis_secondary.tsv"
    for hypothesis in ("AS1", "AS2", "AS3"):
        query = f"hypothesis == '{hypothesis}'"
        for column, name, transform in (
            ("n", "n", "int"), ("diff", "diff", "r3"),
            ("ci_low", "low", "r3"), ("ci_high", "high", "r3"), ("p", "p", "p"),
        ):
            book.add(f"{hypothesis}_{name}", sec, query, column, transform, "preregistered")
    pog = "results/stage4/confirm/secondary_overall.tsv"
    for method, tag in (("BASE-Z", "base"), ("SA-Z", "sa")):
        query = f"method == '{method}'"
        book.add(f"pog_{tag}_top1", pog, query, "top1", "r3", "preregistered", "H2 arm")
        book.add(f"pog_{tag}_host", pog, query, "host_rate", "r3", "preregistered", "H1 arm")
        book.add(f"pog_{tag}_native", pog, query, "native_truth_top1", "r3", "preregistered", "H3 arm")
        book.add(f"pog_{tag}_n", pog, query, "n", "int", "preregistered")
        book.add(f"pog_{tag}_nrisk", pog, query, "n_at_risk", "int", "preregistered")
        book.add(f"pog_{tag}_nnative", pog, query, "n_native_truth", "int", "preregistered")
    ext = "results/stage8/external/tables/metrics.tsv"
    for cohort in ("MET500", "POG570", "aux_rnaseq"):
        for method, tag in (("BASE-Z", "base"), ("SA-Z", "sa"), ("SCOPE", "scope"), ("CUP-AI-Dx", "cup")):
            query = f"subset == 'common_label' and analysis_cohort == '{cohort}' and method == '{method}'"
            prefix = f"ext_{cohort}_{tag}"
            status = "pre-specified descriptive"
            book.add(f"{prefix}_top1", ext, query, "top1", "r3", status)
            book.add(f"{prefix}_host", ext, query, "host_rate", "r3", status)
            book.add(f"{prefix}_n", ext, query, "n", "int", status)
            book.add(f"{prefix}_nrisk", ext, query, "n_at_risk", "int", status)
    abl = "results/stage9/ablation/metrics.tsv"
    for cohort in ("MET500", "POG570", "aux_rnaseq"):
        for method, tag in (
            ("BASE-Z", "base"), ("PURE-Zs", "pure03"), ("PURE-Zs-w", "pure15"),
            ("MIX-Z0", "mix"), ("SA-Z", "sa"),
        ):
            query = f"cohort == '{cohort}' and method == '{method}'"
            prefix = f"abl_{cohort}_{tag}"
            book.add(f"{prefix}_top1", abl, query, "top1", "r3", "post hoc")
            book.add(f"{prefix}_host", abl, query, "host_rate", "r3", "post hoc")
            book.add(f"{prefix}_native", abl, query, "native_truth_top1", "r3", "post hoc")
            book.add(f"{prefix}_n", abl, query, "n", "int", "post hoc")
            book.add(f"{prefix}_nrisk", abl, query, "n_at_risk", "int", "post hoc")
            book.add(f"{prefix}_nnative", abl, query, "n_native_truth", "int", "post hoc")
    liver = "tissue == 'Liver' and rho == 0.6"
    sim = "results/stage5/sim_ext/full.tsv"
    book.add("loho_liver", sim, liver + " and method == 'SA-LOHO-Liver'", "host_rate", "r3", "exploratory")
    book.add("base_liver", sim, liver + " and method == 'BASE-Z'", "host_rate", "r3", "exploratory")
    return book


def documents() -> list[tuple[str, str]]:
    out = []
    for name in ("manuscript/bmc/A_main.md", "manuscript/bmc/cover_letter.md", "manuscript/bmc/Additional_file_1.md"):
        out.append((name, (ROOT / name).read_text()))
    return out


def exempt(file: str, line: str, token: str, start: int) -> str | None:
    if "## References" in file and False:
        return None
    window = line[max(0, start - 24):start + len(token) + 24]
    if re.search(r"top-1", window) and token in {"1", "-1"}:
        return "top-1"
    if "95%" in window and token == "95":
        return "95% CI"
    if re.search(r"\b(H|AH|AS)\d", window) and token in {"1", "2", "3"} and re.search(rf"(H|AH|AS){token}", window):
        return "hypothesis number"
    if re.search(r"(Fig\.|Figure|Table|Additional file)\s+" + re.escape(token.lstrip("−-")), window):
        return "figure or table number"
    if re.search(r"\[\d", line) and token.isdigit() and int(token) <= 67:
        # citation numbers live in brackets; a bare small integer elsewhere is not exempt
        before = line[:start]
        if before.rfind("[") > before.rfind("]"):
            return "citation"
    if re.search(r"Python 3|NumPy|SciPy|pandas|scikit-learn|matplotlib|TensorFlow", line):
        return "software version"
    if re.search(r"20\d\d", token) and re.search(r"October|seed|Accessed", line):
        return "date or seed"
    return None


def pair_table_line(file: str, line: str, book: FactBook) -> list[dict] | None:
    """Pair a markdown table row when its label cells name the source row."""
    if not line.startswith("|") or line.startswith("|---") or "---|" in line:
        return None
    cells = [cell.strip() for cell in line.strip("|").split("|")]
    if not cells or not cells[0]:
        return None
    claims = []

    def take(cell: str, fact_id: str) -> None:
        fact = book.get(fact_id)
        for match in TOKEN.finditer(cell):
            shown = match.group(1)
            expected = fact["displayed"]
            if shown.lstrip("−-") == expected.lstrip("-") or shown == expected.replace("-", "−"):
                claims.append({"fact_id": fact_id, "displayed": shown, "ok": True})
            else:
                claims.append({"fact_id": fact_id, "displayed": shown, "ok": False, "expected": expected})

    head = cells[0]
    if head.startswith("H1 ") or head.startswith("H2 ") or head.startswith("H3 "):
        key = head.split()[0]
        # The row contains n, two rates, the difference and its interval, and p or the bound.
        return None
    if head in {"MET500", "POG570", "Auxiliary RNA-seq"} and len(cells) >= 8 and cells[1] in {
        "Baseline", "Scaled pure, C = 0.03", "Scaled pure, C = 0.15", "Unscaled mixtures", "HostMix-TOO",
    }:
        cohort = {"MET500": "MET500", "POG570": "POG570", "Auxiliary RNA-seq": "aux_rnaseq"}[head]
        tag = {
            "Baseline": "base", "Scaled pure, C = 0.03": "pure03", "Scaled pure, C = 0.15": "pure15",
            "Unscaled mixtures": "mix", "HostMix-TOO": "sa",
        }[cells[1]]
        prefix = f"abl_{cohort}_{tag}"
        for cell, fact_id in (
            (cells[2], f"{prefix}_n"), (cells[3], f"{prefix}_top1"), (cells[4], f"{prefix}_nrisk"),
            (cells[5], f"{prefix}_host"), (cells[6], f"{prefix}_nnative"), (cells[7], f"{prefix}_native"),
        ):
            take(cell, fact_id)
        return claims
    return None


def main() -> None:
    book = build_facts()
    CHECKS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{k: row[k] for k in ("fact_id", "source_file", "selector", "column", "transform", "status", "note")} for row in book.rows]).to_csv(
        CHECKS / "facts.tsv", sep="\t", index=False)
    claims = []
    unpaired = []
    mismatches = []
    in_refs = False
    for file, text in documents():
        in_refs = False
        for lineno, line in enumerate(text.splitlines(), start=1):
            if line.startswith("## References"):
                in_refs = True
            elif line.startswith("## ") and in_refs:
                in_refs = False
            paired = None if in_refs else pair_table_line(file, line, book)
            if paired is not None:
                for item in paired:
                    claims.append({
                        "file": file, "line": lineno, "displayed": item["displayed"],
                        "fact_id": item["fact_id"], "reason": "table",
                        "context": line.strip()[:120],
                    })
                    if not item["ok"]:
                        mismatches.append(item | {"file": file, "line": lineno})
                continue
            for match in TOKEN.finditer(line):
                token = match.group(1)
                kind = "reference list" if in_refs else exempt(file, line, token, match.start())
                if kind:
                    claims.append({
                        "file": file, "line": lineno, "displayed": token, "fact_id": "exempt",
                        "reason": kind, "context": line[max(0, match.start() - 40):match.end() + 40],
                    })
                else:
                    unpaired.append((file, lineno, token, line.strip()[:160]))
    pd.DataFrame(claims).to_csv(CHECKS / "claims_map.tsv", sep="\t", index=False)
    print("facts", len(book.rows))
    print("claims", len(claims))
    print("mismatches", len(mismatches))
    print("unpaired", len(unpaired))
    for row in mismatches[:10]:
        print("MISMATCH", row)
    for row in unpaired[:25]:
        print(f"{row[0]}:{row[1]} {row[2]} | {row[3]}")


if __name__ == "__main__":
    main()
