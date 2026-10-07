"""Bring the figure sheets of Additional file 2 in line with the figures drawn by 57_figures_A15.py.

Sheets are replaced by name; every other sheet is left as it is. The drawn values come from
the drawing tables written by 57 in the same pass, and the S6 sheet joins them by cohort key to
the meta-analysis table so that the counts and variances stay beside the drawn intervals.
"""
from __future__ import annotations

import importlib
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
A13 = importlib.import_module("50_figures_A13")

BOOK = ROOT / "manuscript" / "bmc" / "Additional_file_2.xlsx"
SRC = ROOT / "manuscript" / "bmc" / "figure_source"
META = ROOT / "results" / "stage7" / "confirm" / "tables" / "meta_analysis.tsv"
HANGUL = re.compile(r"[\uac00-\ud7a3]")
FORBIDDEN_COLUMNS = {"sample_id", "patient_id", "PATIENT_ID", "SAMPLE_ID"}
METHOD = {"BASE-Z": "Baseline", "SA-Z": "HostMix-TOO", "SA-pool22": "22-tissue model"}
LAYER = {"RNA-seq": "RNA-seq", "\ub9c8\uc774\ud06c\ub85c\uc5b4\ub808\uc774": "Microarray"}

DESCRIPTION = {
    "Fig2ab_means": "Exploratory simulation. Unweighted means across tissues drawn in Figure 2a and 2b.",
    "Fig2c_leave_one_host_out": "Exploratory simulation. Host-attraction rate at rho = 0.6 for each host-pool tissue, drawn in Figure 2c.",
    "Fig6a_cohorts": "Cohort host-attraction rates drawn in Figure 6a, with the group labels of the figure.",
    "S1_heatmap": "Locked-classifier heatmap source. Three-decimal labels in the figure. Combinations absent from the sheet were not computed (n/c).",
    "S3_esophagus_masking": "Post hoc esophagus-call fractions among biopsies whose true organ is not esophagus. The MET500 linear-deconvolution fraction was not computed (n/c in Figure S3a).",
    "S4_contributions": "Post hoc mean contributions to the esophagus minus bladder logit drawn in Figure S4: the ten largest positive and five most negative genes per classifier. Inputs are rank-normal scores for the baseline and standardized scores for HostMix-TOO.",
    "S5_reconciliation": "Figure S5 bar values: esophagus calls per cohort and classifier, including zero counts, split by the nearer histology centroid.",
    "S6_meta_analysis": "Pre-specified descriptive cohort differences with Wald 95% intervals from the recorded variances, and DerSimonian-Laird layer estimates, drawn in Figure S6. b and c count at-risk biopsies with a host-attraction error under the baseline only and under HostMix-TOO only.",
    "Microarray_layer": "Pre-specified descriptive microarray effects. p values are recorded and are not presented as tests.",
}


def drawn(name: str) -> pd.DataFrame:
    return pd.read_csv(SRC / f"{name}.tsv", sep="\t")


def cell(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if math.isnan(float(value)) else float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if value is None or (isinstance(value, str) and value == ""):
        return None
    if isinstance(value, str) and HANGUL.search(value):
        raise SystemExit(f"untranslated cell: {value}")
    return value


def put(book, name: str, frame: pd.DataFrame, after: str | None = None) -> None:
    bad = FORBIDDEN_COLUMNS.intersection(frame.columns)
    if bad:
        raise SystemExit(f"{name}: identifier columns {sorted(bad)}")
    if name in book.sheetnames:
        index = book.sheetnames.index(name)
        book.remove(book[name])
    else:
        index = book.sheetnames.index(after) + 1
    sheet = book.create_sheet(name, index)
    sheet.append([DESCRIPTION[name]])
    sheet.append([str(column) for column in frame.columns])
    for row in frame.itertuples(index=False):
        sheet.append([cell(value) for value in row])
    print(f"{name}: {len(frame)} rows")


def fig2ab() -> pd.DataFrame:
    table = drawn("Fig2_means").rename(columns={"host_pull_rate": "host_attraction_rate"})
    table["display"] = table["method"].map(METHOD)
    if table["display"].isna().any():
        raise SystemExit("Fig2ab: unmapped method")
    table["pool"] = table["pool"].map({"in": "host pool", "out": "outside the host pool"})
    table = table.sort_values(["pool", "method", "rho"], ascending=[True, True, False])
    return table[["pool", "rho", "method", "display", "n_tissues", "host_attraction_rate"]]


def fig2c() -> pd.DataFrame:
    return drawn("Fig2_loho").rename(columns={"hostmix": "hostmix_too"})


def fig6a() -> pd.DataFrame:
    table = drawn("Fig6a_cohorts").rename(columns={"group": "role", "group_display": "group"})
    return table[["group", "role", "cohort", "method", "host_rate", "n_at_risk", "source"]]


def s3(book) -> pd.DataFrame:
    sheet = book["S3_esophagus_masking"]
    rows = list(sheet.iter_rows(min_row=2, values_only=True))
    table = pd.DataFrame(rows[1:], columns=rows[0])
    table = table.loc[~(table["cohort"].eq("MET500") & table["method"].eq("LD-Z"))].reset_index(drop=True)
    table.loc[len(table)] = {"note": "not computed", "analysis": "MET500", "cohort": "MET500", "method": "LD-Z",
                             "n_truth_not_esophagus": None, "esophagus_absorption": None}
    order = {"MET500": 0, "POG570": 1, "aux_rnaseq": 2, "microarray": 3}
    methods = {"BASE-Z": 0, "SA-Z": 1, "LD-Z": 2, "SCOPE": 3, "CUP-AI-Dx": 4}
    table = table.assign(_a=table["analysis"].map(order), _m=table["method"].map(methods))
    if table[["_a", "_m"]].isna().any().any():
        raise SystemExit("S3: unexpected analysis or method")
    return table.sort_values(["_a", "_m"]).drop(columns=["_a", "_m"])


def s4() -> pd.DataFrame:
    table = drawn("S4_contributions")
    return table[["display", "method", "n_samples", "side", "order", "rank", "gene", "mean_contribution", "input"]]


def s6() -> pd.DataFrame:
    raw = pd.read_csv(META, sep="\t")
    shown = drawn("S6_meta")
    rows = []
    for _, record in raw.iterrows():
        layer = LAYER[record["layer"]]
        if record["status"] == "cohort":
            name = A13.COHORT[record["cohort"]]
            hit = shown.loc[shown["layer"].eq(layer) & shown["row"].eq("cohort") & shown["cohort"].eq(name)]
            if len(hit) != 1 or abs(float(hit["difference"].iloc[0]) - float(record["d"])) > 1e-9 \
                    or int(hit["n_at_risk"].iloc[0]) != int(record["n"]):
                raise SystemExit(f"S6: drawn row does not match {record['cohort']}")
            rows.append({"layer": layer, "row": "cohort", "cohort": name, "cohort_id": record["cohort"],
                         "n_at_risk": int(record["n"]), "b": int(record["b"]), "c": int(record["c"]),
                         "difference": float(record["d"]), "variance": float(record["v"]),
                         "ci_low": float(hit["ci_low"].iloc[0]), "ci_high": float(hit["ci_high"].iloc[0]),
                         "interval": "Wald, recorded variance", "k": None, "I2": None, "tau2": None})
        else:
            hit = shown.loc[shown["layer"].eq(layer) & shown["row"].eq("DerSimonian-Laird")]
            if len(hit) != 1 or abs(float(hit["difference"].iloc[0]) - float(record["estimate"])) > 1e-9 \
                    or int(hit["k"].iloc[0]) != int(record["k"]):
                raise SystemExit(f"S6: drawn pooled row does not match {layer}")
            rows.append({"layer": layer, "row": "DerSimonian-Laird", "cohort": None, "cohort_id": None,
                         "n_at_risk": None, "b": None, "c": None, "difference": float(record["estimate"]),
                         "variance": None, "ci_low": float(record["ci_low"]), "ci_high": float(record["ci_high"]),
                         "interval": "DerSimonian-Laird random effects", "k": int(record["k"]),
                         "I2": float(record["I2"]), "tau2": float(record["tau2"])})
    table = pd.DataFrame(rows)
    if len(table) != len(shown):
        raise SystemExit("S6: row count differs from the drawn rows")
    table["_layer"] = table["layer"].map({"RNA-seq": 0, "Microarray": 1})
    table["_row"] = table["row"].eq("DerSimonian-Laird").astype(int)
    return table.sort_values(["_layer", "_row"], kind="stable").drop(columns=["_layer", "_row"])


def readme(book) -> None:
    old = book["README"]
    current = {row[0]: row[1] for row in old.iter_rows(min_row=2, values_only=True) if row[0]}
    rows = []
    for name in book.sheetnames:
        if name == "README":
            continue
        text = DESCRIPTION.get(name) or current.get(name) or book[name].cell(1, 1).value
        rows.append((name, text))
    index = book.sheetnames.index("README")
    book.remove(old)
    sheet = book.create_sheet("README", index)
    sheet.append(["Sheet", "Status"])
    for row in rows:
        sheet.append(list(row))
    print(f"README: {len(rows)} rows")


def main() -> None:
    book = load_workbook(BOOK)
    put(book, "Fig2ab_means", fig2ab(), after="Fig2_simulation")
    put(book, "Fig2c_leave_one_host_out", fig2c(), after="Fig2ab_means")
    put(book, "Fig6a_cohorts", fig6a())
    put(book, "S1_heatmap", drawn("S1_heatmap"))
    put(book, "S3_esophagus_masking", s3(book))
    put(book, "S4_contributions", s4())
    put(book, "S5_reconciliation", drawn("S5_reconciliation"))
    put(book, "S6_meta_analysis", s6())
    book["Microarray_layer"].cell(1, 1).value = DESCRIPTION["Microarray_layer"]
    readme(book)
    for sheet in book.worksheets:
        for row in sheet.iter_rows(values_only=True):
            for value in row:
                if isinstance(value, str) and HANGUL.search(value):
                    raise SystemExit(f"{sheet.title}: untranslated cell {value}")
    book.save(BOOK)
    print(BOOK, BOOK.stat().st_size, len(book.sheetnames), "sheets")


if __name__ == "__main__":
    main()
