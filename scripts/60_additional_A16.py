"""Bring the ablation sheets of Additional file 2 in line with the A16 revision.

Figure 4 now has three panels and no shares, so Fig4_ablation is replaced by the drawn
values written by 57_figures_A15.py and the RNA-seq shares move to a sheet named for
Table S8. The Table S8 model names follow the main text, and S9_ablation_microarray
holds the stage-10 microarray and missing-gene results. Sheets are replaced by name with
put() and the README of 58_additional2_figures_A15.py; every other sheet is left as it is.
"""
from __future__ import annotations

import importlib
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
B58 = importlib.import_module("58_additional2_figures_A15")

BOOK = B58.BOOK
SCALE_WORDS = re.compile(r"\b(scaled|unscaled|scaling)\b", re.IGNORECASE)

B58.DESCRIPTION.update({
    "Fig4_ablation": "Post hoc ablation. Values drawn in Figure 4: a host-attraction rate in the at-risk sets and "
                     "b top-1 accuracy in the evaluation sets of MET500, POG570, the six auxiliary RNA-seq cohorts "
                     "and the five microarray cohorts (one sample per patient across the microarray layer); c mean "
                     "top-1 loss in MET500 and POG570 under random removal of as many classifier genes as each of "
                     "four microarray platforms lacks (n = number of platforms, ten draws each). Shares are in "
                     "S8_shares and S9_ablation_microarray.",
    "S8_shares": "Post hoc ablation. Shares of the baseline minus HostMix-TOO host-attraction gap attributable to "
                 "the mixtures in the RNA-seq cohorts (Additional file 1: Table S8); reading rule: every share at "
                 "least 0.5. Not drawn in Figure 4.",
})


def fig4(book) -> None:
    frame = B58.drawn("Fig4_ablation")
    frame = frame[["panel", "cohort", "display", "method", "label", "metric", "value", "n", "source"]]
    B58.put(book, "Fig4_ablation", frame)
    if "Fig4_shares" in book.sheetnames:
        book["Fig4_shares"].title = "S8_shares"
    book["S8_shares"].cell(1, 1).value = B58.DESCRIPTION["S8_shares"]


S10 = B58.ROOT / "results" / "stage10" / "ablation_microarray"
SHARE_MODEL = {"mix": "MIX-Z0", "std": "PURE-Zs", "stdw": "PURE-Zs-w"}
GAP = {"accuracy": "gap, baseline minus HostMix-TOO top-1 in the evaluation set",
       "random_removal": "gap, HostMix-TOO minus baseline mean loss under random removal"}
S8_NAMES = {
    "Scaled pure, C = 0.03": "Standardized, pure (C = 0.03)",
    "Scaled pure, C = 0.15": "Standardized, pure (C = 0.15)",
    "Unscaled mixtures": "Unstandardized, mixtures (C = 0.1)",
}

B58.DESCRIPTION["S9_ablation_microarray"] = (
    "Post hoc ablation in the microarray layer and the missing-gene simulation (Additional file 1: Table S9). "
    "part 'microarray': top-1, esophagus calls, host-attraction rate and native-truth top-1 of the five models "
    "per microarray cohort and for the layer ('all'), one sample per patient across the layer, so cohort n add "
    "up to the layer n. part 'missing genes': change in top-1 accuracy of POG570 and MET500 per platform when "
    "the classifier genes absent from that platform were removed (restriction_drop) and over ten random removals "
    "of as many genes (random_mean, random_min, random_max). part 'share' and 'reading': the fraction of "
    "HostMix-TOO's loss relative to the baseline reproduced by each model, the gap it divides, and the category "
    "of the pre-specified reading rule. Models were not retrained.")


def s9() -> pd.DataFrame:
    names = {**B58.A13.COHORT, "all": "Microarray layer", "MET500": "MET500", "POG570": "POG570"}
    metrics = pd.read_csv(S10 / "metrics.tsv", sep="\t")
    if metrics["cohort"].map(names).isna().any():
        raise SystemExit("S9: unmapped cohort")
    for (model, set_name, metric), part in metrics.groupby(["model", "set", "metric"]):
        layer = part.loc[part["cohort"].eq("all"), "n"]
        if len(layer) != 1 or int(part.loc[part["cohort"].ne("all"), "n"].sum()) != int(layer.iloc[0]):
            raise SystemExit(f"S9: cohort n do not add up for {model} {set_name} {metric}")
    rows = [{"part": "microarray", "contrast": None, "model": r.model, "display": B58.A13.label(r.model),
             "set": r.set, "cohort": names[r.cohort], "platform": None, "metric": r.metric,
             "value": float(r.value), "n": int(r.n)} for r in metrics.itertuples()]
    restriction = pd.read_csv(S10 / "restriction.tsv", sep="\t")
    for r in restriction.itertuples():
        for metric in ("restriction_drop", "random_mean", "random_min", "random_max"):
            rows.append({"part": "missing genes", "contrast": None, "model": r.model,
                         "display": B58.A13.label(r.model), "set": "evaluation", "cohort": names[r.cohort],
                         "platform": r.platform, "metric": metric, "value": float(getattr(r, metric)),
                         "n": 10 if metric.startswith("random") else None})
    shares = pd.read_csv(S10 / "shares.tsv", sep="\t")
    for r in shares.itertuples():
        model = SHARE_MODEL[r.name.split("_", 1)[1]]
        rows.append({"part": "share", "contrast": r.contrast, "model": model, "display": B58.A13.label(model),
                     "set": None, "cohort": None, "platform": None, "metric": "share", "value": float(r.value),
                     "n": None})
        rows.append({"part": "share", "contrast": r.contrast, "model": model, "display": B58.A13.label(model),
                     "set": None, "cohort": None, "platform": None,
                     "metric": GAP[r.contrast], "value": float(r.gap), "n": None})
    category = pd.read_csv(S10 / "category.tsv", sep="\t")
    for r in category.itertuples():
        rows.append({"part": "reading", "contrast": r.contrast, "model": None, "display": None, "set": None,
                     "cohort": None, "platform": None, "metric": "category", "value": r.category, "n": None})
    table = pd.DataFrame(rows)
    if table["display"].eq("").any():
        raise SystemExit("S9: unlabeled model")
    return table


def s8_names(book) -> int:
    changed = 0
    for row in book["S8_ablation"].iter_rows():
        for item in row:
            if item.value in S8_NAMES:
                item.value = S8_NAMES[item.value]
                changed += 1
    return changed


POSTHOC_PLANS = {
    "111eee8ff9b29abccf79faadf0a04ea4ae847b5b": "2026-10-02 19:14:53 +0900",
    "f9b2d82c60e41e07256725a3ad379cac8616bfbe": "2026-10-06 19:31:37 +0900",
}


def history() -> None:
    """Append the two post hoc plan rows to the history sheet of Additional file 3 in place."""
    P44 = importlib.import_module("44_release_pack")
    path = B58.ROOT / "manuscript" / "bmc" / "Additional_file_3.xlsx"
    book = load_workbook(path)
    sheet = book["history"]
    present = {row[0] for row in sheet.iter_rows(min_row=3, values_only=True)}
    notes = {full: (stage, role, note) for full, stage, role, note in P44.HISTORY}
    for full, expected in POSTHOC_PLANS.items():
        line = subprocess.check_output(["git", "log", "-1", "--format=%H\t%ci\t%s", full], cwd=B58.ROOT,
                                       text=True).strip()
        commit, when, subject = line.split("\t", 2)
        if commit != full or when != expected:
            raise SystemExit(f"git gives {commit} {when}, expected {full} {expected}")
        stage, role, note = notes[full]
        if full not in present:
            sheet.append([commit, when, stage, role, subject, note])
        print("history", commit, when, role, subject)
    book.save(path)


def check(book) -> None:
    for sheet in book.worksheets:
        for row in sheet.iter_rows(values_only=True):
            for value in row:
                if isinstance(value, str) and B58.HANGUL.search(value):
                    raise SystemExit(f"{sheet.title}: untranslated cell {value}")


def main() -> None:
    book = load_workbook(BOOK)
    fig4(book)
    print("S8_ablation: renamed", s8_names(book), "cells")
    if not B58.DESCRIPTION["S5_reconciliation"].startswith("Post hoc."):
        B58.DESCRIPTION["S5_reconciliation"] = "Post hoc. " + B58.DESCRIPTION["S5_reconciliation"]
    book["S5_reconciliation"].cell(1, 1).value = B58.DESCRIPTION["S5_reconciliation"]
    B58.put(book, "S9_ablation_microarray", s9(), after="S8_ablation")
    B58.readme(book)
    check(book)
    book.save(BOOK)
    left = [(sheet.title, value) for sheet in book.worksheets for row in sheet.iter_rows(values_only=True)
            for value in row if isinstance(value, str) and SCALE_WORDS.search(value)]
    print(BOOK, BOOK.stat().st_size, len(book.sheetnames), "sheets;", len(left), "cells with scaled/unscaled/scaling")
    history()


if __name__ == "__main__":
    main()
