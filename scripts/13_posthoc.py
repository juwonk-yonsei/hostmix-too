#!/usr/bin/env python3
"""Post-hoc tables. Runs after the preregistered confirmation hash commit.

Does not refit models and does not rewrite results/stage4/confirm.
"""
from __future__ import annotations

import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"
os.environ["NUMBA_NUM_THREADS"] = "16"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from decimal import Decimal

import numpy as np
import pandas as pd

from loader import read_pog570_table
from mapping_rules import map_site
from stage3_rules import native_organs
from util import load_paths

GI = ("Stomach", "Esophagus", "Pancreas", "Biliary", "Colorectal")
GI_SET = set(GI)


def load_eval(paths) -> pd.DataFrame:
    labels = pd.read_csv(paths["config"] / "pog570_eval_labels.tsv", sep="\t", dtype=str)
    s1 = read_pog570_table("s1")[["PATIENT_ID", "BIOPSY_SITE"]]
    s1["PATIENT_ID"] = s1["PATIENT_ID"].astype(str)
    s1["standard_site"] = [map_site("pog570_biopsy_site", value) for value in s1["BIOPSY_SITE"]]
    pred = pd.read_parquet(paths["results"] / "stage3" / "pog570_predictions.parquet")
    pred["patient_id"] = pred["patient_id"].astype(str)
    merged = labels.merge(s1, on="PATIENT_ID", how="left")
    merged = merged.merge(pred, left_on="PATIENT_ID", right_on="patient_id", how="left")
    merged = merged.loc[~merged["organ"].isin(["exclude", "NA"])].reset_index(drop=True)
    merged["native"] = [native_organs(site) for site in merged["standard_site"]]
    return merged


def pred_of(frame: pd.DataFrame, method: str) -> np.ndarray:
    return frame[f"{method}__pred"].to_numpy(dtype=object)


def off_diagonal(frame: pd.DataFrame, method: str, k: int = 15) -> pd.DataFrame:
    truth = frame["organ"].to_numpy(dtype=object)
    pred = pred_of(frame, method)
    both = pd.DataFrame({"truth": truth, "pred": pred})
    both = both.loc[both["truth"] != both["pred"]]
    counts = both.value_counts().reset_index(name="n").head(k)
    counts.insert(0, "method", method)
    return counts


def shift_cells(frame: pd.DataFrame, k: int = 15) -> pd.DataFrame:
    both = pd.DataFrame({
        "truth": frame["organ"].to_numpy(dtype=object),
        "base": pred_of(frame, "BASE-Z"),
        "new": pred_of(frame, "SA-Z"),
    })
    changed = both.loc[both["base"] != both["new"]]
    return changed.value_counts().reset_index(name="n").head(k)


def accuracy_by_organ(frame: pd.DataFrame, methods: list[str]) -> pd.DataFrame:
    rows = []
    for organ, sub in frame.groupby("organ", sort=True):
        row = {"organ": organ, "n": int(len(sub))}
        for method in methods:
            pred = pred_of(sub, method)
            row[f"{method}_top1"] = float(np.mean(pred == sub["organ"].to_numpy()))
            row[f"{method}_n_correct"] = int(np.sum(pred == sub["organ"].to_numpy()))
        rows.append(row)
    return pd.DataFrame(rows)


def liver_tables(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    liver = frame.loc[frame["standard_site"] == "liver"].reset_index(drop=True)
    by_truth = accuracy_by_organ(liver, ["BASE-Z", "SA-Z"])
    freq_rows = []
    for method in ("BASE-Z", "SA-Z"):
        counts = pd.Series(pred_of(liver, method)).value_counts()
        for label, count in counts.items():
            freq_rows.append({"method": method, "predicted_organ": label, "n": int(count), "n_liver_eval": int(len(liver))})
    return by_truth, pd.DataFrame(freq_rows)


def gi_tables(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    truth = frame["organ"].to_numpy(dtype=object)
    in_gi = np.array([label in GI_SET for label in truth])
    rows = []
    stomach_rows = []
    for method in ("BASE-Z", "SA-Z"):
        pred = pred_of(frame, method)
        wrong = pred != truth
        inside = in_gi & wrong & np.array([label in GI_SET for label in pred])
        outside = in_gi & wrong & np.array([label not in GI_SET for label in pred])
        rows.append({
            "method": method,
            "n_truth_in_gi": int(in_gi.sum()),
            "n_wrong_inside_gi": int(inside.sum()),
            "n_wrong_outside_gi": int(outside.sum()),
            "n_correct_in_gi": int((in_gi & ~wrong).sum()),
        })
        stomach_pred = pred == "Stomach"
        stomach_rows.append({
            "method": method,
            "n_pred_stomach": int(stomach_pred.sum()),
            "n_pred_stomach_and_truth_stomach": int((stomach_pred & (truth == "Stomach")).sum()),
            "scope": "eval",
        })
    return pd.DataFrame(rows), pd.DataFrame(stomach_rows)


def error_decomposition(frame: pd.DataFrame) -> pd.DataFrame:
    truth = frame["organ"].to_numpy(dtype=object)
    rows = []
    for method in ("BASE-Z", "SA-Z"):
        pred = pred_of(frame, method)
        host = np.array([
            pred[i] in frame["native"].iloc[i] and truth[i] not in frame["native"].iloc[i]
            for i in range(len(frame))
        ])
        wrong = pred != truth
        other = wrong & ~host
        gi_internal = other & np.array([truth[i] in GI_SET and pred[i] in GI_SET for i in range(len(frame))])
        rest = other & ~gi_internal
        rows.append({
            "method": method,
            "n": int(len(frame)),
            "n_error": int(wrong.sum()),
            "n_host": int((wrong & host).sum()),
            "n_other": int(other.sum()),
            "n_other_gi_internal": int(gi_internal.sum()),
            "n_other_rest": int(rest.sum()),
        })
    return pd.DataFrame(rows)


def cohort_compare(paths) -> pd.DataFrame:
    met = pd.read_csv(paths["results"] / "stage2" / "met500_overall.tsv", sep="\t", dtype=str)
    paired = pd.read_csv(paths["results"] / "stage2" / "met500_paired.tsv", sep="\t", dtype=str)
    primary = pd.read_csv(paths["results"] / "stage4" / "confirm" / "primary.tsv", sep="\t", dtype=str)
    overall = pd.read_csv(paths["results"] / "stage4" / "confirm" / "secondary_overall.tsv", sep="\t", dtype=str)

    def pick(frame, method, representation=None):
        hit = frame["method"] == method
        if representation is not None:
            hit = hit & (frame["representation"] == representation)
        return frame.loc[hit].iloc[0]

    base_m, sa_m = pick(met, "BASE", "Z"), pick(met, "SA", "Z")
    pair = paired.loc[
        (paired["method"] == "SA") & (paired["representation"] == "Z") & (paired["baseline"] == "BASE-Z")
    ].iloc[0]
    base_p, sa_p = pick(overall, "BASE-Z"), pick(overall, "SA-Z")
    h1 = primary.loc[primary["hypothesis"] == "H1"].iloc[0]
    h2 = primary.loc[primary["hypothesis"] == "H2"].iloc[0]
    rows = [
        {
            "cohort": "MET500",
            "n_eval": base_m["n"],
            "n_at_risk": base_m["n_at_risk"],
            "n_native_truth": base_m["n_native_truth"],
            "BASE_Z_top1": base_m["top1"],
            "SA_Z_top1": sa_m["top1"],
            "BASE_Z_host_rate": base_m["host_rate"],
            "SA_Z_host_rate": sa_m["host_rate"],
            "BASE_Z_native_top1": base_m["native_truth_top1"],
            "SA_Z_native_top1": sa_m["native_truth_top1"],
            "H1_diff": str(Decimal(sa_m["host_rate"]) - Decimal(base_m["host_rate"])),
            "H2_diff": pair["top1_diff"],
            "H1_diff_source": "recorded host_rate difference",
            "H2_diff_source": "met500_paired.tsv top1_diff",
        },
        {
            "cohort": "POG570",
            "n_eval": base_p["n"],
            "n_at_risk": base_p["n_at_risk"],
            "n_native_truth": base_p["n_native_truth"],
            "BASE_Z_top1": base_p["top1"],
            "SA_Z_top1": sa_p["top1"],
            "BASE_Z_host_rate": base_p["host_rate"],
            "SA_Z_host_rate": sa_p["host_rate"],
            "BASE_Z_native_top1": base_p["native_truth_top1"],
            "SA_Z_native_top1": sa_p["native_truth_top1"],
            "H1_diff": h1["diff"],
            "H2_diff": h2["diff"],
            "H1_diff_source": "primary.tsv",
            "H2_diff_source": "primary.tsv",
        },
    ]
    return pd.DataFrame(rows)


def main() -> None:
    paths = load_paths()
    if not (paths["logs"] / "confirm_output_A4.txt").exists():
        raise SystemExit("confirm output hashes are not recorded")
    frame = load_eval(paths)
    out = paths["results"] / "stage4" / "posthoc"
    out.mkdir(parents=True, exist_ok=True)
    pd.concat([
        off_diagonal(frame, "BASE-Z"),
        off_diagonal(frame, "SA-Z"),
    ], ignore_index=True).to_csv(out / "confusion_offdiag.tsv", sep="\t", index=False)
    shift_cells(frame).to_csv(out / "pred_shift.tsv", sep="\t", index=False)
    accuracy_by_organ(frame, ["BASE-Z", "SA-Z"]).to_csv(out / "accuracy_by_organ_Z.tsv", sep="\t", index=False)
    accuracy_by_organ(frame, ["BASE-K", "SA-K"]).to_csv(out / "accuracy_by_organ_K.tsv", sep="\t", index=False)
    by_truth, freq = liver_tables(frame)
    by_truth.to_csv(out / "liver_by_truth.tsv", sep="\t", index=False)
    freq.to_csv(out / "liver_pred_freq.tsv", sep="\t", index=False)
    gi, stomach = gi_tables(frame)
    gi.to_csv(out / "gi_errors.tsv", sep="\t", index=False)
    stomach.to_csv(out / "stomach_pred.tsv", sep="\t", index=False)
    error_decomposition(frame).to_csv(out / "error_decomposition.tsv", sep="\t", index=False)
    cohort_compare(paths).to_csv(out / "cohort_compare.tsv", sep="\t", index=False)
    print("posthoc", len(frame), flush=True)


if __name__ == "__main__":
    main()
