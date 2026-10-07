#!/usr/bin/env python3
"""Confirmation tables. The label path is the only truth input.

This stage runs it on a fake label file. Point --labels at the real label file
to repeat the same tables. Fake organs are drawn from the stage-1 cohort
marginal counts mapped by the basic organ table. Seed is 1. The real per-sample
organ column is not an input to that draw.
"""
from __future__ import annotations

import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"
os.environ["NUMBA_NUM_THREADS"] = "16"

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import yaml
from scipy.stats import binomtest

from loader import read_pog570_table
from metrics import ORGANS, decompose, macro_f1, percentile_ci, topk_hit
from stage3_rules import FAKE_ORGAN_WEIGHTS, native_organs
from mapping_rules import map_site
from util import SEED, load_paths

N_BOOT = 2000
COMPARISONS = [
    ("SA-Z", "BASE-Z"),
    ("SC-Z", "BASE-Z"),
    ("LD-Z", "BASE-Z"),
    ("NC-Z", "BASE-Z"),
    ("M1-Z", "BASE-Z"),
    ("SA-K", "BASE-K"),
    ("V0-K", "BASE-K"),
    ("BASE-K", "BASE-Z"),
]
PRIMARY = ("SA-Z", "BASE-Z")
SITE_KEEP = ("liver", "lymph_node", "lung", "soft_tissue")


def write_fake(path: Path) -> None:
    s1 = read_pog570_table("s1")
    organs = sorted(FAKE_ORGAN_WEIGHTS)
    weights = np.array([FAKE_ORGAN_WEIGHTS[organ] for organ in organs], dtype=np.float64)
    weights = weights / weights.sum()
    rng = np.random.default_rng(1)
    drawn = rng.choice(np.array(organs, dtype=object), size=len(s1), p=weights)
    frame = pd.DataFrame({
        "PATIENT_ID": s1["PATIENT_ID"].astype(str).to_numpy(),
        "organ": drawn,
        "rule": "fake",
    })
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)
    side = {
        "seed": 1,
        "n": int(len(frame)),
        "source": "stage-1 ANALYSIS_COHORT marginal counts mapped by the basic organ table",
        "weights": FAKE_ORGAN_WEIGHTS,
        "order": "Table S1 row order",
    }
    path.with_suffix(".json").write_text(json.dumps(side, indent=2))


def site_group(site: str) -> str:
    return site if site in SITE_KEEP else "remainder"


def conformal_group(site: str) -> str:
    if site == "liver":
        return "liver"
    if site == "lymph_node":
        return "lymph_node"
    return "other"


def beta_bin(value: float, cut: float) -> str:
    if value != value:
        return "NA"
    if value < 0:
        return "negative"
    if value == 0:
        return "0"
    if value < cut:
        return f"(0, {cut})"
    return f">={cut}"


def tc_bin(value: float, bins: list[dict]) -> str:
    for row in bins:
        left = float(row["left"])
        right = float(row["right"])
        closed = row["closed"]
        if closed == "right":
            ok = (value > left) and (value <= right)
        elif closed == "left":
            ok = (value >= left) and (value < right)
        elif closed == "both":
            ok = (value >= left) and (value <= right)
        else:
            ok = (value > left) and (value < right)
        if ok:
            return row["tertile"]
    return "NA"


def load_predictions(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path)


def method_arrays(frame: pd.DataFrame, method: str) -> tuple[np.ndarray, np.ndarray]:
    pred = frame[f"{method}__pred"].to_numpy(dtype=object)
    cols = [f"{method}__p_{organ}" for organ in ORGANS]
    proba = frame[cols].to_numpy(dtype=np.float64)
    return pred, proba


def metric_row(truth, pred, proba, natives) -> dict:
    stats = decompose(truth, pred, natives)
    native_truth = np.array([bool(natives[i]) and truth[i] in natives[i] for i in range(len(truth))])
    valid = pred != "NA"
    top3 = float(np.mean(topk_hit(proba, ORGANS, truth, 3, valid))) if len(truth) else None
    native_top1 = float(np.mean(pred[native_truth] == truth[native_truth])) if native_truth.any() else None
    return {
        "n": int(len(truth)),
        "top1": float(np.mean(pred == truth)) if len(truth) else None,
        "top3": top3,
        "macro_f1": macro_f1(truth, pred) if len(truth) else None,
        "n_at_risk": stats["n_at_risk"],
        "host_rate": stats["host_rate"],
        "n_native_truth": int(native_truth.sum()),
        "native_truth_top1": native_top1,
    }


def subset_metrics(df, method, mask) -> dict:
    part = df.loc[mask]
    pred, proba = method_arrays(part, method)
    truth = part["organ"].to_numpy(dtype=object)
    natives = part["native"].tolist()
    return metric_row(truth, pred, proba, natives)


def one_sided_greater(n_favor: int, n_against: int) -> float:
    n_disc = n_favor + n_against
    if n_disc == 0:
        return 1.0
    return float(binomtest(n_favor, n_disc, 0.5, alternative="greater").pvalue)


def bootstrap_diff(left: np.ndarray, right: np.ndarray, rng: np.random.Generator):
    """Paired bootstrap of left.mean - right.mean. One Generator call of N_BOOT x n."""
    if len(left) == 0:
        return None, None, None, None
    point = float(left.mean() - right.mean())
    draws = rng.integers(0, len(left), size=(N_BOOT, len(left)))
    diffs = left[draws].mean(axis=1) - right[draws].mean(axis=1)
    low, high = percentile_ci(diffs.tolist())
    return point, float(low), float(high), float(np.percentile(diffs, 5))


def conformal_table(df, thresholds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group = df["conformal_group"].to_numpy()
    for method in ("BASE-Z", "SA-Z"):
        _pred, proba = method_arrays(df, method)
        truth = df["organ"].to_numpy(dtype=object)
        for variant in ("global", "mondrian"):
            for alpha in (0.1, 0.2):
                selected = np.zeros((len(df), len(ORGANS)), dtype=bool)
                if variant == "global":
                    thr = thresholds.loc[
                        (thresholds["method"] == method) & (thresholds["variant"] == "global")
                        & (thresholds["alpha"] == alpha) & (thresholds["level"] == "all"),
                        "threshold",
                    ].iloc[0]
                    value = math.inf if thr == "inf" else float(thr)
                    selected[:] = (1.0 - proba) <= value
                else:
                    for level in ("liver", "lymph_node", "other"):
                        thr = thresholds.loc[
                            (thresholds["method"] == method) & (thresholds["variant"] == "mondrian")
                            & (thresholds["alpha"] == alpha) & (thresholds["level"] == level),
                            "threshold",
                        ].iloc[0]
                        value = math.inf if thr == "inf" else float(thr)
                        mask = group == level
                        selected[mask] = (1.0 - proba[mask]) <= value
                sizes = selected.sum(axis=1)
                covered = np.array([
                    bool(selected[i, ORGANS.index(truth[i])]) if truth[i] in ORGANS else False
                    for i in range(len(truth))
                ])
                single = np.empty(len(truth), dtype=object)
                for i in range(len(truth)):
                    if sizes[i] == 1:
                        single[i] = ORGANS[int(np.flatnonzero(selected[i])[0])]

                def emit(level, mask):
                    n = int(mask.sum())
                    sing = (sizes[mask] == 1) if n else np.array([], dtype=bool)
                    acc = None
                    if n and sing.any():
                        acc = float(np.mean(single[mask][sing] == truth[mask][sing]))
                    rows.append({
                        "method": method, "variant": variant, "alpha": alpha, "level": level, "n": n,
                        "coverage": float(np.mean(covered[mask])) if n else None,
                        "mean_set_size": float(np.mean(sizes[mask])) if n else None,
                        "singleton_rate": float(np.mean(sizes[mask] == 1)) if n else None,
                        "singleton_accuracy": acc,
                    })

                emit("all", np.ones(len(truth), dtype=bool))
                for level in ("liver", "lymph_node", "other"):
                    emit(level, group == level)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--make-fake", type=Path, default=None)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.make_fake is not None:
        write_fake(args.make_fake)
    paths = load_paths()
    plan = yaml.safe_load((paths["config"] / "analysis_plan_A3.yaml").read_text())
    cut = float(plan["beta_cut"])
    tc_bins = plan["tc_tertiles"]
    labels = pd.read_csv(args.labels, sep="\t", dtype=str)
    s1 = read_pog570_table("s1")[["PATIENT_ID", "BIOPSY_SITE", "METASTATIC_OR_RECURRENCE", "TUMOUR_CONTENT"]]
    s1["PATIENT_ID"] = s1["PATIENT_ID"].astype(str)
    s1["standard_site"] = [map_site("pog570_biopsy_site", value) for value in s1["BIOPSY_SITE"]]
    s1["TUMOUR_CONTENT_num"] = pd.to_numeric(s1["TUMOUR_CONTENT"], errors="coerce")
    pred = load_predictions(paths["results"] / "stage3" / "pog570_predictions.parquet")
    pred["patient_id"] = pred["patient_id"].astype(str)
    merged = labels.merge(s1, left_on="PATIENT_ID", right_on="PATIENT_ID", how="left")
    merged = merged.merge(pred, left_on="PATIENT_ID", right_on="patient_id", how="left")
    if int(merged["beta"].isna().sum()) == len(merged) and "beta" not in pred.columns:
        raise SystemExit("predictions did not join")
    if merged["patient_id"].isna().any():
        raise SystemExit("a label id is missing from the prediction file")
    merged["native"] = [native_organs(site) for site in merged["standard_site"]]
    merged["in_eval"] = ~merged["organ"].isin(["exclude", "NA"])
    eval_df = merged.loc[merged["in_eval"]].reset_index(drop=True)
    eval_df["at_risk"] = [
        bool(native) and organ not in native
        for organ, native in zip(eval_df["organ"], eval_df["native"])
    ]
    eval_df["in_native"] = [
        bool(native) and organ in native
        for organ, native in zip(eval_df["organ"], eval_df["native"])
    ]
    eval_df["site_group"] = [site_group(site) for site in eval_df["standard_site"]]
    eval_df["conformal_group"] = [conformal_group(site) for site in eval_df["standard_site"]]
    eval_df["tc_bin"] = [tc_bin(value, tc_bins) for value in eval_df["TUMOUR_CONTENT_num"]]
    eval_df["beta_bin"] = [beta_bin(value, cut) for value in eval_df["beta"].to_numpy(dtype=float)]
    args.out.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(SEED)
    risk = eval_df.loc[eval_df["at_risk"]].reset_index(drop=True)
    native_df = eval_df.loc[eval_df["in_native"]].reset_index(drop=True)
    base_h, _ = method_arrays(risk, "BASE-Z")
    sa_h, _ = method_arrays(risk, "SA-Z")
    base_host = np.array([pred in native for pred, native in zip(base_h, risk["native"])])
    sa_host = np.array([pred in native for pred, native in zip(sa_h, risk["native"])])
    n_base_only_h = int(np.sum(base_host & ~sa_host))
    n_sa_only_h = int(np.sum(sa_host & ~base_host))
    p_h1 = one_sided_greater(n_base_only_h, n_sa_only_h)
    h1_diff, h1_low, h1_high, _h1_low5 = bootstrap_diff(sa_host.astype(float), base_host.astype(float), rng)

    base_ok = method_arrays(eval_df, "BASE-Z")[0] == eval_df["organ"].to_numpy()
    sa_ok = method_arrays(eval_df, "SA-Z")[0] == eval_df["organ"].to_numpy()
    n_base_only_a = int(np.sum(base_ok & ~sa_ok))
    n_sa_only_a = int(np.sum(~base_ok & sa_ok))
    p_h2 = one_sided_greater(n_sa_only_a, n_base_only_a)
    h2_diff, h2_low, h2_high, _h2_low5 = bootstrap_diff(sa_ok.astype(float), base_ok.astype(float), rng)

    base_n = method_arrays(native_df, "BASE-Z")[0] == native_df["organ"].to_numpy()
    sa_n = method_arrays(native_df, "SA-Z")[0] == native_df["organ"].to_numpy()
    h3_diff, h3_low, h3_high, h3_onesided_low = bootstrap_diff(sa_n.astype(float), base_n.astype(float), rng)
    h1_sig = bool(p_h1 < 0.05)
    h2_sig = bool(p_h2 < 0.05)
    rows = [
        {
            "hypothesis": "H1", "comparison": "SA-Z - BASE-Z", "set": "at_risk", "n": int(len(risk)),
            "diff": h1_diff, "ci_low": h1_low, "ci_high": h1_high,
            "n_favor": n_base_only_h, "n_against": n_sa_only_h, "p_one_sided": p_h1,
            "tested": True, "significant": h1_sig, "onesided_low": None,
        },
        {
            "hypothesis": "H2", "comparison": "SA-Z - BASE-Z", "set": "eval", "n": int(len(eval_df)),
            "diff": h2_diff, "ci_low": h2_low, "ci_high": h2_high,
            "n_favor": n_sa_only_a, "n_against": n_base_only_a, "p_one_sided": p_h2,
            "tested": h1_sig, "significant": bool(h1_sig and h2_sig), "onesided_low": None,
        },
        {
            "hypothesis": "H3", "comparison": "SA-Z - BASE-Z", "set": "native_truth", "n": int(len(native_df)),
            "diff": h3_diff, "ci_low": h3_low, "ci_high": h3_high,
            "n_favor": None, "n_against": None, "p_one_sided": None,
            "tested": bool(h1_sig and h2_sig),
            "significant": bool(h1_sig and h2_sig and h3_onesided_low is not None and h3_onesided_low > -0.10),
            "onesided_low": h3_onesided_low,
        },
    ]
    pd.DataFrame(rows).to_csv(args.out / "primary.tsv", sep="\t", index=False)

    overall = []
    for method, _base in COMPARISONS:
        row = subset_metrics(eval_df, method, np.ones(len(eval_df), dtype=bool))
        row["method"] = method
        overall.append(row)
    # BASE-Z is the reference and is not in the left column of every pair. Add it once.
    base_row = subset_metrics(eval_df, "BASE-Z", np.ones(len(eval_df), dtype=bool))
    base_row["method"] = "BASE-Z"
    base_k = subset_metrics(eval_df, "BASE-K", np.ones(len(eval_df), dtype=bool))
    base_k["method"] = "BASE-K"
    pd.DataFrame([base_row, base_k, *overall]).drop_duplicates("method").to_csv(
        args.out / "secondary_overall.tsv", sep="\t", index=False
    )

    site_rows = []
    methods = ["BASE-Z", "BASE-K"] + [method for method, _base in COMPARISONS]
    methods = list(dict.fromkeys(methods))
    for method in methods:
        for level, sub in eval_df.groupby("site_group"):
            row = subset_metrics(eval_df, method, eval_df["site_group"].to_numpy() == level)
            row["method"] = method
            row["site_group"] = level
            site_rows.append(row)
    pd.DataFrame(site_rows).to_csv(args.out / "secondary_by_site.tsv", sep="\t", index=False)

    met = eval_df["METASTATIC_OR_RECURRENCE"].to_numpy() == "Metastatic"
    met_rows = []
    for method in methods:
        row = subset_metrics(eval_df, method, met)
        row["method"] = method
        met_rows.append(row)
    pd.DataFrame(met_rows).to_csv(args.out / "secondary_metastatic.tsv", sep="\t", index=False)

    tc_rows = []
    for method in methods:
        for level in [row["tertile"] for row in tc_bins] + ["NA"]:
            mask = eval_df["tc_bin"].to_numpy() == level
            if not mask.any() and level == "NA":
                continue
            row = subset_metrics(eval_df, method, mask)
            row["method"] = method
            row["tc_bin"] = level
            tc_rows.append(row)
    pd.DataFrame(tc_rows).to_csv(args.out / "secondary_tc.tsv", sep="\t", index=False)

    beta_rows = []
    for method in ("BASE-Z", "SA-Z"):
        for level in ("0", f"(0, {cut})", f">={cut}", "NA", "negative"):
            mask = eval_df["beta_bin"].to_numpy() == level
            if not mask.any() and level == "negative":
                continue
            row = subset_metrics(eval_df, method, mask)
            row["method"] = method
            row["beta_bin"] = level
            beta_rows.append(row)
    pd.DataFrame(beta_rows).to_csv(args.out / "secondary_beta.tsv", sep="\t", index=False)

    thresholds = pd.read_csv(paths["config"] / "conformal_thresholds_A3.tsv", sep="\t", dtype=str)
    thresholds["alpha"] = thresholds["alpha"].astype(float)
    conformal_table(eval_df, thresholds).to_csv(args.out / "conformal.tsv", sep="\t", index=False)
    print("confirm", "n_eval", len(eval_df), "n_risk", int(eval_df["at_risk"].sum()), "n_native", int(eval_df["in_native"].sum()), flush=True)


if __name__ == "__main__":
    main()
