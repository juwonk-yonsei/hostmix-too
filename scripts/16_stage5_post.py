#!/usr/bin/env python3
"""Stage-5 coefficient, selection-rule, and post-hoc diagnosis tables.

Does not set a thread cap. Does not rank models.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import StratifiedKFold

from mapping_rules import PROJECT_TO_ORGAN
from metrics import ORGANS
from stage3_rules import native_organs
from stage5_lib import (
    GI,
    LABEL_HOST,
    POOL10,
    POOL12,
    SEED,
    effective_coef,
    error_parts,
    load_b0_symbols,
    load_paths,
    metric_row,
    paired_mean,
    stage5_dir,
    variant_names,
)

TAUS = [0.0, 0.005, 0.01, 0.02, 0.0327]


def load_coef_models(paths, root):
    out = {
        "BASE-Z": joblib.load(paths["results"] / "stage2" / "models" / "BASE_Z.joblib"),
        "SA-Z": joblib.load(paths["results"] / "stage2" / "models" / "SA_Z.joblib"),
        "SA-m0": joblib.load(root / "variants" / "models" / "SA-m0.joblib"),
    }
    mats = {}
    classes = {}
    for name, bundle in out.items():
        mats[name], classes[name] = effective_coef(bundle)
    return mats, classes


def host_directions(root: Path):
    tissues = json_tissues(root)
    means = np.load(root / "variants" / "gtex_ref_mean_Z.npy")
    train_mean = np.load(root / "variants" / "train_mean_Z.npy")
    return {tissue: means[i].astype(np.float64) - train_mean.astype(np.float64) for i, tissue in enumerate(tissues)}


def json_tissues(root: Path):
    import json
    return json.loads((root / "variants" / "gtex_ref_mean_tissues.json").read_text())


def alignment(root, mats, classes, directions) -> None:
    rows = []
    tissues = POOL10 + POOL12
    for tissue in tissues:
        s = directions[tissue]
        for model in ("BASE-Z", "SA-m0", "SA-Z"):
            scores = mats[model] @ s
            order = np.argsort(-scores, kind="mergesort")
            top = []
            for rank, idx in enumerate(order[:5], start=1):
                top.append((classes[model][int(idx)], float(scores[int(idx)])))
                rows.append({
                    "tissue": tissue, "model": model, "rank": rank,
                    "learning_label": classes[model][int(idx)], "alignment": float(scores[int(idx)]),
                    "L_h": float(scores.max() - np.median(scores)),
                })
    pd.DataFrame(rows).to_csv(root / "mechanism" / "alignment.tsv", sep="\t", index=False)


def coef_drop(root, mats, classes, directions) -> None:
    symbols = np.array(load_b0_symbols(load_paths()))
    tissues = POOL10 + POOL12
    summary = []
    genes = []
    for label, hosts in LABEL_HOST.items():
        for host in hosts:
            for contrast, left in (("SA-m0 - SA-Z", "SA-m0"), ("BASE-Z - SA-Z", "BASE-Z")):
                d = class_row(mats[left], classes[left], label) - class_row(mats["SA-Z"], classes["SA-Z"], label)
                values = []
                for tissue in tissues:
                    rho = float(spearmanr(d, directions[tissue]).statistic)
                    values.append((tissue, rho))
                paired = dict(values)[host]
                higher = sum(rho > paired for _tissue, rho in values)
                arr = np.array([rho for _tissue, rho in values])
                summary.append({
                    "learning_label": label, "host_tissue": host, "contrast": contrast,
                    "spearman": paired, "rank_of_22": int(1 + higher),
                    "spearman_median_22": float(np.median(arr)),
                    "spearman_max_22": float(np.max(arr)),
                })
                order = np.lexsort((symbols, -d))
                for rank, idx in enumerate(order[:20], start=1):
                    genes.append({
                        "learning_label": label, "host_tissue": host, "contrast": contrast,
                        "rank": rank, "symbol": symbols[int(idx)], "d_g": float(d[int(idx)]),
                        "s_h_g": float(directions[host][int(idx)]),
                    })
    pd.DataFrame(summary).to_csv(root / "mechanism" / "coef_spearman.tsv", sep="\t", index=False)
    pd.DataFrame(genes).to_csv(root / "mechanism" / "coef_top_genes.tsv", sep="\t", index=False)


def class_row(matrix, classes, label) -> np.ndarray:
    if label not in classes:
        raise SystemExit(f"learning label {label} is absent; classes={classes}")
    return matrix[classes.index(label)]


def distances(root) -> None:
    z = np.load(root / "variants" / "train_Z_b0.npy")
    rows = pd.read_csv(root / "variants" / "train_rows.tsv", sep="\t", dtype=str)
    organ = rows["learning_label"].map(PROJECT_TO_ORGAN).to_numpy()
    organs = [name for name in GI]
    centers = {}
    for name in sorted(set(organ.tolist())):
        centers[name] = z[organ == name].mean(axis=0)
    table = []
    for left in organs:
        for right in organs:
            dist = float(np.linalg.norm(centers[left] - centers[right]))
            table.append({"organ_a": left, "organ_b": right, "euclidean": dist, "n_a": int((organ == left).sum()), "n_b": int((organ == right).sum())})
    out = pd.DataFrame(table)
    out.to_csv(root / "mechanism" / "gi_distance.tsv", sep="\t", index=False)
    means = []
    for left in organs:
        vals = [float(np.linalg.norm(centers[left] - centers[right])) for right in organs if right != left]
        means.append({"organ": left, "mean_distance_to_other_four": float(np.mean(vals)), "n": int((organ == left).sum())})
    pd.DataFrame(means).to_csv(root / "mechanism" / "gi_distance_mean.tsv", sep="\t", index=False)


def stomach_and_sarcoma(root) -> None:
    for cohort in ("met500", "pog"):
        frame = pd.read_csv(root / "posthoc" / f"{cohort}_samples.tsv", sep="\t", dtype=str)
        truth = frame["truth"].to_numpy(dtype=object)
        summary = []
        by_truth = []
        for method in STOMACH_MODELS_ALL():
            pred = frame[f"pred__{method}"].to_numpy(dtype=object)
            stomach = pred == "Stomach"
            summary.append({
                "cohort": cohort, "method": method,
                "n_pred_stomach": int(stomach.sum()),
                "n_pred_stomach_true": int((stomach & (truth == "Stomach")).sum()),
            })
            for organ in sorted(set(truth.tolist())):
                by_truth.append({
                    "cohort": cohort, "method": method, "truth": organ,
                    "n_truth": int((truth == organ).sum()),
                    "n_pred_stomach": int(((truth == organ) & stomach).sum()),
                })
        pd.DataFrame(summary).to_csv(root / "mechanism" / f"{cohort}_stomach_summary.tsv", sep="\t", index=False)
        pd.DataFrame(by_truth).to_csv(root / "mechanism" / f"{cohort}_stomach_by_truth.tsv", sep="\t", index=False)
        sarc_rows = []
        for method in SARCOMA_MODELS_ALL():
            pred = frame[f"pred__{method}"].to_numpy(dtype=object)
            mask = truth == "Sarcoma"
            counts = pd.Series(pred[mask]).value_counts()
            for label, count in counts.items():
                sarc_rows.append({
                    "cohort": cohort, "method": method, "n_sarcoma": int(mask.sum()),
                    "predicted_organ": label, "n": int(count),
                })
        pd.DataFrame(sarc_rows).to_csv(root / "mechanism" / f"{cohort}_sarcoma_pred.tsv", sep="\t", index=False)


def STOMACH_MODELS_ALL():
    return ["BASE-Z", "SA-m0", "SA-Z", "SA-r30", "SA-r50", "SA-pool3", "SA-pool22", "SA-MLP"]


def SARCOMA_MODELS_ALL():
    return [
        "BASE-Z", "SA-Z", "SA-LOHO-Brain - Cortex",
        "SA-LOHO-Adipose - Subcutaneous", "SA-LOHO-Muscle - Skeletal",
    ]


def discordant(root) -> None:
    frame = pd.read_csv(root / "posthoc" / "pog_samples.tsv", sep="\t", dtype=str)
    truth = frame["truth"].to_numpy(dtype=object)
    natives = [native_organs(site) for site in frame["site"]]
    native_truth = np.array([bool(native) and truth[i] in native for i, native in enumerate(natives)])
    base_ok = frame["pred__BASE-Z"].to_numpy(dtype=object) == truth
    sa_ok = frame["pred__SA-Z"].to_numpy(dtype=object) == truth
    keep = native_truth & (base_ok != sa_ok)
    sub = frame.loc[keep, [
        "id", "site", "truth", "pred__BASE-Z", "pmax__BASE-Z", "pred__SA-Z", "pmax__SA-Z",
        "beta", "TUMOUR_CONTENT", "METASTATIC_OR_RECURRENCE",
    ]].copy()
    sub["which_correct"] = np.where(base_ok[keep], "BASE-Z", "SA-Z")
    sub.to_csv(root / "posthoc" / "native_discordant.tsv", sep="\t", index=False)


def conformal_tables(root, paths) -> None:
    rows = []
    rows.extend(conformal_met500(root))
    rows.extend(conformal_pog(root, paths))
    pd.DataFrame(rows).to_csv(root / "posthoc" / "conformal_errors.tsv", sep="\t", index=False)


def conformal_met500(root) -> list[dict]:
    frame = pd.read_csv(root / "posthoc" / "met500_samples.tsv", sep="\t", dtype=str)
    truth = frame["truth"].to_numpy(dtype=object)
    sites = frame["site"].to_numpy(dtype=object)
    group = np.array(["liver" if site == "liver" else "lymph_node" if site == "lymph_node" else "other" for site in sites])
    natives = [native_organs(site) for site in sites]
    splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    rows = []
    for method in ("BASE-Z", "SA-Z"):
        proba = np.load(root / "posthoc" / f"met500_{method}_proba.npy")
        pred = frame[f"pred__{method}"].to_numpy(dtype=object)
        selected = np.zeros((len(truth), len(ORGANS)), dtype=bool)
        for train_idx, test_idx in splitter.split(np.zeros(len(truth)), group):
            scores = 1.0 - proba[train_idx, [ORGANS.index(label) for label in truth[train_idx]]]
            threshold = conformal_threshold(scores, 0.1)
            selected[test_idx] = (1.0 - proba[test_idx]) <= threshold
        rows.extend(error_set_rows("MET500", method, "crossfit", truth, pred, natives, selected))
    return rows


def conformal_pog(root, paths) -> list[dict]:
    frame = pd.read_csv(root / "posthoc" / "pog_samples.tsv", sep="\t", dtype=str)
    truth = frame["truth"].to_numpy(dtype=object)
    natives = [native_organs(site) for site in frame["site"]]
    thresholds = pd.read_csv(paths["config"] / "conformal_thresholds_A3.tsv", sep="\t", dtype=str)
    rows = []
    for method in ("BASE-Z", "SA-Z"):
        proba = np.load(root / "posthoc" / f"pog_{method}_proba.npy")
        pred = frame[f"pred__{method}"].to_numpy(dtype=object)
        hit = thresholds.loc[
            (thresholds["method"] == method) & (thresholds["variant"] == "global")
            & (thresholds["alpha"] == "0.1") & (thresholds["level"] == "all")
        ]
        threshold = float(hit["threshold"].iloc[0])
        selected = (1.0 - proba) <= threshold
        rows.extend(error_set_rows("POG570", method, "frozen_global", truth, pred, natives, selected))
    return rows


def conformal_threshold(scores: np.ndarray, alpha: float) -> float:
    n_cal = len(scores)
    k = math.ceil((n_cal + 1) * (1.0 - alpha))
    if k > n_cal:
        return math.inf
    return float(np.sort(scores)[k - 1])


def error_set_rows(cohort, method, source, truth, pred, natives, selected) -> list[dict]:
    parts = error_parts(truth, pred, natives)
    truth_in = selected[np.arange(len(truth)), [ORGANS.index(label) for label in truth]]
    sizes = selected.sum(axis=1)
    rows = []
    for kind, mask in (("host", parts["host"] & parts["wrong"]), ("gi_internal", parts["gi"]), ("rest", parts["rest"])):
        n = int(mask.sum())
        rows.append({
            "cohort": cohort, "method": method, "threshold_source": source, "error_type": kind, "n": n,
            "truth_in_set_rate": float(np.mean(truth_in[mask])) if n else None,
            "mean_set_size": float(np.mean(sizes[mask])) if n else None,
        })
    return rows


def rules(root) -> None:
    for cohort, cluster_col in (("met500", "cluster"), ("pog", "cluster")):
        frame = pd.read_csv(root / "posthoc" / f"{cohort}_samples.tsv", sep="\t")
        truth = frame["truth"].to_numpy(dtype=object)
        natives = [native_organs(site) for site in frame["site"].astype(str)]
        clusters = frame[cluster_col].astype(str).to_numpy()
        base = np.load(root / "posthoc" / f"{cohort}_BASE-Z_proba.npy")
        sa = np.load(root / "posthoc" / f"{cohort}_SA-Z_proba.npy")
        base_pred = frame["pred__BASE-Z"].to_numpy(dtype=object)
        sa_pred = frame["pred__SA-Z"].to_numpy(dtype=object)
        beta = pd.to_numeric(frame["beta"], errors="coerce").to_numpy(dtype=float)
        made = {}
        rows = []
        for tau in TAUS:
            name = f"G-beta({tau})"
            use_base = np.isnan(beta) | (beta <= tau)
            proba, pred = switch(base, sa, base_pred, sa_pred, use_base)
            made[name] = pred
            rows.append(rule_metrics(name, truth, pred, proba, natives))
        use_base = np.isnan(beta)
        proba, pred = switch(base, sa, base_pred, sa_pred, use_base)
        made["G-site"] = pred
        rows.append(rule_metrics("G-site", truth, pred, proba, natives))
        avg = (base.astype(np.float64) + sa.astype(np.float64)) / 2.0
        avg_pred = np.asarray(ORGANS, dtype=object)[np.argmax(avg, axis=1)]
        made["AVG"] = avg_pred
        rows.append(rule_metrics("AVG", truth, avg_pred, avg, natives))
        for method in ("BASE-Z", "SA-Z"):
            src = base if method == "BASE-Z" else sa
            src_pred = base_pred if method == "BASE-Z" else sa_pred
            rows.append(rule_metrics(method, truth, src_pred, src, natives))
        pd.DataFrame(rows).to_csv(root / "posthoc" / f"{cohort}_rules.tsv", sep="\t", index=False)
        paired = []
        risk = np.array([bool(native) and truth[i] not in native for i, native in enumerate(natives)])
        for name, pred in made.items():
            for base_name, base_p in (("SA-Z", sa_pred), ("BASE-Z", base_pred)):
                top = paired_mean(pred == truth, base_p == truth, clusters, SEED)
                host_new = error_parts(truth, pred, natives)["host"]
                host_base = error_parts(truth, base_p, natives)["host"]
                host = paired_mean(host_new[risk].astype(float), host_base[risk].astype(float), clusters[risk], SEED)
                paired.append({
                    "comparison": f"{name} - {base_name}", "n": len(truth), "n_at_risk": int(risk.sum()),
                    "top1_diff": top["diff"], "top1_ci_low": top["ci_low"], "top1_ci_high": top["ci_high"],
                    "top1_n_new_only": top["n_new_only"], "top1_n_base_only": top["n_base_only"],
                    "host_diff": host["diff"], "host_ci_low": host["ci_low"], "host_ci_high": host["ci_high"],
                    "host_n_new_only": host["n_new_only"], "host_n_base_only": host["n_base_only"],
                })
        pd.DataFrame(paired).to_csv(root / "posthoc" / f"{cohort}_rules_paired.tsv", sep="\t", index=False)


def switch(base, sa, base_pred, sa_pred, use_base):
    proba = np.where(use_base[:, None], base, sa)
    pred = np.where(use_base, base_pred, sa_pred).astype(object)
    return proba, pred


def rule_metrics(name, truth, pred, proba, natives) -> dict:
    row = metric_row(truth, pred, proba, natives)
    row["rule"] = name
    return row


def main() -> None:
    paths = load_paths()
    root = stage5_dir(paths)
    mats, classes = load_coef_models(paths, root)
    directions = host_directions(root)
    alignment(root, mats, classes, directions)
    coef_drop(root, mats, classes, directions)
    distances(root)
    stomach_and_sarcoma(root)
    discordant(root)
    conformal_tables(root, paths)
    rules(root)
    print("post done", flush=True)


if __name__ == "__main__":
    main()
