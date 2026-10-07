"""Metric definitions locked for stage 1. No cohort labels are inferred here."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import f1_score

from mapping_rules import PROJECT_TO_ORGAN

ORGANS = sorted(set(PROJECT_TO_ORGAN.values()))


def organ_probability(project_proba: np.ndarray, project_classes: list[str]) -> np.ndarray:
    out = np.zeros((project_proba.shape[0], len(ORGANS)), dtype=np.float64)
    for j, project in enumerate(project_classes):
        out[:, ORGANS.index(PROJECT_TO_ORGAN[project])] += project_proba[:, j]
    return out


def predict_from_proba(proba: np.ndarray, classes: list[str], valid: np.ndarray | None = None) -> np.ndarray:
    pred = np.asarray(classes, dtype=object)[np.argmax(proba, axis=1)]
    if valid is not None:
        pred = pred.copy()
        pred[~valid] = "NA"
    return pred


def topk_hit(proba: np.ndarray, classes: list[str], truth: np.ndarray, k: int, valid: np.ndarray | None = None) -> np.ndarray:
    order = np.argsort(-proba, axis=1)[:, :k]
    names = np.asarray(classes, dtype=object)[order]
    hits = np.array([label in row for label, row in zip(truth, names)], dtype=bool)
    if valid is not None:
        hits = hits & valid
    return hits


def macro_f1(truth: np.ndarray, pred: np.ndarray) -> float:
    labels = sorted(set(truth.tolist()))
    if not labels:
        return float("nan")
    return float(f1_score(truth, pred, average="macro", labels=labels, zero_division=0))


def per_class_recall(truth: np.ndarray, pred: np.ndarray) -> list[dict]:
    rows = []
    for label in sorted(set(truth.tolist())):
        mask = truth == label
        rows.append({"label": label, "n": int(mask.sum()), "recall": float(np.mean(pred[mask] == label))})
    return rows


def mask_native(organ_proba: np.ndarray, native_organs: list[list[str]]) -> tuple[np.ndarray, np.ndarray]:
    """Zero native-organ probabilities and renormalize. A zero remainder is an NA prediction."""
    out = organ_proba.copy()
    for i, banned in enumerate(native_organs):
        for organ in banned:
            if organ in ORGANS:
                out[i, ORGANS.index(organ)] = 0.0
    totals = out.sum(axis=1)
    valid = totals > 0
    out[valid] = out[valid] / totals[valid, None]
    pred = predict_from_proba(out, ORGANS, valid)
    return out, pred


def decompose(truth: np.ndarray, pred: np.ndarray, native_sets: list[set[str]]) -> dict:
    truth = np.asarray(truth, dtype=object)
    pred = np.asarray(pred, dtype=object)
    n = len(truth)
    at_risk = np.zeros(n, dtype=bool)
    host = np.zeros(n, dtype=bool)
    for i in range(n):
        native = native_sets[i]
        if not native:
            continue
        if truth[i] not in native:
            at_risk[i] = True
            if pred[i] in native:
                host[i] = True
    error = pred != truth
    n_error = int(error.sum())
    n_host = int((error & host).sum())
    n_risk = int(at_risk.sum())
    return {
        "n": n,
        "n_error": n_error,
        "n_host": n_host,
        "H": (n_host / n_error) if n_error else None,
        "n_at_risk": n_risk,
        "host_rate": (n_host / n_risk) if n_risk else None,
        "n_biopsy_native_truth": int((~at_risk & np.array([bool(s) for s in native_sets])).sum()),
    }


def accuracy_block(truth: np.ndarray, proba: np.ndarray, classes: list[str], pred: np.ndarray, valid: np.ndarray | None = None) -> dict:
    top1 = pred == truth
    top3 = topk_hit(proba, classes, truth, 3, valid)
    return {
        "n": int(len(truth)),
        "top1": float(np.mean(top1)) if len(truth) else None,
        "top3": float(np.mean(top3)) if len(truth) else None,
        "macro_f1": macro_f1(truth, pred) if len(truth) else None,
    }


def percentile_ci(values: list[float]) -> list[float | None]:
    arr = np.asarray([v for v in values if v is not None and np.isfinite(v)], dtype=float)
    if arr.size == 0:
        return [None, None]
    return [float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))]


def bootstrap_metrics(truth, pred, proba, classes, native_sets, n_boot: int = 2000, seed: int = 20261001) -> dict:
    rng = np.random.default_rng(seed)
    n = len(truth)
    draws = rng.integers(0, n, size=(n_boot, n))
    top1, top3, h_values, rate_values = [], [], [], []
    n_h_undefined = 0
    point_top1 = float(np.mean(pred == truth))
    point_top3 = float(np.mean(topk_hit(proba, classes, truth, 3, pred != "NA")))
    point = decompose(truth, pred, native_sets)
    for draw in draws:
        sub_truth = truth[draw]
        sub_pred = pred[draw]
        sub_native = [native_sets[i] for i in draw]
        top1.append(float(np.mean(sub_pred == sub_truth)))
        top3.append(float(np.mean(topk_hit(proba[draw], classes, sub_truth, 3, sub_pred != "NA"))))
        stats = decompose(sub_truth, sub_pred, sub_native)
        if stats["H"] is None:
            n_h_undefined += 1
        else:
            h_values.append(stats["H"])
        if stats["host_rate"] is not None:
            rate_values.append(stats["host_rate"])
    point["top1"] = point_top1
    point["top3"] = point_top3
    point["top1_ci"] = percentile_ci(top1)
    point["top3_ci"] = percentile_ci(top3)
    point["H_ci"] = percentile_ci(h_values)
    point["host_rate_ci"] = percentile_ci(rate_values)
    point["n_bootstrap"] = n_boot
    point["n_bootstrap_H_undefined"] = n_h_undefined
    return point
