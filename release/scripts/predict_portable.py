"""Score a B0 rank matrix with the portable weights. Does not load joblib."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.special import softmax

ORGANS = None


def load_model(directory: Path, name: str) -> dict:
    weights = np.load(directory / f"{name}.npz")
    meta = json.loads((directory / f"{name}.json").read_text())
    return {"weights": weights, "meta": meta}


def predict_proba(model: dict, x_b0: np.ndarray) -> np.ndarray:
    x = np.array(x_b0, dtype=np.float32, copy=True)
    weights = model["weights"]
    if "scaler_mean" in weights.files:
        x -= weights["scaler_mean"]
        x /= weights["scaler_scale"]
    decision = x.astype(np.float64) @ weights["coef"].T + weights["intercept"]
    proba = softmax(decision, axis=1)
    classes = model["meta"]["classes"]
    organ_of = model["meta"]["organ_of_class"]
    organ_names = model["meta"]["organs"]
    keep = [i for i, label in enumerate(classes) if not label.startswith("N_")]
    sub = proba[:, keep]
    totals = sub.sum(axis=1)
    valid = totals > 0
    sub = sub.copy()
    sub[valid] = sub[valid] / totals[valid, None]
    organ = np.zeros((len(x), len(organ_names)), dtype=np.float64)
    for column, index in enumerate(keep):
        label = classes[index]
        organ[:, organ_names.index(organ_of[label])] += sub[:, column]
    return organ


def predict(model: dict, x_b0: np.ndarray) -> np.ndarray:
    organ = predict_proba(model, x_b0)
    names = np.asarray(model["meta"]["organs"], dtype=object)
    return names[np.argmax(organ, axis=1)]
