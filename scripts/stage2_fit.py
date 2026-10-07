"""Fitting helpers for stage 2. C selection uses a scaler inside each fold."""
from __future__ import annotations

import numpy as np
from scipy.optimize import nnls
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

from metrics import ORGANS, macro_f1, organ_probability, predict_from_proba
from stage2_common import C_GRID, SEED, n_drop
from util import SEED as _SEED

assert SEED == _SEED


def fit_logit(X: np.ndarray, y: np.ndarray, C: float, scale: bool) -> dict:
    scaler = None
    X_use = X
    if scale:
        scaler = StandardScaler()
        X_use = scaler.fit_transform(X)
    clf = LogisticRegression(C=float(C), solver="lbfgs", max_iter=5000, class_weight="balanced")
    clf.fit(X_use, y)
    n_iter = int(np.max(clf.n_iter_))
    return {
        "clf": clf,
        "scaler": scaler,
        "C": float(C),
        "n_iter": n_iter,
        "max_iter_hit": bool(n_iter >= 5000),
        "classes": [str(c) for c in clf.classes_],
    }


def cv_select_C(X: np.ndarray, y: np.ndarray, groups: np.ndarray) -> dict:
    splitter = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=SEED)
    rows = []
    for fold, (tr, va) in enumerate(splitter.split(X, y, groups)):
        yva = y[va]
        labels = sorted(set(yva.tolist()))
        for C in C_GRID:
            bundle = fit_logit(X[tr], y[tr], C, scale=True)
            Xva = bundle["scaler"].transform(X[va])
            pred = bundle["clf"].predict(Xva)
            score = macro_f1(yva, pred)
            # macro_f1 uses labels present in y_true. Pass the same yva.
            rows.append({
                "C": float(C),
                "fold": int(fold),
                "macro_f1": float(score),
                "n_iter": bundle["n_iter"],
                "n_valid": int(len(va)),
                "labels_in_valid": labels,
            })
            print("cv", "C", C, "fold", fold, "f1", score, "n_iter", bundle["n_iter"], flush=True)
    means = {}
    for C in C_GRID:
        vals = [row["macro_f1"] for row in rows if row["C"] == float(C)]
        means[float(C)] = float(np.mean(vals))
    best = max(means.values())
    chosen = min(C for C, value in means.items() if value >= best - 1e-12)
    return {"folds": rows, "mean_macro_f1": {str(k): v for k, v in means.items()}, "C": float(chosen)}


def assemble_aug(pure: np.ndarray, mixes: np.ndarray) -> np.ndarray:
    """pure is (n, p). mixes is tumor-major (n * 4, p): tumor0 mix0..3, tumor1 mix0..3, ..."""
    n, p = pure.shape
    out = np.empty((n * 5, p), dtype=np.float32)
    out[0::5] = pure
    for j in range(4):
        out[1 + j::5] = mixes[j::4]
    return out


def repeat_labels(values: np.ndarray) -> np.ndarray:
    return np.repeat(values, 5)


def sensitivity(pure: np.ndarray, mixed: np.ndarray, sd: np.ndarray) -> np.ndarray:
    med = np.median(np.abs(mixed - pure), axis=0)
    out = np.full(med.shape, np.inf, dtype=np.float64)
    ok = sd > 0
    out[ok] = med[ok] / sd[ok]
    return out


def rank_desc(scores: np.ndarray) -> np.ndarray:
    """Higher score first. Ties break toward the smaller index."""
    filled = np.where(np.isfinite(scores), scores, np.inf)
    return np.lexsort((np.arange(scores.size), -filled))


def select_k(scores: np.ndarray, q: float) -> tuple[np.ndarray, np.ndarray]:
    order = rank_desc(scores)
    n_remove = n_drop(scores.size, q)
    dropped = order[:n_remove]
    kept = np.sort(order[n_remove:])
    return kept, dropped


def select_z(scores: np.ndarray, variance: np.ndarray, q: float, n_keep: int = 5000) -> tuple[np.ndarray, np.ndarray]:
    order = rank_desc(scores)
    n_remove = n_drop(scores.size, q)
    dropped = order[:n_remove]
    remain = order[n_remove:]
    var = variance[remain]
    pick = np.lexsort((remain, -var))
    chosen = remain[pick[:n_keep]]
    return np.sort(chosen), dropped


def tumor_organ_from_proba(proba: np.ndarray, classes: list[str]):
    """Drop classes whose names start with N_, then sum to organs. Zero remainder -> NA."""
    keep = np.array([not str(c).startswith("N_") for c in classes])
    tumor_classes = [str(c) for c, flag in zip(classes, keep) if flag]
    sub = proba[:, keep]
    totals = sub.sum(axis=1)
    valid = totals > 0
    sub = sub.copy()
    sub[valid] = sub[valid] / totals[valid, None]
    organ = organ_probability(sub, tumor_classes)
    pred = predict_from_proba(organ, ORGANS, valid)
    raw_argmax = np.asarray(classes, dtype=object)[np.argmax(proba, axis=1)]
    return organ, pred, raw_argmax


def predict_bundle(bundle: dict, X: np.ndarray):
    X_use = X if bundle["scaler"] is None else bundle["scaler"].transform(X)
    proba = bundle["clf"].predict_proba(X_use)
    classes = [str(c) for c in bundle["clf"].classes_]
    organ, pred, raw = tumor_organ_from_proba(proba, classes)
    return organ, pred, raw, proba, classes


def ld_references(tumor_norm: np.ndarray, labels: np.ndarray, host_norm: np.ndarray):
    """Return mu (n_labels, genes) and h (genes,), each renormalized to sum 1e6."""
    from stage2_common import sum_1e6

    label_names = sorted(set(labels.tolist()))
    mus = []
    for label in label_names:
        block = tumor_norm[labels == label].mean(axis=0, keepdims=True)
        mus.append(sum_1e6(block)[0])
    h = sum_1e6(host_norm.mean(axis=0, keepdims=True))[0]
    return np.vstack(mus).astype(np.float32), h.astype(np.float32), label_names


def ld_apply(x: np.ndarray, h: np.ndarray, mu: np.ndarray, keep: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Weighted NNLS. x, h are sum-1e6 on G. mu is (n_labels, genes). keep drops MT- genes."""
    from stage2_common import sum_1e6

    A = np.column_stack([h.astype(np.float64), mu.T.astype(np.float64)])
    A_k = A[keep]
    weight = 1.0 / (A_k.mean(axis=1) + 1.0)
    sw = np.sqrt(weight)
    A_w = A_k * sw[:, None]
    n = x.shape[0]
    beta = np.full(n, np.nan, dtype=np.float64)
    cleaned = np.empty((n, x.shape[1]), dtype=np.float32)
    for i in range(n):
        b_w = x[i, keep].astype(np.float64) * sw
        coef, _ = nnls(A_w, b_w)
        c_h = float(coef[0])
        denom = float(coef.sum())
        beta[i] = c_h / denom if denom > 0 else np.nan
        residual = np.maximum(x[i].astype(np.float64) - c_h * h.astype(np.float64), 0.0)
        cleaned[i] = sum_1e6(residual[None, :])[0]
        if (i + 1) % 500 == 0:
            print("ld", i + 1, "/", n, flush=True)
    return cleaned, beta
