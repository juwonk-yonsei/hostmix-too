#!/usr/bin/env python3
"""Freeze POG570 predictions and MET500 conformal thresholds.

Run only after config/prereg_A3.md and config/analysis_plan_A3.yaml exist.
Reads POG570 expression. Does not read POG570 organ labels or histology.
"""
from __future__ import annotations

import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"
os.environ["NUMBA_NUM_THREADS"] = "16"

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import joblib
import numpy as np
import pandas as pd

from genes import probemap_symbol
from loader import load_pog570_expression, read_pog570_table, verify_pog570_lock_hashes
from mapping_rules import MET500_COHORT_TO_ORGAN, map_site
from metrics import ORGANS
from stage2_common import SITE_MODEL_ORDER, load_gene_pack, load_tpm, sum_1e6
from util import load_paths, sha256_file

KEEP_METHODS = [
    ("BASE", "Z"),
    ("SA", "Z"),
    ("SC", "Z"),
    ("LD", "Z"),
    ("NC", "Z"),
    ("M1", "Z"),
    ("BASE", "K"),
    ("SA", "K"),
    ("V0", "K"),
]


def method_id(method: str, rep: str) -> str:
    return f"{method}-{rep}"


def ensembl_base_to_symbol() -> dict[str, str]:
    """Same unversioned map as genes.build_common_genes. Conflicting bases are dropped."""
    probe = probemap_symbol()
    by_base: dict[str, set[str]] = {}
    for ens, sym in zip(probe["id"].astype(str), probe["gene"].astype(str)):
        by_base.setdefault(ens.split(".")[0], set()).add(sym)
    return {base: next(iter(symbols)) for base, symbols in by_base.items() if len(symbols) == 1}


def collapse_to_g(expr: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray, list[str], dict]:
    mapping = ensembl_base_to_symbol()
    bases = pd.Index(expr.index.astype(str)).str.split(".").str[0]
    symbol = np.asarray(bases.map(mapping), dtype=object)
    keep = pd.notna(symbol)
    sub = expr.loc[keep].astype(np.float32)
    grouped = sub.groupby(symbol[keep], sort=False).sum()
    g_index = pd.Index(genes)
    present = grouped.index.intersection(g_index)
    missing = g_index.difference(grouped.index)
    aligned = grouped.reindex(g_index).fillna(0).astype(np.float32)
    # grouped columns are samples. Transpose to samples x genes.
    matrix = np.ascontiguousarray(aligned.to_numpy(dtype=np.float32).T)
    audit = {
        "n_ensembl_rows": int(expr.shape[0]),
        "n_ensembl_mapped": int(keep.sum()),
        "n_symbols_summed": int(grouped.shape[0]),
        "n_G": len(genes),
        "n_G_present": int(len(present)),
        "n_G_filled_0": int(len(missing)),
        "n_samples": int(expr.shape[1]),
    }
    return matrix, [str(c) for c in expr.columns], audit


def load_bundles(model_dir: Path) -> dict:
    names = ["BASE_Z", "BASE_K", "SA_Z", "SA_K", "NC_Z"]
    names += [f"SC_{site}_Z" for site in SITE_MODEL_ORDER]
    # predict_all also asks for K site models and IF models. Load the full set it expects.
    from importlib import import_module
    ev = import_module("09_eval")
    bundles = {}
    for rep in ("Z", "K"):
        for name in ("BASE", "SA", "NC", "IF20"):
            bundles[f"{name}_{rep}"] = joblib.load(model_dir / f"{name}_{rep}.joblib")
        for site in SITE_MODEL_ORDER:
            for prefix in ("SC", "SCIF20", "SCIF40"):
                bundles[f"{prefix}_{site}_{rep}"] = joblib.load(model_dir / f"{prefix}_{site}_{rep}.joblib")
    return bundles, ev


def load_references(path: Path):
    packed = np.load(path, allow_pickle=True)
    h = {key[3:]: packed[key] for key in packed.files if key.startswith("h__")}
    v0 = {key[4:]: packed[key] for key in packed.files if key.startswith("v0__")}
    return packed["mu"], [str(x) for x in packed["mu_labels"].tolist()], h, v0, packed["keep_mt"], [str(x) for x in packed["genes"].tolist()]


def top3_strings(proba: np.ndarray, pred: np.ndarray) -> list[str]:
    order = np.argsort(-proba, axis=1, kind="mergesort")[:, :3]
    names = np.asarray(ORGANS, dtype=object)[order]
    out = []
    for row, label in zip(names, pred):
        if label == "NA":
            out.append("NA")
        else:
            out.append("|".join(str(x) for x in row))
    return out


def conformal_threshold(scores: np.ndarray, alpha: float) -> tuple[float, int]:
    n_cal = len(scores)
    k = math.ceil((n_cal + 1) * (1.0 - alpha))
    if k > n_cal:
        return math.inf, k
    return float(np.sort(scores)[k - 1]), k


def main() -> None:
    paths = load_paths()
    if not (paths["config"] / "prereg_A3.md").exists() or not (paths["config"] / "analysis_plan_A3.yaml").exists():
        raise SystemExit("prereg_A3.md or analysis_plan_A3.yaml is missing. Section 5 was not started.")
    verify_pog570_lock_hashes()
    out = paths["results"] / "stage3"
    genes, b0_idx, offsets, indices, set_names = load_gene_pack(paths)
    mu, _labels, h_by_site, v0, keep, ref_genes = load_references(out / "references.npz")
    if ref_genes != genes:
        raise SystemExit("reference gene order does not match G_symbols.txt")
    print("load expression", flush=True)
    expr = load_pog570_expression(unlock=True)
    matrix, patient_ids, audit = collapse_to_g(expr, genes)
    del expr
    matrix = sum_1e6(matrix)
    (out / "pog_feature_audit.json").write_text(json.dumps(audit, indent=2))
    print("feature", audit, flush=True)

    s1 = read_pog570_table("s1")[["PATIENT_ID", "BIOPSY_SITE"]]
    site_of = dict(zip(s1["PATIENT_ID"].astype(str), s1["BIOPSY_SITE"]))
    missing = [pid for pid in patient_ids if pid not in site_of]
    if missing:
        raise SystemExit(f"expression ids missing from Table S1: {len(missing)}")
    sites = [map_site("pog570_biopsy_site", site_of[pid]) for pid in patient_ids]

    bundles, ev = load_bundles(paths["results"] / "stage2" / "models")
    print("predict POG", flush=True)
    preds, beta, _natives, _raw = ev.predict_all(
        matrix, sites, bundles, b0_idx, offsets, indices, mu, h_by_site, keep, v0
    )
    frame = {"patient_id": patient_ids, "beta": beta}
    freq_rows = []
    for method, rep in KEEP_METHODS:
        proba, pred = preds[(method, rep)]
        mid = method_id(method, rep)
        frame[f"{mid}__pred"] = pred
        frame[f"{mid}__top3"] = top3_strings(proba, pred)
        for j, organ in enumerate(ORGANS):
            frame[f"{mid}__p_{organ}"] = proba[:, j]
        counts = pd.Series(pred).value_counts()
        for label, count in counts.items():
            freq_rows.append({"method": mid, "predicted_organ": label, "n": int(count)})
        print("pred", mid, "n", len(pred), flush=True)
    pred_path = out / "pog570_predictions.parquet"
    pd.DataFrame(frame).to_parquet(pred_path, index=False)
    pd.DataFrame(freq_rows).to_csv(out / "pred_organ_freq.tsv", sep="\t", index=False)

    print("conformal thresholds", flush=True)
    meta = pd.read_csv(paths["config"] / "met500_eval_samples.tsv", sep="\t", dtype=str)
    meta["organ"] = meta["cohort"].map(MET500_COHORT_TO_ORGAN)
    meta["standard_site"] = [map_site("met500_biopsy_tissue", value) for value in meta["biopsy_tissue"]]
    eval_df = meta.loc[meta["organ"].notna() & (meta["organ"] != "unmapped")].reset_index(drop=True)
    met_samples, met_tpm = load_tpm(paths["data_processed"] / "met500/tpm_G.h5")
    met_row = {sample: i for i, sample in enumerate(met_samples)}
    block = sum_1e6(met_tpm[[met_row[sample] for sample in eval_df["Sample_id"]]])
    met_preds, _beta, _n, _raw = ev.predict_all(
        block, eval_df["standard_site"].tolist(), bundles, b0_idx, offsets, indices, mu, h_by_site, keep, v0
    )
    group = np.array([
        "liver" if site == "liver" else "lymph_node" if site == "lymph_node" else "other"
        for site in eval_df["standard_site"]
    ])
    truth = eval_df["organ"].to_numpy(dtype=object)
    rows = []
    for method, rep in (("BASE", "Z"), ("SA", "Z")):
        proba, _pred = met_preds[(method, rep)]
        scores = 1.0 - proba[np.arange(len(truth)), [ORGANS.index(label) for label in truth]]
        for alpha in (0.1, 0.2):
            threshold, k = conformal_threshold(scores, alpha)
            rows.append({
                "method": method_id(method, rep), "variant": "global", "alpha": alpha,
                "level": "all", "n_cal": int(len(scores)), "k": int(k),
                "threshold": "inf" if math.isinf(threshold) else repr(threshold),
            })
            for level in ("liver", "lymph_node", "other"):
                part = scores[group == level]
                threshold, k = conformal_threshold(part, alpha)
                rows.append({
                    "method": method_id(method, rep), "variant": "mondrian", "alpha": alpha,
                    "level": level, "n_cal": int(len(part)), "k": int(k),
                    "threshold": "inf" if math.isinf(threshold) else repr(threshold),
                })
    thr_path = paths["config"] / "conformal_thresholds_A3.tsv"
    pd.DataFrame(rows).to_csv(thr_path, sep="\t", index=False)
    hash_rows = [
        {"path": "results/stage3/pog570_predictions.parquet", "size_bytes": pred_path.stat().st_size, "sha256": sha256_file(pred_path)},
        {"path": "config/conformal_thresholds_A3.tsv", "size_bytes": thr_path.stat().st_size, "sha256": sha256_file(thr_path)},
    ]
    pd.DataFrame(hash_rows).to_csv(paths["config"] / "frozen_predictions_A3.tsv", sep="\t", index=False)
    print("frozen", hash_rows, flush=True)


if __name__ == "__main__":
    main()
