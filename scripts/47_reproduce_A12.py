#!/usr/bin/env python3
"""Rebuild MET500 and brca Anders scores from raw files and compare stored predictions."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "release" / "scripts"))

from features_from_raw import brca_z, met500_z  # noqa: E402
from predict_portable import load_model, predict, predict_proba  # noqa: E402


def compare(name: str, z: np.ndarray, stored_pred: np.ndarray, stored_proba: np.ndarray, organs: list[str]) -> dict:
    model = load_model(ROOT / "release" / "models", name)
    got_proba = predict_proba(model, z)
    got_pred = predict(model, z)
    if stored_proba.shape != got_proba.shape:
        raise SystemExit(f"{name} probability shape {stored_proba.shape} != {got_proba.shape}")
    order = list(model["meta"]["organs"])
    if organs != order:
        index = [organs.index(organ) for organ in order]
        stored_proba = stored_proba[:, index]
    return {
        "n": int(len(stored_pred)),
        "top1_matches": int(np.sum(got_pred == stored_pred)),
        "max_abs_proba": float(np.max(np.abs(got_proba - stored_proba))),
    }


def main() -> None:
    meta = json.loads((ROOT / "release" / "models" / "SA_Z.json").read_text())
    organs = list(meta["organs"])
    samples = pd.read_csv(ROOT / "results/stage5/posthoc/met500_samples.tsv", sep="\t", dtype=str)
    ids = samples["id"].tolist()
    z = met500_z(ids)
    report = {"MET500": {}}
    for name, column in (("BASE_Z", "BASE-Z"), ("SA_Z", "SA-Z")):
        stored = np.load(ROOT / "results/stage5/posthoc" / f"met500_{column}_proba.npy")
        report["MET500"][name] = compare(name, z, samples[f"pred__{column}"].to_numpy(), stored, organs)
    brca, brca_ids = brca_z()
    stored_z = np.load(ROOT / "data/processed/aux/brca_iatlas_anders_2022/Z.npy")
    report["brca_anders"] = {"n_features": int(len(brca_ids)), "z_max_abs": float(np.max(np.abs(brca - stored_z)))}
    frame = pd.read_parquet(ROOT / "results/stage7/confirm/predictions.parquet")
    part = frame.loc[frame["cohort"] == "brca_iatlas_anders_2022"].copy()
    part = part.set_index("sample_id").loc[brca_ids].reset_index()
    for name, column in (("BASE_Z", "BASE-Z"), ("SA_Z", "SA-Z")):
        stored = np.column_stack([part[f"{column}__p_{organ}"].to_numpy(dtype=np.float64) for organ in organs])
        report["brca_anders"][name] = compare(name, brca, part[f"{column}__pred"].to_numpy(), stored, organs)
    out = ROOT / "manuscript" / "checks" / "reproduction_A12.json"
    out.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
