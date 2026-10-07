#!/usr/bin/env python3
"""Add the post hoc unstandardized mixture model to release/models and check it.

The stage-9 joblib file is copied byte for byte; the portable NumPy archive and JSON
are written with export_one() of 46_portable_A12.py. MET500 predictions rebuilt from
the stored coefficients, through both the joblib bundle and the portable files, are
compared with the stage-9 predictions. Nothing is refitted.
"""
from __future__ import annotations

import hashlib
import importlib
import json
import shutil
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "release" / "scripts"))
P46 = importlib.import_module("46_portable_A12")

from predict_portable import load_model, predict, predict_proba  # noqa: E402
from stage2_common import load_gene_pack, load_tpm, sum_1e6  # noqa: E402
from stage5_lib import predict_organ, rank_b0  # noqa: E402

SOURCE = ROOT / "results" / "stage9" / "ablation" / "models" / "MIX-Z0.joblib"
PREDICTIONS = ROOT / "results" / "stage9" / "ablation" / "predictions.tsv"
OUT = P46.OUT
NAME = "MIX_Z0_posthoc"
STATUS = "post hoc ablation model; not evaluated in a preregistered confirmation"
CHECK = ROOT / "manuscript" / "checks" / "posthoc_model_check.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    bundle = joblib.load(SOURCE)
    clf = bundle["clf"]
    if bundle["scaler"] is not None or float(clf.C) != 0.1:
        raise SystemExit(f"{SOURCE.name}: scaler {bundle['scaler']}, C {clf.C}")
    symbols_before = sha256(OUT / "G_symbols.txt")
    shutil.copyfile(SOURCE, OUT / f"{NAME}.joblib")
    P46.export_one(NAME, SOURCE)
    if sha256(OUT / "G_symbols.txt") != symbols_before:
        raise SystemExit("G_symbols.txt changed")
    meta = json.loads((OUT / f"{NAME}.json").read_text())
    meta["status"] = STATUS
    meta["training"] = ("pure TCGA tumors and the HostMix-TOO mixtures (same mixture records as SA_Z), "
                        "no feature standardization, C = 0.1")
    (OUT / f"{NAME}.json").write_text(json.dumps(meta))
    if "scaler_mean" in np.load(OUT / f"{NAME}.npz", allow_pickle=False).files:
        raise SystemExit("portable archive has a scaler")

    stored = pd.read_csv(PREDICTIONS, sep="\t")
    stored = stored.loc[stored["cohort"].eq("MET500") & stored["method"].eq("MIX-Z0")].reset_index(drop=True)
    paths = P46.paths()
    _genes, b0_idx, *_rest = load_gene_pack(paths)
    samples, tpm = load_tpm(paths["data_processed"] / "met500" / "tpm_G.h5")
    row = {sample: i for i, sample in enumerate(samples)}
    matrix = rank_b0(sum_1e6(tpm[[row[sample] for sample in stored["sample_id"]]]), b0_idx)
    organ, pred = predict_organ(joblib.load(OUT / f"{NAME}.joblib"), matrix)
    portable = load_model(OUT, NAME)
    got = predict_proba(portable, matrix)
    top = predict(portable, matrix)
    truth = stored["pred"].to_numpy(dtype=object)
    report = {
        "model": NAME,
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": sha256(SOURCE),
        "files": {f"{NAME}{suffix}": sha256(OUT / f"{NAME}{suffix}") for suffix in (".joblib", ".npz", ".json")},
        "n_met500": int(len(stored)),
        "joblib_top1_mismatch_vs_stage9": int(np.sum(pred != truth)),
        "portable_top1_mismatch_vs_stage9": int(np.sum(top != truth)),
        "max_abs_proba_joblib_vs_portable": float(np.max(np.abs(organ - got))),
    }
    report["pass"] = bool(report["joblib_top1_mismatch_vs_stage9"] == 0
                          and report["portable_top1_mismatch_vs_stage9"] == 0
                          and report["max_abs_proba_joblib_vs_portable"] <= 1e-6)
    CHECK.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["pass"]:
        raise SystemExit("reproduction check failed")


if __name__ == "__main__":
    main()
