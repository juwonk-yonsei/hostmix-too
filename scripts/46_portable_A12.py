#!/usr/bin/env python3
"""Export portable weights and check them against the stored joblib models."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import joblib
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "release" / "scripts"))

from mapping_rules import PROJECT_TO_ORGAN  # noqa: E402
from metrics import ORGANS  # noqa: E402
from predict_portable import load_model, predict, predict_proba  # noqa: E402
from stage2_common import load_gene_pack, load_tpm, sum_1e6  # noqa: E402
from stage5_lib import predict_organ, rank_b0  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "release" / "models"


def paths() -> dict:
    out = {"root": ROOT}
    for line in (ROOT / "config" / "paths.yaml").read_text().splitlines():
        text = line.split("#", 1)[0].strip()
        if not text or ":" not in text:
            continue
        key, value = text.split(":", 1)
        key = key.strip()
        if key == "n_threads":
            continue
        value = value.strip()
        if key != "seed":
            out[key] = ROOT / value
    return out


def export_one(name: str, source: Path) -> None:
    bundle = joblib.load(source)
    clf = bundle["clf"]
    classes = [str(c) for c in clf.classes_]
    payload = {
        "coef": clf.coef_.astype(np.float64),
        "intercept": clf.intercept_.astype(np.float64),
    }
    if bundle["scaler"] is not None:
        payload["scaler_mean"] = bundle["scaler"].mean_.astype(np.float64)
        payload["scaler_scale"] = bundle["scaler"].scale_.astype(np.float64)
    np.savez(OUT / f"{name}.npz", **payload)
    genes, b0_idx, *_rest = load_gene_pack(paths())
    meta = {
        "name": name,
        "C": float(clf.C),
        "class_weight": str(clf.class_weight),
        "solver": str(clf.solver),
        "max_iter": int(clf.max_iter),
        "classes": classes,
        "organs": list(ORGANS),
        "organ_of_class": {label: PROJECT_TO_ORGAN[label] for label in classes},
        "b0_genes": [genes[int(i)] for i in b0_idx],
        "z_rule": "within-sample average ranks, then norm.ppf((rank-0.5)/n_genes) on all G columns, then B0",
        "n_G": len(genes),
    }
    (OUT / f"{name}.json").write_text(json.dumps(meta))
    (OUT / "G_symbols.txt").write_text("\n".join(genes) + "\n")


def main() -> None:
    cfg = paths()
    export_one("BASE_Z", ROOT / "results/stage2/models/BASE_Z.joblib")
    export_one("SA_Z", ROOT / "results/stage2/models/SA_Z.joblib")
    samples = [line for line in (ROOT / "results/stage5/posthoc/met500_samples.tsv").read_text().splitlines()[1:]]
    ids = [line.split("\t")[0] for line in samples]
    raw_samples, tpm = load_tpm(cfg["data_processed"] / "met500" / "tpm_G.h5")
    row = {sample: i for i, sample in enumerate(raw_samples)}
    genes, b0_idx, *_rest = load_gene_pack(cfg)
    matrix = rank_b0(sum_1e6(tpm[[row[sample] for sample in ids]]), b0_idx)
    report = {}
    for name in ("BASE_Z", "SA_Z"):
        bundle = joblib.load(ROOT / "results/stage2/models" / f"{name}.joblib")
        organ, pred = predict_organ(bundle, matrix)
        portable = load_model(OUT, name)
        got = predict_proba(portable, matrix)
        top = predict(portable, matrix)
        gap = float(np.max(np.abs(organ - got)))
        mismatch = int(np.sum(pred != top))
        report[name] = {"max_abs_proba": gap, "top1_mismatch": mismatch, "n": int(len(ids))}
        if gap > 1e-10 or mismatch:
            raise SystemExit(report)
    (ROOT / "manuscript" / "checks" / "portable_equivalence.json").write_text(json.dumps(report, indent=2))
    print(report)


if __name__ == "__main__":
    main()
