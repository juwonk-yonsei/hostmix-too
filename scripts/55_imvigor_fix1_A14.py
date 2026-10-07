"""Recalculate IMvigor210 esophagus-minus-bladder contributions.

Baseline input is the rank-normal score. HostMix-TOO input is that score
after the frozen standard scaler. Does not retrain and does not overwrite
results/stage9/imvigor_contribution.tsv.
"""
from __future__ import annotations

from pathlib import Path

import joblib
import numpy as np
import pandas as pd

ROOT = Path("/home/kangjw/WORK/Cisplatin/2026_NewProj/proj_A_host_too")
OUT = ROOT / "results/stage9/imvigor_contribution_fix1.tsv"
ORIGINAL = ROOT / "results/stage9/imvigor_contribution.tsv"


def class_index(bundle, label: str) -> int:
    classes = [str(item) for item in bundle["clf"].classes_]
    if label not in classes:
        raise SystemExit(f"missing class {label}")
    return classes.index(label)


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"refusing to overwrite {OUT}")
    genes = [gene for gene in (ROOT / "results/B0_genes.txt").read_text().splitlines() if gene and gene != "symbol"]
    confirm = pd.read_parquet(
        ROOT / "results/stage7/confirm/predictions.parquet",
        columns=["cohort", "sample_id", "selected_native", "BASE-Z__pred", "SA-Z__pred"],
    )
    confirm["selected_native"] = confirm["selected_native"].astype(str).str.lower().isin(["true", "1"])
    native = confirm.loc[confirm["cohort"].eq("blca_iatlas_imvigor210_2017") & confirm["selected_native"]].copy()
    ids = (ROOT / "data/processed/aux/blca_iatlas_imvigor210_2017/sample_ids.txt").read_text().splitlines()
    z = np.load(ROOT / "data/processed/aux/blca_iatlas_imvigor210_2017/Z.npy")
    index = {sample: i for i, sample in enumerate(ids)}
    native["row"] = native["sample_id"].astype(str).map(index)
    if native["row"].isna().any():
        raise SystemExit("IMvigor Z id missing")
    rows = []
    for method, path, n_expected, scaled in (
        ("BASE-Z", "BASE_Z.joblib", 145, False),
        ("SA-Z", "SA_Z.joblib", 65, True),
    ):
        bundle = joblib.load(ROOT / "results/stage2/models" / path)
        coef = np.asarray(bundle["clf"].coef_, dtype=np.float64)
        delta = coef[class_index(bundle, "ESCA")] - coef[class_index(bundle, "BLCA")]
        called = native.loc[native[f"{method}__pred"].eq("Esophagus")]
        if len(called) != n_expected:
            raise SystemExit(f"{method} esophagus {len(called)}")
        block = z[[int(i) for i in called["row"]]].astype(np.float64)
        if scaled:
            scaler = bundle["scaler"]
            mu = np.asarray(scaler.mean_, dtype=np.float64)
            sigma = np.asarray(scaler.scale_, dtype=np.float64)
            if mu.shape != delta.shape or sigma.shape != delta.shape:
                raise SystemExit(f"scaler shape {mu.shape} coef {delta.shape}")
            if np.any(sigma == 0):
                raise SystemExit("zero scaler scale")
            model_input = (block - mu) / sigma
        else:
            model_input = block
        mean_contrib = (delta * model_input).mean(axis=0)
        order = np.argsort(-np.abs(mean_contrib))
        for rank, gene_i in enumerate(order, start=1):
            rows.append({
                "method": method,
                "n_samples": int(len(called)),
                "rank": rank,
                "gene": genes[gene_i],
                "mean_contribution": float(mean_contrib[gene_i]),
                "input": "standardized rank-normal score" if scaled else "rank-normal score",
            })
    frame = pd.DataFrame(rows)
    original = pd.read_csv(ORIGINAL, sep="\t")
    base_new = frame.loc[frame["method"].eq("BASE-Z")].set_index("gene")["mean_contribution"]
    base_old = original.loc[original["method"].eq("BASE-Z")].set_index("gene")["mean_contribution"]
    shared = base_old.index.intersection(base_new.index)
    gap = (base_new.loc[shared] - base_old.loc[shared].astype(float)).abs().max()
    if gap > 1e-8:
        raise SystemExit(f"baseline contribution disagrees with the original file by {gap}")
    top = frame.loc[frame["method"].eq("BASE-Z")].nsmallest(5, "rank")["gene"].tolist()
    expected = ["NPIPB5", "NPIPB4", "REG3A", "IRX2", "REG1A"]
    if top != expected:
        raise SystemExit(f"baseline top genes {top}")
    frame.to_csv(OUT, sep="\t", index=False)
    print("wrote", OUT, "rows", len(frame))
    print("baseline top", top)
    print("hostmix top", frame.loc[frame["method"].eq("SA-Z")].nsmallest(5, "rank")["gene"].tolist())


if __name__ == "__main__":
    main()
