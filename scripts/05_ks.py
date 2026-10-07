#!/usr/bin/env python3
import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import json

import h5py
import numpy as np
import pandas as pd

from ks_fast import ks_matrix, original_getes, right_and_tie
from mapping_rules import SITE_NATIVE
from util import SEED, load_paths


def decode_samples(raw) -> list[str]:
    out = []
    for item in raw:
        if isinstance(item, bytes):
            out.append(item.decode())
        else:
            out.append(str(item))
    return out


def load_tpm(path: Path):
    with h5py.File(path, "r") as handle:
        samples = decode_samples(handle["samples"][:])
        tpm = handle["tpm"][:]
    return samples, tpm


def validate(met_tpm, tcga_tpm, offsets, indices, n_sets: int) -> dict:
    rng = np.random.default_rng(SEED)
    zero_frac = (met_tpm == 0).mean(axis=1)
    high = np.argsort(zero_frac)[-30:]
    pairs = []
    # 400 pairs: high-zero MET500. 300 pairs: random MET500. 300 pairs: TCGA.
    for _ in range(400):
        pairs.append(("met", int(rng.choice(high)), int(rng.integers(0, n_sets))))
    for _ in range(300):
        pairs.append(("met", int(rng.integers(0, met_tpm.shape[0])), int(rng.integers(0, n_sets))))
    for _ in range(300):
        pairs.append(("tcga", int(rng.integers(0, tcga_tpm.shape[0])), int(rng.integers(0, n_sets))))
    # Score the paired matrices once with the fast method, then compare 1000 entries.
    blocks = {}
    for label, matrix in (("met", met_tpm), ("tcga", tcga_tpm)):
        # Validate on a sample subset to avoid scoring every row when only some are used.
        used = sorted({i for src, i, _ in pairs if src == label})
        sub = matrix[used].astype(np.float32, copy=False)
        right, tie = right_and_tie(sub)
        scores = ks_matrix(right, tie, offsets, indices)
        blocks[label] = (used, scores)
    worst = 0.0
    sign_mismatch = 0
    n = 0
    for src, sample_i, set_i in pairs:
        used, scores = blocks[src]
        row = used.index(sample_i)
        fast = float(scores[row, set_i])
        members = indices[offsets[set_i]:offsets[set_i + 1]]
        matrix = met_tpm if src == "met" else tcga_tpm
        slow = original_getes(matrix[sample_i].astype(np.float64), members)
        worst = max(worst, abs(fast - slow))
        if fast * slow < 0:
            sign_mismatch += 1
        n += 1
    return {
        "n_pairs": n,
        "max_abs_diff": worst,
        "n_sign_mismatch": sign_mismatch,
        "passed": bool(worst <= 1e-6 and sign_mismatch == 0),
        "n_high_zero_met500_samples": int(len(high)),
        "high_zero_fraction_min": float(zero_frac[high].min()),
    }


def score_samples(tpm: np.ndarray, offsets, indices) -> np.ndarray:
    parts = []
    batch = 2000
    for start in range(0, tpm.shape[0], batch):
        block = tpm[start:start + batch].astype(np.float32, copy=False)
        right, tie = right_and_tie(block)
        parts.append(ks_matrix(right, tie, offsets, indices))
        print("scored", start + block.shape[0], "/", tpm.shape[0], flush=True)
    return np.vstack(parts)


def save_parquet(path: Path, samples: list[str], scores: np.ndarray, names: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(scores, index=pd.Index(samples, name="sample"), columns=names)
    frame.to_parquet(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compute", action="store_true")
    args = parser.parse_args()
    paths = load_paths()
    genes = pd.read_csv(paths["data_processed"] / "genes/G_symbols.txt", header=None)[0].astype(str).tolist()
    packed = np.load(paths["data_processed"] / "genes/gene_sets.npz", allow_pickle=True)
    offsets = packed["offsets"]
    indices = packed["indices"]
    names = [str(x) for x in packed["names"].tolist()]
    print("G", len(genes), "sets", len(names))

    met_samples, met_tpm = load_tpm(paths["data_processed"] / "met500/tpm_G.h5")
    toil_samples, toil_tpm = load_tpm(paths["data_processed"] / "toil/tpm_G.h5")
    assert met_tpm.shape[1] == len(genes)
    assert toil_tpm.shape[1] == len(genes)
    # Validation uses a TCGA slice of 400 expression rows, not a label join.
    split = pd.read_csv(paths["config"] / "split_tcga.tsv", sep="\t", dtype=str)
    tcga_ids = split.loc[split["split"].isin(["TCGA-train", "TCGA-test", "TCGA-met"]), "sample"]
    tcga_index = {sample: i for i, sample in enumerate(toil_samples)}
    tcga_rows = [tcga_index[s] for s in tcga_ids if s in tcga_index]
    report = validate(met_tpm, toil_tpm[tcga_rows], offsets, indices, len(names))
    report["implementation"] = "fast" if report["passed"] else "REJECTED"
    out = paths["results"] / "ks_validation.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit("Fast KS failed validation. Full scoring was not run.")
    if not args.compute:
        print("validation only")
        return

    feature_dir = paths["data_processed"] / "features"
    # TCGA learning + metastatic rows.
    save_parquet(feature_dir / "ks_tcga.parquet", [toil_samples[i] for i in tcga_rows], score_samples(toil_tpm[tcga_rows], offsets, indices), names)
    # GTEx-ref samples in tissues named by site_native, plus the three simulation tissues' sim split is not scored here.
    gtex = pd.read_csv(paths["config"] / "split_gtex.tsv", sep="\t", dtype=str)
    needed = set()
    for tissues, _, _ in SITE_NATIVE.values():
        needed.update(tissues)
    ref = gtex.loc[(gtex["split"] == "GTEx-ref") & (gtex["tissue"].isin(needed))]
    ref_rows = [tcga_index[s] for s in ref["sample"] if s in tcga_index]
    missing = int((~ref["sample"].isin(tcga_index)).sum())
    print("gtex ref rows", len(ref_rows), "missing", missing)
    save_parquet(
        feature_dir / "ks_gtex_ref.parquet",
        [toil_samples[i] for i in ref_rows],
        score_samples(toil_tpm[ref_rows], offsets, indices),
        names,
    )
    save_parquet(feature_dir / "ks_met500.parquet", met_samples, score_samples(met_tpm, offsets, indices), names)
    print("ks parquet written")


if __name__ == "__main__":
    main()
