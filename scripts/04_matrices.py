#!/usr/bin/env python3
import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import json

import h5py
import numpy as np
import pandas as pd

from util import load_paths, now_iso


def histogram_median(counts: np.ndarray, edges: np.ndarray):
    total = int(counts.sum())
    if total == 0:
        return None
    half = total / 2.0
    index = int(np.searchsorted(np.cumsum(counts), half, side="left"))
    index = min(index, len(counts) - 1)
    return float(0.5 * (edges[index] + edges[index + 1]))


def build_one(path: Path, id_to_col: dict[str, int], n_genes: int, out_path: Path, scale_path: Path) -> None:
    if out_path.exists() and scale_path.exists():
        print("skip existing", out_path.name)
        return
    # dtype cannot include the gene-id column; pandas applies it before index_col.
    reader = pd.read_csv(path, sep="\t", index_col=0, chunksize=400)
    first = next(reader)
    samples = [str(c) for c in first.columns]
    n_samples = len(samples)
    matrix = np.zeros((n_samples, n_genes), dtype=np.float32)
    first_values = first.astype(np.float32).to_numpy()
    is_log = bool(np.nanmin(first_values) < 0)
    edges = np.linspace(-15.0, 30.0, 45001) if is_log else None
    counts = np.zeros(len(edges) - 1, dtype=np.int64) if is_log else None
    n_outside = 0
    raw_min = np.inf
    raw_max = -np.inf
    n_values = 0
    n_negative = 0
    n_nonfinite = 0
    n_rows = 0
    n_kept_rows = 0
    linear_values = []

    def consume(chunk: pd.DataFrame) -> None:
        nonlocal raw_min, raw_max, n_values, n_negative, n_nonfinite, n_rows, n_kept_rows, n_outside
        values = chunk.astype(np.float32).to_numpy()
        finite = np.isfinite(values)
        n_nonfinite += int((~finite).sum())
        if finite.any():
            raw_min = min(raw_min, float(values[finite].min()))
            raw_max = max(raw_max, float(values[finite].max()))
            n_negative += int(np.sum(values[finite] < 0))
            if is_log:
                finite_values = values[finite]
                n_outside += int(np.sum((finite_values < edges[0]) | (finite_values > edges[-1])))
                counts[:] += np.histogram(finite_values, bins=edges)[0]
            else:
                linear_values.append(values[finite].ravel().astype(np.float32, copy=True))
        n_values += int(values.size)
        n_rows += int(values.shape[0])
        keep_rows = []
        keep_cols = []
        for i, gene_id in enumerate(chunk.index.astype(str)):
            col = id_to_col.get(gene_id)
            if col is not None:
                keep_rows.append(i)
                keep_cols.append(col)
        if not keep_rows:
            return
        block = values[keep_rows].astype(np.float32, copy=True)
        if is_log:
            block = np.maximum(np.exp2(block) - np.float32(0.001), np.float32(0))
        block = np.nan_to_num(block, nan=0.0, posinf=0.0, neginf=0.0)
        for row_i, col in enumerate(keep_cols):
            matrix[:, col] += block[row_i]
        n_kept_rows += len(keep_rows)

    consume(first)
    for i, chunk in enumerate(reader, start=1):
        consume(chunk)
        if i % 20 == 0:
            print(path.name, "chunks", i, "rows", n_rows, flush=True)
    if (not is_log) and raw_min < 0:
        raise SystemExit(f"{path.name} looked linear in the first chunk but later values are negative")
    if is_log:
        median = histogram_median(counts, edges)
        median_method = "histogram median, bin width 0.001, range [-15, 30]"
    else:
        all_linear = np.concatenate(linear_values) if linear_values else np.array([])
        median = float(np.median(all_linear)) if all_linear.size else None
        median_method = "exact median of all finite raw values"
        del all_linear
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_path, "w") as handle:
        handle.create_dataset("tpm", data=matrix)
        handle.create_dataset("samples", data=np.array(samples, dtype="S80"))
        handle.attrs["is_log_input"] = is_log
        handle.attrs["n_kept_rows"] = n_kept_rows
    scale = {
        "file": str(path.name),
        "is_log_input": is_log,
        "raw_min": None if raw_min == np.inf else raw_min,
        "raw_max": None if raw_max == -np.inf else raw_max,
        "raw_median": median,
        "median_method": median_method,
        "n_values": n_values,
        "n_negative": n_negative,
        "n_nonfinite": n_nonfinite,
        "n_outside_histogram_range": n_outside if is_log else 0,
        "n_rows": n_rows,
        "n_samples": n_samples,
        "n_kept_ensembl_rows": n_kept_rows,
        "datetime": now_iso(),
    }
    scale_path.write_text(json.dumps(scale, indent=2))
    print("wrote", out_path.name, "log" if is_log else "linear", "min", scale["raw_min"], "max", scale["raw_max"])


def main() -> None:
    paths = load_paths()
    gene_dir = paths["data_processed"] / "genes"
    symbols = pd.read_csv(gene_dir / "G_symbols.txt", header=None)[0].astype(str).tolist()
    symbol_to_col = {symbol: i for i, symbol in enumerate(symbols)}
    # Rebuild ensembl -> column from the probemap using the same rule as genes.py.
    from genes import probemap_symbol

    probe = probemap_symbol()
    id_to_col = {}
    for ens, sym in zip(probe["id"], probe["gene"]):
        col = symbol_to_col.get(sym)
        if col is not None:
            id_to_col[ens] = col
    with h5py.File(gene_dir / "G_index.h5", "w") as handle:
        handle.create_dataset("genes", data=np.array(symbols, dtype="S32"))
    out_dir = paths["data_processed"]
    build_one(
        paths["data_raw"] / "toil/TcgaTargetGtex_rsem_gene_tpm.gz",
        id_to_col,
        len(symbols),
        out_dir / "toil/tpm_G.h5",
        paths["results"] / "audit/toil_scale.json",
    )
    build_one(
        paths["data_raw"] / "met500/M.mx.txt.gz",
        id_to_col,
        len(symbols),
        out_dir / "met500/fpkm_sum_G.h5",
        paths["results"] / "audit/met500_fpkm_scale.json",
    )
    # FPKM on G -> TPM = FPKM / sum(FPKM) * 1e6. Zeros-sum samples stay zero.
    with h5py.File(out_dir / "met500/fpkm_sum_G.h5", "r") as src, h5py.File(out_dir / "met500/tpm_G.h5", "w") as dst:
        fpkm = src["tpm"][:]
        totals = fpkm.sum(axis=1, keepdims=True)
        totals[totals == 0] = np.nan
        tpm = fpkm / totals * np.float32(1e6)
        tpm = np.nan_to_num(tpm, nan=0.0).astype(np.float32)
        dst.create_dataset("tpm", data=tpm)
        dst.create_dataset("samples", data=src["samples"][:])
        dst.attrs["zero_sum_samples"] = int(np.sum(fpkm.sum(axis=1) == 0))
    print("met500 tpm done", tpm.shape)


if __name__ == "__main__":
    main()
