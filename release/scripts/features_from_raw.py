"""Gene map and within-sample rank scores from received raw files.

Does not load joblib and does not set a thread cap.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

ROOT = Path(__file__).resolve().parents[2]
MODELS = Path(__file__).resolve().parents[1] / "models"


def genes_and_b0() -> tuple[list[str], list[str]]:
    genes = (MODELS / "G_symbols.txt").read_text().splitlines()
    meta = json.loads((MODELS / "SA_Z.json").read_text())
    return genes, list(meta["b0_genes"])


def rank_normal(values: np.ndarray) -> np.ndarray:
    out = np.empty(values.shape, dtype=np.float64)
    n_genes = values.shape[1]
    for start in range(0, values.shape[0], 500):
        block = values[start:start + 500]
        ranks = rankdata(block, method="average", axis=1)
        out[start:start + block.shape[0]] = norm.ppf((ranks - 0.5) / n_genes)
    return out


def take_b0(full: np.ndarray, genes: list[str], b0: list[str]) -> np.ndarray:
    index = {gene: i for i, gene in enumerate(genes)}
    columns = np.array([index[gene] for gene in b0], dtype=np.int32)
    return full[:, columns].astype(np.float32)


def ensembl_to_column(genes: list[str]) -> dict[str, int]:
    symbol_to_col = {symbol: i for i, symbol in enumerate(genes)}
    hgnc = pd.read_csv(ROOT / "data/raw/hgnc/hgnc_complete_set.txt", sep="\t", dtype=str, low_memory=False)
    symbols = set(hgnc.loc[hgnc["locus_group"] == "protein-coding gene", "symbol"].dropna().astype(str))
    probe = pd.read_csv(ROOT / "data/raw/toil/gencode.v23.annotation.gene.probemap", sep="\t", dtype=str)
    probe = probe.loc[probe["gene"].isin(symbols), ["id", "gene"]].drop_duplicates()
    mapping = {}
    for ens, sym in zip(probe["id"].astype(str), probe["gene"].astype(str)):
        col = symbol_to_col.get(sym)
        if col is not None:
            mapping[ens] = col
    return mapping


def met500_tpm(sample_ids: list[str]) -> np.ndarray:
    """Versioned Ensembl rows summed onto G, then FPKM converted to TPM."""
    genes, _b0 = genes_and_b0()
    id_to_col = ensembl_to_column(genes)
    path = ROOT / "data/raw/met500/M.mx.txt.gz"
    reader = pd.read_csv(path, sep="\t", index_col=0, chunksize=400)
    first = next(reader)
    samples = [str(column) for column in first.columns]
    wanted = {sample: i for i, sample in enumerate(samples)}
    missing = [sample for sample in sample_ids if sample not in wanted]
    if missing:
        raise SystemExit(f"MET500 raw file lacks {len(missing)} requested samples")
    matrix = np.zeros((len(sample_ids), len(genes)), dtype=np.float32)
    row_of = {sample: i for i, sample in enumerate(sample_ids)}
    columns = [wanted[sample] for sample in sample_ids]

    def consume(chunk: pd.DataFrame) -> None:
        values = chunk.to_numpy(dtype=np.float32)
        keep_rows = []
        keep_cols = []
        for i, gene_id in enumerate(chunk.index.astype(str)):
            col = id_to_col.get(gene_id)
            if col is not None:
                keep_rows.append(i)
                keep_cols.append(col)
        if not keep_rows:
            return
        block = np.nan_to_num(values[keep_rows][:, columns], nan=0.0, posinf=0.0, neginf=0.0)
        for row_i, col in enumerate(keep_cols):
            matrix[:, col] += block[row_i]

    consume(first)
    for chunk in reader:
        consume(chunk)
    totals = matrix.sum(axis=1, keepdims=True).astype(np.float64)
    totals[totals == 0] = np.nan
    tpm = np.nan_to_num(matrix.astype(np.float64) / totals * 1e6, nan=0.0)
    return tpm.astype(np.float32)


def met500_z(sample_ids: list[str]) -> np.ndarray:
    genes, b0 = genes_and_b0()
    tpm = met500_tpm(sample_ids)
    totals = tpm.sum(axis=1, keepdims=True).astype(np.float64)
    totals[totals == 0] = np.nan
    scaled = np.nan_to_num(tpm.astype(np.float64) / totals * 1e6, nan=0.0).astype(np.float32)
    return take_b0(rank_normal(scaled), genes, b0)


def brca_z() -> tuple[np.ndarray, list[str]]:
    """cBioPortal Hugo symbols, column sum 1e6, ranks on genes present in G."""
    genes, b0 = genes_and_b0()
    path = ROOT / "data/processed/aux/downloads/brca_iatlas_anders_2022__data_mrna_seq_tpm.txt"
    frame = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    symbol_col = "Hugo_Symbol" if "Hugo_Symbol" in frame.columns else frame.columns[0]
    frame = frame.loc[frame[symbol_col].notna() & (frame[symbol_col] != "") & (frame[symbol_col] != "NA")]
    ids = (ROOT / "data/processed/aux/brca_iatlas_anders_2022/sample_ids.txt").read_text().splitlines()
    ids = [item for item in ids if item]
    sample_cols = [col for col in frame.columns if col not in {symbol_col, "Entrez_Gene_Id"}]
    values = frame[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    values.index = frame[symbol_col].astype(str).to_numpy()
    grouped = values.groupby(level=0).sum()
    missing = [sample for sample in ids if sample not in grouped.columns]
    if missing:
        raise SystemExit(f"brca raw file lacks {len(missing)} stored samples")
    grouped = grouped[ids]
    if float(np.nanmin(grouped.to_numpy(dtype=float))) < 0:
        raise SystemExit("brca file has negative values")
    totals = grouped.sum(axis=0).replace(0, np.nan)
    tpm = grouped.div(totals, axis=1) * 1e6
    present = [gene for gene in genes if gene in tpm.index]
    block = np.nan_to_num(tpm.loc[present].to_numpy(dtype=np.float64).T, nan=0.0)
    scored = rank_normal(block)
    present_index = {gene: i for i, gene in enumerate(present)}
    z = np.zeros((len(ids), len(b0)), dtype=np.float32)
    for j, symbol in enumerate(b0):
        if symbol in present_index:
            z[:, j] = scored[:, present_index[symbol]].astype(np.float32)
    return z, ids
