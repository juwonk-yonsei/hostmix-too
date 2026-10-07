"""Shared stage-2 settings. Rules here are fixed before section 5 runs.

Stage 5 removes the thread caps that this module used to set on import.
Stage 1-2 scripts that still set those variables do so in their own files.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import h5py
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.stats import rankdata
from sklearn.model_selection import train_test_split

from ks_fast import ks_matrix, right_and_tie
from mapping_rules import PROJECT_TO_ORGAN, SITE_NATIVE
from util import SEED, load_paths

CANDIDATE_TISSUES = [
    "Spleen",
    "Whole Blood",
    "Cells - Ebv-Transformed Lymphocytes",
    "Small Intestine - Terminal Ileum",
]
# Tie-break for proxy correlation is this list order (earlier wins).
HPA_TARGETS = ("lymph node", "bone marrow")

# Nine site models. bone and bone_marrow share bone_marrow.
# skin and subcutaneous share skin. omentum and peritoneum share omentum.
SITE_MODEL_ORDER = [
    "adrenal",
    "bone_marrow",
    "brain",
    "liver",
    "lung",
    "lymph_node",
    "omentum",
    "skin",
    "soft_tissue",
]
SITE_TO_MODEL = {
    "liver": "liver",
    "lung": "lung",
    "lymph_node": "lymph_node",
    "bone_marrow": "bone_marrow",
    "bone": "bone_marrow",
    "soft_tissue": "soft_tissue",
    "skin": "skin",
    "subcutaneous": "skin",
    "brain": "brain",
    "adrenal": "adrenal",
    "omentum": "omentum",
    "peritoneum": "omentum",
}
# Fixed GTEx tissues. Proxy tissues are filled after section 2.
FIXED_SITE_TISSUES = {
    "liver": ["Liver"],
    "lung": ["Lung"],
    "soft_tissue": ["Adipose - Subcutaneous", "Muscle - Skeletal"],
    "skin": ["Skin - Not Sun Exposed (Suprapubic)"],
    "brain": ["Brain - Cortex"],
    "adrenal": ["Adrenal Gland"],
    "omentum": ["Adipose - Visceral (Omentum)"],
}
SIM_TISSUES_FIXED = ["Liver", "Lung", "Adipose - Subcutaneous"]
SIM_RHOS = [1.0, 0.8, 0.6, 0.4, 0.2]
SIM_TISSUE_TO_SITE = {
    "Liver": "liver",
    "Lung": "lung",
    "Adipose - Subcutaneous": "soft_tissue",
}
C_GRID = [0.03, 0.1, 0.3]
N_MIX = 4
RHO_LOW = 0.15
RHO_HIGH = 1.0
IF_N = 1000
IF_RHO = 0.6
IF_Q = {"IF20": 0.20, "IF40": 0.40}
BASE_C = {"Z": 0.1, "K": 0.1}


def model_names(ln_proxy: str, bm_proxy: str) -> list[str]:
    names = ["IF_pool", "SA"]
    for site in SITE_MODEL_ORDER:
        names.append(f"SC_{site}")
        names.append(f"IF_SC_{site}")
    return sorted(names)


def augmentation_seeds(names: list[str]) -> dict[str, int]:
    """One Generator(SEED), names in sorted order, each seed is one integers() draw."""
    rng = np.random.default_rng(SEED)
    ordered = sorted(names)
    return {name: int(rng.integers(0, 2**31 - 1)) for name in ordered}


def n_drop(n_features: int, q: float) -> int:
    """Number of features removed: floor(n * q)."""
    return int(np.floor(n_features * q))


def decode_samples(raw) -> list[str]:
    return [item.decode() if isinstance(item, bytes) else str(item) for item in raw]


def load_tpm(path: Path):
    with h5py.File(path, "r") as handle:
        samples = decode_samples(handle["samples"][:])
        tpm = handle["tpm"][:]
    return samples, tpm


def sum_1e6(tpm: np.ndarray) -> np.ndarray:
    totals = tpm.sum(axis=1, keepdims=True).astype(np.float64)
    totals[totals == 0] = np.nan
    out = np.nan_to_num(tpm.astype(np.float64) / totals * 1e6, nan=0.0)
    return out.astype(np.float32)


def rank_z(values: np.ndarray) -> np.ndarray:
    """Within-sample rank-normal scores on the columns that are passed in."""
    out = np.empty(values.shape, dtype=np.float64)
    n_genes = values.shape[1]
    for start in range(0, values.shape[0], 500):
        block = values[start:start + 500]
        ranks = rankdata(block, method="average", axis=1)
        out[start:start + block.shape[0]] = norm.ppf((ranks - 0.5) / n_genes)
    return out


def rank_z_take(values: np.ndarray, column_sets: list[np.ndarray]) -> list[np.ndarray]:
    """Rank on all columns, then keep each column set. Ranks use |G| = values.shape[1]."""
    n = values.shape[0]
    n_genes = values.shape[1]
    outs = [np.empty((n, len(idx)), dtype=np.float32) for idx in column_sets]
    for start in range(0, n, 400):
        block = values[start:start + 400]
        ranks = rankdata(block, method="average", axis=1)
        z = norm.ppf((ranks - 0.5) / n_genes).astype(np.float32)
        stop = start + block.shape[0]
        for out, idx in zip(outs, column_sets):
            out[start:stop] = z[:, idx]
    return outs


def score_ks(tpm: np.ndarray, offsets, indices) -> np.ndarray:
    parts = []
    for start in range(0, tpm.shape[0], 2000):
        block = np.ascontiguousarray(tpm[start:start + 2000], dtype=np.float32)
        right, tie = right_and_tie(block)
        parts.append(ks_matrix(right, tie, offsets, indices))
        print("ks", start + block.shape[0], "/", tpm.shape[0], flush=True)
    return np.vstack(parts) if parts else np.empty((0, int(offsets.shape[0] - 1)), dtype=np.float32)


def load_gene_pack(paths):
    genes = pd.read_csv(paths["data_processed"] / "genes/G_symbols.txt", header=None)[0].astype(str).tolist()
    packed = np.load(paths["data_processed"] / "genes/gene_sets.npz", allow_pickle=True)
    offsets = packed["offsets"]
    indices = packed["indices"]
    names = [str(x) for x in packed["names"].tolist()]
    b0 = pd.read_csv(paths["results"] / "B0_genes.txt")["symbol"].astype(str).tolist()
    gene_index = {gene: i for i, gene in enumerate(genes)}
    b0_idx = np.array([gene_index[gene] for gene in b0], dtype=np.int32)
    return genes, b0_idx, offsets, indices, names


def train_frame(paths) -> pd.DataFrame:
    split = pd.read_csv(paths["config"] / "split_tcga.tsv", sep="\t", dtype=str)
    train = split.loc[split["split"] == "TCGA-train"].copy()
    return train.reset_index(drop=True)


def test_frame(paths) -> pd.DataFrame:
    split = pd.read_csv(paths["config"] / "split_tcga.tsv", sep="\t", dtype=str)
    return split.loc[split["split"] == "TCGA-test"].reset_index(drop=True)


def simulation_tumor_index(test: pd.DataFrame) -> np.ndarray:
    """Same 1,000 TCGA-test rows as stage 1: stratified split, then sorted positions."""
    labels = test["learning_label"].to_numpy()
    take, _ = train_test_split(
        np.arange(len(test)), train_size=min(1000, len(test)), stratify=labels, random_state=SEED
    )
    return np.sort(take)


def if_tumor_index(train: pd.DataFrame) -> np.ndarray:
    labels = train["learning_label"].to_numpy()
    take, _ = train_test_split(
        np.arange(len(train)), train_size=min(IF_N, len(train)), stratify=labels, random_state=SEED
    )
    return np.sort(take)


def gtex_ref_ids(paths, tissues: list[str]) -> dict[str, list[str]]:
    gtex = pd.read_csv(paths["config"] / "split_gtex.tsv", sep="\t", dtype=str)
    ref = gtex.loc[gtex["split"] == "GTEx-ref"]
    out = {}
    for tissue in tissues:
        out[tissue] = ref.loc[ref["tissue"] == tissue, "sample"].tolist()
    return out


def gtex_sim_ids(paths, tissues: list[str]) -> dict[str, list[str]]:
    gtex = pd.read_csv(paths["config"] / "split_gtex.tsv", sep="\t", dtype=str)
    sim = gtex.loc[gtex["split"] == "GTEx-sim"]
    return {tissue: sim.loc[sim["tissue"] == tissue, "sample"].tolist() for tissue in tissues}


def site_native_organs(site: str) -> set[str]:
    return set(SITE_NATIVE.get(site, SITE_NATIVE["other"])[2])


def learning_labels_sorted(train: pd.DataFrame) -> list[str]:
    return sorted(train["learning_label"].unique())


def organ_from_project_proba(proba: np.ndarray, classes: list[str]):
    from metrics import ORGANS, organ_probability, predict_from_proba

    organ = organ_probability(proba, classes)
    pred = predict_from_proba(organ, ORGANS)
    return organ, pred
