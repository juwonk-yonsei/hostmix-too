"""Stage-5 helpers. This module does not set a thread cap."""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.metrics import log_loss

from mapping_rules import PROJECT_TO_ORGAN
from metrics import ORGANS, macro_f1, percentile_ci, topk_hit
from stage2_common import (
    N_MIX,
    RHO_HIGH,
    RHO_LOW,
    SIM_RHOS,
    SITE_TO_MODEL,
    gtex_ref_ids,
    gtex_sim_ids,
    load_gene_pack,
    load_tpm,
    rank_z,
    simulation_tumor_index,
    site_native_organs,
    sum_1e6,
    test_frame,
    train_frame,
)
from stage2_fit import predict_bundle
from util import SEED, load_paths

STAGE5_SEED = 20261005
GI = ("Stomach", "Esophagus", "Pancreas", "Biliary", "Colorectal")
GI_SET = set(GI)
POOL10 = [
    "Adipose - Subcutaneous",
    "Adipose - Visceral (Omentum)",
    "Adrenal Gland",
    "Brain - Cortex",
    "Liver",
    "Lung",
    "Muscle - Skeletal",
    "Skin - Not Sun Exposed (Suprapubic)",
    "Spleen",
    "Whole Blood",
]
POOL12 = [
    "Stomach",
    "Colon - Transverse",
    "Pancreas",
    "Kidney - Cortex",
    "Breast - Mammary Tissue",
    "Ovary",
    "Prostate",
    "Thyroid",
    "Esophagus - Mucosa",
    "Minor Salivary Gland",
    "Cervix - Ectocervix",
    "Bladder",
]
# Overrides from instruction 5 §2.2. Other tissues follow the single site_native row.
TISSUE_TO_SITE = {
    "Adipose - Subcutaneous": "soft_tissue",
    "Muscle - Skeletal": "soft_tissue",
    "Adipose - Visceral (Omentum)": "omentum",
    "Skin - Not Sun Exposed (Suprapubic)": "skin",
    "Spleen": "lymph_node",
    "Whole Blood": "bone_marrow",
    "Liver": "liver",
    "Lung": "lung",
    "Adrenal Gland": "adrenal",
    "Brain - Cortex": "brain",
    "Stomach": "stomach",
    "Colon - Transverse": "colon_rectum",
    "Pancreas": "pancreas",
    "Kidney - Cortex": "kidney",
    "Breast - Mammary Tissue": "breast",
    "Ovary": "ovary",
    "Prostate": "prostate",
    "Thyroid": "thyroid",
    "Esophagus - Mucosa": "esophagus",
    "Minor Salivary Gland": "head_neck",
    "Cervix - Ectocervix": "cervix",
    "Bladder": "bladder",
}
STAGE2_SIM_TISSUES = ["Liver", "Lung", "Adipose - Subcutaneous", "Spleen", "Whole Blood"]
NAMED_MET_SITES = {
    "Regional Lymph Node",
    "Distant Metastasis",
    "Regional Cutaneous or Subcutaneous Tissue (includes satellite and in-transit metastasis)",
}
LABEL_HOST = {
    "LIHC": ["Liver"],
    "DLBC": ["Spleen"],
    "LUAD": ["Lung"],
    "LUSC": ["Lung"],
    "GBM": ["Brain - Cortex"],
    "LGG": ["Brain - Cortex"],
    "SARC": ["Adipose - Subcutaneous", "Muscle - Skeletal"],
    "LAML": ["Whole Blood"],
    "SKCM": ["Skin - Not Sun Exposed (Suprapubic)"],
    "ACC": ["Adrenal Gland"],
    "PCPG": ["Adrenal Gland"],
}


def variant_names() -> list[str]:
    names = [
        "SA-m0", "SA-m1", "SA-m8", "SA-C0.01", "SA-C0.1", "SA-r30", "SA-r50",
        "SA-pool3", "SA-pool22", "BASE-MLP", "SA-MLP",
    ]
    names.extend(f"SA-LOHO-{tissue}" for tissue in POOL10)
    return sorted(names)


def variant_seeds(names: list[str] | None = None) -> dict[str, int]:
    ordered = list(names) if names is not None else variant_names()
    if ordered != sorted(ordered):
        raise SystemExit("variant seeds require sorted names")
    rng = np.random.default_rng(STAGE5_SEED)
    return {name: int(rng.integers(0, 2**31 - 1)) for name in ordered}


def c_for(name: str) -> float | None:
    if name in ("BASE-MLP", "SA-MLP"):
        return None
    if name == "SA-m0":
        return 0.15
    if name == "SA-m1":
        return 0.075
    if name == "SA-m8":
        return 0.03 * 5.0 / 9.0
    if name == "SA-C0.01":
        return 0.01
    if name == "SA-C0.1":
        return 0.1
    return 0.03


def _host(rng, tissues, ids_by_tissue, pool, pooled: bool):
    if pooled:
        return pool[int(rng.integers(0, len(pool)))]
    tissue = tissues[int(rng.integers(0, len(tissues)))]
    group = ids_by_tissue[tissue]
    return tissue, group[int(rng.integers(0, len(group)))]


def draw_mixes(tumors, tissues, ids_by_tissue, rng, pooled: bool, rho_fixed, n_mix: int = N_MIX, rho_low: float = RHO_LOW):
    """Same draw as scripts/08_train.py draw_mixes when rho_low is 0.15 and n_mix is 4."""
    records = []
    pool = [(tissue, sample) for tissue in tissues for sample in ids_by_tissue[tissue]]
    for tumor in tumors:
        if rho_fixed is None:
            for _ in range(n_mix):
                tissue, host = _host(rng, tissues, ids_by_tissue, pool, pooled)
                records.append({
                    "tumor": tumor,
                    "host": host,
                    "tissue": tissue,
                    "rho": float(rng.uniform(rho_low, RHO_HIGH)),
                })
        elif pooled:
            tissue, host = _host(rng, tissues, ids_by_tissue, pool, True)
            records.append({"tumor": tumor, "host": host, "tissue": tissue, "rho": float(rho_fixed)})
        else:
            for tissue in tissues:
                group = ids_by_tissue[tissue]
                host = group[int(rng.integers(0, len(group)))]
                records.append({"tumor": tumor, "host": host, "tissue": tissue, "rho": float(rho_fixed)})
    return records


def materialize(tumor_norm, train_pos, bank, records) -> np.ndarray:
    out = np.empty((len(records), tumor_norm.shape[1]), dtype=np.float32)
    for start in range(0, len(records), 2000):
        chunk = records[start:start + 2000]
        rho = np.array([row["rho"] for row in chunk], dtype=np.float32)[:, None]
        t_idx = np.array([train_pos[row["tumor"]] for row in chunk])
        host = np.stack([bank[row["host"]] for row in chunk])
        out[start:start + len(chunk)] = rho * tumor_norm[t_idx] + (np.float32(1.0) - rho) * host
    return out


def rank_b0(values: np.ndarray, b0_idx: np.ndarray, chunk: int = 400) -> np.ndarray:
    n = values.shape[0]
    out = np.empty((n, len(b0_idx)), dtype=np.float32)
    for start in range(0, n, chunk):
        block = values[start:start + chunk]
        scored = rank_z(block)
        out[start:start + block.shape[0]] = scored[:, b0_idx].astype(np.float32)
    return out


def load_references(path: Path):
    packed = np.load(path, allow_pickle=True)
    h = {key[3:]: packed[key] for key in packed.files if key.startswith("h__")}
    genes = [str(x) for x in packed["genes"].tolist()]
    return packed["mu"], h, packed["keep_mt"], genes


def predict_organ(bundle: dict, x_b0: np.ndarray):
    organ, pred, _raw, _proba, _classes = predict_bundle(bundle, x_b0)
    return organ, pred


def organ_of_labels(labels: np.ndarray) -> np.ndarray:
    return np.array([PROJECT_TO_ORGAN[label] for label in labels], dtype=object)


def substitute_site(site: str) -> bool:
    return site not in SITE_TO_MODEL


def error_parts(truth: np.ndarray, pred: np.ndarray, natives: list[set[str]]) -> dict:
    n = len(truth)
    host = np.zeros(n, dtype=bool)
    for i, native in enumerate(natives):
        if native and truth[i] not in native and pred[i] in native:
            host[i] = True
    wrong = pred != truth
    gi = np.array([
        wrong[i] and not host[i] and truth[i] in GI_SET and pred[i] in GI_SET
        for i in range(n)
    ])
    rest = wrong & ~host & ~gi
    risk = np.array([bool(native) and truth[i] not in native for i, native in enumerate(natives)])
    native_truth = np.array([bool(native) and truth[i] in native for i, native in enumerate(natives)])
    return {
        "host": host,
        "wrong": wrong,
        "gi": gi,
        "rest": rest,
        "risk": risk,
        "native_truth": native_truth,
        "n_gi": int(gi.sum()),
        "n_host": int((wrong & host).sum()),
        "n_rest": int(rest.sum()),
        "n_error": int(wrong.sum()),
    }


def metric_row(truth, pred, proba, natives) -> dict:
    parts = error_parts(truth, pred, natives)
    risk = parts["risk"]
    native_truth = parts["native_truth"]
    pulled = risk & np.array([pred[i] in natives[i] for i in range(len(pred))])
    stomach = pred == "Stomach"
    return {
        "n": int(len(truth)),
        "top1": float(np.mean(pred == truth)) if len(truth) else None,
        "top3": float(np.mean(topk_hit(proba, ORGANS, truth, 3, pred != "NA"))) if len(truth) else None,
        "macro_f1": macro_f1(truth, pred) if len(truth) else None,
        "n_at_risk": int(risk.sum()),
        "host_rate": float(pulled.sum() / risk.sum()) if risk.any() else None,
        "n_native_truth": int(native_truth.sum()),
        "native_truth_top1": float(np.mean(pred[native_truth] == truth[native_truth])) if native_truth.any() else None,
        "n_gi_internal": parts["n_gi"],
        "n_pred_stomach": int(stomach.sum()),
        "n_pred_stomach_true": int((stomach & (truth == "Stomach")).sum()),
        "n_sarc_to_brain": int(((truth == "Sarcoma") & (pred == "Brain")).sum()),
    }


def cluster_indices(cluster_ids: np.ndarray, n_boot: int, seed: int) -> list[np.ndarray]:
    codes, inverse = np.unique(cluster_ids.astype(str), return_inverse=True)
    groups = [np.flatnonzero(inverse == i) for i in range(len(codes))]
    rng = np.random.default_rng(seed)
    draws = []
    n_groups = len(groups)
    for _ in range(n_boot):
        picked = rng.integers(0, n_groups, size=n_groups)
        draws.append(np.concatenate([groups[i] for i in picked]))
    return draws


def paired_mean(new_flag: np.ndarray, base_flag: np.ndarray, cluster_ids: np.ndarray, seed: int = SEED) -> dict:
    new_flag = np.asarray(new_flag, dtype=float)
    base_flag = np.asarray(base_flag, dtype=float)
    point = float(new_flag.mean() - base_flag.mean()) if len(new_flag) else None
    n_new_only = int(np.sum((new_flag == 1) & (base_flag == 0)))
    n_base_only = int(np.sum((new_flag == 0) & (base_flag == 1)))
    diffs = []
    undefined = 0
    if len(new_flag):
        for draw in cluster_indices(cluster_ids, 2000, seed):
            if draw.size == 0:
                undefined += 1
                continue
            diffs.append(float(new_flag[draw].mean() - base_flag[draw].mean()))
    low, high = percentile_ci(diffs)
    return {
        "diff": point,
        "ci_low": low,
        "ci_high": high,
        "n_new_only": n_new_only,
        "n_base_only": n_base_only,
        "n_boot_used": len(diffs),
        "n_boot_undefined": undefined,
    }


def mcnemar_two_sided(base_ok: np.ndarray, new_ok: np.ndarray) -> dict:
    n_base_only = int(np.sum(base_ok & ~new_ok))
    n_new_only = int(np.sum(~base_ok & new_ok))
    n_disc = n_base_only + n_new_only
    pvalue = None if n_disc == 0 else float(binomtest(n_new_only, n_disc, 0.5, alternative="two-sided").pvalue)
    return {"n_base_only": n_base_only, "n_new_only": n_new_only, "p_two_sided": pvalue}


def holdout_log_loss(clf, x_hold, y_hold, labels) -> float:
    proba = clf.predict_proba(x_hold)
    class_order = [str(c) for c in clf.classes_]
    columns = np.zeros((len(y_hold), len(labels)), dtype=np.float64)
    for j, label in enumerate(labels):
        if label in class_order:
            columns[:, j] = proba[:, class_order.index(label)]
    return float(log_loss(y_hold, columns, labels=labels))


def write_disk(path: Path) -> None:
    import subprocess
    text = subprocess.check_output(["df", "-h"], text=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def stage5_dir(paths) -> Path:
    root = paths["results"] / "stage5"
    for name in ("tcga_test", "tcga_met", "sim_ext", "variants", "mechanism", "posthoc", "aux"):
        (root / name).mkdir(parents=True, exist_ok=True)
    return root


def dump_json(path: Path, payload) -> None:
    path.write_text(json.dumps(payload, indent=2))


def load_b0_symbols(paths) -> list[str]:
    return pd.read_csv(paths["results"] / "B0_genes.txt")["symbol"].astype(str).tolist()


def effective_coef(bundle: dict) -> tuple[np.ndarray, list[str]]:
    coef = np.asarray(bundle["clf"].coef_, dtype=np.float64)
    classes = [str(c) for c in bundle["clf"].classes_]
    if bundle["scaler"] is None:
        scale = np.ones(coef.shape[1], dtype=np.float64)
    else:
        scale = np.asarray(bundle["scaler"].scale_, dtype=np.float64)
    centered = coef / scale
    centered = centered - centered.mean(axis=0, keepdims=True)
    return centered, classes
