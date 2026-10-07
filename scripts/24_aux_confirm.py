#!/usr/bin/env python3
"""Auxiliary confirmation. No arguments and --real without unlock read no data files."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "results" / "stage6" / "opened_paths.txt"
SYN = ROOT / "results" / "stage6" / "synthetic"


def fail_locked() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("")
    raise SystemExit(2)


argv = sys.argv[1:]
PLANTED_DIR = None
if argv == ["--sets-only"]:
    MODE = "sets-only"
elif len(argv) == 2 and argv[0] == "--planted":
    MODE = "planted"
    PLANTED_DIR = Path(argv[1])
elif argv in (["--synthetic"], ["--real"]):
    MODE = argv[0][2:]
else:
    fail_locked()
if MODE == "real" and not (ROOT / "config" / "unlock_A7.md").exists():
    fail_locked()

import hashlib
import json
import math
import subprocess

import joblib
import numpy as np
import pandas as pd
import yaml
from scipy.stats import binom
from sklearn.model_selection import train_test_split

sys.path.insert(0, str(ROOT / "scripts"))

from mapping_rules import PROJECT_TO_ORGAN, SITE_NATIVE
from metrics import ORGANS, macro_f1, mask_native, percentile_ci, topk_hit
from stage2_common import (
    SEED,
    SITE_TO_MODEL,
    gtex_sim_ids,
    load_gene_pack,
    load_tpm,
    rank_z,
    score_ks,
    simulation_tumor_index,
    site_native_organs,
    sum_1e6,
    test_frame,
)
from stage2_fit import ld_apply, predict_bundle
from stage5_lib import GI_SET, load_references, rank_b0

POOL_OUT = {
    "kidney", "pancreas", "bladder", "thyroid", "ovary", "breast", "stomach",
    "colon_rectum", "prostate", "esophagus", "head_neck", "cervix",
}
MIX_SITES = ["liver", "lung", "lymph_node", "bone", "kidney", "soft_tissue"]
HOST_TISSUE = {
    "liver": "Liver",
    "lung": "Lung",
    "lymph_node": "Spleen",
    "bone": "Whole Blood",
    "kidney": "Kidney - Cortex",
    "soft_tissue": "Adipose - Subcutaneous",
}
PRIMARY_OVERRIDE = {
    "Melanoma": "skin",
    "Lymphoid": "lymph_node",
    "Myeloid": "bone_marrow",
    "Mesothelioma": "pleura",
}
METHODS = [
    "BASE-Z", "SA-Z", "SA-G", "SA-pool22", "SC-Z", "LD-Z", "NC-Z", "M1-Z",
    "BASE-K", "SA-K", "V0-K", "SA-MLP",
]
N_BOOT = 2000
SA_G_TAU = 0.02


def project_paths() -> dict:
    raw = yaml.safe_load((ROOT / "config" / "paths.yaml").read_text())
    out = {"root": ROOT}
    for key, value in raw.items():
        if key in ("seed", "n_threads"):
            continue
        out[key] = ROOT / value
    return out


def reset_log() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("")


def is_real_aux(path: Path) -> bool:
    path = path if path.is_absolute() else ROOT / path
    try:
        rel = path.resolve().relative_to(ROOT.resolve())
    except ValueError:
        return False
    parts = rel.parts
    if parts[:2] == ("data", "aux"):
        return True
    return len(parts) >= 3 and parts[:3] == ("data", "processed", "aux")


def note(path: Path, purpose: str) -> None:
    path = Path(path)
    if MODE == "synthetic" and purpose != "hash" and is_real_aux(path):
        with LOG.open("a") as handle:
            handle.write(f"blocked\t{path}\n")
        raise SystemExit(f"synthetic run tried to open a real auxiliary cohort file: {path}")
    with LOG.open("a") as handle:
        handle.write(f"{purpose}\t{path}\n")


def sha256_file(path: Path) -> str:
    note(path, "hash")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(1 << 20)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def check_freeze() -> None:
    for name in ("frozen_models_A6.tsv", "frozen_inputs_A6.tsv", "frozen_code_A7.tsv"):
        table_path = ROOT / "config" / name
        note(table_path, "hash")
        table = pd.read_csv(table_path, sep="\t", dtype=str)
        for row in table.itertuples(index=False):
            path = ROOT / row.path
            digest = sha256_file(path)
            if digest != row.sha256:
                raise SystemExit(f"freeze hash mismatch {row.path}")


def check_unlock() -> None:
    unlock = ROOT / "config" / "unlock_A7.md"
    note(unlock, "config")
    text = unlock.read_text()
    freeze = subprocess.check_output(
        ["git", "log", "-1", "--format=%H", "--", "config/frozen_models_A6.tsv"],
        cwd=ROOT, text=True,
    ).strip()
    if freeze not in text:
        raise SystemExit(2)


def pack_offsets(member_lists):
    offsets = np.zeros(len(member_lists) + 1, dtype=np.int64)
    chunks = []
    for i, members in enumerate(member_lists):
        chunks.append(np.asarray(members, dtype=np.int32))
        offsets[i + 1] = offsets[i] + chunks[-1].size
    indices = np.concatenate(chunks) if chunks else np.empty(0, dtype=np.int32)
    return offsets, indices


def z_and_k(block, present, genes, b0, offsets, indices):
    present_index = {gene: i for i, gene in enumerate(present)}
    z_all = rank_z(block)
    z = np.zeros((block.shape[0], len(b0)), dtype=np.float32)
    for j, symbol in enumerate(b0):
        if symbol in present_index:
            z[:, j] = z_all[:, present_index[symbol]].astype(np.float32)
    lists = []
    set_ok = []
    n_sets = int(offsets.shape[0] - 1)
    for g in range(n_sets):
        members = []
        for gene_i in indices[offsets[g]:offsets[g + 1]]:
            symbol = genes[int(gene_i)]
            if symbol in present_index:
                members.append(present_index[symbol])
        if len(members) < 5:
            set_ok.append(False)
            lists.append(np.array([0], dtype=np.int32))
        else:
            set_ok.append(True)
            lists.append(np.array(members, dtype=np.int32))
    valid = [members for members, ok in zip(lists, set_ok) if ok]
    packed_off, packed_idx = pack_offsets(valid)
    scored = score_ks(block.astype(np.float32), packed_off, packed_idx)
    k = np.zeros((block.shape[0], n_sets), dtype=np.float32)
    cursor = 0
    for g, ok in enumerate(set_ok):
        if ok:
            k[:, g] = scored[:, cursor]
            cursor += 1
    return z, k


def primary_sites() -> tuple[dict[str, str], list[dict]]:
    found: dict[str, list[str]] = {}
    for site, (_tissues, _projects, organs) in SITE_NATIVE.items():
        for organ in organs:
            found.setdefault(organ, []).append(site)
    chosen = {}
    skipped = []
    for organ in sorted(set(PROJECT_TO_ORGAN.values())):
        sites = found.get(organ, [])
        if organ in PRIMARY_OVERRIDE:
            chosen[organ] = PRIMARY_OVERRIDE[organ]
        elif len(sites) == 1:
            chosen[organ] = sites[0]
        else:
            skipped.append({"organ": organ, "sites": sites})
    return chosen, skipped


def native_strings(path: Path) -> dict[str, str]:
    note(path, "config")
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return dict(zip(table["standard_site"], table["native_organs"]))


def select_patients(labels: pd.DataFrame) -> pd.DataFrame:
    labels = labels.copy()
    risk_n = labels.loc[labels["in_risk"]].groupby("cohort").size()
    confirm = {}
    for cohort, layer in labels.groupby("cohort")["layer"].first().items():
        enough = int(risk_n.get(cohort, 0)) >= 10
        confirm[cohort] = layer == "RNA-seq" and enough
    labels["include_confirm"] = [bool(confirm[c]) and (not ex) for c, ex in zip(labels["cohort"], labels["excluded"])]
    labels["layer_test"] = [int(risk_n.get(c, 0)) >= 10 and (not ex) for c, ex in zip(labels["cohort"], labels["excluded"])]
    for column, mask_col in (
        ("selected_eval", "in_eval"),
        ("selected_risk", "in_risk"),
        ("selected_native", "in_native"),
        ("selected_pool_out", "in_pool_out_risk"),
    ):
        labels[column] = False
        part = labels.loc[labels[mask_col]].copy()
        part["library_rank"] = np.where(part["library"].eq("polyA"), 0, 1)
        part = part.sort_values(["cohort", "patient_id", "library_rank", "sample_id"])
        keep = part.drop_duplicates(["cohort", "patient_id"]).index
        labels.loc[keep, column] = True
    return labels


def build_synthetic(paths) -> None:
    genes, b0_idx, offsets, indices, set_names = load_gene_pack(paths)
    b0 = [genes[int(i)] for i in b0_idx]
    note(paths["data_processed"] / "toil/tpm_G.h5", "synthetic-source")
    samples, tpm = load_tpm(paths["data_processed"] / "toil/tpm_G.h5")
    row = {sample: i for i, sample in enumerate(samples)}
    test = test_frame(paths)
    take = simulation_tumor_index(test)
    tumors = test.iloc[take].reset_index(drop=True)
    if len(tumors) != 1000:
        raise SystemExit(f"TCGA-test draw is {len(tumors)}")
    organs = [PROJECT_TO_ORGAN[label] for label in tumors["learning_label"]]
    tumor_ids = tumors["sample"].astype(str).tolist()
    tumor_mat = sum_1e6(tpm[[row[sample] for sample in tumor_ids]])
    host_ids = gtex_sim_ids(paths, list(HOST_TISSUE.values()))
    host_bank = {}
    for tissue, ids in host_ids.items():
        if not ids:
            raise SystemExit(f"GTEx-sim has no samples for {tissue}")
        host_bank[tissue] = sum_1e6(tpm[[row[sample] for sample in ids]])
    del tpm
    primary, skipped = primary_sites()
    natives = native_strings(ROOT / "config" / "mappings" / "site_native.tsv")
    rng = np.random.default_rng(20261006)
    order = rng.permutation(1000)
    splits = {
        "syn_rna_1": order[:334],
        "syn_rna_2": order[334:667],
        "syn_array": order[667:],
    }
    n_drop = int(round(0.3 * len(genes)))
    dropped = np.sort(rng.choice(len(genes), size=n_drop, replace=False))
    keep_gene = np.ones(len(genes), dtype=bool)
    keep_gene[dropped] = False
    layers = {"syn_rna_1": "RNA-seq", "syn_rna_2": "RNA-seq", "syn_array": "마이크로어레이"}
    rows = []
    matrices = {name: [] for name in splits}
    ids_by = {name: [] for name in splits}
    overlap_from = [tumor_ids[int(i)] for i in splits["syn_rna_1"][:40]]
    for cohort, indexes in splits.items():
        gene_mask = keep_gene if cohort == "syn_array" else np.ones(len(genes), dtype=bool)
        present = [gene for gene, flag in zip(genes, gene_mask) if flag]
        present_idx = np.flatnonzero(gene_mask)
        for local_i, tumor_i in enumerate(indexes):
            tumor_i = int(tumor_i)
            organ = organs[tumor_i]
            patient = tumor_ids[tumor_i]
            if cohort == "syn_rna_2" and local_i < 40:
                patient = overlap_from[local_i]
            pure = tumor_mat[tumor_i]
            for site in MIX_SITES:
                tissue = HOST_TISSUE[site]
                host_i = int(rng.integers(0, len(host_ids[tissue])))
                rho = float(rng.uniform(0.3, 0.9))
                mixed = np.float32(rho) * pure + np.float32(1.0 - rho) * host_bank[tissue][host_i]
                mixed = sum_1e6(mixed[None, present_idx])[0]
                sample_id = f"{tumor_ids[tumor_i]}|{site}|mix"
                matrices[cohort].append(mixed)
                ids_by[cohort].append(sample_id)
                rows.append(label_row(cohort, layers[cohort], sample_id, patient, organ, site, natives))
            site = primary.get(organ)
            if site is None:
                continue
            pure_part = sum_1e6(pure[None, present_idx])[0]
            sample_id = f"{tumor_ids[tumor_i]}|{site}|primary"
            matrices[cohort].append(pure_part)
            ids_by[cohort].append(sample_id)
            rows.append(label_row(cohort, layers[cohort], sample_id, patient, organ, site, natives))
        block = np.vstack(matrices[cohort]).astype(np.float32)
        save_cohort(SYN, cohort, block, present, ids_by[cohort], genes, b0, offsets, indices)
        del block, matrices[cohort]
    labels = select_patients(pd.DataFrame(rows))
    label_path = SYN / "config" / "aux_eval_labels.tsv"
    label_path.parent.mkdir(parents=True, exist_ok=True)
    labels.to_csv(label_path, sep="\t", index=False)
    design = {
        "seed": 20261006,
        "n_tumors": 1000,
        "cohorts": {name: int(len(idx)) for name, idx in splits.items()},
        "n_genes_dropped_array": int(n_drop),
        "host_tissue": HOST_TISSUE,
        "primary_override": PRIMARY_OVERRIDE,
        "skipped_primary_organs": skipped,
        "patient_overlap": "first 40 syn_rna_2 tumors use the first 40 syn_rna_1 patient ids",
        "soft_tissue_host": "Adipose - Subcutaneous",
    }
    (SYN / "design.json").write_text(json.dumps(design, indent=2, ensure_ascii=False))


def label_row(cohort, layer, sample_id, patient, organ, site, natives) -> dict:
    native = site_native_organs(site)
    in_eval = True
    in_risk = bool(native) and organ not in native
    in_native = bool(native) and organ in native
    return {
        "cohort": cohort,
        "layer": layer,
        "sample_id": sample_id,
        "patient_id": patient,
        "library": "",
        "organ": organ,
        "organ_rule": "synthetic",
        "raw_site": site,
        "standard_site": site,
        "native_organs": natives.get(site, ""),
        "in_eval": in_eval,
        "in_risk": in_risk,
        "in_native": in_native,
        "pool_out_site": site in POOL_OUT,
        "in_pool_out_risk": site in POOL_OUT and in_risk,
        "excluded": False,
        "exclude_reason": "",
    }


def save_cohort(base, cohort, block, present, ids, genes, b0, offsets, indices) -> None:
    feat = base / "data" / "processed" / "aux" / cohort
    feat.mkdir(parents=True, exist_ok=True)
    z, k = z_and_k(block, present, genes, b0, offsets, indices)
    np.save(feat / "Z.npy", z)
    np.save(feat / "K.npy", k)
    (feat / "sample_ids.txt").write_text("\n".join(ids) + "\n")
    frame = pd.DataFrame(block, columns=present)
    frame.insert(0, "sample_id", ids)
    ld_path = base / "data" / "aux" / "ld_input" / f"{cohort}.parquet"
    ld_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(ld_path, index=False)


def feature_cohort(cohort: str, library: str) -> str:
    if cohort == "prad_su2c_2019":
        if library == "polyA":
            return "prad_su2c_2019__fpkm_polya"
        if library == "capture":
            return "prad_su2c_2019__fpkm_capture"
        raise SystemExit(f"SU2C library is {library}")
    return cohort


def load_matrix(base: Path, cohort: str, library: str, kind: str) -> tuple[np.ndarray, list[str]]:
    folder = feature_cohort(cohort, library)
    path = base / "data" / "processed" / "aux" / folder / f"{kind}.npy"
    ids_path = base / "data" / "processed" / "aux" / folder / "sample_ids.txt"
    note(path, "feature")
    note(ids_path, "feature")
    values = np.load(path)
    ids = [line for line in ids_path.read_text().splitlines() if line]
    return values, ids


def load_ld(base: Path, cohort: str) -> pd.DataFrame:
    name = "prad_su2c_2019" if cohort == "prad_su2c_2019" else cohort
    path = base / "data" / "aux" / "ld_input" / f"{name}.parquet"
    note(path, "ld")
    return pd.read_parquet(path)


def predict_methods(z, k, sites, ld_z, beta, models, v0_means):
    base_p, base_pred, _raw, _proba, _classes = predict_bundle(models["BASE_Z"], z)
    sa_p, sa_pred, *_rest = predict_bundle(models["SA_Z"], z)
    nc_p, nc_pred, *_rest = predict_bundle(models["NC_Z"], z)
    pool_p, pool_pred, *_rest = predict_bundle(models["SA-pool22"], z)
    mlp_p, mlp_pred, *_rest = predict_bundle(models["SA-MLP"], z)
    base_k_p, base_k_pred, *_rest = predict_bundle(models["BASE_K"], k)
    sa_k_p, sa_k_pred, *_rest = predict_bundle(models["SA_K"], k)
    natives = [sorted(site_native_organs(site)) for site in sites]
    m1_p, m1_pred = mask_native(base_p, natives)
    use_base = np.isnan(beta) | (beta <= SA_G_TAU)
    sag_pred = np.where(use_base, base_pred, sa_pred)
    sag_p = np.where(use_base[:, None], base_p, sa_p)
    sc_p = np.zeros_like(base_p)
    sc_pred = np.empty(len(sites), dtype=object)
    ld_p = np.zeros_like(base_p)
    ld_pred = np.empty(len(sites), dtype=object)
    k_v0 = k.astype(np.float64).copy()
    for site in sorted(set(sites)):
        idx = np.array([i for i, value in enumerate(sites) if value == site])
        model_site = SITE_TO_MODEL.get(site)
        if model_site is None:
            sc_p[idx] = sa_p[idx]
            sc_pred[idx] = sa_pred[idx]
            ld_p[idx] = base_p[idx]
            ld_pred[idx] = base_pred[idx]
        else:
            part_p, part_pred, *_rest = predict_bundle(models[f"SC_{model_site}_Z"], z[idx])
            sc_p[idx] = part_p
            sc_pred[idx] = part_pred
            part_p, part_pred, *_rest = predict_bundle(models["BASE_Z"], ld_z[idx])
            ld_p[idx] = part_p
            ld_pred[idx] = part_pred
        host = v0_means.get(site)
        if host is not None:
            k_v0[idx] = (k[idx] - 0.3 * host) / 0.7
    v0_p, v0_pred, *_rest = predict_bundle(models["BASE_K"], k_v0.astype(np.float32))
    return {
        "BASE-Z": (base_p, base_pred),
        "SA-Z": (sa_p, sa_pred),
        "SA-G": (sag_p, sag_pred),
        "SA-pool22": (pool_p, pool_pred),
        "SC-Z": (sc_p, sc_pred),
        "LD-Z": (ld_p, ld_pred),
        "NC-Z": (nc_p, nc_pred),
        "M1-Z": (m1_p, m1_pred),
        "BASE-K": (base_k_p, base_k_pred),
        "SA-K": (sa_k_p, sa_k_pred),
        "V0-K": (v0_p, v0_pred),
        "SA-MLP": (mlp_p, mlp_pred),
    }


def ld_for_rows(frame: pd.DataFrame, genes, h_all, mu, keep_mt, b0_idx) -> tuple[np.ndarray, np.ndarray]:
    present = [col for col in frame.columns if col != "sample_id"]
    present_index = {gene: i for i, gene in enumerate(present)}
    gene_index = {gene: i for i, gene in enumerate(genes)}
    use = []
    for gene, flag in zip(genes, keep_mt):
        if flag and gene in present_index:
            use.append(gene)
    if not use:
        beta = np.full(len(frame), np.nan)
        return np.zeros((len(frame), len(b0_idx)), dtype=np.float32), beta
    idx = np.array([gene_index[gene] for gene in use], dtype=np.int32)
    cols = [present_index[gene] for gene in use]
    values = frame[present].to_numpy(dtype=np.float32)
    x = sum_1e6(values[:, cols])
    h = sum_1e6(h_all[idx][None, :])[0]
    mu_sub = sum_1e6(mu[:, idx])
    cleaned, beta = ld_apply(x, h, mu_sub, np.ones(len(use), dtype=bool))
    full = np.zeros((len(frame), len(genes)), dtype=np.float32)
    full[:, idx] = cleaned
    return rank_b0(full, b0_idx), beta


def score_base(base: Path, labels: pd.DataFrame, paths) -> pd.DataFrame:
    genes, b0_idx, _offsets, _indices, set_names = load_gene_pack(paths)
    ref_path = paths["results"] / "stage3" / "references.npz"
    note(ref_path, "model")
    mu, h_by_site, keep_mt, ref_genes = load_references(ref_path)
    if ref_genes != genes:
        raise SystemExit("reference gene order differs from G")
    models = {}
    model_dir = paths["results"] / "stage2" / "models"
    for name in ("BASE_Z", "SA_Z", "NC_Z", "BASE_K", "SA_K"):
        path = model_dir / f"{name}.joblib"
        note(path, "model")
        models[name] = joblib.load(path)
    for site in sorted(set(SITE_TO_MODEL.values())):
        path = model_dir / f"SC_{site}_Z.joblib"
        note(path, "model")
        models[f"SC_{site}_Z"] = joblib.load(path)
    var = paths["results"] / "stage5" / "variants" / "models"
    for name in ("SA-pool22", "SA-MLP"):
        path = var / f"{name}.joblib"
        note(path, "model")
        models[name] = joblib.load(path)
    v0_means = v0_site_means(paths, set_names)
    # Empty label fields are NaN under read_csv even with dtype=str. Excluded rows have an empty site.
    for column in ("library", "standard_site", "organ", "native_organs", "patient_id", "sample_id"):
        if column in labels.columns:
            labels[column] = labels[column].fillna("").astype(str)
    pieces = []
    for (cohort, library), part in labels.groupby(["cohort", "library"], sort=False, dropna=False):
        z, z_ids = load_matrix(base, cohort, library, "Z")
        k, k_ids = load_matrix(base, cohort, library, "K")
        if z_ids != k_ids:
            raise SystemExit(f"{cohort} Z and K ids differ")
        id_row = {sample: i for i, sample in enumerate(z_ids)}
        order = [id_row[sample] for sample in part["sample_id"]]
        ld = load_ld(base, cohort)
        ld = ld.set_index("sample_id").loc[part["sample_id"].tolist()].reset_index()
        sites = part["standard_site"].tolist()
        beta = np.full(len(part), np.nan, dtype=np.float64)
        ld_z = z[order]
        grouped = {}
        for site in sorted(set(sites)):
            model_site = SITE_TO_MODEL.get(site)
            idx = np.array([i for i, value in enumerate(sites) if value == site])
            if model_site is None:
                continue
            sub = ld.iloc[idx]
            cleaned, part_beta = ld_for_rows(sub, genes, h_by_site[model_site], mu, keep_mt, b0_idx)
            ld_z_site = np.zeros_like(z[order])
            ld_z_site[idx] = cleaned
            ld_z = ld_z.copy()
            ld_z[idx] = cleaned
            beta[idx] = part_beta
            grouped[site] = True
        pred = predict_methods(z[order], k[order], sites, ld_z, beta, models, v0_means)
        extra = {"beta": beta}
        for method, (proba, labels_pred) in pred.items():
            extra[f"{method}__pred"] = labels_pred
            for j, organ in enumerate(ORGANS):
                extra[f"{method}__p_{organ}"] = proba[:, j]
        pieces.append(pd.concat([part.reset_index(drop=True), pd.DataFrame(extra)], axis=1))
        print("scored", cohort, library, len(part), flush=True)
    out = pd.concat(pieces, ignore_index=True)
    if len(out) != len(labels):
        raise SystemExit(f"prediction rows {len(out)} != label rows {len(labels)}")
    return out


def v0_site_means(paths, set_names) -> dict[str, np.ndarray]:
    split = paths["config"] / "split_gtex.tsv"
    ks_path = paths["data_processed"] / "features/ks_gtex_ref.parquet"
    note(split, "v0")
    note(ks_path, "v0")
    gtex = pd.read_csv(split, sep="\t", dtype=str)
    ref = gtex.loc[gtex["split"] == "GTEx-ref"]
    ks = pd.read_parquet(ks_path)
    tissue_mean = {}
    for tissue, sub in ref.groupby("tissue"):
        ids = [sample for sample in sub["sample"] if sample in ks.index]
        if ids:
            tissue_mean[tissue] = ks.loc[ids, set_names].to_numpy(dtype=np.float64).mean(axis=0)
    site_mean = {}
    for site, (tissues, _projects, _organs) in SITE_NATIVE.items():
        present = [tissue_mean[tissue] for tissue in tissues if tissue in tissue_mean]
        if present:
            site_mean[site] = np.mean(np.vstack(present), axis=0)
    return site_mean


def as_bool(frame: pd.DataFrame, column: str) -> pd.Series:
    if frame[column].dtype == bool:
        return frame[column]
    return frame[column].astype(str).eq("True")


def native_list(frame: pd.DataFrame) -> list[set[str]]:
    return [set(filter(None, str(text).split("|"))) for text in frame["native_organs"]]


def method_pred(frame: pd.DataFrame, method: str) -> np.ndarray:
    return frame[f"{method}__pred"].to_numpy(dtype=object)


def method_proba(frame: pd.DataFrame, method: str) -> np.ndarray:
    return np.column_stack([frame[f"{method}__p_{organ}"].to_numpy(dtype=np.float64) for organ in ORGANS])


def one_sided(n_favor: int, n_against: int) -> float:
    m = n_favor + n_against
    if m == 0:
        return 1.0
    return float(binom.sf(n_favor - 1, m, 0.5))


def bootstrap_diff(left: np.ndarray, right: np.ndarray, rng: np.random.Generator):
    if len(left) == 0:
        return None, None, None, None
    point = float(left.mean() - right.mean())
    draws = rng.integers(0, len(left), size=(N_BOOT, len(left)))
    diffs = left[draws].mean(axis=1) - right[draws].mean(axis=1)
    low, high = percentile_ci(diffs.tolist())
    return point, low, high, float(np.percentile(diffs, 5))


def holm(pvalues: dict[str, float], alpha: float = 0.05) -> dict[str, bool]:
    ordered = sorted(pvalues, key=lambda name: (pvalues[name], name))
    rejected = {name: False for name in pvalues}
    for i, name in enumerate(ordered):
        if pvalues[name] > alpha / (len(ordered) - i):
            break
        rejected[name] = True
    return rejected


def host_flags(frame: pd.DataFrame, method: str) -> np.ndarray:
    pred = method_pred(frame, method)
    truth = frame["organ"].to_numpy(dtype=object)
    natives = native_list(frame)
    return np.array([
        bool(natives[i]) and truth[i] not in natives[i] and pred[i] in natives[i]
        for i in range(len(frame))
    ])


def dedupe_patients(frame: pd.DataFrame) -> pd.DataFrame:
    """One row per patient_id. Within a cohort the label columns already do this.
    The same patient_id in two cohorts is kept once: polyA, then sample_id.
    A frame that already has unique patient ids keeps its row order.
    """
    if frame.empty or frame["patient_id"].is_unique:
        return frame.reset_index(drop=True)
    part = frame.copy()
    library = part["library"].fillna("")
    part["_library_rank"] = np.where(library.eq("polyA"), 0, 1)
    part = part.sort_values(["_library_rank", "sample_id", "cohort"], kind="mergesort")
    part = part.drop_duplicates("patient_id", keep="first").drop(columns="_library_rank")
    return part.reset_index(drop=True)


def correct_flags(frame: pd.DataFrame, method: str) -> np.ndarray:
    return method_pred(frame, method) == frame["organ"].to_numpy(dtype=object)


def gi_count(frame: pd.DataFrame, method: str) -> int:
    pred = method_pred(frame, method)
    truth = frame["organ"].to_numpy(dtype=object)
    natives = native_list(frame)
    n = 0
    for i in range(len(frame)):
        host = bool(natives[i]) and truth[i] not in natives[i] and pred[i] in natives[i]
        if pred[i] != truth[i] and not host and truth[i] in GI_SET and pred[i] in GI_SET:
            n += 1
    return n


def metric_dict(frame: pd.DataFrame, method: str) -> dict:
    if frame.empty:
        return {"n": 0}
    pred = method_pred(frame, method)
    proba = method_proba(frame, method)
    truth = frame["organ"].to_numpy(dtype=object)
    natives = native_list(frame)
    risk = np.array([bool(native) and truth[i] not in native for i, native in enumerate(natives)])
    native_truth = np.array([bool(native) and truth[i] in native for i, native in enumerate(natives)])
    pulled = risk & host_flags(frame, method)
    valid = pred != "NA"
    return {
        "n": int(len(frame)),
        "top1": float(np.mean(pred == truth)),
        "top3": float(np.mean(topk_hit(proba, ORGANS, truth, 3, valid))),
        "macro_f1": macro_f1(truth, pred),
        "n_at_risk": int(risk.sum()),
        "host_rate": float(pulled.sum() / risk.sum()) if risk.any() else None,
        "n_native_truth": int(native_truth.sum()),
        "native_truth_top1": float(np.mean(pred[native_truth] == truth[native_truth])) if native_truth.any() else None,
        "n_gi_internal": gi_count(frame, method),
        "n_stomach_pred": int((pred == "Stomach").sum()),
        "n_stomach_truth": int((truth == "Stomach").sum()),
    }


def hypothesis_rows(frame: pd.DataFrame, rng) -> tuple[list[dict], list[dict]]:
    risk = dedupe_patients(frame.loc[as_bool(frame, "selected_risk")])
    eval_df = dedupe_patients(frame.loc[as_bool(frame, "selected_eval")])
    native = dedupe_patients(frame.loc[as_bool(frame, "selected_native")])
    pool = dedupe_patients(frame.loc[as_bool(frame, "selected_pool_out")])
    primary = []
    pairs = {
        "AH1": (risk, host_flags(risk, "BASE-Z") & ~host_flags(risk, "SA-Z"), host_flags(risk, "SA-Z") & ~host_flags(risk, "BASE-Z"),
                host_flags(risk, "SA-Z").astype(float), host_flags(risk, "BASE-Z").astype(float)),
        "AH2": (eval_df, correct_flags(eval_df, "SA-Z") & ~correct_flags(eval_df, "BASE-Z"), correct_flags(eval_df, "BASE-Z") & ~correct_flags(eval_df, "SA-Z"),
                correct_flags(eval_df, "SA-Z").astype(float), correct_flags(eval_df, "BASE-Z").astype(float)),
        "AH3": (native, correct_flags(native, "SA-Z") & ~correct_flags(native, "BASE-Z"), correct_flags(native, "BASE-Z") & ~correct_flags(native, "SA-Z"),
                correct_flags(native, "SA-Z").astype(float), correct_flags(native, "BASE-Z").astype(float)),
    }
    computed = {}
    for name in ("AH1", "AH2", "AH3"):
        part, favor, against, left, right = pairs[name]
        n_favor = int(favor.sum()) if len(part) else 0
        n_against = int(against.sum()) if len(part) else 0
        pvalue = None if len(part) == 0 else one_sided(n_favor, n_against)
        diff, low, high, low5 = bootstrap_diff(left, right, rng) if len(part) else (None, None, None, None)
        computed[name] = {
            "hypothesis": name, "n": int(len(part)), "n_favor": n_favor, "n_against": n_against,
            "p": pvalue, "diff_sa_minus_comparator": diff, "ci_low": low, "ci_high": high, "onesided_low": low5,
            "status": "검정 불가" if len(part) == 0 else "계산",
        }
    ah1_sig = computed["AH1"]["status"] != "검정 불가" and computed["AH1"]["p"] <= 0.05
    ah2_sig = ah1_sig and computed["AH2"]["status"] != "검정 불가" and computed["AH2"]["p"] <= 0.05
    computed["AH1"]["tested"] = computed["AH1"]["status"] != "검정 불가"
    computed["AH1"]["reject"] = bool(ah1_sig)
    computed["AH2"]["tested"] = bool(ah1_sig and computed["AH2"]["status"] != "검정 불가")
    computed["AH2"]["reject"] = bool(ah2_sig)
    computed["AH3"]["tested"] = bool(ah2_sig and computed["AH3"]["status"] != "검정 불가")
    computed["AH3"]["reject"] = bool(computed["AH3"]["tested"] and computed["AH3"]["onesided_low"] is not None and computed["AH3"]["onesided_low"] > -0.10)
    primary = [computed[name] for name in ("AH1", "AH2", "AH3")]
    secondary_pairs = {
        "AS1": (native, correct_flags(native, "SA-G"), correct_flags(native, "SA-Z")),
        "AS2": (pool, host_flags(pool, "SA-Z"), host_flags(pool, "SA-pool22")),
        "AS3": (eval_df, correct_flags(eval_df, "SA-G"), correct_flags(eval_df, "SA-Z")),
    }
    # AS2 favor is SA-Z host and not pool22 host, so left indicator for the rate difference is pool22 - SA.
    raw = {}
    for name, (part, favor_mask_a, favor_mask_b) in secondary_pairs.items():
        if name == "AS2":
            favor = favor_mask_a & ~favor_mask_b
            against = favor_mask_b & ~favor_mask_a
            left = favor_mask_b.astype(float)
            right = favor_mask_a.astype(float)
        else:
            favor = favor_mask_a & ~favor_mask_b
            against = favor_mask_b & ~favor_mask_a
            left = favor_mask_a.astype(float)
            right = favor_mask_b.astype(float)
        n_favor = int(favor.sum()) if len(part) else 0
        n_against = int(against.sum()) if len(part) else 0
        pvalue = None if len(part) == 0 else one_sided(n_favor, n_against)
        diff, low, high, low5 = bootstrap_diff(left, right, rng) if len(part) else (None, None, None, None)
        raw[name] = {
            "hypothesis": name, "n": int(len(part)), "n_favor": n_favor, "n_against": n_against,
            "p": pvalue if pvalue is not None else 1.0, "diff": diff, "ci_low": low, "ci_high": high,
            "status": "검정 불가" if len(part) == 0 else "계산",
        }
    usable = {name: row["p"] for name, row in raw.items() if row["status"] != "검정 불가"}
    rejected = holm(usable) if usable else {}
    secondary = []
    for name in ("AS1", "AS2", "AS3"):
        row = raw[name]
        row["tested"] = row["status"] != "검정 불가"
        row["reject"] = bool(rejected.get(name, False))
        secondary.append(row)
    return primary, secondary


def conformal_table(frame: pd.DataFrame, thresholds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    thresholds = thresholds.copy()
    thresholds["alpha"] = pd.to_numeric(thresholds["alpha"])
    group = np.array(["liver" if site == "liver" else "lymph_node" if site == "lymph_node" else "other" for site in frame["standard_site"]])
    for layer, layer_df in frame.groupby("layer", sort=False):
        index = layer_df.index.to_numpy()
        for method in ("BASE-Z", "SA-Z"):
            proba = method_proba(layer_df, method)
            truth = layer_df["organ"].to_numpy(dtype=object)
            layer_group = group[frame.index.get_indexer(layer_df.index)]
            for variant in ("global", "mondrian"):
                for alpha in (0.1, 0.2):
                    selected = np.zeros((len(layer_df), len(ORGANS)), dtype=bool)
                    if variant == "global":
                        thr = thresholds.loc[
                            (thresholds["method"] == method) & (thresholds["variant"] == "global")
                            & (thresholds["alpha"] == alpha) & (thresholds["level"] == "all"),
                            "threshold",
                        ].iloc[0]
                        value = math.inf if thr == "inf" else float(thr)
                        selected[:] = (1.0 - proba) <= value
                    else:
                        for level in ("liver", "lymph_node", "other"):
                            thr = thresholds.loc[
                                (thresholds["method"] == method) & (thresholds["variant"] == "mondrian")
                                & (thresholds["alpha"] == alpha) & (thresholds["level"] == level),
                                "threshold",
                            ].iloc[0]
                            value = math.inf if thr == "inf" else float(thr)
                            mask = layer_group == level
                            selected[mask] = (1.0 - proba[mask]) <= value
                    sizes = selected.sum(axis=1)
                    covered = np.array([
                        bool(selected[i, ORGANS.index(truth[i])]) if truth[i] in ORGANS else False
                        for i in range(len(truth))
                    ])
                    single = np.empty(len(truth), dtype=object)
                    for i in range(len(truth)):
                        if sizes[i] == 1:
                            single[i] = ORGANS[int(np.flatnonzero(selected[i])[0])]

                    def emit(level, mask, layer=layer, method=method, variant=variant, alpha=alpha):
                        n = int(mask.sum())
                        sing = (sizes[mask] == 1) if n else np.array([], dtype=bool)
                        acc = None
                        if n and sing.any():
                            acc = float(np.mean(single[mask][sing] == truth[mask][sing]))
                        rows.append({
                            "layer": layer, "method": method, "variant": variant, "alpha": alpha, "level": level,
                            "n": n, "coverage": float(np.mean(covered[mask])) if n else None,
                            "mean_set_size": float(np.mean(sizes[mask])) if n else None,
                            "singleton_rate": float(np.mean(sizes[mask] == 1)) if n else None,
                            "singleton_accuracy": acc,
                        })

                    emit("all", np.ones(len(truth), dtype=bool))
                    for level in ("liver", "lymph_node", "other"):
                        emit(level, layer_group == level)
    return pd.DataFrame(rows)


def dl_row(records: list[dict]) -> dict:
    if len(records) < 2:
        return {"status": "하지 않음", "k": len(records)}
    d = np.array([row["d"] for row in records], dtype=float)
    v = np.array([row["v"] for row in records], dtype=float)
    w = 1.0 / v
    d_fe = float(np.sum(w * d) / np.sum(w))
    q = float(np.sum(w * (d - d_fe) ** 2))
    df = len(records) - 1
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    w_re = 1.0 / (v + tau2)
    estimate = float(np.sum(w_re * d) / np.sum(w_re))
    se = math.sqrt(1.0 / float(np.sum(w_re)))
    i2 = 0.0 if q <= 0 else max(0.0, (q - df) / q)
    return {
        "status": "계산", "k": len(records), "estimate": estimate,
        "ci_low": estimate - 1.96 * se, "ci_high": estimate + 1.96 * se, "I2": i2, "tau2": tau2,
    }


def meta_for(frame: pd.DataFrame) -> list[dict]:
    rows = []
    for layer, layer_df in frame.groupby("layer", sort=False):
        records = []
        for cohort, part in layer_df.groupby("cohort", sort=False):
            risk = part.loc[as_bool(part, "selected_risk")]
            n = len(risk)
            if n < 10:
                continue
            b = int((host_flags(risk, "BASE-Z") & ~host_flags(risk, "SA-Z")).sum())
            c = int((host_flags(risk, "SA-Z") & ~host_flags(risk, "BASE-Z")).sum())
            d = (c - b) / n
            b_v, c_v = b, c
            if b == 0 or c == 0:
                b_v += 0.5
                c_v += 0.5
            v = ((b_v + c_v) - (b_v - c_v) ** 2 / n) / n ** 2
            records.append({"cohort": cohort, "n": n, "b": b, "c": c, "d": d, "v": v})
        summary = dl_row(records)
        summary["layer"] = layer
        rows.append(summary)
        for record in records:
            record["layer"] = layer
            record["status"] = "cohort"
            rows.append(record)
    return rows


def write_tables(pred: pd.DataFrame, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261001)
    confirm = pred.loc[as_bool(pred, "include_confirm") & pred["layer"].eq("RNA-seq")].copy()
    primary, secondary = hypothesis_rows(confirm, rng)
    pd.DataFrame(primary).to_csv(out / "hypothesis_primary.tsv", sep="\t", index=False)
    pd.DataFrame(secondary).to_csv(out / "hypothesis_secondary.tsv", sep="\t", index=False)
    metric_rows = []
    for layer, part in pred.groupby("layer", sort=False):
        use = part.loc[as_bool(part, "selected_eval")]
        for method in METHODS:
            row = metric_dict(use, method)
            row.update({"layer": layer, "cohort": "all", "standard_site": "all", "method": method, "marked_descriptive": layer != "RNA-seq"})
            metric_rows.append(row)
    for (layer, cohort), part in pred.groupby(["layer", "cohort"], sort=False):
        use = part.loc[as_bool(part, "selected_eval")]
        descriptive = not bool(as_bool(part, "include_confirm").any())
        for method in METHODS:
            row = metric_dict(use, method)
            row.update({"layer": layer, "cohort": cohort, "standard_site": "all", "method": method, "marked_descriptive": descriptive})
            metric_rows.append(row)
    for (layer, site), part in pred.groupby(["layer", "standard_site"], sort=False):
        use = part.loc[as_bool(part, "selected_eval")]
        for method in METHODS:
            row = metric_dict(use, method)
            row.update({"layer": layer, "cohort": "all", "standard_site": site, "method": method, "marked_descriptive": layer != "RNA-seq"})
            metric_rows.append(row)
    pd.DataFrame(metric_rows).to_csv(out / "metrics.tsv", sep="\t", index=False)
    array = pred.loc[pred["layer"].eq("마이크로어레이")].copy()
    array_primary, array_secondary = hypothesis_rows(array.loc[as_bool(array, "layer_test")], rng)
    for row in array_primary + array_secondary:
        row["family"] = "microarray_descriptive"
        row["reject"] = False
        row["tested"] = False
    pd.DataFrame(array_primary + array_secondary).to_csv(out / "microarray_effects.tsv", sep="\t", index=False)
    pd.DataFrame(meta_for(pred.loc[as_bool(pred, "layer_test")])).to_csv(out / "meta_analysis.tsv", sep="\t", index=False)
    thresholds = pd.read_csv(ROOT / "config" / "conformal_thresholds_A3.tsv", sep="\t", dtype=str)
    note(ROOT / "config" / "conformal_thresholds_A3.tsv", "config")
    conformal_rows = pred.loc[as_bool(pred, "in_eval") & ~as_bool(pred, "excluded")].reset_index(drop=True)
    conformal_table(conformal_rows, thresholds).to_csv(out / "conformal.tsv", sep="\t", index=False)
    sensitivity(pred, rng).to_csv(out / "sensitivity.tsv", sep="\t", index=False)


def sensitivity(pred: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    confirm = pred.loc[as_bool(pred, "include_confirm") & pred["layer"].eq("RNA-seq")]
    all_samples = confirm.loc[as_bool(confirm, "in_eval") & ~as_bool(confirm, "excluded")]
    if len(all_samples):
        patients = all_samples["patient_id"].to_numpy()
        unique = pd.unique(patients)
        groups = [np.flatnonzero(patients == patient) for patient in unique]
        left_flag = host_flags(all_samples, "SA-Z")
        right_flag = host_flags(all_samples, "BASE-Z")
        left = left_flag.astype(float)
        right = right_flag.astype(float)
        point = float(left.mean() - right.mean())
        n_favor = int((right_flag & ~left_flag).sum())
        n_against = int((left_flag & ~right_flag).sum())
        draws = rng.integers(0, len(unique), size=(N_BOOT, len(unique)))
        diffs = np.empty(N_BOOT)
        for i, draw in enumerate(draws):
            idx = np.concatenate([groups[j] for j in draw])
            diffs[i] = left[idx].mean() - right[idx].mean()
        low, high = percentile_ci(diffs.tolist())
        rows.append({
            "analysis": "all_samples_cluster_AH1", "n": int(len(all_samples)),
            "n_favor": n_favor, "n_against": n_against,
            "diff": point, "ci_low": low, "ci_high": high,
        })
    capture = confirm.loc[confirm["library"].ne("capture")]
    capture_risk = dedupe_patients(capture.loc[as_bool(capture, "selected_risk")])
    if len(capture_risk):
        diff = float(host_flags(capture_risk, "SA-Z").mean() - host_flags(capture_risk, "BASE-Z").mean())
    else:
        diff = None
    rows.append({"analysis": "drop_SU2C_capture_AH1", "n": int(len(capture_risk)), "diff": diff, "ci_low": None, "ci_high": None})
    if "GSE209998" in set(confirm["cohort"]):
        kept = confirm.loc[confirm["cohort"].ne("GSE209998")]
        risk = dedupe_patients(kept.loc[as_bool(kept, "selected_risk")])
        rows.append({
            "analysis": "drop_GSE209998_AH1", "n": int(len(risk)),
            "diff": float(host_flags(risk, "SA-Z").mean() - host_flags(risk, "BASE-Z").mean()) if len(risk) else None,
            "ci_low": None, "ci_high": None,
        })
    else:
        rows.append({"analysis": "drop_GSE209998_AH1", "n": None, "diff": None, "ci_low": None, "ci_high": None, "note": "확인 검정에 없음"})
    for cohort in sorted(confirm["cohort"].unique()):
        kept = confirm.loc[confirm["cohort"].ne(cohort)]
        for name, mask_col, left_m, right_m, left_flag, right_flag in (
            ("AH1", "selected_risk", "SA-Z", "BASE-Z", host_flags, host_flags),
            ("AH2", "selected_eval", "SA-Z", "BASE-Z", correct_flags, correct_flags),
        ):
            part = dedupe_patients(kept.loc[as_bool(kept, mask_col)])
            if len(part) == 0:
                diff = None
            else:
                diff = float(left_flag(part, left_m).mean() - right_flag(part, right_m).mean())
            rows.append({"analysis": f"leave_out_{cohort}_{name}", "n": int(len(part)), "diff": diff, "ci_low": None, "ci_high": None})
    return pd.DataFrame(rows)


def data_base() -> Path:
    if MODE == "synthetic":
        return SYN
    return ROOT


def label_path() -> Path:
    if MODE == "synthetic":
        return SYN / "config" / "aux_eval_labels.tsv"
    return ROOT / "config" / "aux_eval_labels.tsv"


BOOL_COLS = (
    "in_eval", "in_risk", "in_native", "in_pool_out_risk", "excluded", "include_confirm",
    "layer_test", "selected_eval", "selected_risk", "selected_native", "selected_pool_out", "pool_out_site",
)


def read_labels(path: Path) -> pd.DataFrame:
    note(path, "labels")
    labels = pd.read_csv(path, sep="\t", dtype=str)
    for column in BOOL_COLS:
        if column in labels.columns:
            labels[column] = labels[column].eq("True")
    if "library" in labels.columns:
        labels["library"] = labels["library"].fillna("")
    return labels


def sets_only() -> None:
    """Label-file set sizes. Opens the label file and nothing else."""
    log_path = ROOT / "results" / "stage7" / "sets_only" / "opened_paths.txt"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    label = ROOT / "config" / "aux_eval_labels.tsv"
    allowed = {label.resolve()}
    log_path.write_text("")

    def guarded_note(path: Path, purpose: str) -> None:
        resolved = Path(path).resolve()
        with log_path.open("a") as handle:
            handle.write(f"{purpose}\t{path}\n")
        if resolved not in allowed:
            raise SystemExit(f"sets-only opened a file outside the label file: {path}")

    global note
    note = guarded_note
    labels = read_labels(label)
    membership = labels.drop(columns=[c for c in ("selected_eval", "selected_risk", "selected_native", "selected_pool_out", "include_confirm", "layer_test") if c in labels.columns])
    # select_patients recomputes include_confirm and the selected columns from membership.
    fresh = select_patients(membership)
    rows = []
    for column in ("selected_eval", "selected_risk", "selected_native", "selected_pool_out"):
        left = labels[column].to_numpy()
        right = fresh[column].to_numpy()
        if not np.array_equal(left, right):
            raise SystemExit(f"recomputed {column} does not match the label file")
    confirm = labels.loc[labels["include_confirm"] & labels["layer"].eq("RNA-seq")]
    summary = {
        "confirm_eval": int(confirm["selected_eval"].sum()),
        "confirm_risk": int(confirm["selected_risk"].sum()),
        "confirm_native": int(confirm["selected_native"].sum()),
        "confirm_pool": int(confirm["selected_pool_out"].sum()),
    }
    expected = {"confirm_eval": 729, "confirm_risk": 427, "confirm_native": 315, "confirm_pool": 71}
    if summary != expected:
        raise SystemExit(f"confirm set sizes {summary} != {expected}")
    gse = labels.loc[labels["cohort"].eq("GSE209998")]
    gse_summary = {
        "eval": int(gse["selected_eval"].sum()),
        "risk": int(gse["selected_risk"].sum()),
        "native": int(gse["selected_native"].sum()),
        "pool": int(gse["selected_pool_out"].sum()),
        "confirm": int(gse["include_confirm"].sum()),
    }
    if gse_summary != {"eval": 123, "risk": 64, "native": 44, "pool": 3, "confirm": 0}:
        raise SystemExit(f"GSE209998 sizes {gse_summary}")
    cohort_expected = {
        "blca_iatlas_imvigor210_2017": (276, 82, 194, 67),
        "brca_iatlas_anders_2022": (30, 17, 13, 1),
        "mel_dfci_2019": (121, 35, 83, 3),
        "paad_iatlas_prince_2022": (58, 56, 0, 0),
        "GSE50760": (18, 18, 18, 0),
        "prad_su2c_2019": (226, 219, 7, 0),
    }
    for cohort, expect in cohort_expected.items():
        part = confirm.loc[confirm["cohort"].eq(cohort)]
        got = (
            int(part["selected_eval"].sum()), int(part["selected_risk"].sum()),
            int(part["selected_native"].sum()), int(part["selected_pool_out"].sum()),
        )
        if got != expect:
            raise SystemExit(f"{cohort} sizes {got} != {expect}")
    out = ROOT / "results" / "stage7" / "sets_only" / "set_sizes.tsv"
    pd.DataFrame([
        {"scope": "confirm", **summary},
        {"scope": "GSE209998", "confirm_eval": gse_summary["eval"], "confirm_risk": gse_summary["risk"],
         "confirm_native": gse_summary["native"], "confirm_pool": gse_summary["pool"]},
    ]).to_csv(out, sep="\t", index=False)
    print("sets-only", summary, "GSE209998", gse_summary, flush=True)


def _assign_prediction(frame: pd.DataFrame, idx, method: str, label: str) -> None:
    frame.loc[idx, f"{method}__pred"] = label
    for organ in ORGANS:
        frame.loc[idx, f"{method}__p_{organ}"] = 1.0 if organ == label else 0.0


def _host_label(native_text: str, truth: str) -> str:
    organs = [item for item in str(native_text).split("|") if item and item != truth]
    if not organs:
        raise SystemExit(f"no host class for {truth} / {native_text}")
    return organs[0]


def _other_label(truth: str) -> str:
    return "Stomach" if truth != "Stomach" else "Lung"


def build_planted_predictions(labels: pd.DataFrame) -> pd.DataFrame:
    """Known predictions on the synthetic labels plus single-sample fixture patients."""
    base = labels.copy()
    rng = np.random.default_rng(20261013)
    rna = base["include_confirm"] & base["layer"].eq("RNA-seq")
    d1_pool = base.loc[rna & base["in_risk"] & ~base["selected_risk"]].sort_values(["cohort", "sample_id"])
    if len(d1_pool) < 10:
        raise SystemExit(f"D1 candidates {len(d1_pool)}")
    d1_idx = rng.choice(d1_pool.index.to_numpy(), size=10, replace=False)
    array = base["layer"].eq("마이크로어레이") & base["in_risk"]
    d2_pool = base.loc[array].sort_values(["cohort", "sample_id"])
    if len(d2_pool) < 20:
        raise SystemExit(f"D2 candidates {len(d2_pool)}")
    d2_idx = rng.choice(d2_pool.index.to_numpy(), size=20, replace=False)
    fixtures = []
    groups = {
        "in": (41, "liver", "Breast", "Liver|Biliary", False, True),
        "out": (13, "kidney", "Breast", "Kidney", True, True),
        "native": (11, "liver", "Liver", "Liver|Biliary", False, False),
    }
    # counts: in = P1 30 + P2 5 + P3 6; out = P4 12 + P5 1; native = P6 3 + P7 8
    made = {}
    for name, (count, site, organ, native, pool, risk) in groups.items():
        ids = [f"PLANT-{name}-{i:03d}" for i in range(count)]
        made[name] = [ids[i] for i in rng.permutation(count)]
        for patient in ids:
            fixtures.append({
                "cohort": "syn_rna_1", "layer": "RNA-seq", "sample_id": f"{patient}|one",
                "patient_id": patient, "library": "", "organ": organ, "organ_rule": "planted",
                "raw_site": site, "standard_site": site, "native_organs": native,
                "in_eval": True, "in_risk": risk, "in_native": not risk,
                "pool_out_site": pool, "in_pool_out_risk": pool and risk,
                "excluded": False, "exclude_reason": "",
            })
    extra = pd.DataFrame(fixtures)
    combined = pd.concat([base, extra], ignore_index=True)
    combined = select_patients(combined)
    # Fixture patients have one sample, so selection keeps that row.
    single = combined.groupby(["layer", "patient_id"]).size()
    rna_single = set(single[single == 1].index)
    for patient in made["in"] + made["out"] + made["native"]:
        if ("RNA-seq", patient) not in {(layer, pid) for layer, pid in rna_single}:
            raise SystemExit(f"{patient} is not a single-sample RNA-seq patient")
    frame = combined.copy()
    n = len(frame)
    truth = frame["organ"].to_numpy(dtype=object)
    frame["beta"] = np.where(frame["in_native"].to_numpy(), np.nan, 0.5)
    organ_col = {organ: i for i, organ in enumerate(ORGANS)}
    base_p = np.zeros((n, len(ORGANS)), dtype=np.float64)
    for i, organ in enumerate(truth):
        j = organ_col.get(organ)
        if j is not None:
            base_p[i, j] = 1.0
    for method in METHODS:
        frame[f"{method}__pred"] = truth
        for j, organ in enumerate(ORGANS):
            frame[f"{method}__p_{organ}"] = base_p[:, j]
    # The permutation order is the planting order within each group.
    in_ids = made["in"]
    out_ids = made["out"]
    native_ids = made["native"]
    plan = [
        ("P1", in_ids[:30], {"BASE-Z": "host"}),
        ("P2", in_ids[30:35], {"SA-Z": "host"}),
        ("P3", in_ids[35:41], {"BASE-Z": "host", "SA-Z": "host"}),
        ("P4", out_ids[:12], {"BASE-Z": "host", "SA-Z": "host"}),
        ("P5", out_ids[12:13], {"SA-pool22": "host"}),
        ("P6", native_ids[:3], {"BASE-Z": "other"}),
        ("P7", native_ids[3:11], {"SA-Z": "other", "SA-pool22": "other"}),
    ]
    for _name, patients, edits in plan:
        rows = frame.loc[frame["patient_id"].isin(patients)]
        for idx, row in rows.iterrows():
            for method, kind in edits.items():
                label = _host_label(row["native_organs"], row["organ"]) if kind == "host" else _other_label(row["organ"])
                _assign_prediction(frame, idx, method, label)
    for idx in d1_idx:
        row = frame.loc[idx]
        _assign_prediction(frame, idx, "BASE-Z", _host_label(row["native_organs"], row["organ"]))
    for idx in d2_idx:
        row = frame.loc[idx]
        _assign_prediction(frame, idx, "BASE-Z", _host_label(row["native_organs"], row["organ"]))
    use_base = frame["beta"].isna() | (frame["beta"] <= SA_G_TAU)
    for organ in ORGANS:
        frame[f"SA-G__p_{organ}"] = np.where(use_base, frame[f"BASE-Z__p_{organ}"], frame[f"SA-Z__p_{organ}"])
    frame["SA-G__pred"] = np.where(use_base, frame["BASE-Z__pred"], frame["SA-Z__pred"])
    if len(frame) != n:
        raise SystemExit("planted frame length changed")
    return frame


def planted() -> None:
    directory = PLANTED_DIR
    if directory is None:
        raise SystemExit(2)
    directory = directory if directory.is_absolute() else ROOT / directory
    real = (ROOT / "config" / "aux_eval_labels.tsv").resolve()
    source = (SYN / "config" / "aux_eval_labels.tsv").resolve()
    if source == real or "synthetic" not in source.parts:
        raise SystemExit("planted mode refuses the real label file")
    reset_log()
    labels = read_labels(source)
    frame = build_planted_predictions(labels)
    directory.mkdir(parents=True, exist_ok=True)
    label_out = directory / "labels.tsv"
    frame.drop(columns=[c for c in frame.columns if "__" in c or c == "beta"]).to_csv(label_out, sep="\t", index=False)
    pred_path = directory / "predictions.parquet"
    frame.to_parquet(pred_path, index=False)
    write_tables(frame, directory / "tables")
    print("planted", directory, flush=True)


def main() -> None:
    if MODE == "sets-only":
        sets_only()
        return
    if MODE == "planted":
        planted()
        return
    reset_log()
    if MODE == "real":
        check_unlock()
    check_freeze()
    paths = project_paths()
    if MODE == "synthetic":
        build_synthetic(paths)
    note(label_path(), "labels")
    labels = pd.read_csv(label_path(), sep="\t", dtype=str)
    for column in ("in_eval", "in_risk", "in_native", "in_pool_out_risk", "excluded", "include_confirm", "layer_test", "selected_eval", "selected_risk", "selected_native", "selected_pool_out", "pool_out_site"):
        if column in labels.columns:
            labels[column] = labels[column].eq("True")
    pred = score_base(data_base(), labels, paths)
    out = ROOT / "results" / "stage7" / "confirm" if MODE == "real" else ROOT / "results" / "stage6" / MODE
    out.mkdir(parents=True, exist_ok=True)
    pred_path = out / "predictions.parquet"
    pred.to_parquet(pred_path, index=False)
    digest = hashlib.sha256(pred_path.read_bytes()).hexdigest()
    (out / "predictions.sha256").write_text(digest + "\n")
    write_tables(pred, out / "tables")
    print("predictions", digest, flush=True)


if __name__ == "__main__":
    main()
