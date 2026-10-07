#!/usr/bin/env python3
import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"
os.environ["NUMBA_NUM_THREADS"] = "16"

"""Train stage-2 models. MET500 is not read. GTEx-sim is not read.

Augmentation seeds come from augmentation_seeds(). SA chooses a host tissue
uniformly from the sorted host-pool names, then a GTEx-ref sample of that tissue.
A site model with more than one tissue draws uniformly from the pooled samples.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import joblib
import numpy as np
import pandas as pd

from stage2_common import (
    BASE_C,
    FIXED_SITE_TISSUES,
    IF_RHO,
    N_MIX,
    RHO_HIGH,
    RHO_LOW,
    SITE_MODEL_ORDER,
    augmentation_seeds,
    gtex_ref_ids,
    if_tumor_index,
    load_gene_pack,
    load_tpm,
    model_names,
    rank_z,
    score_ks,
    sum_1e6,
    train_frame,
)
from stage2_fit import (
    assemble_aug,
    cv_select_C,
    fit_logit,
    repeat_labels,
    select_k,
    select_z,
    sensitivity,
)
from util import load_paths


def pool_tissues(site_tissues: dict[str, list[str]]) -> list[str]:
    names = []
    for site in SITE_MODEL_ORDER:
        names.extend(site_tissues[site])
    return sorted(set(names))


def _host(rng, tissues, ids_by_tissue, pool, pooled: bool):
    if pooled:
        return pool[int(rng.integers(0, len(pool)))]
    tissue = tissues[int(rng.integers(0, len(tissues)))]
    group = ids_by_tissue[tissue]
    return tissue, group[int(rng.integers(0, len(group)))]


def draw_mixes(tumors: list[str], tissues: list[str], ids_by_tissue: dict, rng, pooled: bool, rho_fixed: float | None):
    """Training draws 4 mixes with Uniform rho.

    Fixed rho and pooled=False draws one host from each listed tissue (pool IF).
    Fixed rho and pooled=True draws one host from the pooled samples (site IF).
    """
    records = []
    pool = [(tissue, sample) for tissue in tissues for sample in ids_by_tissue[tissue]]
    for tumor in tumors:
        if rho_fixed is None:
            for _ in range(N_MIX):
                tissue, host = _host(rng, tissues, ids_by_tissue, pool, pooled)
                records.append({
                    "tumor": tumor,
                    "host": host,
                    "tissue": tissue,
                    "rho": float(rng.uniform(RHO_LOW, RHO_HIGH)),
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
        out[start:start + len(chunk)] = rho * tumor_norm[t_idx] + (1.0 - rho) * host
    return out


def save_bundle(path: Path, bundle: dict, extra: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(bundle)
    payload.update(extra)
    joblib.dump(payload, path)
    print("saved", path.name, "C", bundle["C"], "n_iter", bundle["n_iter"], "hit", bundle["max_iter_hit"], flush=True)


def top_names(scores: np.ndarray, names: list[str], k: int = 20) -> list[dict]:
    order = np.lexsort((np.arange(scores.size), -np.where(np.isfinite(scores), scores, np.inf)))
    rows = []
    for rank, idx in enumerate(order[:k], start=1):
        rows.append({"rank": rank, "name": names[int(idx)], "s": None if not np.isfinite(scores[idx]) else float(scores[idx])})
    return rows


def main() -> None:
    paths = load_paths()
    if not (paths["config"] / "analysis_plan_A2.yaml").exists():
        raise SystemExit("config/analysis_plan_A2.yaml is missing. Training was not started.")
    choice = json.loads((paths["results"] / "stage2/proxy_choice.json").read_text())
    ln_proxy = choice["ln_proxy"]
    bm_proxy = choice["bm_proxy"]
    site_tissues = dict(FIXED_SITE_TISSUES)
    site_tissues["lymph_node"] = [ln_proxy]
    site_tissues["bone_marrow"] = [bm_proxy]
    tissues = pool_tissues(site_tissues)
    names = model_names(ln_proxy, bm_proxy)
    seeds = augmentation_seeds(names)
    out = paths["results"] / "stage2"
    out.mkdir(parents=True, exist_ok=True)
    (out / "augmentation_seeds.json").write_text(json.dumps({
        "rule": "sorted model names, one numpy Generator(20261001), integers(0, 2**31-1) in that order",
        "names_sorted": sorted(names),
        "seeds": seeds,
        "host_pool_sorted": tissues,
        "site_tissues": site_tissues,
        "sa_host_draw": "uniform tissue from host_pool_sorted, then uniform GTEx-ref sample",
        "site_host_draw": "uniform among pooled GTEx-ref samples of that site",
        "rho_mix": "numpy Generator.uniform(0.15, 1.0), half-open",
    }, indent=2))

    genes, b0_idx, offsets, indices, set_names = load_gene_pack(paths)
    train = train_frame(paths)
    tumors = train["sample"].tolist()
    train_pos = {sample: i for i, sample in enumerate(tumors)}
    print("load tpm", flush=True)
    samples, tpm = load_tpm(paths["data_processed"] / "toil/tpm_G.h5")
    row = {sample: i for i, sample in enumerate(samples)}
    missing = [sample for sample in tumors if sample not in row]
    if missing:
        raise SystemExit(f"train samples missing {len(missing)}")
    tumor_norm = sum_1e6(tpm[[row[sample] for sample in tumors]])
    ids_by_tissue = gtex_ref_ids(paths, tissues)
    bank = {}
    host_norm = {}
    for tissue, ids in ids_by_tissue.items():
        if not ids:
            raise SystemExit(f"no GTEx-ref for {tissue}")
        host_norm[tissue] = sum_1e6(tpm[[row[sample] for sample in ids]])
        for i, sample in enumerate(ids):
            bank[sample] = host_norm[tissue][i]
    del tpm
    print("pure features", flush=True)
    z_pure = rank_z(tumor_norm).astype(np.float32)
    ks = pd.read_parquet(paths["data_processed"] / "features/ks_tcga.parquet")
    k_pure = ks.loc[tumors, set_names].to_numpy(dtype=np.float32)
    del ks
    var_z = z_pure.var(axis=0, ddof=1).astype(np.float64)
    sd_z = z_pure.std(axis=0, ddof=1).astype(np.float64)
    sd_k = k_pure.std(axis=0, ddof=1).astype(np.float64)
    y = train["learning_label"].to_numpy()
    model_dir = out / "models"
    print("BASE", flush=True)
    base_z = fit_logit(z_pure[:, b0_idx], y, BASE_C["Z"], scale=False)
    save_bundle(model_dir / "BASE_Z.joblib", base_z, {"features": "B0_genes", "representation": "Z"})
    base_k = fit_logit(k_pure, y, BASE_C["K"], scale=True)
    save_bundle(model_dir / "BASE_K.joblib", base_k, {"features": "all_sets", "representation": "K"})

    def fit_pair(tag: str, Xz, Xk, y_fit, C_z, C_k, extra_z, extra_k):
        bundle_z = fit_logit(Xz, y_fit, C_z, scale=True)
        save_bundle(model_dir / f"{tag}_Z.joblib", bundle_z, extra_z)
        bundle_k = fit_logit(Xk, y_fit, C_k, scale=True)
        save_bundle(model_dir / f"{tag}_K.joblib", bundle_k, extra_k)
        return bundle_z, bundle_k

    print("SA augment", flush=True)
    sa_rng = np.random.default_rng(seeds["SA"])
    sa_records = draw_mixes(tumors, tissues, ids_by_tissue, sa_rng, pooled=False, rho_fixed=None)
    pd.DataFrame(sa_records).to_parquet(out / "mixes_SA.parquet", index=False)
    mixed = materialize(tumor_norm, train_pos, bank, sa_records)
    z_mix = rank_z(mixed).astype(np.float32)[:, b0_idx]
    k_mix = score_ks(mixed, offsets, indices)
    del mixed
    Xz = assemble_aug(z_pure[:, b0_idx], z_mix)
    Xk = assemble_aug(k_pure, k_mix)
    del z_mix, k_mix
    groups = repeat_labels(train["patient"].to_numpy())
    y_aug = repeat_labels(y)
    print("SA CV Z", Xz.shape, flush=True)
    cv_z = cv_select_C(Xz, y_aug, groups)
    print("SA CV K", Xk.shape, flush=True)
    cv_k = cv_select_C(Xk, y_aug, groups)
    (out / "cv_SA.json").write_text(json.dumps({"Z": cv_z, "K": cv_k}, indent=2))
    C_z, C_k = cv_z["C"], cv_k["C"]
    print("chosen C", C_z, C_k, flush=True)
    fit_pair("SA", Xz, Xk, y_aug, C_z, C_k, {"representation": "Z"}, {"representation": "K"})
    del Xz, Xk

    if_pos = if_tumor_index(train)
    if_tumors = [tumors[i] for i in if_pos]
    (out / "if_tumor_ids.txt").write_text("\n".join(if_tumors) + "\n")

    print("IF pool", flush=True)
    if_rng = np.random.default_rng(seeds["IF_pool"])
    if_records = draw_mixes(if_tumors, tissues, ids_by_tissue, if_rng, pooled=False, rho_fixed=IF_RHO)
    pd.DataFrame(if_records).to_parquet(out / "mixes_IF_pool.parquet", index=False)
    mixed = materialize(tumor_norm, train_pos, bank, if_records)
    z_mix = rank_z(mixed).astype(np.float32)
    k_mix = score_ks(mixed, offsets, indices)
    pure_rows = np.array([train_pos[row["tumor"]] for row in if_records])
    s_z = sensitivity(z_pure[pure_rows], z_mix, sd_z)
    s_k = sensitivity(k_pure[pure_rows], k_mix, sd_k)
    del mixed, z_mix, k_mix
    kept_z, dropped_z = select_z(s_z, var_z, 0.20)
    kept_k, dropped_k = select_k(s_k, 0.20)
    pd.DataFrame({"gene": [genes[i] for i in dropped_z]}).to_csv(out / "IF_pool_Z_dropped.txt", index=False)
    pd.DataFrame({"set": [set_names[i] for i in dropped_k]}).to_csv(out / "IF_pool_K_dropped.txt", index=False)
    (out / "IF_pool_top20.json").write_text(json.dumps({
        "Z": top_names(s_z, genes),
        "K": top_names(s_k, set_names),
        "n_dropped_Z": int(len(dropped_z)),
        "n_dropped_K": int(len(dropped_k)),
        "n_kept_Z": int(len(kept_z)),
        "n_zero_sd_Z": int(np.sum(sd_z == 0)),
        "n_zero_sd_K": int(np.sum(sd_k == 0)),
    }, indent=2))
    fit_pair(
        "IF20",
        z_pure[:, kept_z],
        k_pure[:, kept_k],
        y,
        C_z,
        C_k,
        {"representation": "Z", "feature_index": kept_z.astype(int).tolist()},
        {"representation": "K", "feature_index": kept_k.astype(int).tolist()},
    )

    summaries = []
    for site in SITE_MODEL_ORDER:
        print("site", site, flush=True)
        site_list = site_tissues[site]
        pooled = len(site_list) > 1
        if_rng = np.random.default_rng(seeds[f"IF_SC_{site}"])
        if_records = draw_mixes(if_tumors, site_list, ids_by_tissue, if_rng, pooled=pooled, rho_fixed=IF_RHO)
        pd.DataFrame(if_records).to_parquet(out / f"mixes_IF_SC_{site}.parquet", index=False)
        mixed = materialize(tumor_norm, train_pos, bank, if_records)
        z_mix = rank_z(mixed).astype(np.float32)
        k_mix = score_ks(mixed, offsets, indices)
        pure_rows = np.array([train_pos[row["tumor"]] for row in if_records])
        s_z = sensitivity(z_pure[pure_rows], z_mix, sd_z)
        s_k = sensitivity(k_pure[pure_rows], k_mix, sd_k)
        del mixed, z_mix, k_mix
        z20, z20_drop = select_z(s_z, var_z, 0.20)
        z40, z40_drop = select_z(s_z, var_z, 0.40)
        k20, k20_drop = select_k(s_k, 0.20)
        k40, k40_drop = select_k(s_k, 0.40)
        (out / f"IF_SC_{site}_top20.json").write_text(json.dumps({
            "Z": top_names(s_z, genes),
            "K": top_names(s_k, set_names),
            "n_dropped_Z20": int(len(z20_drop)),
            "n_dropped_Z40": int(len(z40_drop)),
            "n_dropped_K20": int(len(k20_drop)),
            "n_dropped_K40": int(len(k40_drop)),
        }, indent=2))
        sc_rng = np.random.default_rng(seeds[f"SC_{site}"])
        sc_records = draw_mixes(tumors, site_list, ids_by_tissue, sc_rng, pooled=pooled, rho_fixed=None)
        pd.DataFrame(sc_records).to_parquet(out / f"mixes_SC_{site}.parquet", index=False)
        mixed = materialize(tumor_norm, train_pos, bank, sc_records)
        # One rank, three gene sets.
        from stage2_common import rank_z_take
        z_b0, z_if20, z_if40 = rank_z_take(mixed, [b0_idx, z20, z40])
        k_mix = score_ks(mixed, offsets, indices)
        del mixed
        y_aug = repeat_labels(y)
        fit_pair(
            f"SC_{site}",
            assemble_aug(z_pure[:, b0_idx], z_b0),
            assemble_aug(k_pure, k_mix),
            y_aug, C_z, C_k,
            {"representation": "Z", "site": site},
            {"representation": "K", "site": site},
        )
        fit_pair(
            f"SCIF20_{site}",
            assemble_aug(z_pure[:, z20], z_if20),
            assemble_aug(k_pure[:, k20], k_mix[:, k20]),
            y_aug, C_z, C_k,
            {"representation": "Z", "site": site, "feature_index": z20.astype(int).tolist()},
            {"representation": "K", "site": site, "feature_index": k20.astype(int).tolist()},
        )
        fit_pair(
            f"SCIF40_{site}",
            assemble_aug(z_pure[:, z40], z_if40),
            assemble_aug(k_pure[:, k40], k_mix[:, k40]),
            y_aug, C_z, C_k,
            {"representation": "Z", "site": site, "feature_index": z40.astype(int).tolist()},
            {"representation": "K", "site": site, "feature_index": k40.astype(int).tolist()},
        )
        del z_b0, z_if20, z_if40, k_mix
        summaries.append(site)

    print("NC", flush=True)
    normal_z = []
    normal_k = []
    normal_y = []
    for tissue in tissues:
        block = host_norm[tissue]
        normal_z.append(rank_z(block).astype(np.float32)[:, b0_idx])
        normal_k.append(score_ks(block, offsets, indices))
        normal_y.append(np.array([f"N_{tissue}"] * block.shape[0], dtype=object))
    Xz = np.vstack([z_pure[:, b0_idx], *normal_z])
    Xk = np.vstack([k_pure, *normal_k])
    y_nc = np.concatenate([y.astype(object), *normal_y])
    counts = {label: int(np.sum(y_nc == label)) for label in sorted(set(y_nc.tolist())) if str(label).startswith("N_")}
    (out / "NC_class_counts.json").write_text(json.dumps(counts, indent=2))
    fit_pair("NC", Xz, Xk, y_nc, C_z, C_k, {"representation": "Z"}, {"representation": "K"})
    print("train done", summaries, flush=True)


if __name__ == "__main__":
    main()
