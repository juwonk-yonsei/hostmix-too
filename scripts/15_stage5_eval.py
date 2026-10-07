#!/usr/bin/env python3
"""Stage-5 evaluation on TCGA, the extended simulation, MET500, and POG570.

POG570 numbers from this script are post hoc. Does not set a thread cap.
Does not modify frozen predictions, models, thresholds, or labels.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import joblib
import numpy as np
import pandas as pd

from genes import probemap_symbol
from loader import load_pog570_expression, read_pog570_table, read_skcm_clinical
from mapping_rules import MET500_COHORT_TO_ORGAN, map_site
from metrics import ORGANS, mask_native
from stage2_fit import ld_apply
from stage3_rules import native_organs
from stage5_lib import (
    GI_SET,
    NAMED_MET_SITES,
    POOL10,
    SIM_RHOS,
    TISSUE_TO_SITE,
    error_parts,
    load_gene_pack,
    load_paths,
    load_references,
    load_tpm,
    metric_row,
    mcnemar_two_sided,
    organ_of_labels,
    paired_mean,
    predict_organ,
    rank_b0,
    site_native_organs,
    stage5_dir,
    substitute_site,
    sum_1e6,
    test_frame,
    train_frame,
    variant_names,
    SITE_TO_MODEL,
)
from util import SEED

DELTA_RHOS = {0.8, 0.6, 0.4}
STOMACH_MODELS = ["BASE-Z", "SA-m0", "SA-Z", "SA-r30", "SA-r50", "SA-pool3", "SA-pool22", "SA-MLP"]
SARCOMA_MODELS = [
    "BASE-Z", "SA-Z", "SA-LOHO-Brain - Cortex",
    "SA-LOHO-Adipose - Subcutaneous", "SA-LOHO-Muscle - Skeletal",
]
SIM_BASE_MODELS = ["BASE-Z", "SA-Z", "SC-Z", "LD-Z", "NC-Z", "M1-Z"]


def load_z_models(paths, root: Path) -> dict:
    out = {}
    model_dir = paths["results"] / "stage2" / "models"
    for name in ("BASE_Z", "SA_Z", "NC_Z", "BASE_K", "SA_K"):
        out[name] = joblib.load(model_dir / f"{name}.joblib")
    for site in sorted(set(SITE_TO_MODEL.values())):
        out[f"SC_{site}_Z"] = joblib.load(model_dir / f"SC_{site}_Z.joblib")
    var_dir = root / "variants" / "models"
    missing = [name for name in variant_names() if not (var_dir / f"{name}.joblib").exists()]
    if missing:
        raise SystemExit(f"missing variant models: {missing}")
    for name in variant_names():
        out[name] = joblib.load(var_dir / f"{name}.joblib")
    return out


def single_ids(models) -> list[str]:
    return ["BASE-Z", "SA-Z", *variant_names()]


def collapse_to_g(expr: pd.DataFrame, genes: list[str]):
    mapping = ensembl_base_to_symbol()
    bases = pd.Index(expr.index.astype(str)).str.split(".").str[0]
    symbol = np.asarray(bases.map(mapping), dtype=object)
    keep = pd.notna(symbol)
    sub = expr.loc[keep].astype(np.float32)
    grouped = sub.groupby(symbol[keep], sort=False).sum()
    g_index = pd.Index(genes)
    present = grouped.index.intersection(g_index)
    missing = g_index.difference(grouped.index)
    aligned = grouped.reindex(g_index).fillna(0).astype(np.float32)
    matrix = np.ascontiguousarray(aligned.to_numpy(dtype=np.float32).T)
    audit = {
        "n_ensembl_rows": int(expr.shape[0]),
        "n_ensembl_mapped": int(keep.sum()),
        "n_symbols_summed": int(grouped.shape[0]),
        "n_G": len(genes),
        "n_G_present": int(len(present)),
        "n_G_filled_0": int(len(missing)),
        "n_samples": int(expr.shape[1]),
    }
    return matrix, [str(c) for c in expr.columns], audit


def ensembl_base_to_symbol() -> dict[str, str]:
    probe = probemap_symbol()
    by_base: dict[str, set[str]] = {}
    for ens, sym in zip(probe["id"].astype(str), probe["gene"].astype(str)):
        by_base.setdefault(ens.split(".")[0], set()).add(sym)
    return {base: next(iter(symbols)) for base, symbols in by_base.items() if len(symbols) == 1}


def k_matrix(paths, samples: list[str], set_names: list[str]) -> np.ndarray:
    ks = pd.read_parquet(paths["data_processed"] / "features/ks_tcga.parquet")
    return ks.loc[samples, set_names].to_numpy(dtype=np.float32)


def tcga_tables(paths, root, models, b0_idx, set_names) -> None:
    split = pd.read_csv(paths["config"] / "split_tcga.tsv", sep="\t", dtype=str)
    split["exclude_from_eval"] = split["exclude_from_eval"].astype(str).str.lower().eq("true")
    test = split.loc[split["split"] == "TCGA-test"].reset_index(drop=True)
    met = split.loc[(split["split"] == "TCGA-met") & (~split["exclude_from_eval"])].reset_index(drop=True)
    if len(test) != 1872 or len(met) != 373:
        raise SystemExit(f"unexpected TCGA sizes test={len(test)} met={len(met)}")
    samples, tpm = load_tpm(paths["data_processed"] / "toil/tpm_G.h5")
    index = {sample: i for i, sample in enumerate(samples)}
    frames = {"TCGA-test": test, "TCGA-met": met}
    out_dirs = {"TCGA-test": root / "tcga_test", "TCGA-met": root / "tcga_met"}
    for name, frame in frames.items():
        block = sum_1e6(tpm[[index[sample] for sample in frame["sample"]]])
        z = rank_b0(block, b0_idx)
        k = k_matrix(paths, frame["sample"].tolist(), set_names)
        truth = organ_of_labels(frame["learning_label"].to_numpy())
        cluster = frame["patient"].to_numpy()
        preds = {}
        rows = []
        spec = {
            "BASE-Z": ("BASE_Z", z), "SA-Z": ("SA_Z", z), "NC-Z": ("NC_Z", z),
            "BASE-K": ("BASE_K", k), "SA-K": ("SA_K", k),
        }
        for method, (key, matrix) in spec.items():
            proba, pred = predict_organ(models[key], matrix)
            preds[method] = (proba, pred)
            metrics = metric_row(truth, pred, proba, [set() for _ in truth])
            metrics["method"] = method
            rows.append(metrics)
        pd.DataFrame(rows).to_csv(out_dirs[name] / "overall.tsv", sep="\t", index=False)
        paired = []
        for new, base in (("SA-Z", "BASE-Z"), ("SA-K", "BASE-K")):
            diff = paired_mean(preds[new][1] == truth, preds[base][1] == truth, cluster, SEED)
            mc = mcnemar_two_sided(preds[base][1] == truth, preds[new][1] == truth)
            paired.append({"comparison": f"{new} - {base}", "n": len(truth), **diff, **mc})
        pd.DataFrame(paired).to_csv(out_dirs[name] / "paired.tsv", sep="\t", index=False)
        if name == "TCGA-test":
            organ_rows = []
            for organ in sorted(set(truth.tolist())):
                mask = truth == organ
                organ_rows.append({
                    "organ": organ, "n": int(mask.sum()),
                    "BASE-Z_top1": float(np.mean(preds["BASE-Z"][1][mask] == organ)),
                    "SA-Z_top1": float(np.mean(preds["SA-Z"][1][mask] == organ)),
                })
            pd.DataFrame(organ_rows).to_csv(out_dirs[name] / "by_organ.tsv", sep="\t", index=False)
        else:
            write_met_types(frame, preds, out_dirs[name] / "by_sample_type.tsv")
        print("tcga", name, flush=True)
    del tpm
    score_tcga_variants(paths, root, models, b0_idx)


def score_tcga_variants(paths, root, models, b0_idx) -> None:
    split = pd.read_csv(paths["config"] / "split_tcga.tsv", sep="\t", dtype=str)
    test = split.loc[split["split"] == "TCGA-test"].reset_index(drop=True)
    if len(test) != 1872:
        raise SystemExit(f"unexpected TCGA-test size {len(test)}")
    samples, tpm = load_tpm(paths["data_processed"] / "toil/tpm_G.h5")
    index = {sample: i for i, sample in enumerate(samples)}
    block = sum_1e6(tpm[[index[sample] for sample in test["sample"]]])
    del tpm
    z = rank_b0(block, b0_idx)
    truth = organ_of_labels(test["learning_label"].to_numpy())
    rows = []
    for method in ["BASE-Z", "SA-Z", *variant_names()]:
        key = {"BASE-Z": "BASE_Z", "SA-Z": "SA_Z"}.get(method, method)
        proba, pred = predict_organ(models[key], z)
        metrics = metric_row(truth, pred, proba, [set() for _ in truth])
        rows.append({"method": method, "n": metrics["n"], "top1": metrics["top1"], "macro_f1": metrics["macro_f1"]})
    pd.DataFrame(rows).to_csv(root / "tcga_test" / "variant_overall.tsv", sep="\t", index=False)
    print("tcga variants", len(rows), flush=True)


def write_met_types(frame, preds, path: Path) -> None:
    clinical = read_skcm_clinical()[["sampleID", "tumor_tissue_site"]].drop_duplicates("sampleID")
    merged = frame[["sample", "learning_label"]].merge(clinical, left_on="sample", right_on="sampleID", how="left")
    raw = merged["tumor_tissue_site"].where(merged["learning_label"].eq("SKCM"), other="non-SKCM")
    group = np.where(raw.isin(NAMED_MET_SITES), raw, "기타")
    rows = []
    for label in [
        "Regional Lymph Node", "Distant Metastasis",
        "Regional Cutaneous or Subcutaneous Tissue (includes satellite and in-transit metastasis)",
        "기타",
    ]:
        mask = group == label
        row = {"sample_type": label, "n": int(mask.sum())}
        for method in ("BASE-Z", "SA-Z", "NC-Z", "BASE-K", "SA-K"):
            pred = preds[method][1]
            row[f"{method}_top1"] = float(np.mean(pred[mask] == organ_of_labels(frame["learning_label"].to_numpy())[mask])) if mask.any() else None
        rows.append(row)
    detail = raw.fillna("<NA>").value_counts()
    pd.DataFrame({"raw_value": detail.index.astype(str), "n": detail.to_numpy()}).to_csv(
        path.with_name("sample_type_raw_counts.tsv"), sep="\t", index=False
    )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def sim_tables(paths, root, models, b0_idx, mu, h_by_site, keep) -> None:
    mixes = pd.read_parquet(root / "sim_ext" / "mixes.parquet")
    tumors = pd.read_csv(root / "sim_ext" / "tumors.tsv", sep="\t", dtype=str)["tumor"].tolist()
    test = test_frame(paths)
    take_ids = test["sample"].tolist()
    truth_of = dict(zip(test["sample"], organ_of_labels(test["learning_label"].to_numpy())))
    truth = np.array([truth_of[sample] for sample in tumors], dtype=object)
    samples, tpm = load_tpm(paths["data_processed"] / "toil/tpm_G.h5")
    row = {sample: i for i, sample in enumerate(samples)}
    tumor_block = sum_1e6(tpm[[row[sample] for sample in tumors]])
    pure_z = rank_b0(tumor_block, b0_idx)
    host_ids = mixes["host"].drop_duplicates().tolist()
    host_mat = sum_1e6(tpm[[row[sample] for sample in host_ids]])
    host = {sample: host_mat[i] for i, sample in enumerate(host_ids)}
    del tpm, host_mat
    pure_methods = ["BASE-Z", "SA-Z", "SA-m0"] + [f"SA-LOHO-{tissue}" for tissue in POOL10]
    pure_pred = {}
    for method in pure_methods:
        key = {"BASE-Z": "BASE_Z", "SA-Z": "SA_Z"}.get(method, method)
        pure_pred[method] = predict_organ(models[key], pure_z)[0]
    rows = []
    delta_rows = []
    stomach_rows = []
    sarcoma_rows = []
    methods = SIM_BASE_MODELS + variant_names()
    for tissue, rho in mixes[["tissue", "rho"]].drop_duplicates().itertuples(index=False):
        sub = mixes.loc[(mixes["tissue"] == tissue) & np.isclose(mixes["rho"], rho)]
        if sub["tumor"].tolist() != tumors:
            raise SystemExit(f"mix order {tissue} {rho}")
        site = TISSUE_TO_SITE[tissue]
        native = site_native_organs(site)
        natives = [native for _ in tumors]
        hosts = np.stack([host[sample] for sample in sub["host"]])
        rho_v = np.float32(float(rho))
        mixed = rho_v * tumor_block + (np.float32(1.0) - rho_v) * hosts
        z = rank_b0(mixed, b0_idx)
        pred_map = predict_sim_methods(mixed, z, site, native, models, b0_idx, mu, h_by_site, keep)
        flag = substitute_site(site)
        for method in methods:
            proba, pred = pred_map[method]
            stats = metric_row(truth, pred, proba, natives)
            stats.update({
                "tissue": tissue, "site": site, "rho": float(rho), "method": method,
                "host_pull_rate": stats["host_rate"], "substitute": bool(flag and method in ("SC-Z", "LD-Z")),
            })
            rows.append(stats)
        if float(rho) in DELTA_RHOS:
            delta_rows.extend(delta_block(tissue, site, float(rho), truth, native, pred_map, pure_pred))
        if tissue == "Liver":
            stomach_rows.extend(rate_block("stomach", tissue, float(rho), truth, pred_map, STOMACH_MODELS))
        if tissue in ("Adipose - Subcutaneous", "Muscle - Skeletal"):
            sarcoma_rows.extend(rate_block("sarcoma", tissue, float(rho), truth, pred_map, SARCOMA_MODELS))
        print("sim", tissue, rho, flush=True)
    pd.DataFrame(rows).to_csv(root / "sim_ext" / "full.tsv", sep="\t", index=False)
    pd.DataFrame(delta_rows).to_csv(root / "mechanism" / "delta_values.tsv", sep="\t", index=False)
    summarize_delta(pd.DataFrame(delta_rows)).to_csv(root / "mechanism" / "delta_summary.tsv", sep="\t", index=False)
    pd.DataFrame(stomach_rows).to_csv(root / "mechanism" / "stomach_sim.tsv", sep="\t", index=False)
    pd.DataFrame(sarcoma_rows).to_csv(root / "mechanism" / "sarcoma_sim.tsv", sep="\t", index=False)


def predict_sim_methods(mixed, z, site, native, models, b0_idx, mu, h_by_site, keep):
    out = {}
    base_p, base_pred = predict_organ(models["BASE_Z"], z)
    sa_p, sa_pred = predict_organ(models["SA_Z"], z)
    out["BASE-Z"] = (base_p, base_pred)
    out["SA-Z"] = (sa_p, sa_pred)
    out["NC-Z"] = predict_organ(models["NC_Z"], z)
    masked, m1_pred = mask_native(base_p, [sorted(native) for _ in range(len(z))])
    out["M1-Z"] = (masked, m1_pred)
    model_site = SITE_TO_MODEL.get(site)
    if model_site is None:
        out["SC-Z"] = (sa_p, sa_pred)
        out["LD-Z"] = (base_p, base_pred)
    else:
        out["SC-Z"] = predict_organ(models[f"SC_{model_site}_Z"], z)
        cleaned, _beta = ld_apply(mixed, h_by_site[model_site], mu, keep)
        out["LD-Z"] = predict_organ(models["BASE_Z"], rank_b0(cleaned, b0_idx))
    for name in variant_names():
        out[name] = predict_organ(models[name], z)
    return out


def delta_block(tissue, site, rho, truth, native, pred_map, pure_pred):
    outside = np.array([label not in native for label in truth])
    methods = ["BASE-Z", "SA-m0", "SA-Z"]
    loho = f"SA-LOHO-{tissue}"
    if tissue in POOL10:
        methods.append(loho)
    native_idx = [ORGANS.index(organ) for organ in native if organ in ORGANS]
    rows = []
    for method in methods:
        mixed = pred_map[method][0]
        pure = pure_pred[method]
        for i in np.flatnonzero(outside):
            rows.append({
                "tissue": tissue, "site": site, "rho": rho, "method": method,
                "delta_native": log_mass(mixed[i], native_idx) - log_mass(pure[i], native_idx),
                "delta_true": log_mass(mixed[i], [ORGANS.index(truth[i])]) - log_mass(pure[i], [ORGANS.index(truth[i])]),
            })
    return rows


def log_mass(proba, indices) -> float:
    if not indices:
        total = 0.0
    else:
        total = float(np.sum(proba[indices]))
    return float(np.log(max(total, 1e-12)))


def summarize_delta(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if frame.empty:
        return pd.DataFrame(rows)
    for keys, sub in frame.groupby(["tissue", "site", "rho", "method"], sort=False):
        row = dict(zip(["tissue", "site", "rho", "method"], keys))
        row["n"] = int(len(sub))
        for column in ("delta_native", "delta_true"):
            values = sub[column].to_numpy(dtype=float)
            q25, q50, q75 = np.percentile(values, [25, 50, 75])
            row[f"{column}_q25"] = float(q25)
            row[f"{column}_q50"] = float(q50)
            row[f"{column}_q75"] = float(q75)
        rows.append(row)
    return pd.DataFrame(rows)


def rate_block(kind, tissue, rho, truth, pred_map, methods):
    if kind == "stomach":
        mask = np.array([label in GI_SET for label in truth])
    else:
        mask = truth == "Sarcoma"
    rows = []
    n = int(mask.sum())
    for method in methods:
        pred = pred_map[method][1][mask]
        sub_truth = truth[mask]
        if kind == "stomach":
            parts = error_parts(sub_truth, pred, [site_native_organs("liver") for _ in pred])
            gi_rate = float(parts["gi"].mean()) if n else None
            stomach_rate = float(np.mean(pred == "Stomach")) if n else None
            rows.append({
                "tissue": tissue, "rho": rho, "method": method, "n": n,
                "gi_internal_error_rate": gi_rate, "stomach_pred_rate": stomach_rate,
            })
        else:
            rows.append({
                "tissue": tissue, "rho": rho, "method": method, "n": n,
                "brain_pred_rate": float(np.mean(pred == "Brain")) if n else None,
            })
    return rows


def cohort_predictions(paths, root, models, b0_idx) -> None:
    genes, _b0, _o, _i, _s = load_gene_pack(paths)
    write_met500(paths, root, models, b0_idx)
    write_pog(paths, root, models, b0_idx, genes)


def write_met500(paths, root, models, b0_idx) -> None:
    meta = pd.read_csv(paths["config"] / "met500_eval_samples.tsv", sep="\t", dtype=str)
    beta = pd.read_csv(paths["results"] / "stage2" / "ld_beta_samples.tsv", sep="\t", dtype=str)
    meta = meta.merge(beta[["Sample_id", "beta"]], on="Sample_id", how="left")
    meta["organ"] = meta["cohort"].map(MET500_COHORT_TO_ORGAN)
    meta["standard_site"] = [map_site("met500_biopsy_tissue", value) for value in meta["biopsy_tissue"]]
    frame = meta.loc[meta["organ"].notna() & (meta["organ"] != "unmapped")].reset_index(drop=True)
    if len(frame) != 437:
        raise SystemExit(f"MET500 eval n={len(frame)}")
    samples, tpm = load_tpm(paths["data_processed"] / "met500/tpm_G.h5")
    row = {sample: i for i, sample in enumerate(samples)}
    block = sum_1e6(tpm[[row[sample] for sample in frame["Sample_id"]]])
    del tpm
    z = rank_b0(block, b0_idx)
    truth = frame["organ"].to_numpy(dtype=object)
    natives = [native_organs(site) for site in frame["standard_site"]]
    preds, pmax, proba = predict_variants(models, z)
    save_cohort(root / "posthoc" / "met500_samples.tsv", frame["Sample_id"], frame["sample_source"], truth, frame["standard_site"], frame["beta"], preds, pmax, None)
    np.save(root / "posthoc" / "met500_BASE-Z_proba.npy", proba["BASE-Z"])
    np.save(root / "posthoc" / "met500_SA-Z_proba.npy", proba["SA-Z"])
    write_metric_tables(root / "posthoc", "met500", truth, natives, frame["sample_source"].to_numpy(), frame["standard_site"].to_numpy(), preds, proba)
    print("met500", len(frame), flush=True)


def write_pog(paths, root, models, b0_idx, genes) -> None:
    print("load POG expression", flush=True)
    expr = load_pog570_expression(unlock=True)
    matrix, patient_ids, audit = collapse_to_g(expr, genes)
    del expr
    matrix = sum_1e6(matrix)
    (root / "posthoc" / "pog_feature_audit.json").write_text(pd.Series(audit).to_json())
    z = rank_b0(matrix, b0_idx)
    preds, pmax, proba = predict_variants(models, z)
    frozen = pd.read_parquet(paths["results"] / "stage3" / "pog570_predictions.parquet")
    frozen["patient_id"] = frozen["patient_id"].astype(str)
    check = check_frozen(patient_ids, proba, preds, frozen)
    (root / "posthoc" / "pog_frozen_check.json").write_text(pd.Series(check).to_json())
    print("pog frozen", check, flush=True)
    if check["max_abs_proba"] > 1e-6 or check["top1_mismatch"]:
        raise SystemExit("recomputed POG probabilities do not match the frozen predictions")
    labels = pd.read_csv(paths["config"] / "pog570_eval_labels.tsv", sep="\t", dtype=str)
    s1 = read_pog570_table("s1")
    s1["PATIENT_ID"] = s1["PATIENT_ID"].astype(str)
    order = pd.DataFrame({"patient_id": patient_ids})
    merged = order.merge(labels, left_on="patient_id", right_on="PATIENT_ID", how="left")
    merged = merged.merge(
        s1[["PATIENT_ID", "BIOPSY_SITE", "METASTATIC_OR_RECURRENCE", "TUMOUR_CONTENT"]],
        on="PATIENT_ID", how="left",
    )
    merged = merged.merge(frozen[["patient_id", "beta"]], on="patient_id", how="left")
    merged["standard_site"] = [map_site("pog570_biopsy_site", value) for value in merged["BIOPSY_SITE"]]
    keep = ~merged["organ"].isin(["exclude", "NA"])
    if int(keep.sum()) != 512:
        raise SystemExit(f"POG eval n={int(keep.sum())}")
    idx = np.flatnonzero(keep.to_numpy())
    sub_preds = {method: pred[idx] for method, pred in preds.items()}
    sub_pmax = {method: value[idx] for method, value in pmax.items()}
    sub_proba = {method: value[idx] for method, value in proba.items()}
    part = merged.loc[keep].reset_index(drop=True)
    truth = part["organ"].to_numpy(dtype=object)
    natives = [native_organs(site) for site in part["standard_site"]]
    extra = part[["METASTATIC_OR_RECURRENCE", "TUMOUR_CONTENT"]].reset_index(drop=True)
    save_cohort(root / "posthoc" / "pog_samples.tsv", part["patient_id"], part["patient_id"], truth, part["standard_site"], part["beta"], sub_preds, sub_pmax, extra)
    np.save(root / "posthoc" / "pog_BASE-Z_proba.npy", sub_proba["BASE-Z"])
    np.save(root / "posthoc" / "pog_SA-Z_proba.npy", sub_proba["SA-Z"])
    write_metric_tables(root / "posthoc", "pog", truth, natives, part["patient_id"].to_numpy(), part["standard_site"].to_numpy(), sub_preds, sub_proba)
    print("pog", len(part), flush=True)


def predict_variants(models, z):
    preds, pmax, proba = {}, {}, {}
    for method in single_ids(models):
        key = {"BASE-Z": "BASE_Z", "SA-Z": "SA_Z"}.get(method, method)
        organ, pred = predict_organ(models[key], z)
        preds[method] = pred
        proba[method] = organ.astype(np.float32)
        pmax[method] = organ.max(axis=1).astype(np.float32)
    return preds, pmax, proba


def check_frozen(patient_ids, proba, preds, frozen) -> dict:
    frozen = frozen.set_index("patient_id").loc[patient_ids]
    max_abs = 0.0
    top_mismatch = 0
    for method in ("BASE-Z", "SA-Z"):
        for organ in ORGANS:
            stored = frozen[f"{method}__p_{organ}"].to_numpy(dtype=np.float64)
            got = proba[method][:, ORGANS.index(organ)].astype(np.float64)
            max_abs = max(max_abs, float(np.max(np.abs(stored - got))))
        top_mismatch += int(np.sum(frozen[f"{method}__pred"].to_numpy(dtype=object) != preds[method]))
    return {"max_abs_proba": max_abs, "top1_mismatch": top_mismatch, "n": len(patient_ids)}


def save_cohort(path, ids, clusters, truth, sites, beta, preds, pmax, extra) -> None:
    frame = pd.DataFrame({
        "id": ids.to_numpy() if hasattr(ids, "to_numpy") else np.asarray(ids),
        "cluster": np.asarray(clusters),
        "truth": truth,
        "site": np.asarray(sites),
        "beta": np.asarray(beta),
    })
    for method, pred in preds.items():
        frame[f"pred__{method}"] = pred
        frame[f"pmax__{method}"] = pmax[method]
    if extra is not None:
        frame = pd.concat([frame, extra], axis=1)
    frame.to_csv(path, sep="\t", index=False)


def write_metric_tables(folder, prefix, truth, natives, clusters, sites, preds, proba) -> None:
    rows = []
    for method, pred in preds.items():
        row = metric_row(truth, pred, proba[method], natives)
        row["method"] = method
        rows.append(row)
    pd.DataFrame(rows).to_csv(folder / f"{prefix}_overall.tsv", sep="\t", index=False)
    site_rows = []
    risk = np.array([bool(native) and truth[i] not in native for i, native in enumerate(natives)])
    for site in sorted(set(sites.tolist())):
        mask = (sites == site) & risk
        if not mask.any():
            continue
        for method, pred in preds.items():
            pulled = np.array([pred[i] in natives[i] for i in np.flatnonzero(mask)])
            site_rows.append({
                "site": site, "method": method, "n_at_risk": int(mask.sum()),
                "host_rate": float(np.mean(pulled)),
            })
    pd.DataFrame(site_rows).to_csv(folder / f"{prefix}_by_site.tsv", sep="\t", index=False)
    paired = []
    for method in variant_names():
        for base in ("SA-Z", "BASE-Z"):
            top = paired_mean(preds[method] == truth, preds[base] == truth, clusters, SEED)
            host_new = error_parts(truth, preds[method], natives)["host"]
            host_base = error_parts(truth, preds[base], natives)["host"]
            risk_ids = clusters[risk]
            host = paired_mean(host_new[risk].astype(float), host_base[risk].astype(float), risk_ids, SEED)
            paired.append({
                "comparison": f"{method} - {base}", "n": len(truth), "n_at_risk": int(risk.sum()),
                "top1_diff": top["diff"], "top1_ci_low": top["ci_low"], "top1_ci_high": top["ci_high"],
                "top1_n_new_only": top["n_new_only"], "top1_n_base_only": top["n_base_only"],
                "host_diff": host["diff"], "host_ci_low": host["ci_low"], "host_ci_high": host["ci_high"],
                "host_n_new_only": host["n_new_only"], "host_n_base_only": host["n_base_only"],
            })
    pd.DataFrame(paired).to_csv(folder / f"{prefix}_paired.tsv", sep="\t", index=False)


def main() -> None:
    paths = load_paths()
    root = stage5_dir(paths)
    genes, b0_idx, _offsets, _indices, set_names = load_gene_pack(paths)
    mu, h_by_site, keep, ref_genes = load_references(paths["results"] / "stage3" / "references.npz")
    if ref_genes != genes:
        raise SystemExit("reference gene order does not match G")
    models = load_z_models(paths, root)
    print("TCGA", flush=True)
    tcga_tables(paths, root, models, b0_idx, set_names)
    print("simulation", flush=True)
    sim_tables(paths, root, models, b0_idx, mu, h_by_site, keep)
    print("cohorts", flush=True)
    cohort_predictions(paths, root, models, b0_idx)
    print("eval done", flush=True)


def tcga_variants_only() -> None:
    paths = load_paths()
    root = stage5_dir(paths)
    _genes, b0_idx, _offsets, _indices, _set_names = load_gene_pack(paths)
    models = {}
    model_dir = paths["results"] / "stage2" / "models"
    models["BASE_Z"] = joblib.load(model_dir / "BASE_Z.joblib")
    models["SA_Z"] = joblib.load(model_dir / "SA_Z.joblib")
    var_dir = root / "variants" / "models"
    for name in variant_names():
        models[name] = joblib.load(var_dir / f"{name}.joblib")
    score_tcga_variants(paths, root, models, b0_idx)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--tcga-variants-only":
        tcga_variants_only()
    else:
        main()
