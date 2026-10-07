#!/usr/bin/env python3
"""Post hoc ablation. The plan file must already be committed.

Does not set a thread cap and does not read n_threads.
Does not import loader.py or genes.py.
Mixture-row helpers match scripts/08_train.py; that module is not imported
because it sets thread caps on import.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from stage2_common import (  # noqa: E402
    N_MIX,
    RHO_HIGH,
    RHO_LOW,
    gtex_ref_ids,
    load_gene_pack,
    load_tpm,
    rank_z,
    sum_1e6,
    train_frame,
)
from stage2_fit import assemble_aug, fit_logit, repeat_labels  # noqa: E402
from stage5_lib import paired_mean, predict_organ, rank_b0  # noqa: E402
import importlib  # noqa: E402

score = importlib.import_module("31_score_A8")

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage9" / "ablation"
PLAN = ROOT / "config" / "analysis_plan_A12.yaml"


def project_paths() -> dict:
    root = ROOT
    out = {"root": root}
    for line in (root / "config" / "paths.yaml").read_text().splitlines():
        text = line.split("#", 1)[0].strip()
        if not text or ":" not in text:
            continue
        key, value = text.split(":", 1)
        key = key.strip()
        if key == "n_threads":
            continue
        value = value.strip().strip('"').strip("'")
        if key == "seed":
            out["seed"] = int(value)
        else:
            out[key] = root / value
    return out


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(1 << 20)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def draw_mixes(tumors, tissues, ids_by_tissue, rng, pooled: bool):
    records = []
    pool = [(tissue, sample) for tissue in tissues for sample in ids_by_tissue[tissue]]
    for tumor in tumors:
        for _ in range(N_MIX):
            if pooled:
                tissue, host = pool[int(rng.integers(0, len(pool)))]
            else:
                tissue = tissues[int(rng.integers(0, len(tissues)))]
                group = ids_by_tissue[tissue]
                host = group[int(rng.integers(0, len(group)))]
            records.append({
                "tumor": tumor,
                "host": host,
                "tissue": tissue,
                "rho": float(rng.uniform(RHO_LOW, RHO_HIGH)),
            })
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


def pipeline_table() -> pd.DataFrame:
    specs = [
        ("BASE-Z", "BASE_Z.joblib", 7486, "fixed at 0.1; not selected by cross-validation", "none", "none"),
        ("SA-Z", "SA_Z.joblib", 37430, "3-fold StratifiedGroupKFold macro-F1 on {0.03, 0.1, 0.3}; smallest C within 1e-12 of the best; report A2 chose 0.03", "Uniform(0.15, 1.0), 4 mixes", "StandardScaler"),
        ("NC-Z", "NC_Z.joblib", 8589, "same C as SA-Z", "none; normal classes added", "StandardScaler"),
        ("SC-Z liver", "SC_liver_Z.joblib", 37430, "same C as SA-Z", "liver tissue only, same rho rule", "StandardScaler"),
        ("IF20-Z", "IF20_Z.joblib", 7486, "same C as SA-Z", "none; sensitivity-filtered genes", "StandardScaler"),
    ]
    rows = []
    folder = ROOT / "results" / "stage2" / "models"
    for name, filename, n_rows, how, mixture, scale_note in specs:
        bundle = joblib.load(folder / filename)
        clf = bundle["clf"]
        scaler = bundle["scaler"]
        rows.append({
            "model": name,
            "preprocess": "none" if scaler is None else type(scaler).__name__,
            "C": float(clf.C),
            "class_weight": str(clf.class_weight),
            "solver": str(clf.solver),
            "max_iter": int(clf.max_iter),
            "multi_class": str(getattr(clf, "multi_class", "")),
            "n_features": int(clf.n_features_in_),
            "n_classes": int(len(clf.classes_)),
            "n_train_rows": int(n_rows),
            "C_selection": how,
            "mixture": mixture,
            "scaling_note": scale_note,
        })
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "pipeline_table.tsv", sep="\t", index=False)
    return frame


def coef_gap(stored, fitted) -> float:
    stored_classes = [str(c) for c in stored["clf"].classes_]
    fitted_classes = [str(c) for c in fitted["clf"].classes_]
    if stored_classes != fitted_classes:
        raise SystemExit(f"class order {stored_classes} vs {fitted_classes}")
    gap = np.max(np.abs(stored["clf"].coef_ - fitted["clf"].coef_))
    gap = max(gap, np.max(np.abs(stored["clf"].intercept_ - fitted["clf"].intercept_)))
    if stored["scaler"] is not None:
        gap = max(gap, np.max(np.abs(stored["scaler"].mean_ - fitted["scaler"].mean_)))
        gap = max(gap, np.max(np.abs(stored["scaler"].scale_ - fitted["scaler"].scale_)))
    return float(gap)


def feature_cohort(cohort: str, library: str) -> str:
    if cohort == "prad_su2c_2019":
        if library == "polyA":
            return "prad_su2c_2019__fpkm_polya"
        if library == "capture":
            return "prad_su2c_2019__fpkm_capture"
        raise SystemExit(f"SU2C library {library}")
    return cohort


def ensembl_base_to_symbol(paths) -> dict[str, str]:
    hgnc = pd.read_csv(paths["data_raw"] / "hgnc" / "hgnc_complete_set.txt", sep="\t", dtype=str, low_memory=False)
    symbols = set(hgnc.loc[hgnc["locus_group"] == "protein-coding gene", "symbol"].dropna().astype(str))
    probe = pd.read_csv(paths["data_raw"] / "toil" / "gencode.v23.annotation.gene.probemap", sep="\t", dtype=str)
    probe = probe.loc[probe["gene"].isin(symbols), ["id", "gene"]].drop_duplicates()
    by_base: dict[str, set[str]] = {}
    for ens, sym in zip(probe["id"].astype(str), probe["gene"].astype(str)):
        by_base.setdefault(ens.split(".")[0], set()).add(sym)
    return {base: next(iter(group)) for base, group in by_base.items() if len(group) == 1}


def collapse_to_g(expr: pd.DataFrame, genes: list[str], mapping: dict[str, str]):
    bases = pd.Index(expr.index.astype(str)).str.split(".").str[0]
    symbol = np.asarray(bases.map(mapping), dtype=object)
    keep = pd.notna(symbol)
    sub = expr.loc[keep].astype(np.float32)
    grouped = sub.groupby(symbol[keep], sort=False).sum()
    aligned = grouped.reindex(pd.Index(genes)).fillna(0).astype(np.float32)
    matrix = np.ascontiguousarray(aligned.to_numpy(dtype=np.float32).T)
    return matrix, [str(c) for c in expr.columns]


def flags(frame: pd.DataFrame, column: str, kind: str) -> tuple[np.ndarray, np.ndarray]:
    if kind == "host":
        part = frame.loc[frame["at_risk"]]
        pred = part[column].to_numpy(dtype=object)
        values = np.array([
            score.host_flag(truth, pred_i, native)
            for truth, pred_i, native in zip(part["organ"], pred, part["native"])
        ], dtype=float)
    elif kind == "top1":
        part = frame.loc[frame["in_top1"]] if "in_top1" in frame.columns else frame
        values = (part[column].to_numpy(dtype=object) == part["organ"].to_numpy(dtype=object)).astype(float)
    elif kind == "native":
        part = frame.loc[frame["in_native"]]
        values = (part[column].to_numpy(dtype=object) == part["organ"].to_numpy(dtype=object)).astype(float)
    else:
        raise SystemExit(kind)
    return values, part["patient_id"].to_numpy(dtype=object)


def score_frame(frame: pd.DataFrame, column: str) -> dict:
    row = score.metric_row(frame, column, column, column + "_top3", None)
    native = frame.loc[frame["in_native"]]
    row["native_truth_top1_set"] = float(np.mean(native[column].to_numpy(dtype=object) == native["organ"].to_numpy(dtype=object)))
    row["n_native_set"] = int(len(native))
    return row


def main() -> None:
    if not PLAN.exists():
        raise SystemExit("config/analysis_plan_A12.yaml is missing. Fitting was not started.")
    OUT.mkdir(parents=True, exist_ok=True)
    paths = project_paths()
    table = pipeline_table()
    print(table.to_string(index=False), flush=True)
    genes, b0_idx, _offsets, _indices, _names = load_gene_pack(paths)
    train = train_frame(paths)
    if len(train) != 7486:
        raise SystemExit(f"TCGA-train {len(train)}")
    tumors = train["sample"].tolist()
    train_pos = {sample: i for i, sample in enumerate(tumors)}
    seed_info = json.loads((paths["results"] / "stage2" / "augmentation_seeds.json").read_text())
    tissues = seed_info["host_pool_sorted"]
    print("load toil", flush=True)
    samples, tpm = load_tpm(paths["data_processed"] / "toil" / "tpm_G.h5")
    row = {sample: i for i, sample in enumerate(samples)}
    tumor_norm = sum_1e6(tpm[[row[sample] for sample in tumors]])
    ids_by_tissue = gtex_ref_ids(paths, tissues)
    bank = {}
    for tissue, ids in ids_by_tissue.items():
        block = sum_1e6(tpm[[row[sample] for sample in ids]])
        for i, sample in enumerate(ids):
            bank[sample] = block[i]
    del tpm
    print("pure Z", flush=True)
    z_pure = rank_z(tumor_norm).astype(np.float32)[:, b0_idx]
    y = train["learning_label"].to_numpy()
    print("regenerate SA mixes", flush=True)
    rng = np.random.default_rng(int(seed_info["seeds"]["SA"]))
    records = draw_mixes(tumors, tissues, ids_by_tissue, rng, pooled=False)
    stored_mix = pd.read_parquet(paths["results"] / "stage2" / "mixes_SA.parquet")
    made = pd.DataFrame(records)
    same_rows = (
        len(made) == len(stored_mix)
        and made["tumor"].tolist() == stored_mix["tumor"].astype(str).tolist()
        and made["host"].tolist() == stored_mix["host"].astype(str).tolist()
        and made["tissue"].tolist() == stored_mix["tissue"].astype(str).tolist()
        and np.allclose(made["rho"].to_numpy(), stored_mix["rho"].to_numpy(dtype=float), rtol=0, atol=0)
    )
    if not same_rows:
        (OUT / "regeneration_check.json").write_text(json.dumps({"mix_records_match": False}, indent=2))
        raise SystemExit("regenerated SA mixture records do not match the stage-2 parquet; MIX-Z0 was not fit")
    print("materialize", len(records), flush=True)
    mixed = materialize(tumor_norm, train_pos, bank, records)
    del tumor_norm, bank
    z_mix = rank_z(mixed).astype(np.float32)[:, b0_idx]
    del mixed
    x_aug = assemble_aug(z_pure, z_mix)
    del z_mix
    y_aug = repeat_labels(y)
    print("refit SA", x_aug.shape, flush=True)
    refit = fit_logit(x_aug, y_aug, 0.03, scale=True)
    stored = joblib.load(paths["results"] / "stage2" / "models" / "SA_Z.joblib")
    gap = coef_gap(stored, refit)
    print("coef gap", gap, flush=True)
    met = score.met500()
    print("MET500 features", flush=True)
    met_samples, met_tpm = load_tpm(paths["data_processed"] / "met500" / "tpm_G.h5")
    met_row = {sample: i for i, sample in enumerate(met_samples)}
    missing = [sample for sample in met["sample_id"] if sample not in met_row]
    if missing:
        raise SystemExit(f"MET500 ids missing {len(missing)}")
    met_z = rank_b0(sum_1e6(met_tpm[[met_row[sample] for sample in met["sample_id"]]]), b0_idx)
    del met_tpm
    _organ, met_pred = predict_organ(refit, met_z)
    mismatch = int(np.sum(met_pred != met["SA-Z"].to_numpy(dtype=object)))
    check = {"mix_records_match": True, "coef_max_abs": gap, "met500_top1_mismatch": mismatch, "pass": bool(gap <= 1e-6 and mismatch == 0)}
    (OUT / "regeneration_check.json").write_text(json.dumps(check, indent=2))
    print(check, flush=True)
    if not check["pass"]:
        raise SystemExit("regeneration check failed; MIX-Z0 was not fit")
    print("fit ablation", flush=True)
    models = {
        "PURE-Zs": fit_logit(z_pure, y, 0.03, scale=True),
        "PURE-Zs-w": fit_logit(z_pure, y, 0.15, scale=True),
        "MIX-Z0": fit_logit(x_aug, y_aug, 0.1, scale=False),
        "SA-Z": stored,
        "BASE-Z": joblib.load(paths["results"] / "stage2" / "models" / "BASE_Z.joblib"),
    }
    (OUT / "models").mkdir(parents=True, exist_ok=True)
    for name in ("PURE-Zs", "PURE-Zs-w", "MIX-Z0"):
        joblib.dump(models[name], OUT / "models" / f"{name}.joblib")
    del x_aug, z_pure
    print("POG features", flush=True)
    expr = pd.read_csv(paths["data_raw"] / "POG570" / "POG570_TPM_expression.txt.gz", sep="\t", index_col=0)
    pog_matrix, pog_ids = collapse_to_g(expr, genes, ensembl_base_to_symbol(paths))
    del expr
    pog_z = rank_b0(sum_1e6(pog_matrix), b0_idx)
    del pog_matrix
    pog = score.pog570()
    pog_pos = {sample: i for i, sample in enumerate(pog_ids)}
    if any(sample not in pog_pos for sample in pog["sample_id"]):
        raise SystemExit("POG sample missing from expression")
    pog_take = np.array([pog_pos[sample] for sample in pog["sample_id"]])
    print("aux features", flush=True)
    aux = score.auxiliary()
    lib = pd.read_parquet(
        paths["results"] / "stage7" / "confirm" / "predictions.parquet",
        columns=["cohort", "sample_id", "patient_id", "library"],
    )
    lib["sample_id"] = lib["sample_id"].astype(str)
    lib["patient_id"] = lib["patient_id"].astype(str)
    aux = aux.merge(lib.drop_duplicates(["cohort", "sample_id", "patient_id"]), on=["cohort", "sample_id", "patient_id"], how="left")
    if aux["library"].isna().any():
        raise SystemExit("aux library join failed")
    cache = {}
    rows_z = []
    for cohort, library, sample in zip(aux["cohort"], aux["library"], aux["sample_id"]):
        folder = feature_cohort(str(cohort), str(library))
        if folder not in cache:
            ids = [line for line in (paths["data_processed"] / "aux" / folder / "sample_ids.txt").read_text().splitlines() if line]
            cache[folder] = (np.load(paths["data_processed"] / "aux" / folder / "Z.npy"), {item: i for i, item in enumerate(ids)})
        matrix, index = cache[folder]
        if sample not in index:
            raise SystemExit(f"{folder} missing {sample}")
        rows_z.append(matrix[index[sample]])
    aux_z = np.vstack(rows_z).astype(np.float32)
    cohorts = {"MET500": (met, met_z), "POG570": (pog, pog_z[pog_take]), "aux_rnaseq": (aux, aux_z)}
    pred_rows = []
    metric_rows = []
    for cohort, (frame, matrix) in cohorts.items():
        for name, bundle in models.items():
            organ, pred = predict_organ(bundle, matrix)
            frame[name] = pred
            frame[name + "_top3"] = score.organ_top3(organ)
            part = frame[["sample_id", "patient_id"]].copy()
            part.insert(0, "cohort", cohort)
            part["method"] = name
            part["pred"] = pred
            pred_rows.append(part)
            row = score_frame(frame, name)
            row.update({"cohort": cohort, "method": name})
            metric_rows.append(row)
            print(cohort, name, row["top1"], row["host_rate"], row["native_truth_top1_set"], flush=True)
    predictions = pd.concat(pred_rows, ignore_index=True)
    pred_path = OUT / "predictions.tsv"
    predictions.to_csv(pred_path, sep="\t", index=False)
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "metrics.tsv", sep="\t", index=False)
    contrasts = [
        ("PURE-Zs", "BASE-Z", "PURE-Zs_minus_BASE-Z"),
        ("SA-Z", "PURE-Zs", "SA-Z_minus_PURE-Zs"),
        ("PURE-Zs-w", "BASE-Z", "PURE-Zs-w_minus_BASE-Z"),
        ("SA-Z", "PURE-Zs-w", "SA-Z_minus_PURE-Zs-w"),
        ("MIX-Z0", "BASE-Z", "MIX-Z0_minus_BASE-Z"),
        ("SA-Z", "MIX-Z0", "SA-Z_minus_MIX-Z0"),
    ]
    diff_rows = []
    for cohort, (frame, _matrix) in cohorts.items():
        for left, right, label in contrasts:
            for kind in ("host", "top1", "native"):
                a, patients = flags(frame, left, kind)
                b, patients_b = flags(frame, right, kind)
                if not np.array_equal(patients, patients_b):
                    raise SystemExit("patient alignment")
                stats = paired_mean(a, b, patients, seed=20261001)
                stats.update({"cohort": cohort, "contrast": label, "metric": kind, "n": int(len(a))})
                diff_rows.append(stats)
    diffs = pd.DataFrame(diff_rows)
    diffs.to_csv(OUT / "differences.tsv", sep="\t", index=False)
    share_rows = []
    passed = True
    for cohort, part in metrics.groupby("cohort"):
        rates = {row.method: float(row.host_rate) for row in part.itertuples(index=False)}
        gap_host = rates["BASE-Z"] - rates["SA-Z"]
        shares = {
            "PURE-Zs": (rates["PURE-Zs"] - rates["SA-Z"]) / gap_host,
            "PURE-Zs-w": (rates["PURE-Zs-w"] - rates["SA-Z"]) / gap_host,
            "MIX-Z0": (rates["BASE-Z"] - rates["MIX-Z0"]) / gap_host,
        }
        for name, value in shares.items():
            ok = bool(np.isfinite(value) and value >= 0.5)
            passed = passed and ok
            share_rows.append({"cohort": cohort, "share": name, "g": gap_host, "value": value, "at_least_0.5": ok, **{k: rates[k] for k in rates}})
    shares_frame = pd.DataFrame(share_rows)
    shares_frame.to_csv(OUT / "shares.tsv", sep="\t", index=False)
    (OUT / "reading.json").write_text(json.dumps({
        "rule": "all three shares >= 0.5 in MET500, POG570, and aux_rnaseq",
        "pass": bool(passed),
    }, indent=2))
    hash_path = ROOT / "config" / "ablation_A12_sha256.tsv"
    hash_path.write_text("path\tsha256\n" + f"results/stage9/ablation/predictions.tsv\t{sha256_file(pred_path)}\n")
    print("reading pass", passed, flush=True)


if __name__ == "__main__":
    main()
