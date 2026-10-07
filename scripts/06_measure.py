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
from joblib import Parallel, delayed
from scipy.stats import norm, rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from ks_fast import ks_matrix, pack_sets, right_and_tie
from loader import read_pog570_table, read_skcm_clinical
from mapping_rules import (
    MET500_COHORT_TO_ORGAN,
    POG570_COHORT_TO_ORGAN,
    PROJECT_TO_ORGAN,
    SIM_HOST_NATIVE_ORGANS,
    SITE_NATIVE,
    map_site,
)
from metrics import (
    ORGANS,
    accuracy_block,
    bootstrap_metrics,
    decompose,
    macro_f1,
    mask_native,
    organ_probability,
    per_class_recall,
    predict_from_proba,
    topk_hit,
)
from util import SEED, assert_disjoint, load_paths

C_GRID = (0.01, 0.1, 1.0, 10.0)
RHOS = (1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2)
SIM_TISSUES = ("Liver", "Lung", "Adipose - Subcutaneous")


def decode_samples(raw) -> list[str]:
    return [item.decode() if isinstance(item, bytes) else str(item) for item in raw]


def rank_z(values: np.ndarray) -> np.ndarray:
    """Within-sample rank-normal scores. Average ranks for ties."""
    out = np.empty(values.shape, dtype=np.float64)
    n_genes = values.shape[1]
    for start in range(0, values.shape[0], 500):
        block = values[start:start + 500]
        ranks = rankdata(block, method="average", axis=1)
        out[start:start + block.shape[0]] = norm.ppf((ranks - 0.5) / n_genes)
    return out


def select_genes(z_train: np.ndarray, n_keep: int = 5000) -> np.ndarray:
    # Higher variance first. Ties break toward the smaller gene index.
    variance = z_train.var(axis=0, ddof=1)
    order = np.lexsort((np.arange(variance.size), -variance))
    return order[:n_keep]


def fit_logit(X, y, C: float):
    # Final fit keeps the process-wide 16-thread BLAS limit.
    clf = LogisticRegression(C=C, solver="lbfgs", max_iter=5000, class_weight="balanced")
    clf.fit(X, y)
    return clf


def cv_select_C(X: np.ndarray, y: np.ndarray, scale: bool) -> dict:
    splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    folds = list(splitter.split(X, y))

    def one(C, fold, train_idx, valid_idx):
        with threadpool_limits(limits=1):
            Xtr, Xva = X[train_idx], X[valid_idx]
            ytr, yva = y[train_idx], y[valid_idx]
            if scale:
                scaler = StandardScaler()
                Xtr = scaler.fit_transform(Xtr)
                Xva = scaler.transform(Xva)
            clf = LogisticRegression(C=float(C), solver="lbfgs", max_iter=5000, class_weight="balanced")
            clf.fit(Xtr, ytr)
            pred = clf.predict(Xva)
            labels = sorted(set(yva.tolist()))
            score = float(f1_score_safe(yva, pred, labels))
            return {"C": float(C), "fold": int(fold), "macro_f1": score, "n_iter": int(np.max(clf.n_iter_))}

    jobs = [(C, fold, tr, va) for fold, (tr, va) in enumerate(folds) for C in C_GRID]
    rows = Parallel(n_jobs=8)(delayed(one)(*job) for job in jobs)
    frame = pd.DataFrame(rows)
    means = frame.groupby("C")["macro_f1"].mean().sort_index()
    best = float(means.max())
    chosen = float(means[means >= best - 1e-12].index.min())
    return {"folds": rows, "mean_macro_f1": {str(k): float(v) for k, v in means.items()}, "C": chosen}


def f1_score_safe(y_true, y_pred, labels) -> float:
    from sklearn.metrics import f1_score

    return float(f1_score(y_true, y_pred, average="macro", labels=labels, zero_division=0))


def native_sets(sites: list[str]) -> list[set[str]]:
    out = []
    for site in sites:
        tissues, projects, organs = SITE_NATIVE.get(site, SITE_NATIVE["other"])
        out.append(set(organs))
    return out


def confusion_top(truth: np.ndarray, pred: np.ndarray, k: int = 15) -> pd.DataFrame:
    frame = pd.DataFrame({"truth": truth, "pred": pred})
    frame = frame.loc[frame["truth"] != frame["pred"]]
    counts = frame.value_counts().reset_index(name="n").head(k)
    return counts


def main() -> None:
    paths = load_paths()
    plan = paths["config"] / "analysis_plan_A1.yaml"
    if not plan.exists():
        raise SystemExit("config/analysis_plan_A1.yaml is missing. Section 10 was not started.")
    results = paths["results"]
    results.mkdir(parents=True, exist_ok=True)
    split = pd.read_csv(paths["config"] / "split_tcga.tsv", sep="\t", dtype=str)
    split["exclude_from_eval"] = split["exclude_from_eval"].astype(str).str.lower().eq("true")
    assert_disjoint(
        split.loc[split["split"] == "TCGA-train", "patient"],
        split.loc[split["split"] == "TCGA-test", "patient"],
        "TCGA patients",
    )
    genes = pd.read_csv(paths["data_processed"] / "genes/G_symbols.txt", header=None)[0].astype(str).tolist()
    with h5py.File(paths["data_processed"] / "toil/tpm_G.h5", "r") as handle:
        toil_samples = decode_samples(handle["samples"][:])
        toil_tpm = handle["tpm"][:]
    sample_row = {sample: i for i, sample in enumerate(toil_samples)}
    train = split.loc[split["split"] == "TCGA-train"].copy()
    test = split.loc[split["split"] == "TCGA-test"].copy()
    met = split.loc[(split["split"] == "TCGA-met") & (~split["exclude_from_eval"])].copy()
    for frame, name in ((train, "train"), (test, "test"), (met, "met")):
        missing = [s for s in frame["sample"] if s not in sample_row]
        if missing:
            raise SystemExit(f"{name} samples missing from the Toil matrix: {len(missing)}")
    print("rank-normal TCGA train", len(train), flush=True)
    z_train = rank_z(toil_tpm[[sample_row[s] for s in train["sample"]]])
    gene_idx = select_genes(z_train, 5000)
    pd.Series([genes[i] for i in gene_idx], name="symbol").to_csv(results / "B0_genes.txt", index=False)
    X_train = z_train[:, gene_idx]
    y_train = train["learning_label"].to_numpy()
    print("CV B0", flush=True)
    cv_b0 = cv_select_C(X_train, y_train, scale=False)
    clf_b0 = fit_logit(X_train, y_train, cv_b0["C"])
    print("B0 C", cv_b0["C"], "n_iter", int(np.max(clf_b0.n_iter_)), flush=True)

    ks_tcga = pd.read_parquet(paths["data_processed"] / "features/ks_tcga.parquet")
    ks_train = ks_tcga.loc[train["sample"]].to_numpy(dtype=np.float64)
    set_names = list(ks_tcga.columns)
    print("CV B1", ks_train.shape, flush=True)
    cv_b1 = cv_select_C(ks_train, y_train, scale=True)
    scaler = StandardScaler().fit(ks_train)
    clf_b1 = fit_logit(scaler.transform(ks_train), y_train, cv_b1["C"])
    print("B1 C", cv_b1["C"], "n_iter", int(np.max(clf_b1.n_iter_)), flush=True)
    (results / "cv_summary.json").write_text(json.dumps({"B0": cv_b0, "B1": cv_b1}, indent=2))

    def project_and_organ(X_b0, ks_rows):
        proba_b0 = clf_b0.predict_proba(X_b0)
        organ_b0 = organ_probability(proba_b0, list(clf_b0.classes_))
        proba_b1 = clf_b1.predict_proba(scaler.transform(ks_rows))
        organ_b1 = organ_probability(proba_b1, list(clf_b1.classes_))
        return proba_b0, organ_b0, proba_b1, organ_b1

    def eval_split(frame: pd.DataFrame, level_labels: np.ndarray) -> dict:
        z = rank_z(toil_tpm[[sample_row[s] for s in frame["sample"]]])[:, gene_idx]
        ks_rows = ks_tcga.loc[frame["sample"]].to_numpy(dtype=np.float64)
        proba_b0, organ_b0, proba_b1, organ_b1 = project_and_organ(z, ks_rows)
        pred_b0 = clf_b0.predict(z)
        pred_b1 = clf_b1.predict(scaler.transform(ks_rows))
        organ_pred_b0 = predict_from_proba(organ_b0, ORGANS)
        organ_pred_b1 = predict_from_proba(organ_b1, ORGANS)
        organ_truth = np.array([PROJECT_TO_ORGAN[label] for label in level_labels], dtype=object)
        return {
            "project_B0": {**accuracy_block(level_labels, proba_b0, list(clf_b0.classes_), pred_b0), "recall": per_class_recall(level_labels, pred_b0)},
            "project_B1": {**accuracy_block(level_labels, proba_b1, list(clf_b1.classes_), pred_b1), "recall": per_class_recall(level_labels, pred_b1)},
            "organ_B0": {**accuracy_block(organ_truth, organ_b0, ORGANS, organ_pred_b0), "recall": per_class_recall(organ_truth, organ_pred_b0)},
            "organ_B1": {**accuracy_block(organ_truth, organ_b1, ORGANS, organ_pred_b1), "recall": per_class_recall(organ_truth, organ_pred_b1)},
            "organ_truth": organ_truth,
            "organ_pred_b0": organ_pred_b0,
            "organ_pred_b1": organ_pred_b1,
            "project_pred_b0": pred_b0,
            "project_pred_b1": pred_b1,
        }

    print("internal test", flush=True)
    test_eval = eval_split(test, test["learning_label"].to_numpy())
    print("internal met", flush=True)
    met_eval = eval_split(met, met["learning_label"].to_numpy())
    internal = {
        "TCGA-test": {k: v for k, v in test_eval.items() if k in {"project_B0", "project_B1", "organ_B0", "organ_B1"}},
        "TCGA-met": {k: v for k, v in met_eval.items() if k in {"project_B0", "project_B1", "organ_B0", "organ_B1"}},
        "TCGA-met_excluded_train_overlap": int(((split["split"] == "TCGA-met") & split["exclude_from_eval"]).sum()),
    }
    # Per-project accuracy on TCGA-met.
    met_project_rows = []
    for model, pred in (("B0", met_eval["project_pred_b0"]), ("B1", met_eval["project_pred_b1"])):
        truth = met["learning_label"].to_numpy()
        for project in sorted(set(truth.tolist())):
            mask = truth == project
            met_project_rows.append({"model": model, "project": project, "n": int(mask.sum()), "top1": float(np.mean(pred[mask] == truth[mask]))})
    pd.DataFrame(met_project_rows).to_csv(results / "tcga_met_by_project.tsv", sep="\t", index=False)

    skcm = read_skcm_clinical()
    skcm = skcm[["sampleID", "tumor_tissue_site"]].drop_duplicates("sampleID")
    skcm_rows = []
    for split_name, frame, pred_b0, pred_b1 in (
        ("TCGA-test", test, test_eval["project_pred_b0"], test_eval["project_pred_b1"]),
        ("TCGA-met", met, met_eval["project_pred_b0"], met_eval["project_pred_b1"]),
    ):
        sub = frame.loc[frame["learning_label"] == "SKCM", ["sample", "learning_label"]].copy()
        if sub.empty:
            continue
        merged = sub.merge(skcm, left_on="sample", right_on="sampleID", how="left")
        # Predictions are aligned to `frame`, not to the SKCM subset. Reindex.
        positions = [frame["sample"].tolist().index(s) for s in merged["sample"]]
        merged["pred_b0"] = pred_b0[positions]
        merged["pred_b1"] = pred_b1[positions]
        for site, part in merged.groupby(merged["tumor_tissue_site"].fillna("<NA>")):
            skcm_rows.append({
                "split": split_name,
                "tumor_tissue_site": site,
                "n": int(len(part)),
                "B0_top1": float(np.mean(part["pred_b0"].to_numpy() == "SKCM")),
                "B1_top1": float(np.mean(part["pred_b1"].to_numpy() == "SKCM")),
            })
    pd.DataFrame(skcm_rows).to_csv(results / "skcm_by_tumor_tissue_site.tsv", sep="\t", index=False)
    internal["skcm_clinical_join_note"] = "Joined on the 15-character sample id to SKCM_clinicalMatrix.sampleID. Field used: tumor_tissue_site."
    (results / "internal_eval.json").write_text(json.dumps(internal, indent=2))

    print("MET500", flush=True)
    met500 = pd.read_csv(paths["config"] / "met500_eval_samples.tsv", sep="\t", dtype=str)
    met500["organ"] = met500["cohort"].map(MET500_COHORT_TO_ORGAN)
    met500["standard_site"] = [map_site("met500_biopsy_tissue", v) for v in met500["biopsy_tissue"]]
    met500["library"] = met500["Sample_id"].str.extract(r"-(capt|poly)-", expand=False)
    met500["tc_num"] = pd.to_numeric(met500["tc"], errors="coerce")
    eval_df = met500.loc[met500["organ"].notna() & (met500["organ"] != "unmapped")].copy()
    ks_met = pd.read_parquet(paths["data_processed"] / "features/ks_met500.parquet")
    with h5py.File(paths["data_processed"] / "met500/tpm_G.h5", "r") as handle:
        met_samples = decode_samples(handle["samples"][:])
        met_tpm = handle["tpm"][:]
    met_row = {sample: i for i, sample in enumerate(met_samples)}
    missing = [s for s in eval_df["Sample_id"] if s not in ks_met.index or s not in met_row]
    if missing:
        raise SystemExit(f"MET500 eval ids missing from features: {len(missing)}")
    z_met = rank_z(met_tpm[[met_row[s] for s in eval_df["Sample_id"]]])[:, gene_idx]
    ks_met_rows = ks_met.loc[eval_df["Sample_id"]].to_numpy(dtype=np.float64)
    proba_b0, organ_b0, proba_b1, organ_b1 = project_and_organ(z_met, ks_met_rows)
    sites = eval_df["standard_site"].tolist()
    natives = native_sets(sites)
    truth = eval_df["organ"].to_numpy(dtype=object)
    # V0 host means from GTEx-ref.
    gtex = pd.read_csv(paths["config"] / "split_gtex.tsv", sep="\t", dtype=str)
    ks_gtex = pd.read_parquet(paths["data_processed"] / "features/ks_gtex_ref.parquet")
    tissue_mean = {}
    for tissue, sub in gtex.loc[gtex["split"] == "GTEx-ref"].groupby("tissue"):
        ids = [s for s in sub["sample"] if s in ks_gtex.index]
        if ids:
            tissue_mean[tissue] = ks_gtex.loc[ids].to_numpy(dtype=np.float64).mean(axis=0)
    site_mean = {}
    for site, (tissues, _, _) in SITE_NATIVE.items():
        present = [tissue_mean[t] for t in tissues if t in tissue_mean]
        if present:
            site_mean[site] = np.mean(np.vstack(present), axis=0)
    ks_v0 = ks_met_rows.copy()
    for i, site in enumerate(sites):
        host = site_mean.get(site)
        if host is not None:
            ks_v0[i] = (ks_met_rows[i] - 0.3 * host) / 0.7
    organ_v0 = organ_probability(clf_b1.predict_proba(scaler.transform(ks_v0)), list(clf_b1.classes_))
    pred = {
        "B0": predict_from_proba(organ_b0, ORGANS),
        "B1": predict_from_proba(organ_b1, ORGANS),
        "B1+V0": predict_from_proba(organ_v0, ORGANS),
    }
    masked_b0, pred["B0+M1"] = mask_native(organ_b0, [sorted(s) for s in natives])
    masked_b1, pred["B1+M1"] = mask_native(organ_b1, [sorted(s) for s in natives])
    proba = {"B0": organ_b0, "B0+M1": masked_b0, "B1": organ_b1, "B1+M1": masked_b1, "B1+V0": organ_v0}
    met_summary = {}
    for name in ("B0", "B0+M1", "B1", "B1+M1", "B1+V0"):
        print("bootstrap", name, flush=True)
        block = bootstrap_metrics(truth, pred[name], proba[name], ORGANS, natives)
        block["macro_f1"] = macro_f1(truth, pred[name])
        met_summary[name] = block
    (results / "met500_summary.json").write_text(json.dumps(met_summary, indent=2))

    # Stratified error tables for B0 and B1.
    eval_df = eval_df.reset_index(drop=True)
    eval_df["tc_tertile"] = pd.qcut(eval_df["tc_num"], 3, labels=["T1", "T2", "T3"], duplicates="drop").astype(str)
    strat_rows = []
    for model in ("B0", "B1"):
        for column in ("standard_site", "tc_tertile", "library"):
            for key, index in eval_df.groupby(column).groups.items():
                idx = np.asarray(list(index))
                stats = decompose(truth[idx], pred[model][idx], [natives[i] for i in idx])
                stats.update({"model": model, "factor": column, "level": key})
                strat_rows.append(stats)
    pd.DataFrame(strat_rows).to_csv(results / "met500_error_by_factor.tsv", sep="\t", index=False)

    # Biopsy-site-native truth accuracy.
    native_truth_rows = []
    is_native_truth = np.array([truth[i] in natives[i] and len(natives[i]) > 0 for i in range(len(truth))])
    for name in proba:
        if is_native_truth.any():
            native_truth_rows.append({
                "model": name,
                "n": int(is_native_truth.sum()),
                "top1": float(np.mean(pred[name][is_native_truth] == truth[is_native_truth])),
            })
        else:
            native_truth_rows.append({"model": name, "n": 0, "top1": None})
    pd.DataFrame(native_truth_rows).to_csv(results / "met500_native_truth_accuracy.tsv", sep="\t", index=False)

    confusion_top(truth, pred["B1"], 15).to_csv(results / "b1_confusion_top15.tsv", sep="\t", index=False)
    full_counts = pd.DataFrame({"truth": truth, "pred": pred["B1"]}).value_counts().reset_index(name="n")
    full_counts.to_csv(results / "b1_confusion_all.tsv", sep="\t", index=False)
    organ_rows = []
    for organ in sorted(set(truth.tolist())):
        mask = truth == organ
        organ_rows.append({
            "organ": organ,
            "n": int(mask.sum()),
            "B0_top1": float(np.mean(pred["B0"][mask] == truth[mask])),
            "B1_top1": float(np.mean(pred["B1"][mask] == truth[mask])),
        })
    pd.DataFrame(organ_rows).to_csv(results / "met500_by_organ.tsv", sep="\t", index=False)

    for model in ("B0", "B1", "B0+M1", "B1+M1", "B1+V0"):
        host_rows = []
        stats_host = pred[model]
        for i in range(len(truth)):
            if not natives[i] or truth[i] in natives[i] or stats_host[i] not in natives[i]:
                continue
            order = np.argsort(-proba[model][i])[:3]
            top3 = [f"{ORGANS[j]}:{proba[model][i, j]:.4f}" for j in order]
            host_rows.append({
                "sample_id": eval_df.loc[i, "Sample_id"],
                "sample_source": eval_df.loc[i, "sample_source"],
                "standard_site": sites[i],
                "tc": eval_df.loc[i, "tc"],
                "truth": truth[i],
                "pred": stats_host[i],
                "top3": "|".join(top3),
            })
        pd.DataFrame(host_rows).to_csv(results / f"host_errors_{model.replace('+','')}.tsv", sep="\t", index=False)
    host_b1 = pd.read_csv(results / "host_errors_B1.tsv", sep="\t") if (results / "host_errors_B1.tsv").exists() else pd.DataFrame()
    if len(host_b1):
        host_b1["standard_site"].value_counts().rename_axis("standard_site").reset_index(name="n").to_csv(
            results / "host_errors_B1_by_site.tsv", sep="\t", index=False
        )
    (results / "met500_eval_n.json").write_text(json.dumps({
        "n_eval_unit": int(len(met500)),
        "n_mapped": int(len(eval_df)),
        "n_unmapped": int((met500["organ"] == "unmapped").sum()),
        "tc_tertile_edges": [float(x) for x in pd.qcut(eval_df["tc_num"], 3, retbins=True, duplicates="drop")[1]],
    }, indent=2))

    print("simulation", flush=True)
    sim = run_simulation(
        test, toil_tpm, sample_row, genes, gene_idx, gtex, ks_tcga, set_names,
        clf_b0, clf_b1, scaler, paths,
    )
    sim.to_csv(results / "simulation.tsv", sep="\t", index=False)

    print("POG risk", flush=True)
    risk = pog_risk()
    risk.to_csv(results / "pog570_risk_by_site.tsv", sep="\t", index=False)
    print("measure done")


def run_simulation(test, toil_tpm, sample_row, genes, gene_idx, gtex, ks_tcga, set_names, clf_b0, clf_b1, scaler, paths):
    labels = test["learning_label"].to_numpy()
    organs = np.array([PROJECT_TO_ORGAN[x] for x in labels], dtype=object)
    n_take = min(1000, len(test))
    take, _ = train_test_split(np.arange(len(test)), train_size=n_take, stratify=labels, random_state=SEED)
    take = np.sort(take)
    tumor_ids = test["sample"].to_numpy()[take]
    tumor_organ = organs[take]
    tumor_tpm = toil_tpm[[sample_row[s] for s in tumor_ids]].astype(np.float64)
    tumor_sum = tumor_tpm.sum(axis=1, keepdims=True)
    tumor_sum[tumor_sum == 0] = np.nan
    tumor_tpm = np.nan_to_num(tumor_tpm / tumor_sum * 1e6, nan=0.0)
    packed = np.load(paths["data_processed"] / "genes/gene_sets.npz", allow_pickle=True)
    offsets, indices = packed["offsets"], packed["indices"]
    gtex_sim = gtex.loc[gtex["split"] == "GTEx-sim"]
    pools = {}
    for tissue in SIM_TISSUES:
        ids = gtex_sim.loc[gtex_sim["tissue"] == tissue, "sample"].tolist()
        rows = [sample_row[s] for s in ids if s in sample_row]
        pool = toil_tpm[rows].astype(np.float64)
        totals = pool.sum(axis=1, keepdims=True)
        totals[totals == 0] = np.nan
        pools[tissue] = np.nan_to_num(pool / totals * 1e6, nan=0.0)
    rng = np.random.default_rng(SEED)
    rows = []
    sim_scores = []
    b0_classes = list(clf_b0.classes_)
    b1_classes = list(clf_b1.classes_)
    for tissue in SIM_TISSUES:
        native = set(SIM_HOST_NATIVE_ORGANS[tissue])
        host_pool = pools[tissue]
        for rho in RHOS:
            host_idx = rng.choice(host_pool.shape[0], size=len(tumor_ids), replace=True)
            mixed = rho * tumor_tpm + (1.0 - rho) * host_pool[host_idx]
            z = rank_z(mixed.astype(np.float32))[:, gene_idx]
            right, tie = right_and_tie(mixed.astype(np.float32))
            ks = ks_matrix(right, tie, offsets, indices).astype(np.float64)
            organ_b0 = organ_probability(clf_b0.predict_proba(z), b0_classes)
            organ_b1 = organ_probability(clf_b1.predict_proba(scaler.transform(ks)), b1_classes)
            pred_b0 = predict_from_proba(organ_b0, ORGANS)
            pred_b1 = predict_from_proba(organ_b1, ORGANS)
            for model, pred in (("B0", pred_b0), ("B1", pred_b1)):
                wrong = pred[pred != tumor_organ]
                top = pd.Series(wrong).value_counts().head(3)
                denom = np.array([organ not in native for organ in tumor_organ])
                pulled = denom & np.array([p in native for p in pred])
                rows.append({
                    "tissue": tissue,
                    "rho": rho,
                    "model": model,
                    "n": int(len(tumor_organ)),
                    "top1": float(np.mean(pred == tumor_organ)),
                    "n_denom_non_native_truth": int(denom.sum()),
                    "host_pull_rate": float(pulled.sum() / denom.sum()) if denom.any() else None,
                    "wrong_top1": "" if top.empty else f"{top.index[0]}:{int(top.iloc[0])}",
                    "wrong_top2": "" if len(top) < 2 else f"{top.index[1]}:{int(top.iloc[1])}",
                    "wrong_top3": "" if len(top) < 3 else f"{top.index[2]}:{int(top.iloc[2])}",
                })
            print("sim", tissue, rho, flush=True)
            index = [f"{sample}|{tissue}|rho={rho}" for sample in tumor_ids]
            sim_scores.append(pd.DataFrame(ks.astype(np.float32), index=index, columns=set_names))
    sim_path = paths["data_processed"] / "features/ks_sim.parquet"
    pd.concat(sim_scores).to_parquet(sim_path)
    print("wrote", sim_path.name, flush=True)
    return pd.DataFrame(rows)


def pog_risk() -> pd.DataFrame:
    s1 = read_pog570_table("s1")
    s1["standard_site"] = [map_site("pog570_biopsy_site", v) for v in s1["BIOPSY_SITE"]]
    s1["organ"] = s1["ANALYSIS_COHORT"].map(POG570_COHORT_TO_ORGAN)
    rows = []
    for site, sub in s1.groupby("standard_site"):
        native = set(SITE_NATIVE.get(site, SITE_NATIVE["other"])[2])
        unmapped = sub["organ"].isna() | sub["organ"].eq("unmapped")
        if not native:
            at_risk = 0
            native_truth = 0
        else:
            at_risk = int(((~unmapped) & ~sub["organ"].isin(native)).sum())
            native_truth = int(((~unmapped) & sub["organ"].isin(native)).sum())
        rows.append({
            "standard_site": site,
            "n": int(len(sub)),
            "n_unmapped_diagnosis": int(unmapped.sum()),
            "n_at_risk": at_risk,
            "n_native_truth": native_truth,
            "native_organs": "|".join(sorted(native)),
        })
    return pd.DataFrame(rows).sort_values("n", ascending=False)


if __name__ == "__main__":
    main()
