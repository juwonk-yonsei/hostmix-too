"""Remaining stage-8 diagnostics that do not need SCOPE output."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path("/home/kangjw/WORK/Cisplatin/2026_NewProj/proj_A_host_too/scripts")))
import importlib
diagnostics = importlib.import_module("32_diagnostics_A8")

ROOT = Path("/home/kangjw/WORK/Cisplatin/2026_NewProj/proj_A_host_too")
OUT = ROOT / "results/stage8/diagnostics"
GI6 = ["Esophagus", "Stomach", "Colorectal", "Pancreas", "Biliary", "Liver"]
MARKERS = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "hepatocyte": ["ALB", "APOA1", "TF", "HP"],
    "stromal": ["COL1A1", "COL1A2", "DCN", "LUM"],
}


def write(frame: pd.DataFrame, name: str) -> None:
    diagnostics.write(frame, name)


def percentiles(matrix: pd.DataFrame) -> pd.DataFrame:
    from scipy.stats import rankdata
    values = matrix.to_numpy(dtype=float)
    ranks = np.empty_like(values, dtype=float)
    for column in range(values.shape[1]):
        ranks[:, column] = rankdata(values[:, column], method="average")
    return pd.DataFrame((ranks - 1) / max(values.shape[0] - 1, 1), index=matrix.index, columns=matrix.columns)


def paad(frame: pd.DataFrame) -> None:
    part = frame.loc[frame["cohort"].eq("paad_iatlas_prince_2022") & frame["include_confirm"] & frame["selected_eval"]].copy()
    rows = []
    for i, row in part.reset_index(drop=True).iterrows():
        record = {"sample_id": row["sample_id"], "patient_id": row["patient_id"]}
        for method in ("BASE-Z", "SA-Z", "LD-Z"):
            prob = np.array([row[f"{method}__p_{organ}"] for organ in diagnostics.ORGANS], dtype=float)
            order = np.argsort(-prob)
            names = [diagnostics.ORGANS[j] for j in order[:3]]
            record[f"{method}_top3"] = "|".join(f"{name}:{prob[diagnostics.ORGANS.index(name)]:.6g}" for name in names)
            pancreas = diagnostics.ORGANS.index("Pancreas")
            record[f"{method}_pancreas_rank"] = int(np.where(order == pancreas)[0][0] + 1)
            record[f"{method}_pancreas_probability"] = float(prob[pancreas])
        rows.append(record)
    detail = pd.DataFrame(rows)
    values = pd.read_csv(
        ROOT / "data/processed/aux/downloads/paad_iatlas_prince_2022__data_mrna_seq_expression.txt",
        sep="\t", usecols=["Hugo_Symbol", *detail["sample_id"].astype(str).tolist()],
    )
    values = (np.power(2.0, values.set_index("Hugo_Symbol").apply(pd.to_numeric, errors="coerce")) - 1).clip(lower=0)
    values = values.groupby(level=0).sum()
    pct = percentiles(values)
    for group, genes in MARKERS.items():
        present = [gene for gene in genes if gene in pct.index]
        detail[f"{group}_mean_percentile"] = pct.loc[present, detail["sample_id"].astype(str)].mean(axis=0).to_numpy()
        detail[f"{group}_n_genes"] = len(present)
    detail["epithelial_tertile"] = pd.qcut(detail["epithelial_mean_percentile"], 3, labels=["T1", "T2", "T3"], duplicates="drop").astype(str)
    write(detail, "paad_ranks.tsv")
    tertile = detail.groupby("epithelial_tertile")["SA-Z_pancreas_rank"].median().rename("median_SA-Z_pancreas_rank").reset_index()
    tertile["n"] = detail.groupby("epithelial_tertile").size().to_numpy()
    write(tertile, "paad_epithelial_tertiles.tsv")


def gi_hub(frame: pd.DataFrame) -> None:
    rows = []

    def add(analysis, cohort, truth, pred, method):
        truth = np.asarray(truth, dtype=object)
        pred = np.asarray(pred, dtype=object)
        for gi_truth, mask in ((True, np.isin(truth, GI6)), (False, ~np.isin(truth, GI6))):
            sub = pred[mask]
            row = {"analysis": analysis, "cohort": cohort, "method": method, "truth_in_gi6": bool(gi_truth), "n": int(mask.sum())}
            for organ in GI6:
                row[f"n_pred_{organ}"] = int(np.sum(sub == organ))
            row["esophagus_absorption"] = float(np.mean(sub == "Esophagus")) if (not gi_truth or True) and len(sub) else np.nan
            rows.append(row)

    aux = frame.loc[frame["include_confirm"] & frame["selected_eval"] & frame["cohort"].isin(diagnostics.RNA)]
    micro = frame.loc[frame["layer"].eq("마이크로어레이") & frame["selected_eval"]]
    for method in ("BASE-Z", "SA-Z", "LD-Z"):
        add("aux_rnaseq", "all", aux["organ"], aux[f"{method}__pred"], method)
        for cohort, sub in aux.groupby("cohort"):
            add("aux_rnaseq", cohort, sub["organ"], sub[f"{method}__pred"], method)
        add("microarray", "all", micro["organ"], micro[f"{method}__pred"], method)
    pog = pd.read_parquet(ROOT / "results/stage3/pog570_predictions.parquet")
    labels = pd.read_csv(ROOT / "config/pog570_eval_labels.tsv", sep="\t", dtype=str)
    labels = labels.loc[~labels["organ"].isin(["exclude", "NA"])]
    merged = labels.merge(pog, left_on="PATIENT_ID", right_on="patient_id")
    met = pd.read_csv(ROOT / "results/stage5/posthoc/met500_samples.tsv", sep="\t", dtype=str)
    for method, column in (("BASE-Z", "BASE-Z__pred"), ("SA-Z", "SA-Z__pred"), ("LD-Z", "LD-Z__pred")):
        add("POG570", "POG570", merged["organ"], merged[column], method)
    for method, column in (("BASE-Z", "pred__BASE-Z"), ("SA-Z", "pred__SA-Z")):
        add("MET500", "MET500", met["truth"], met[column], method)
    tools = pd.read_csv(ROOT / "results/stage8/external/tables/sample_predictions.tsv", sep="\t", dtype=str)
    aux_ids = set(aux["sample_id"].astype(str))
    tool_aux = tools.loc[tools["analysis_cohort"].eq("aux_rnaseq") & tools["sample_id"].isin(aux_ids)]
    for method in ("SCOPE", "CUP-AI-Dx"):
        add("aux_rnaseq", "all", tool_aux["organ"], tool_aux[method], method)
        for cohort, sub in tool_aux.groupby("cohort"):
            add("aux_rnaseq", cohort, sub["organ"], sub[method], method)
        add("MET500", "MET500", tools.loc[tools["analysis_cohort"].eq("MET500"), "organ"], tools.loc[tools["analysis_cohort"].eq("MET500"), method], method)
        add("POG570", "POG570", tools.loc[tools["analysis_cohort"].eq("POG570"), "organ"], tools.loc[tools["analysis_cohort"].eq("POG570"), method], method)
    # esophagus absorption among non-esophagus truth is the gi_truth False row's esophagus rate only if GI6 truth includes esophagus.
    # Recompute the requested rate on truth != Esophagus.
    extra = []
    blocks = {
        ("aux_rnaseq", "all", "BASE-Z"): (aux["organ"], aux["BASE-Z__pred"]),
        ("aux_rnaseq", "all", "SA-Z"): (aux["organ"], aux["SA-Z__pred"]),
        ("aux_rnaseq", "all", "LD-Z"): (aux["organ"], aux["LD-Z__pred"]),
        ("POG570", "POG570", "BASE-Z"): (merged["organ"], merged["BASE-Z__pred"]),
        ("POG570", "POG570", "SA-Z"): (merged["organ"], merged["SA-Z__pred"]),
        ("POG570", "POG570", "LD-Z"): (merged["organ"], merged["LD-Z__pred"]),
        ("MET500", "MET500", "BASE-Z"): (met["truth"], met["pred__BASE-Z"]),
        ("MET500", "MET500", "SA-Z"): (met["truth"], met["pred__SA-Z"]),
        ("microarray", "all", "BASE-Z"): (micro["organ"], micro["BASE-Z__pred"]),
        ("microarray", "all", "SA-Z"): (micro["organ"], micro["SA-Z__pred"]),
        ("microarray", "all", "LD-Z"): (micro["organ"], micro["LD-Z__pred"]),
        ("aux_rnaseq", "all", "SCOPE"): (tool_aux["organ"], tool_aux["SCOPE"]),
        ("aux_rnaseq", "all", "CUP-AI-Dx"): (tool_aux["organ"], tool_aux["CUP-AI-Dx"]),
        ("MET500", "MET500", "SCOPE"): (tools.loc[tools["analysis_cohort"].eq("MET500"), "organ"], tools.loc[tools["analysis_cohort"].eq("MET500"), "SCOPE"]),
        ("MET500", "MET500", "CUP-AI-Dx"): (tools.loc[tools["analysis_cohort"].eq("MET500"), "organ"], tools.loc[tools["analysis_cohort"].eq("MET500"), "CUP-AI-Dx"]),
        ("POG570", "POG570", "SCOPE"): (tools.loc[tools["analysis_cohort"].eq("POG570"), "organ"], tools.loc[tools["analysis_cohort"].eq("POG570"), "SCOPE"]),
        ("POG570", "POG570", "CUP-AI-Dx"): (tools.loc[tools["analysis_cohort"].eq("POG570"), "organ"], tools.loc[tools["analysis_cohort"].eq("POG570"), "CUP-AI-Dx"]),
    }
    for (analysis, cohort, method), (truth, pred) in blocks.items():
        truth = np.asarray(truth, dtype=object)
        pred = np.asarray(pred, dtype=object)
        mask = truth != "Esophagus"
        extra.append({
            "analysis": analysis, "cohort": cohort, "method": method,
            "n_truth_not_esophagus": int(mask.sum()),
            "esophagus_absorption": float(np.mean(pred[mask] == "Esophagus")) if mask.any() else np.nan,
        })
    write(pd.DataFrame(rows), "gi_hub_counts.tsv")
    write(pd.DataFrame(extra), "esophagus_absorption.tsv")


def beta_vs_top1(frame: pd.DataFrame) -> None:
    missing = {}
    path = ROOT / "results/stage6/ld_inputs.json"
    for record in json.loads(path.read_text()):
        missing[record.get("cohort")] = record.get("n_ld_missing")
    rows = []
    use = frame.loc[frame["selected_eval"]]
    for cohort, sub in use.groupby("cohort"):
        base = float((sub["BASE-Z__pred"] == sub["organ"]).mean())
        sa = float((sub["SA-Z__pred"] == sub["organ"]).mean())
        rows.append({
            "cohort": cohort, "n": int(len(sub)), "top1_BASE-Z": base, "top1_SA-Z": sa,
            "diff_sa_minus_base": sa - base, "n_ld_missing": missing.get(cohort),
        })
    write(pd.DataFrame(rows), "beta_missing_vs_top1.tsv")


def tcga_paad_markers() -> None:
    split = pd.read_csv(ROOT / "config/split_tcga.tsv", sep="\t", dtype=str)
    test = split.loc[split["split"].eq("TCGA-test") & split["learning_label"].isin(["PAAD", "COADREAD"])]
    genes = (ROOT / "data/processed/genes/G_symbols.txt").read_text().splitlines()
    gene_index = {gene: i for i, gene in enumerate(genes)}
    with h5py.File(ROOT / "data/processed/toil/tpm_G.h5", "r") as handle:
        samples = [item.decode() if isinstance(item, bytes) else str(item) for item in handle["samples"][:]]
        sample_index = {sample: i for i, sample in enumerate(samples)}
        wanted = [sample for sample in test["sample"] if sample in sample_index]
        order = np.argsort([sample_index[sample] for sample in wanted])
        sorted_samples = [wanted[i] for i in order]
        rows = handle["tpm"][[sample_index[sample] for sample in sorted_samples]]
    matrix = pd.DataFrame(rows, index=sorted_samples, columns=genes).T
    pct = percentiles(matrix)
    labels = test.set_index("sample")["learning_label"]
    out_rows = []
    for label, samples_of in labels.groupby(labels):
        present_samples = [sample for sample in samples_of.index if sample in pct.columns]
        for group, group_genes in {**MARKERS, "colon": ["CDX1", "CDX2", "VIL1", "CDH17", "SATB2", "LGALS4", "CEACAM5", "KRT20"]}.items():
            present = [gene for gene in group_genes if gene in gene_index]
            if not present_samples or not present:
                continue
            out_rows.append({
                "learning_label": label, "marker_group": group, "n_samples": len(present_samples),
                "n_genes": len(present), "denominator": "G genes in tpm_G.h5",
                "mean_percentile": float(pct.loc[present, present_samples].to_numpy().mean()),
            })
    write(pd.DataFrame(out_rows), "tcga_test_marker_percentiles.tsv")


def main() -> None:
    frame = diagnostics.confirm()
    diagnostics.effect_bootstrap(frame)
    paad(frame)
    gi_hub(frame)
    beta_vs_top1(frame)
    tcga_paad_markers()
    print("extra diagnostics", OUT, flush=True)


if __name__ == "__main__":
    main()
