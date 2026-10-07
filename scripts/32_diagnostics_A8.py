"""Stage-8 post-hoc diagnostics. Reads existing predictions and features only."""
from __future__ import annotations

import sys
from pathlib import Path

import h5py
import joblib
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path("/home/kangjw/WORK/Cisplatin/2026_NewProj/proj_A_host_too/scripts")))
from mapping_rules import PROJECT_TO_ORGAN, SITE_NATIVE, map_site  # noqa: E402

ROOT = Path("/home/kangjw/WORK/Cisplatin/2026_NewProj/proj_A_host_too")
OUT = ROOT / "results/stage8/diagnostics"
NOTE = "사후"
ORGANS = sorted(set(PROJECT_TO_ORGAN.values()))
RNA = [
    "blca_iatlas_imvigor210_2017", "brca_iatlas_anders_2022", "mel_dfci_2019",
    "paad_iatlas_prince_2022", "GSE50760", "prad_su2c_2019",
]
GI = ["Esophagus", "Stomach", "Colorectal", "Pancreas", "Biliary", "Liver"]
GROUPS = {
    "real_host_tissue": {"adrenal", "brain", "liver", "lung", "skin", "soft_tissue"},
    "proxy_host": {"bone", "lymph_node"},
    "pool_out": {"kidney", "pancreas", "bladder", "thyroid", "ovary", "breast", "stomach", "colon_rectum", "prostate", "esophagus", "head_neck", "cervix"},
    "other": {"pleura"},
}
MARKERS = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "hepatocyte": ["ALB", "APOA1", "TF", "HP"],
    "stromal": ["COL1A1", "COL1A2", "DCN", "LUM"],
}
COLON = ["CDX1", "CDX2", "VIL1", "CDH17", "SATB2", "LGALS4", "CEACAM5", "KRT20"]


def write(frame: pd.DataFrame, name: str) -> None:
    frame.insert(0, "note", NOTE)
    frame.to_csv(OUT / name, sep="\t", index=False, float_format="%.16g")


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin(["true", "1"])


def natives(series: pd.Series) -> list[list[str]]:
    out = []
    for value in series.fillna("").astype(str):
        out.append([item for item in value.split("|") if item and item != "nan"])
    return out


def host(truth, pred, native) -> bool:
    return bool(native) and truth not in native and pred in native


def confirm() -> pd.DataFrame:
    frame = pd.read_parquet(ROOT / "results/stage7/confirm/predictions.parquet")
    for column in ("include_confirm", "selected_eval", "selected_risk", "selected_native", "in_eval", "in_risk", "in_native"):
        if column in frame.columns:
            frame[column] = as_bool(frame[column])
    return frame


def top_classes(frame: pd.DataFrame, method: str, label: str) -> list[dict]:
    pred = frame[f"{method}__pred"].astype(str)
    prob = frame[[f"{method}__p_{organ}" for organ in ORGANS]].to_numpy(dtype=float)
    pmax = prob.max(axis=1)
    rows = []
    for name, n in pred.value_counts().head(5).items():
        mask = pred.eq(name).to_numpy()
        rows.append({
            "slice": label, "method": method, "predicted": name, "n": int(n),
            "mean_max_probability": float(pmax[mask].mean()), "n_slice": int(len(frame)),
        })
    return rows


def imvigor_confusion(frame: pd.DataFrame) -> None:
    part = frame.loc[frame["cohort"].eq("blca_iatlas_imvigor210_2017") & frame["include_confirm"]]
    rows = []
    for label, mask in (
        ("native_truth", part["selected_native"]),
        ("eval", part["selected_eval"]),
    ):
        sub = part.loc[mask]
        for method in ("BASE-Z", "SA-Z", "LD-Z"):
            rows.extend(top_classes(sub, method, label))
    write(pd.DataFrame(rows), "imvigor_top_classes.tsv")


def assay_table() -> None:
    rows = []
    root = ROOT / "results/stage5/aux/cohorts"
    for path in sorted(root.glob("*/samples.tsv")):
        samples = pd.read_csv(path, sep="\t", dtype=str)
        subtype = [column for column in samples.columns if "subtype" in column.lower() or "lund" in column.lower()]
        preservation = samples["SPECIMEN_PRESERVATION_TYPE"].fillna("").value_counts().to_dict() if "SPECIMEN_PRESERVATION_TYPE" in samples.columns else {}
        rows.append({
            "cohort": path.parent.name,
            "n_rows": int(len(samples)),
            "library_values": "|".join(f"{k}:{v}" for k, v in samples["library"].fillna("").value_counts().items()) if "library" in samples.columns else "",
            "preservation_values": "|".join(f"{k}:{v}" for k, v in preservation.items()),
            "subtype_columns": "|".join(subtype),
        })
    write(pd.DataFrame(rows), "assay_fields.tsv")


def tcga_top1() -> None:
    organs = pd.read_csv(ROOT / "results/stage5/tcga_test/by_organ.tsv", sep="\t")
    keep = organs.loc[organs["organ"].isin(["Bladder", "Pancreas"])].copy()
    keep["class_distribution"] = "not stored in results/stage5/tcga_test"
    write(keep, "tcga_test_top1.tsv")


def su2c_pairs(frame: pd.DataFrame) -> None:
    part = frame.loc[frame["cohort"].eq("prad_su2c_2019")].copy()
    rows = []
    for patient, sub in part.groupby("patient_id"):
        libraries = set(sub["library"].astype(str))
        if not {"polyA", "capture"} <= libraries:
            continue
        poly = sub.loc[sub["library"].eq("polyA")].iloc[0]
        cap = sub.loc[sub["library"].eq("capture")].iloc[0]
        row = {"patient_id": patient}
        for method in ("BASE-Z", "SA-Z"):
            row[f"{method}_polyA"] = poly[f"{method}__pred"]
            row[f"{method}_capture"] = cap[f"{method}__pred"]
            row[f"{method}_agree"] = poly[f"{method}__pred"] == cap[f"{method}__pred"]
            row[f"{method}_polyA_correct"] = poly[f"{method}__pred"] == poly["organ"]
            row[f"{method}_capture_correct"] = cap[f"{method}__pred"] == cap["organ"]
        rows.append(row)
    detail = pd.DataFrame(rows)
    write(detail, "su2c_library_pairs.tsv")
    summary = []
    for method in ("BASE-Z", "SA-Z"):
        summary.append({
            "method": method, "n_pairs": int(len(detail)),
            "top1_polyA": float(detail[f"{method}_polyA_correct"].mean()) if len(detail) else np.nan,
            "top1_capture": float(detail[f"{method}_capture_correct"].mean()) if len(detail) else np.nan,
            "agreement": float(detail[f"{method}_agree"].mean()) if len(detail) else np.nan,
        })
    write(pd.DataFrame(summary), "su2c_library_summary.tsv")


def site_group(site: str) -> str:
    for name, sites in GROUPS.items():
        if site in sites:
            return name
    return "ungrouped"


def native_of(site: str) -> list[str]:
    return list(SITE_NATIVE.get(str(site), ([], [], []))[2])


def redistribution(frame: pd.DataFrame) -> None:
    rows = []

    def record(analysis, cohort, site, truth, base, sa, native):
        if host(truth, base, native) and not host(truth, sa, native):
            rows.append({
                "analysis": analysis, "cohort": cohort, "site_group": site_group(str(site)),
                "truth": truth, "BASE-Z": base, "SA-Z": sa,
                "SA_correct": sa == truth, "SA_other_error": sa != truth,
            })

    aux = frame.loc[frame["include_confirm"] & frame["selected_risk"] & frame["cohort"].isin(RNA)]
    native = natives(aux["native_organs"])
    for cohort, site, truth, base, sa, organs in zip(
        aux["cohort"], aux["standard_site"].fillna(""), aux["organ"], aux["BASE-Z__pred"], aux["SA-Z__pred"], native,
    ):
        record("aux_rnaseq", cohort, site, truth, base, sa, organs)
    met = pd.read_csv(ROOT / "results/stage5/posthoc/met500_samples.tsv", sep="\t", dtype=str)
    for site, truth, base, sa in zip(met["site"], met["truth"], met["pred__BASE-Z"], met["pred__SA-Z"]):
        record("MET500", "MET500", site, truth, base, sa, native_of(site))
    labels = pd.read_csv(ROOT / "config/pog570_eval_labels.tsv", sep="\t", dtype=str)
    labels = labels.loc[~labels["organ"].isin(["exclude", "NA"])]
    demo = pd.read_excel(ROOT / "data/raw/POG570/Table_S1_Demographics.xlsx", dtype=str)
    demo["standard_site"] = [map_site("pog570_biopsy_site", value) for value in demo["BIOPSY_SITE"]]
    pred = pd.read_parquet(ROOT / "results/stage3/pog570_predictions.parquet")
    merged = labels.merge(demo[["PATIENT_ID", "standard_site"]], on="PATIENT_ID").merge(pred, left_on="PATIENT_ID", right_on="patient_id")
    for site, truth, base, sa in zip(merged["standard_site"], merged["organ"], merged["BASE-Z__pred"], merged["SA-Z__pred"]):
        record("POG570", "POG570", site, truth, base, sa, native_of(site))
    detail = pd.DataFrame(rows)
    write(detail, "error_shift_cases.tsv")
    summary = []
    for keys, sub in detail.groupby(["analysis", "cohort", "site_group"], sort=False):
        summary.append({
            "analysis": keys[0], "cohort": keys[1], "site_group": keys[2], "n_favor": int(len(sub)),
            "n_sa_correct": int(sub["SA_correct"].sum()),
            "n_sa_other_error": int((~sub["SA_correct"]).sum()),
            "other_error_classes": "|".join(f"{k}:{v}" for k, v in sub.loc[~sub["SA_correct"], "SA-Z"].value_counts().items()),
        })
    write(pd.DataFrame(summary), "error_shift_summary.tsv")
    totals = detail.groupby("analysis").agg(n_favor=("SA_correct", "size"), n_sa_correct=("SA_correct", "sum")).reset_index()
    write(totals, "error_shift_totals.tsv")


def effect_bootstrap(frame: pd.DataFrame) -> None:
    part = frame.loc[frame["include_confirm"] & frame["cohort"].isin(RNA)].copy()
    part["native"] = natives(part["native_organs"])
    rows = []
    cohorts = ["none"] + list(RNA) + ["drop_blca_and_paad"]
    specs = {
        "AH1": ("selected_risk", "host"),
        "AH2": ("selected_eval", "top1"),
        "AH3": ("selected_native", "top1"),
    }
    for drop in cohorts:
        use = part
        if drop == "none":
            use = part
        elif drop == "drop_blca_and_paad":
            use = part.loc[~part["cohort"].isin(["blca_iatlas_imvigor210_2017", "paad_iatlas_prince_2022"])]
        elif drop in RNA:
            use = part.loc[part["cohort"] != drop]
        for hypothesis, (flag, kind) in specs.items():
            sub = use.loc[use[flag]].copy()
            if kind == "host":
                left = np.array([host(t, p, n) for t, p, n in zip(sub["organ"], sub["SA-Z__pred"], sub["native"])])
                right = np.array([host(t, p, n) for t, p, n in zip(sub["organ"], sub["BASE-Z__pred"], sub["native"])])
            else:
                left = (sub["SA-Z__pred"].to_numpy() == sub["organ"].to_numpy())
                right = (sub["BASE-Z__pred"].to_numpy() == sub["organ"].to_numpy())
            patients = sub["patient_id"].astype(str).to_numpy()
            unique = pd.unique(patients)
            grouped = {patient: np.flatnonzero(patients == patient) for patient in unique}
            keys = list(grouped)
            rng = np.random.default_rng(20261001)
            draws = np.empty(2000)
            for draw in range(2000):
                chosen = rng.choice(len(keys), size=len(keys), replace=True)
                idx = np.concatenate([grouped[keys[i]] for i in chosen]) if len(keys) else np.array([], dtype=int)
                draws[draw] = left[idx].mean() - right[idx].mean() if len(idx) else np.nan
            low, high = np.percentile(draws, [2.5, 97.5])
            rows.append({
                "dropped": drop, "hypothesis": hypothesis, "n": int(len(sub)),
                "diff_sa_minus_base": float(left.mean() - right.mean()) if len(sub) else np.nan,
                "ci_low": float(low), "ci_high": float(high),
                "onesided_p05": float(np.percentile(draws, 5)) if hypothesis == "AH3" else np.nan,
                "favor": int((~left & right).sum()) if kind == "host" else int((left & ~right).sum()),
                "against": int((left & ~right).sum()) if kind == "host" else int((~left & right).sum()),
            })
    write(pd.DataFrame(rows), "leave_one_cohort.tsv")


def training_counts() -> None:
    split = pd.read_csv(ROOT / "config/split_tcga.tsv", sep="\t", dtype=str)
    train = split.loc[split["split"].eq("TCGA-train") & ~split["exclude_from_eval"].astype(str).str.lower().eq("true")]
    rows = []
    for label in ("ESCA", "STAD", "PAAD", "CHOL"):
        sub = train.loc[train["learning_label"].eq(label)]
        diseases = "|".join(f"{k}:{v}" for k, v in sub["primary disease or tissue"].value_counts().items())
        rows.append({"learning_label": label, "n_train": int(len(sub)), "disease_strings": diseases})
    write(pd.DataFrame(rows), "tcga_train_gi_counts.tsv")


def colon_and_models(frame: pd.DataFrame) -> None:
    part = frame.loc[
        frame["cohort"].eq("GSE41258") & frame["standard_site"].eq("colon_rectum") & frame["selected_eval"]
    ].copy()
    rows = []
    for method in ("BASE-Z", "SA-Z", "LD-Z"):
        pred = part[f"{method}__pred"].astype(str)
        for name, n in pred.value_counts().items():
            rows.append({"method": method, "predicted": name, "n": int(n), "n_slice": int(len(part))})
    write(pd.DataFrame(rows), "gse41258_colon_confusion.tsv")
    z = np.load(ROOT / "data/processed/aux/GSE41258/Z.npy")
    ids = (ROOT / "data/processed/aux/GSE41258/sample_ids.txt").read_text().splitlines()
    genes = [gene for gene in (ROOT / "results/B0_genes.txt").read_text().splitlines() if gene and gene != "symbol"]
    index = {sample: i for i, sample in enumerate(ids)}
    gene_index = {gene: i for i, gene in enumerate(genes)}
    use = [sample for sample in part["sample_id"].astype(str) if sample in index]
    block = z[[index[sample] for sample in use]]
    z_rows = []
    for gene in COLON:
        if gene not in gene_index:
            z_rows.append({"gene": gene, "in_B0": False})
            continue
        values = block[:, gene_index[gene]]
        z_rows.append({"gene": gene, "in_B0": True, "n": int(len(values)), "mean_Z": float(values.mean()), "median_Z": float(np.median(values))})
    write(pd.DataFrame(z_rows), "gse41258_colon_marker_Z.tsv")
    # Linear models: learning-label logits. Colorectal is the single class COADREAD.
    mis = part.loc[part["SA-Z__pred"].ne("Colorectal") | part["BASE-Z__pred"].ne("Colorectal"), ["sample_id", "BASE-Z__pred", "SA-Z__pred"]]
    contrib_rows = []
    for method, path in (("BASE-Z", "BASE_Z.joblib"), ("SA-Z", "SA_Z.joblib")):
        bundle = joblib.load(ROOT / "results/stage2/models" / path)
        classes = list(bundle["classes"])
        coef = bundle["clf"].coef_
        features = genes
        if coef.shape[1] != len(genes):
            raise SystemExit(f"{method} coef {coef.shape} genes {len(genes)}")
        order = list(range(len(genes)))
        wrong = part.loc[part[f"{method}__pred"].ne("Colorectal")]
        if wrong.empty or "COADREAD" not in classes:
            continue
        sample_rows = z[[index[sample] for sample in wrong["sample_id"].astype(str)]][:, order]
        crc = classes.index("COADREAD")
        # predicted learning label is the argmax project, recovered from the organ only when one project maps to it
        organ_to_projects = {}
        for project, organ in PROJECT_TO_ORGAN.items():
            organ_to_projects.setdefault(organ, []).append(project)
        chunks = []
        for sample_i, organ in enumerate(wrong[f"{method}__pred"].astype(str)):
            projects = [project for project in organ_to_projects.get(organ, []) if project in classes]
            if len(projects) != 1:
                continue
            pred_i = classes.index(projects[0])
            chunks.append((coef[pred_i] - coef[crc]) * sample_rows[sample_i])
        if not chunks:
            continue
        mean = np.mean(chunks, axis=0)
        top = np.argsort(-np.abs(mean))[:20]
        for rank, feature_i in enumerate(top, start=1):
            contrib_rows.append({
                "method": method, "n_misclassified": int(len(wrong)), "n_single_project_class": int(len(chunks)), "rank": rank,
                "feature": features[feature_i], "mean_contribution": float(mean[feature_i]),
                "model_classes": "32 learning labels", "colorectal_class": "COADREAD",
            })
    write(pd.DataFrame(contrib_rows), "gse41258_logit_contribution.tsv")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frame = confirm()
    imvigor_confusion(frame)
    assay_table()
    tcga_top1()
    su2c_pairs(frame)
    redistribution(frame)
    effect_bootstrap(frame)
    training_counts()
    colon_and_models(frame)
    print("diagnostics written", OUT, flush=True)


if __name__ == "__main__":
    main()
