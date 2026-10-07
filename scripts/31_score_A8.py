"""Score the pre-specified SCOPE and CUP-AI-Dx comparison. No hypothesis test."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.metrics import f1_score

sys.path.insert(0, str(Path("/home/kangjw/WORK/Cisplatin/2026_NewProj/proj_A_host_too/scripts")))
from mapping_rules import PROJECT_TO_ORGAN, SITE_NATIVE, map_site  # noqa: E402

ROOT = Path("/home/kangjw/WORK/Cisplatin/2026_NewProj/proj_A_host_too")
ORGANS = sorted(set(PROJECT_TO_ORGAN.values()))
OUT = ROOT / "results/stage8/external/tables"
NOTE = "사전 지정 서술, 가설 검정 아님"
GI = {"Stomach", "Esophagus", "Pancreas", "Biliary", "Colorectal"}
GROUPS = {
    "real_host_tissue": {"adrenal", "brain", "liver", "lung", "skin", "soft_tissue"},
    "proxy_host": {"bone", "lymph_node"},
    "pool_out": {"kidney", "pancreas", "bladder", "thyroid", "ovary", "breast", "stomach", "colon_rectum", "prostate", "esophagus", "head_neck", "cervix"},
    "other": {"pleura"},
}
RNA = [
    "blca_iatlas_imvigor210_2017", "brca_iatlas_anders_2022", "mel_dfci_2019",
    "paad_iatlas_prince_2022", "GSE50760", "prad_su2c_2019",
]
METHODS = ["BASE-Z", "SA-Z", "SCOPE", "CUP-AI-Dx"]


def native_organs(site: str) -> list[str]:
    return list(SITE_NATIVE.get(str(site), ([], [], []))[2])


def site_group(site: str) -> str:
    for name, sites in GROUPS.items():
        if site in sites:
            return name
    return "ungrouped"


def class_map(tool: str) -> dict[str, str]:
    frame = pd.read_csv(ROOT / "config/external_tool_label_map.tsv", sep="\t", dtype=str)
    part = frame.loc[(frame["row_type"] == "class") & (frame["tool"] == tool)]
    return dict(zip(part["source"], part["target"]))


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin(["true", "1"])


def top_organs(probabilities: pd.DataFrame, mapper: dict[str, str] | None) -> list[list[str]]:
    columns = list(probabilities.columns)
    values = probabilities.to_numpy(dtype=float)
    order = np.argsort(-values, axis=1)
    rows = []
    for ranks in order:
        seen = []
        for index in ranks:
            organ = columns[index] if mapper is None else mapper.get(columns[index], "")
            if organ and organ not in seen:
                seen.append(organ)
            if len(seen) == 3:
                break
        rows.append(seen)
    return rows


def host_flag(truth: str, pred: str, native: list[str]) -> bool:
    return bool(native) and truth not in native and pred in native


def gi_flag(truth: str, pred: str, native: list[str]) -> bool:
    return pred != truth and not host_flag(truth, pred, native) and truth in GI and pred in GI


def correct_primary(pred: str, truth: str, ns: bool) -> bool:
    return (not ns) and pred == truth and pred != ""


def metric_row(frame: pd.DataFrame, method: str, pred_col: str, top_col: str, ns_col: str | None) -> dict:
    scored = frame.loc[frame["in_top1"]] if "in_top1" in frame.columns else frame
    truth = scored["organ"].to_numpy(dtype=object)
    pred = scored[pred_col].fillna("").to_numpy(dtype=object)
    ns = scored[ns_col].to_numpy(dtype=bool) if ns_col else np.zeros(len(scored), dtype=bool)
    primary = np.array([correct_primary(p, t, flag) for p, t, flag in zip(pred, truth, ns)])
    sensitive = pred == truth
    top3 = np.array([t in organs for t, organs in zip(truth, scored[top_col])])
    risk = frame.loc[frame["at_risk"]]
    risk_pred = risk[pred_col].fillna("").to_numpy(dtype=object)
    host = np.array([host_flag(t, p, n) for t, p, n in zip(risk["organ"], risk_pred, risk["native"])])
    gi = np.array([gi_flag(t, p, n) for t, p, n in zip(truth, pred, scored["native"])])
    native_truth = scored["in_native"].to_numpy(dtype=bool)
    labels = sorted(set(truth.tolist()))
    scored_pred = np.where(ns, "NS_INCORRECT", pred)
    return {
        "n": int(len(scored)),
        "top1": float(primary.mean()) if len(scored) else np.nan,
        "top1_ns_as_organ": float(sensitive.mean()) if len(scored) else np.nan,
        "top3": float(top3.mean()) if len(scored) else np.nan,
        "macro_f1": float(f1_score(truth, scored_pred, labels=labels, average="macro", zero_division=0)) if len(scored) else np.nan,
        "n_at_risk": int(len(risk)),
        "host_rate": float(host.mean()) if len(risk) else np.nan,
        "n_native_truth": int(native_truth.sum()),
        "native_truth_top1": float(primary[native_truth].mean()) if native_truth.any() else np.nan,
        "gi_internal_error_count": int(gi.sum()),
        "n_unmapped_top1": int((pred == "").sum()),
        "n_scope_ns_top1": int(ns.sum()),
    }


def paired(frame: pd.DataFrame, left: str, right: str, kind: str) -> dict:
    if kind == "top1":
        frame = frame.loc[frame["in_top1"]].reset_index(drop=True) if "in_top1" in frame.columns else frame
        a = frame[f"{left}_correct"].to_numpy(dtype=bool)
        b = frame[f"{right}_correct"].to_numpy(dtype=bool)
    else:
        frame = frame.loc[frame["at_risk"]].reset_index(drop=True)
        a = frame[f"{left}_host"].to_numpy(dtype=bool)
        b = frame[f"{right}_host"].to_numpy(dtype=bool)
    diff = float(a.mean() - b.mean()) if len(a) else np.nan
    only_left = int((a & ~b).sum())
    only_right = int((~a & b).sum())
    discordant = only_left + only_right
    p_value = 1.0 if discordant == 0 else float(binomtest(only_left, discordant, 0.5, alternative="two-sided").pvalue)
    patients = frame["patient_id"].astype(str).to_numpy()
    unique = pd.unique(patients)
    rng = np.random.default_rng(20261001)
    draws = np.empty(2000, dtype=float)
    grouped = {patient: np.flatnonzero(patients == patient) for patient in unique}
    keys = list(grouped)
    for draw in range(2000):
        chosen = rng.choice(len(keys), size=len(keys), replace=True)
        rows = np.concatenate([grouped[keys[i]] for i in chosen])
        draws[draw] = a[rows].mean() - b[rows].mean()
    low, high = np.percentile(draws, [2.5, 97.5])
    return {
        "n": int(len(frame)), "diff": diff, "ci_low": float(low), "ci_high": float(high),
        "n_left_only": only_left, "n_right_only": only_right, "mcnemar_p": p_value,
    }


def distribution(frame: pd.DataFrame, method: str, column: str, level: str) -> pd.DataFrame:
    counts = frame[column].fillna("").astype(str).value_counts()
    out = counts.rename_axis("predicted").reset_index(name="n")
    out.insert(0, "method", method)
    out.insert(1, "level", level)
    return out


def load_tool(cohort: str, scope_map: dict[str, str], cup_map: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    paths = [ROOT / f"results/stage8/external/scope/{cohort}"]
    cup_paths = [ROOT / f"results/stage8/external/cup/{cohort}.tsv"]
    if cohort == "GSE50760":
        paths.append(ROOT / "results/stage8/external/scope/GSE50760_risk")
        cup_paths.append(ROOT / "results/stage8/external/cup/GSE50760_risk.tsv")
    tops, probs = [], []
    for path in paths:
        top = pd.read_csv(path / "SCOPE_topPredictions.txt", sep="\t", dtype=str)
        top = top.loc[top["rank_pred"] == "1"].copy()
        prob = pd.read_csv(path / "SCOPE_meanProbabilities.tsv", sep="\t", index_col=0)
        prob.index = prob.index.astype(str)
        tops.append(top)
        probs.append(prob)
    top = pd.concat(tops, ignore_index=True)
    probs = pd.concat(probs)
    if top["sample_name"].duplicated().any():
        raise SystemExit(f"duplicate SCOPE top1 {cohort}")
    top["organ"] = top["label"].map(scope_map).fillna("")
    top["ns"] = top["label"].str.endswith("_NS")
    top["top3"] = top_organs(probs.loc[top["sample_name"].tolist()], scope_map)
    scope = top.set_index("sample_name")[["label", "organ", "ns", "top3"]]
    cups = []
    for path in cup_paths:
        cup = pd.read_csv(path, sep="\t", index_col=0)
        cup.index = cup.index.astype(str)
        cups.append(cup)
    cup = pd.concat(cups)
    classes = [column for column in cup.columns if column != "pred"]
    cup["organ"] = cup["pred"].map(cup_map).fillna("")
    cup["top3"] = top_organs(cup[classes], cup_map)
    return scope, cup[["pred", "organ", "top3"]]


def attach(frame: pd.DataFrame, scope: pd.DataFrame, cup: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["sample_id"] = out["sample_id"].astype(str)
    missing_scope = [sample for sample in out["sample_id"] if sample not in scope.index]
    missing_cup = [sample for sample in out["sample_id"] if sample not in cup.index]
    if missing_scope or missing_cup:
        raise SystemExit(f"missing tool rows scope {len(missing_scope)} cup {len(missing_cup)}")
    out["SCOPE_class"] = scope.loc[out["sample_id"], "label"].to_numpy()
    out["SCOPE"] = scope.loc[out["sample_id"], "organ"].to_numpy()
    out["SCOPE_ns"] = scope.loc[out["sample_id"], "ns"].to_numpy(dtype=bool)
    out["SCOPE_top3"] = scope.loc[out["sample_id"], "top3"].tolist()
    out["CUP_class"] = cup.loc[out["sample_id"], "pred"].to_numpy()
    out["CUP-AI-Dx"] = cup.loc[out["sample_id"], "organ"].to_numpy()
    out["CUP-AI-Dx_top3"] = cup.loc[out["sample_id"], "top3"].tolist()
    return out


def organ_top3(proba: np.ndarray) -> list[list[str]]:
    order = np.argsort(-proba, axis=1)[:, :3]
    names = np.asarray(ORGANS, dtype=object)
    return [names[row].tolist() for row in order]


def met500() -> pd.DataFrame:
    samples = pd.read_csv(ROOT / "results/stage5/posthoc/met500_samples.tsv", sep="\t", dtype=str)
    if len(samples) != 437:
        raise SystemExit(f"MET500 {len(samples)}")
    base = np.load(ROOT / "results/stage5/posthoc/met500_BASE-Z_proba.npy")
    sa = np.load(ROOT / "results/stage5/posthoc/met500_SA-Z_proba.npy")
    if base.shape != (437, len(ORGANS)) or sa.shape != base.shape:
        raise SystemExit(f"MET500 proba {base.shape} {sa.shape}")
    native = [native_organs(site) for site in samples["site"]]
    frame = pd.DataFrame({
        "cohort": "MET500",
        "sample_id": samples["id"],
        "patient_id": samples["cluster"],
        "organ": samples["truth"],
        "standard_site": samples["site"],
        "native": native,
        "BASE-Z": samples["pred__BASE-Z"],
        "SA-Z": samples["pred__SA-Z"],
        "BASE-Z_top3": organ_top3(base),
        "SA-Z_top3": organ_top3(sa),
    })
    frame["at_risk"] = [bool(n) and truth not in n for truth, n in zip(frame["organ"], native)]
    frame["in_native"] = [truth in n for truth, n in zip(frame["organ"], native)]
    if int(frame["at_risk"].sum()) != 361:
        raise SystemExit(f"MET500 risk {int(frame['at_risk'].sum())}")
    return frame


def pog570() -> pd.DataFrame:
    labels = pd.read_csv(ROOT / "config/pog570_eval_labels.tsv", sep="\t", dtype=str)
    labels = labels.loc[~labels["organ"].isin(["exclude", "NA"])].copy()
    demo = pd.read_excel(ROOT / "data/raw/POG570/Table_S1_Demographics.xlsx", dtype=str)
    demo["standard_site"] = [map_site("pog570_biopsy_site", value) for value in demo["BIOPSY_SITE"]]
    labels = labels.merge(demo[["PATIENT_ID", "standard_site"]], on="PATIENT_ID", how="left")
    pred = pd.read_parquet(ROOT / "results/stage3/pog570_predictions.parquet")
    pred["patient_id"] = pred["patient_id"].astype(str)
    labels["PATIENT_ID"] = labels["PATIENT_ID"].astype(str)
    merged = labels.merge(pred, left_on="PATIENT_ID", right_on="patient_id", how="left")
    if merged["BASE-Z__pred"].isna().any() or len(merged) != 512:
        raise SystemExit(f"POG join {len(merged)}")
    native = [native_organs(site) for site in merged["standard_site"]]
    base_p = merged[[f"BASE-Z__p_{organ}" for organ in ORGANS]].to_numpy(dtype=float)
    sa_p = merged[[f"SA-Z__p_{organ}" for organ in ORGANS]].to_numpy(dtype=float)
    frame = pd.DataFrame({
        "cohort": "POG570",
        "sample_id": merged["PATIENT_ID"],
        "patient_id": merged["PATIENT_ID"],
        "organ": merged["organ"],
        "standard_site": merged["standard_site"],
        "native": native,
        "BASE-Z": merged["BASE-Z__pred"],
        "SA-Z": merged["SA-Z__pred"],
        "BASE-Z_top3": organ_top3(base_p),
        "SA-Z_top3": organ_top3(sa_p),
    })
    frame["at_risk"] = [bool(n) and truth not in n for truth, n in zip(frame["organ"], native)]
    frame["in_native"] = [truth in n for truth, n in zip(frame["organ"], native)]
    if int(frame["at_risk"].sum()) != 378 or int(frame["in_native"].sum()) != 91:
        raise SystemExit(f"POG sets risk {int(frame['at_risk'].sum())} native {int(frame['in_native'].sum())}")
    return frame


def auxiliary() -> pd.DataFrame:
    columns = ["cohort", "sample_id", "patient_id", "organ", "standard_site", "native_organs",
               "include_confirm", "selected_eval", "selected_risk", "selected_native",
               "BASE-Z__pred", "SA-Z__pred"]
    columns += [f"BASE-Z__p_{organ}" for organ in ORGANS]
    columns += [f"SA-Z__p_{organ}" for organ in ORGANS]
    pred = pd.read_parquet(ROOT / "results/stage7/confirm/predictions.parquet", columns=columns)
    for name in ("include_confirm", "selected_eval", "selected_risk", "selected_native"):
        pred[name] = as_bool(pred[name])
    part = pred.loc[pred["include_confirm"] & (pred["selected_eval"] | pred["selected_risk"]) & pred["cohort"].isin(RNA)].copy()
    if int((part["selected_eval"]).sum()) != 729 or int(part["selected_risk"].sum()) != 427 or int(part["selected_native"].sum()) != 315:
        raise SystemExit(f"aux sets {int(part['selected_eval'].sum())} {int(part['selected_risk'].sum())} {int(part['selected_native'].sum())}")
    native = [str(value).split("|") if str(value) not in {"", "nan"} else [] for value in part["native_organs"]]
    base_p = part[[f"BASE-Z__p_{organ}" for organ in ORGANS]].to_numpy(dtype=float)
    sa_p = part[[f"SA-Z__p_{organ}" for organ in ORGANS]].to_numpy(dtype=float)
    frame = pd.DataFrame({
        "cohort": part["cohort"].to_numpy(),
        "analysis_cohort": "aux_rnaseq",
        "sample_id": part["sample_id"].astype(str),
        "patient_id": part["patient_id"].astype(str),
        "organ": part["organ"],
        "standard_site": part["standard_site"].fillna(""),
        "native": native,
        "BASE-Z": part["BASE-Z__pred"],
        "SA-Z": part["SA-Z__pred"],
        "BASE-Z_top3": organ_top3(base_p),
        "SA-Z_top3": organ_top3(sa_p),
        "at_risk": part["selected_risk"].to_numpy(),
        "in_native": part["selected_native"].to_numpy(),
        "in_top1": part["selected_eval"].to_numpy(),
    })
    if int(frame["at_risk"].sum()) != 427 or int(frame["in_native"].sum()) != 315 or int(frame["in_top1"].sum()) != 729:
        raise SystemExit(f"aux sets {int(frame['at_risk'].sum())} {int(frame['in_native'].sum())}")
    return frame


def score_block(name: str, frame: pd.DataFrame, scope_map, cup_map, rows: dict[str, list]) -> pd.DataFrame:
    pieces = []
    if name == "aux_rnaseq":
        for cohort, part in frame.groupby("cohort", sort=False):
            scope, cup = load_tool(cohort, scope_map, cup_map)
            pieces.append(attach(part, scope, cup))
        use = pd.concat(pieces, ignore_index=True)
    else:
        scope, cup = load_tool(name, scope_map, cup_map)
        use = attach(frame, scope, cup)
    use["analysis_cohort"] = name
    use["BASE-Z_ns"] = False
    use["SA-Z_ns"] = False
    use["CUP-AI-Dx_ns"] = False
    for method, pred_col, top_col, ns_col in (
        ("BASE-Z", "BASE-Z", "BASE-Z_top3", "BASE-Z_ns"),
        ("SA-Z", "SA-Z", "SA-Z_top3", "SA-Z_ns"),
        ("SCOPE", "SCOPE", "SCOPE_top3", "SCOPE_ns"),
        ("CUP-AI-Dx", "CUP-AI-Dx", "CUP-AI-Dx_top3", "CUP-AI-Dx_ns"),
    ):
        use[f"{method}_correct"] = [correct_primary(p, t, flag) for p, t, flag in zip(use[pred_col], use["organ"], use[ns_col])]
        use[f"{method}_host"] = [host_flag(t, p, n) for t, p, n in zip(use["organ"], use[pred_col], use["native"])]
    scope_organs = set(scope_map.values())
    cup_organs = set(cup_map.values())
    common_organs = scope_organs & cup_organs & set(ORGANS)
    use["common_label"] = use["organ"].isin(common_organs)
    for subset, mask in (("full", np.ones(len(use), dtype=bool)), ("common_label", use["common_label"].to_numpy())):
        part = use.loc[mask].reset_index(drop=True)
        scored_rows = part.loc[part["in_top1"]] if "in_top1" in part.columns else part
        for method, pred_col, top_col, ns_col in (
            ("BASE-Z", "BASE-Z", "BASE-Z_top3", "BASE-Z_ns"),
            ("SA-Z", "SA-Z", "SA-Z_top3", "SA-Z_ns"),
            ("SCOPE", "SCOPE", "SCOPE_top3", "SCOPE_ns"),
            ("CUP-AI-Dx", "CUP-AI-Dx", "CUP-AI-Dx_top3", "CUP-AI-Dx_ns"),
        ):
            row = metric_row(part, method, pred_col, top_col, ns_col)
            row.update({"analysis_cohort": name, "subset": subset, "method": method, "note": NOTE})
            rows["metrics"].append(row)
            rows["distribution"].append(distribution(scored_rows, method, pred_col, "organ").assign(analysis_cohort=name, subset=subset, note=NOTE))
            raw_col = "SCOPE_class" if method == "SCOPE" else "CUP_class" if method == "CUP-AI-Dx" else pred_col
            rows["distribution"].append(distribution(scored_rows, method, raw_col, "raw_class").assign(analysis_cohort=name, subset=subset, note=NOTE))
        for left, right in (("SA-Z", "SCOPE"), ("SA-Z", "CUP-AI-Dx"), ("BASE-Z", "SCOPE"), ("BASE-Z", "CUP-AI-Dx")):
            for kind in ("top1", "host_rate"):
                stats = paired(part, left, right, kind)
                stats.update({
                    "analysis_cohort": name, "subset": subset, "contrast": f"{left} minus {right}",
                    "metric": kind, "note": NOTE,
                })
                rows["paired"].append(stats)
        risk = part.loc[part["at_risk"]].copy()
        risk["site_group"] = [site_group(site) for site in risk["standard_site"]]
        for group, sub in risk.groupby("site_group", sort=False):
            for method, pred_col in (("BASE-Z", "BASE-Z"), ("SA-Z", "SA-Z"), ("SCOPE", "SCOPE"), ("CUP-AI-Dx", "CUP-AI-Dx")):
                host = [host_flag(t, p, n) for t, p, n in zip(sub["organ"], sub[pred_col], sub["native"])]
                rows["site"].append({
                    "analysis_cohort": name, "subset": subset, "site_group": group, "method": method,
                    "n_at_risk": int(len(sub)), "host_rate": float(np.mean(host)) if len(sub) else np.nan, "note": NOTE,
                })
    return use


def overlaps(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    meta = pd.read_csv(
        ROOT / "results/stage7/external_tools/CUP-AI-Dx/data/ExternalDataMeta.csv",
        usecols=[0], index_col=0,
    )
    external = set(meta.index.astype(str))
    rows = []
    for name, frame in frames.items():
        ids = set(frame["sample_id"].astype(str)) | set(frame["patient_id"].astype(str))
        rows.append({"analysis_cohort": name, "source": "CUP-AI-Dx ExternalDataMeta", "n_overlap": int(len(ids & external)), "note": NOTE})
    return pd.DataFrame(rows)


def write(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, float_format="%.16g")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    scope_map = class_map("SCOPE")
    cup_map = class_map("CUP-AI-Dx")
    frames = {"MET500": met500(), "POG570": pog570(), "aux_rnaseq": auxiliary()}
    rows = {"metrics": [], "paired": [], "site": [], "distribution": []}
    scored = []
    for name, frame in frames.items():
        scored.append(score_block(name, frame, scope_map, cup_map, rows))
        print(name, "scored", flush=True)
    write(pd.DataFrame(rows["metrics"]), OUT / "metrics.tsv")
    write(pd.DataFrame(rows["paired"]), OUT / "paired.tsv")
    write(pd.DataFrame(rows["site"]), OUT / "host_by_site_group.tsv")
    write(pd.concat(rows["distribution"], ignore_index=True), OUT / "class_distribution.tsv")
    write(overlaps(frames), OUT / "id_overlap.tsv")
    keep = ["analysis_cohort", "cohort", "sample_id", "patient_id", "organ", "standard_site", "at_risk", "in_native", "common_label",
            "BASE-Z", "SA-Z", "SCOPE", "SCOPE_class", "SCOPE_ns", "CUP-AI-Dx", "CUP_class"]
    sample_table = pd.concat(scored, ignore_index=True)
    sample_table["cohort"] = sample_table["cohort"].astype(str)
    write(sample_table[keep], OUT / "sample_predictions.tsv")
    print("tables", OUT, flush=True)


if __name__ == "__main__":
    main()
