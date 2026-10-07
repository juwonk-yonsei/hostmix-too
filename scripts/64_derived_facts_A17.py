"""Post hoc recounts of stored predictions for numbers that have no stored aggregate (A17 §4.7).

Outputs go to results/stage10/derived/ and are post hoc. No model is refit and nothing is rescored:
the counts are read from results/stage7/confirm/predictions.parquet.

microarray_rates.tsv  one sample per patient across the microarray layer (polyA library first, then
                      sample id, then cohort), evaluation and at-risk sets, per cohort and for the layer
aux_six_rates.tsv     the six confirmatory RNA-seq cohorts pooled, evaluation and at-risk sets as stored
gate_decisions.tsv    the six confirmatory RNA-seq cohorts pooled, evaluation and native-truth sets: samples
                      whose recorded LD beta sends the gated model to the baseline (beta <= 0.02 or missing,
                      the rule sa_g_rule of config/analysis_plan_A6.yaml)
screening_decisions.tsv
                      cohort search log, decision counts using the last row of each study identifier;
                      read from results/stage5/aux/screening_log.tsv
screening_criteria.tsv
                      excluded studies by failed criterion, same last-row rule
patient_overlap.tsv   patients of the evaluation and at-risk sets present in both GSE74685 and prad_fhcrc,
                      read from config/aux_eval_labels.tsv
duplicate_removed.tsv samples flagged by the duplicate screen (flag_candidate), unique sample ids per
                      cohort with the expression libraries of one cohort taken together; read from
                      results/stage5/aux/cohorts/*/samples.tsv
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "results" / "stage7" / "confirm" / "predictions.parquet"
OUT = ROOT / "results" / "stage10" / "derived"
COHORTS = ROOT / "results" / "stage5" / "aux" / "cohorts"
METHODS = ["BASE-Z", "SA-Z", "SA-G", "SA-pool22", "SC-Z", "LD-Z", "NC-Z", "M1-Z", "BASE-K", "SA-K", "V0-K", "SA-MLP"]
MICROARRAY = "\ub9c8\uc774\ud06c\ub85c\uc5b4\ub808\uc774"
AUX = ["blca_iatlas_imvigor210_2017", "brca_iatlas_anders_2022", "mel_dfci_2019", "paad_iatlas_prince_2022",
       "GSE50760", "prad_su2c_2019"]


def flag(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin(["true", "1"])


def one_per_patient(frame: pd.DataFrame) -> pd.DataFrame:
    if frame["patient_id"].is_unique:
        return frame
    rank = np.where(frame["library"].fillna("").eq("polyA"), 0, 1)
    ordered = frame.assign(_rank=rank).sort_values(["_rank", "sample_id", "cohort"], kind="mergesort")
    return ordered.drop_duplicates("patient_id", keep="first").drop(columns="_rank")


def attracted(frame: pd.DataFrame, method: str) -> np.ndarray:
    natives = [set(filter(None, str(text).split("|"))) for text in frame["native_organs"]]
    pred = frame[f"{method}__pred"].to_numpy(dtype=object)
    truth = frame["organ"].to_numpy(dtype=object)
    return np.array([bool(natives[i]) and truth[i] not in natives[i] and pred[i] in natives[i]
                     for i in range(len(frame))], dtype=bool)


def rates(evaluation: pd.DataFrame, risk: pd.DataFrame, label: str) -> list[dict]:
    rows = []
    for method in METHODS:
        correct = evaluation[f"{method}__pred"].to_numpy(dtype=object) == evaluation["organ"].to_numpy(dtype=object)
        hits = attracted(risk, method)
        rows.append({"cohort": label, "method": method, "n_evaluation": int(len(evaluation)),
                     "n_correct": int(correct.sum()), "top1": float(correct.mean()) if len(evaluation) else np.nan,
                     "n_at_risk": int(len(risk)), "n_attracted": int(hits.sum()),
                     "host_rate": float(hits.mean()) if len(risk) else np.nan})
    return rows


def main() -> None:
    columns = ["cohort", "patient_id", "sample_id", "library", "organ", "native_organs", "layer", "layer_test",
               "selected_eval", "selected_risk", "selected_native", "beta", *[f"{m}__pred" for m in METHODS]]
    pred = pd.read_parquet(PRED, columns=columns)
    OUT.mkdir(parents=True, exist_ok=True)

    array = pred.loc[pred["layer"].eq(MICROARRAY) & flag(pred["layer_test"])]
    evaluation = one_per_patient(array.loc[flag(array["selected_eval"])])
    risk = one_per_patient(array.loc[flag(array["selected_risk"])])
    rows = rates(evaluation, risk, "all")
    for cohort in sorted(set(array["cohort"])):
        rows += rates(evaluation.loc[evaluation["cohort"].eq(cohort)], risk.loc[risk["cohort"].eq(cohort)], cohort)
    micro = pd.DataFrame(rows)
    micro.to_csv(OUT / "microarray_rates.tsv", sep="\t", index=False)

    part = pred.loc[pred["cohort"].isin(AUX)]
    aux = pd.DataFrame(rates(part.loc[flag(part["selected_eval"])], part.loc[flag(part["selected_risk"])], "aux_six"))
    aux.to_csv(OUT / "aux_six_rates.tsv", sep="\t", index=False)

    tau = float(yaml.safe_load((ROOT / "config" / "analysis_plan_A6.yaml").read_text())["sa_g_tau"])
    gate = []
    for set_name, column in (("evaluation", "selected_eval"), ("native truth", "selected_native")):
        chosen = part.loc[flag(part[column])]
        to_base = chosen["beta"].isna() | (chosen["beta"] <= tau)
        gate.append({"cohort": "aux_six", "set": set_name, "n": int(len(chosen)), "n_baseline": int(to_base.sum()),
                     "n_hostmix": int((~to_base).sum()), "tau": tau})
    gate = pd.DataFrame(gate)
    gate.to_csv(OUT / "gate_decisions.tsv", sep="\t", index=False)

    log = pd.read_csv(ROOT / "results" / "stage5" / "aux" / "screening_log.tsv", sep="\t", dtype=str)
    last = log.groupby("id", sort=False).tail(1)
    screening = (last.groupby(["path", "decision"]).size().rename("n_studies").reset_index())
    screening.to_csv(OUT / "screening_decisions.tsv", sep="\t", index=False)
    criteria = (last.loc[last["decision"].eq("\uc81c\uc678")].groupby("failed_criterion").size()
                .rename("n_studies").reset_index())
    criteria.to_csv(OUT / "screening_criteria.tsv", sep="\t", index=False)

    labels = pd.read_csv(ROOT / "config" / "aux_eval_labels.tsv", sep="\t", dtype=str)
    overlap = []
    for set_name, column in (("evaluation", "in_eval"), ("at risk", "in_risk")):
        ids = {c: set(labels.loc[labels["cohort"].eq(c) & labels[column].eq("True"), "patient_id"])
               for c in ("GSE74685", "prad_fhcrc")}
        overlap.append({"set": set_name, "n_GSE74685": len(ids["GSE74685"]), "n_prad_fhcrc": len(ids["prad_fhcrc"]),
                        "n_both": len(ids["GSE74685"] & ids["prad_fhcrc"])})
    overlap = pd.DataFrame(overlap)
    overlap.to_csv(OUT / "patient_overlap.tsv", sep="\t", index=False)

    flagged: dict = {}
    for path in sorted(COHORTS.glob("*/samples.tsv")):
        cohort = path.parent.name.split("__")[0]
        frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        ids = flagged.setdefault(cohort, {"libraries": 0, "screened": False, "ids": set()})
        ids["libraries"] += 1
        if "flag_candidate" in frame:
            ids["screened"] = True
            ids["ids"] |= set(frame.loc[frame["flag_candidate"].eq("True"), "sample_id"])
    dup = pd.DataFrame([{"cohort": c, "libraries": v["libraries"], "screened": v["screened"],
                         "n_removed": len(v["ids"]) if v["screened"] else None} for c, v in flagged.items()])
    dup["n_removed"] = dup["n_removed"].astype("Int64")
    dup.to_csv(OUT / "duplicate_removed.tsv", sep="\t", index=False)

    meta = {"status": "post hoc", "source": "results/stage7/confirm/predictions.parquet; "
            "results/stage5/aux/cohorts/*/samples.tsv",
            "script": "scripts/64_derived_facts_A17.py",
            "files": {"microarray_rates.tsv": len(micro), "aux_six_rates.tsv": len(aux),
                      "duplicate_removed.tsv": len(dup), "gate_decisions.tsv": len(gate),
                      "screening_decisions.tsv": len(screening), "screening_criteria.tsv": len(criteria),
                      "patient_overlap.tsv": len(overlap)},
            "microarray_layer_n": {"evaluation": int(len(evaluation)), "at_risk": int(len(risk))}}
    (OUT / "README.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
