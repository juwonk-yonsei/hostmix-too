"""Post-hoc split of the stage-8 classifier comparison. Does not rescore."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path("/home/kangjw/WORK/Cisplatin/2026_NewProj/proj_A_host_too/scripts")))
from mapping_rules import SITE_NATIVE  # noqa: E402

ROOT = Path("/home/kangjw/WORK/Cisplatin/2026_NewProj/proj_A_host_too")
OUT = ROOT / "results/stage9"
NOTE = "사후"
METHODS = ["BASE-Z", "SA-Z", "SCOPE", "CUP-AI-Dx"]
RNA = [
    "blca_iatlas_imvigor210_2017",
    "brca_iatlas_anders_2022",
    "mel_dfci_2019",
    "paad_iatlas_prince_2022",
    "GSE50760",
    "prad_su2c_2019",
]
POOL_OUT = {
    "kidney", "pancreas", "bladder", "thyroid", "ovary", "breast", "stomach",
    "colon_rectum", "prostate", "esophagus", "head_neck", "cervix",
}


def native_list(site: str) -> list[str]:
    return list(SITE_NATIVE.get(str(site), ([], [], []))[2])


def rates(frame: pd.DataFrame, method: str) -> dict:
    truth = frame["organ"].astype(str)
    pred = frame[method].astype(str)
    correct = pred.eq(truth)
    host = []
    for organ, call, site in zip(truth, pred, frame["standard_site"].astype(str)):
        native = native_list(site)
        host.append(bool(native) and call in native and call != organ)
    risk = frame["at_risk"].astype(str).str.lower().isin(["true", "1"])
    return {
        "n_top1": int(len(frame)),
        "top1": float(correct.mean()) if len(frame) else float("nan"),
        "n_at_risk": int(risk.sum()),
        "host_rate": float(pd.Series(host)[risk.to_numpy()].mean()) if risk.any() else float("nan"),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pred = pd.read_parquet(
        ROOT / "results/stage7/confirm/predictions.parquet",
        columns=["cohort", "sample_id", "include_confirm", "selected_eval", "selected_risk"],
    )
    for column in ("include_confirm", "selected_eval", "selected_risk"):
        pred[column] = pred[column].astype(str).str.lower().isin(["true", "1"])
    eval_ids = set(pred.loc[pred["include_confirm"] & pred["selected_eval"], "sample_id"].astype(str))
    scored = pd.read_csv(ROOT / "results/stage8/external/tables/sample_predictions.tsv", sep="\t", dtype=str)
    scored["at_risk"] = scored["at_risk"].astype(str).str.lower().isin(["true", "1"])
    scored["in_native"] = scored["in_native"].astype(str).str.lower().isin(["true", "1"])
    rows = []
    for cohort in RNA:
        part = scored.loc[scored["cohort"].eq(cohort) & scored["sample_id"].isin(eval_ids)].copy()
        risk = scored.loc[scored["cohort"].eq(cohort) & scored["at_risk"]].copy()
        for method in METHODS:
            top = rates(part.assign(at_risk=False), method)
            host = rates(risk, method)
            rows.append({
                "note": NOTE, "analysis_cohort": "aux_rnaseq", "cohort": cohort, "method": method,
                "n": top["n_top1"], "top1": top["top1"], "n_at_risk": host["n_at_risk"], "host_rate": host["host_rate"],
            })
    by_cohort = pd.DataFrame(rows)
    if not (by_cohort.groupby("method")["n"].sum() == 729).all():
        raise SystemExit(by_cohort.groupby("method")["n"].sum().to_dict())
    if not (by_cohort.groupby("method")["n_at_risk"].sum() == 427).all():
        raise SystemExit(by_cohort.groupby("method")["n_at_risk"].sum().to_dict())
    by_cohort.to_csv(OUT / "subcohort_metrics.tsv", sep="\t", index=False)

    set_rows = []
    for name, frame in (
        ("MET500", scored.loc[scored["analysis_cohort"].eq("MET500")]),
        ("POG570", scored.loc[scored["analysis_cohort"].eq("POG570")]),
        ("aux_rnaseq", scored.loc[scored["analysis_cohort"].eq("aux_rnaseq") & scored["sample_id"].isin(eval_ids)]),
    ):
        slices = {
            "native_truth": frame.loc[frame["in_native"]],
            "risk": frame.loc[frame["at_risk"]],
            "pool_out_risk": frame.loc[frame["at_risk"] & frame["standard_site"].isin(POOL_OUT)],
        }
        if name == "aux_rnaseq":
            slices["risk"] = scored.loc[scored["analysis_cohort"].eq("aux_rnaseq") & scored["at_risk"]]
            slices["pool_out_risk"] = slices["risk"].loc[slices["risk"]["standard_site"].isin(POOL_OUT)]
        for set_name, part in slices.items():
            for method in METHODS:
                truth = part["organ"].astype(str)
                call = part[method].astype(str)
                host = []
                for organ, pred_organ, site in zip(truth, call, part["standard_site"].astype(str)):
                    native = native_list(site)
                    host.append(bool(native) and pred_organ in native and pred_organ != organ)
                set_rows.append({
                    "note": NOTE, "analysis_cohort": name, "set": set_name, "method": method,
                    "n": int(len(part)), "top1": float(call.eq(truth).mean()) if len(part) else float("nan"),
                    "host_rate": float(sum(host) / len(host)) if host else float("nan"),
                })
    sets = pd.DataFrame(set_rows)
    sets.to_csv(OUT / "set_metrics.tsv", sep="\t", index=False)

    aux_eval = scored.loc[scored["analysis_cohort"].eq("aux_rnaseq") & scored["sample_id"].isin(eval_ids)].copy()
    cup_only = aux_eval.loc[aux_eval["CUP-AI-Dx"].eq(aux_eval["organ"]) & aux_eval["SA-Z"].ne(aux_eval["organ"])]
    pieces = []
    for column in ("cohort", "standard_site", "organ", "SA-Z"):
        counts = cup_only[column].astype(str).value_counts()
        pieces.append(pd.DataFrame({
            "note": NOTE, "field": column, "level": counts.index.astype(str), "n": counts.to_numpy(),
            "n_cup_correct_sa_wrong": int(len(cup_only)), "n_eval": int(len(aux_eval)),
        }))
    pd.concat(pieces, ignore_index=True).to_csv(OUT / "cup_correct_sa_wrong.tsv", sep="\t", index=False)
    print("subcohort", len(by_cohort), "cup_only", len(cup_only), "eval", len(aux_eval), flush=True)
    print(sets.loc[sets["method"].eq("SA-Z"), ["analysis_cohort", "set", "n", "top1", "host_rate"]].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
