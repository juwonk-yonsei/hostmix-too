#!/usr/bin/env python3
"""Independent check of confirmation sets and statistics.

Reads labels and predictions only. Does not import the confirmation script
or its project modules. pandas, numpy, and scipy only.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
BOOL_COLS = (
    "in_eval", "in_risk", "in_native", "in_pool_out_risk", "excluded", "include_confirm",
    "layer_test", "selected_eval", "selected_risk", "selected_native", "selected_pool_out", "pool_out_site",
)
METHODS = (
    "BASE-Z", "SA-Z", "SA-G", "SA-pool22", "SC-Z", "LD-Z", "NC-Z", "M1-Z",
    "BASE-K", "SA-K", "V0-K", "SA-MLP",
)


def read_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    for column in BOOL_COLS:
        if column in frame.columns:
            frame[column] = frame[column].eq("True")
    if "library" in frame.columns:
        frame["library"] = frame["library"].fillna("")
    return frame


def natives_of(text: str) -> set[str]:
    return set(filter(None, str(text).split("|")))


def as_bool(frame: pd.DataFrame, column: str) -> pd.Series:
    if frame[column].dtype == bool:
        return frame[column]
    return frame[column].astype(str).eq("True")


def dedupe_patients(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty or frame["patient_id"].is_unique:
        return frame.reset_index(drop=True)
    part = frame.copy()
    library = part["library"].fillna("")
    part["_library_rank"] = np.where(library.eq("polyA"), 0, 1)
    part = part.sort_values(["_library_rank", "sample_id", "cohort"], kind="mergesort")
    part = part.drop_duplicates("patient_id", keep="first").drop(columns="_library_rank")
    return part.reset_index(drop=True)


def one_sided(n_favor: int, n_against: int) -> float:
    m = n_favor + n_against
    if m == 0:
        return 1.0
    return float(binom.sf(n_favor - 1, m, 0.5))


def holm(pvalues: dict[str, float], alpha: float = 0.05) -> dict[str, bool]:
    ordered = sorted(pvalues, key=lambda name: (pvalues[name], name))
    rejected = {name: False for name in pvalues}
    for i, name in enumerate(ordered):
        if pvalues[name] > alpha / (len(ordered) - i):
            break
        rejected[name] = True
    return rejected


def host_error(frame: pd.DataFrame, method: str) -> np.ndarray:
    pred = frame[f"{method}__pred"].to_numpy(dtype=object)
    truth = frame["organ"].to_numpy(dtype=object)
    native = [natives_of(text) for text in frame["native_organs"]]
    return np.array([
        bool(native[i]) and truth[i] not in native[i] and pred[i] in native[i]
        for i in range(len(frame))
    ])


def correct(frame: pd.DataFrame, method: str) -> np.ndarray:
    return frame[f"{method}__pred"].to_numpy(dtype=object) == frame["organ"].to_numpy(dtype=object)


def hypothesis(frame: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    risk = dedupe_patients(frame.loc[as_bool(frame, "selected_risk")])
    eval_df = dedupe_patients(frame.loc[as_bool(frame, "selected_eval")])
    native = dedupe_patients(frame.loc[as_bool(frame, "selected_native")])
    pool = dedupe_patients(frame.loc[as_bool(frame, "selected_pool_out")])

    def pack(name, part, favor, against, left, right):
        n_favor = int(favor.sum()) if len(part) else 0
        n_against = int(against.sum()) if len(part) else 0
        diff = None if len(part) == 0 else float(left.mean() - right.mean())
        return {
            "hypothesis": name, "n": int(len(part)), "n_favor": n_favor, "n_against": n_against,
            "p": None if len(part) == 0 else one_sided(n_favor, n_against),
            "diff": diff, "status": "검정 불가" if len(part) == 0 else "계산",
        }

    primary_raw = {
        "AH1": pack(
            "AH1", risk,
            host_error(risk, "BASE-Z") & ~host_error(risk, "SA-Z"),
            host_error(risk, "SA-Z") & ~host_error(risk, "BASE-Z"),
            host_error(risk, "SA-Z").astype(float), host_error(risk, "BASE-Z").astype(float),
        ),
        "AH2": pack(
            "AH2", eval_df,
            correct(eval_df, "SA-Z") & ~correct(eval_df, "BASE-Z"),
            correct(eval_df, "BASE-Z") & ~correct(eval_df, "SA-Z"),
            correct(eval_df, "SA-Z").astype(float), correct(eval_df, "BASE-Z").astype(float),
        ),
        "AH3": pack(
            "AH3", native,
            correct(native, "SA-Z") & ~correct(native, "BASE-Z"),
            correct(native, "BASE-Z") & ~correct(native, "SA-Z"),
            correct(native, "SA-Z").astype(float), correct(native, "BASE-Z").astype(float),
        ),
    }
    ah1_sig = primary_raw["AH1"]["status"] != "검정 불가" and primary_raw["AH1"]["p"] <= 0.05
    ah2_sig = ah1_sig and primary_raw["AH2"]["status"] != "검정 불가" and primary_raw["AH2"]["p"] <= 0.05
    primary_raw["AH1"]["tested"] = primary_raw["AH1"]["status"] != "검정 불가"
    primary_raw["AH1"]["reject"] = bool(ah1_sig)
    primary_raw["AH2"]["tested"] = bool(ah1_sig and primary_raw["AH2"]["status"] != "검정 불가")
    primary_raw["AH2"]["reject"] = bool(ah2_sig)
    primary_raw["AH3"]["tested"] = bool(ah2_sig and primary_raw["AH3"]["status"] != "검정 불가")
    secondary_raw = {
        "AS1": pack(
            "AS1", native,
            correct(native, "SA-G") & ~correct(native, "SA-Z"),
            correct(native, "SA-Z") & ~correct(native, "SA-G"),
            correct(native, "SA-G").astype(float), correct(native, "SA-Z").astype(float),
        ),
        "AS2": pack(
            "AS2", pool,
            host_error(pool, "SA-Z") & ~host_error(pool, "SA-pool22"),
            host_error(pool, "SA-pool22") & ~host_error(pool, "SA-Z"),
            host_error(pool, "SA-pool22").astype(float), host_error(pool, "SA-Z").astype(float),
        ),
        "AS3": pack(
            "AS3", eval_df,
            correct(eval_df, "SA-G") & ~correct(eval_df, "SA-Z"),
            correct(eval_df, "SA-Z") & ~correct(eval_df, "SA-G"),
            correct(eval_df, "SA-G").astype(float), correct(eval_df, "SA-Z").astype(float),
        ),
    }
    usable = {name: row["p"] for name, row in secondary_raw.items() if row["status"] != "검정 불가"}
    rejected = holm(usable) if usable else {}
    for name, row in secondary_raw.items():
        row["tested"] = row["status"] != "검정 불가"
        row["reject"] = bool(rejected.get(name, False))
    return [primary_raw[name] for name in ("AH1", "AH2", "AH3")], [secondary_raw[name] for name in ("AS1", "AS2", "AS3")]


def layer_metrics(frame: pd.DataFrame) -> list[dict]:
    rows = []
    use = frame.loc[as_bool(frame, "selected_eval")]
    for method in METHODS:
        if f"{method}__pred" not in frame.columns:
            continue
        pred = use[f"{method}__pred"].to_numpy(dtype=object)
        truth = use["organ"].to_numpy(dtype=object)
        native = [natives_of(text) for text in use["native_organs"]]
        risk = np.array([bool(native[i]) and truth[i] not in native[i] for i in range(len(use))])
        native_truth = np.array([bool(native[i]) and truth[i] in native[i] for i in range(len(use))])
        host = risk & host_error(use, method) if len(use) else risk
        rows.append({
            "method": method,
            "n": int(len(use)),
            "top1": float(np.mean(pred == truth)) if len(use) else None,
            "n_at_risk": int(risk.sum()),
            "host_rate": float(host.sum() / risk.sum()) if risk.any() else None,
            "n_native_truth": int(native_truth.sum()),
            "native_truth_top1": float(np.mean(pred[native_truth] == truth[native_truth])) if native_truth.any() else None,
        })
    return rows


def conformal_n(frame: pd.DataFrame) -> list[dict]:
    use = frame.loc[as_bool(frame, "in_eval") & ~as_bool(frame, "excluded")]
    rows = []
    for layer, part in use.groupby("layer", sort=False):
        level = np.array([
            "liver" if site == "liver" else "lymph_node" if site == "lymph_node" else "other"
            for site in part["standard_site"]
        ])
        rows.append({"layer": layer, "level": "all", "n": int(len(part))})
        for name in ("liver", "lymph_node", "other"):
            rows.append({"layer": layer, "level": name, "n": int((level == name).sum())})
    return rows


def meta_records(frame: pd.DataFrame) -> list[dict]:
    rows = []
    scoped = frame.loc[as_bool(frame, "layer_test")]
    for layer, layer_df in scoped.groupby("layer", sort=False):
        records = []
        for cohort, part in layer_df.groupby("cohort", sort=False):
            risk = part.loc[as_bool(part, "selected_risk")]
            n = len(risk)
            if n < 10:
                continue
            b = int((host_error(risk, "BASE-Z") & ~host_error(risk, "SA-Z")).sum())
            c = int((host_error(risk, "SA-Z") & ~host_error(risk, "BASE-Z")).sum())
            d = (c - b) / n
            b_v, c_v = b, c
            if b == 0 or c == 0:
                b_v += 0.5
                c_v += 0.5
            v = ((b_v + c_v) - (b_v - c_v) ** 2 / n) / n ** 2
            records.append({"layer": layer, "cohort": cohort, "n": n, "b": b, "c": c, "d": d, "v": v})
        if len(records) >= 2:
            dvec = np.array([row["d"] for row in records], dtype=float)
            vvec = np.array([row["v"] for row in records], dtype=float)
            w = 1.0 / vvec
            d_fe = float(np.sum(w * dvec) / np.sum(w))
            q = float(np.sum(w * (dvec - d_fe) ** 2))
            df = len(records) - 1
            csum = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
            tau2 = max(0.0, (q - df) / csum) if csum > 0 else 0.0
            w_re = 1.0 / (vvec + tau2)
            estimate = float(np.sum(w_re * dvec) / np.sum(w_re))
            rows.append({"layer": layer, "kind": "dl", "k": len(records), "estimate": estimate})
        for record in records:
            record["kind"] = "cohort"
            rows.append(record)
    return rows


def sensitivity_sizes(frame: pd.DataFrame) -> list[dict]:
    confirm = frame.loc[as_bool(frame, "include_confirm") & frame["layer"].eq("RNA-seq")]
    all_samples = confirm.loc[as_bool(confirm, "in_eval") & ~as_bool(confirm, "excluded")]
    base = host_error(all_samples, "BASE-Z") if len(all_samples) else np.array([], dtype=bool)
    sa = host_error(all_samples, "SA-Z") if len(all_samples) else np.array([], dtype=bool)
    rows = [{
        "analysis": "all_samples_cluster_AH1", "n": int(len(all_samples)),
        "n_favor": int((base & ~sa).sum()) if len(all_samples) else 0,
        "n_against": int((sa & ~base).sum()) if len(all_samples) else 0,
    }]
    capture = confirm.loc[confirm["library"].ne("capture")]
    rows.append({"analysis": "drop_SU2C_capture_AH1", "n": int(len(dedupe_patients(capture.loc[as_bool(capture, "selected_risk")])))})
    if "GSE209998" in set(confirm["cohort"]):
        kept = confirm.loc[confirm["cohort"].ne("GSE209998")]
        rows.append({"analysis": "drop_GSE209998_AH1", "n": int(len(dedupe_patients(kept.loc[as_bool(kept, "selected_risk")])))})
    else:
        rows.append({"analysis": "drop_GSE209998_AH1", "n": None})
    for cohort in sorted(confirm["cohort"].unique()):
        kept = confirm.loc[confirm["cohort"].ne(cohort)]
        rows.append({"analysis": f"leave_out_{cohort}_AH1", "n": int(dedupe_patients(kept.loc[as_bool(kept, "selected_risk")]).shape[0])})
        rows.append({"analysis": f"leave_out_{cohort}_AH2", "n": int(dedupe_patients(kept.loc[as_bool(kept, "selected_eval")]).shape[0])})
    return rows


def join_labels(labels: pd.DataFrame, pred: pd.DataFrame) -> pd.DataFrame:
    key = ["cohort", "sample_id"]
    merged = labels.merge(pred, on=key, how="outer", suffixes=("", "_pred"), indicator=True)
    if not (merged["_merge"] == "both").all():
        raise SystemExit("label and prediction keys differ")
    if len(merged) != len(labels):
        raise SystemExit(f"joined rows {len(merged)} != label rows {len(labels)}")
    # Prefer prediction copies of shared descriptive columns when present.
    for column in labels.columns:
        other = f"{column}_pred"
        if other in merged.columns and column not in key:
            merged[column] = merged[other].where(merged[other].notna(), merged[column])
            merged = merged.drop(columns=other)
    merged = merged.drop(columns="_merge")
    return merged


def close(left, right) -> bool:
    if left is None or right is None or (isinstance(left, float) and math.isnan(left) and (right is None or (isinstance(right, float) and math.isnan(right)))):
        return True
    if isinstance(left, (bool, np.bool_)) or isinstance(right, (bool, np.bool_)):
        return bool(left) == bool(right)
    try:
        if pd.isna(left) and pd.isna(right):
            return True
    except TypeError:
        pass
    try:
        return abs(float(left) - float(right)) <= 1e-9
    except (TypeError, ValueError):
        return str(left) == str(right)


def compare_hypotheses(computed: list[dict], table: pd.DataFrame, family: str) -> list[dict]:
    rows = []
    lookup = {row["hypothesis"]: row for row in computed}
    for record in table.to_dict("records"):
        name = record["hypothesis"]
        got = lookup[name]
        for field, key in (("n", "n"), ("n_favor", "n_favor"), ("n_against", "n_against"), ("p", "p"), ("tested", "tested")):
            script = record.get(field)
            if field == "p" and (script == "" or script is None):
                script_value = None
            elif field == "tested":
                script_value = str(script) in ("True", "true", "1")
            else:
                script_value = None if script in ("", None) else float(script) if field != "n" else int(float(script))
            if field in ("n", "n_favor", "n_against"):
                script_value = int(float(script))
            match = close(got[key], script_value)
            rows.append({"family": family, "hypothesis": name, "field": field, "independent": got[key], "script": script_value, "match": match})
        diff_field = "diff_sa_minus_comparator" if "diff_sa_minus_comparator" in record and record.get("diff_sa_minus_comparator") not in ("", None) else "diff"
        script_diff = record.get(diff_field)
        script_diff = None if script_diff in ("", None) else float(script_diff)
        rows.append({"family": family, "hypothesis": name, "field": "diff", "independent": got["diff"], "script": script_diff, "match": close(got["diff"], script_diff)})
        if name != "AH3":
            script_reject = str(record.get("reject")) in ("True", "true", "1")
            rows.append({"family": family, "hypothesis": name, "field": "reject", "independent": got["reject"], "script": script_reject, "match": bool(got["reject"]) == script_reject})
    return rows


def select_counts(labels: pd.DataFrame) -> dict:
    """One sample per (cohort, patient), which is how the frozen label columns were built.
    include_confirm is taken from the label file: GSE209998 is descriptive even though its risk count is at least 10.
    """
    labels = labels.copy()
    out = {}
    for column, mask_col in (
        ("selected_eval", "in_eval"), ("selected_risk", "in_risk"),
        ("selected_native", "in_native"), ("selected_pool_out", "in_pool_out_risk"),
    ):
        part = labels.loc[labels[mask_col]].copy()
        part["_library_rank"] = np.where(part["library"].eq("polyA"), 0, 1)
        part = part.sort_values(["cohort", "patient_id", "_library_rank", "sample_id"], kind="mergesort")
        keep = part.drop_duplicates(["cohort", "patient_id"]).index
        chosen = pd.Series(False, index=labels.index)
        chosen.loc[keep] = True
        if column in labels.columns and not np.array_equal(chosen.to_numpy(), labels[column].to_numpy()):
            raise SystemExit(f"recomputed {column} differs from the label file")
        out[column] = chosen
    conf = labels["include_confirm"] & labels["layer"].eq("RNA-seq")
    return {
        "confirm_eval": int(out["selected_eval"][conf].sum()),
        "confirm_risk": int(out["selected_risk"][conf].sum()),
        "confirm_native": int(out["selected_native"][conf].sum()),
        "confirm_pool": int(out["selected_pool_out"][conf].sum()),
    }


def main() -> None:
    mode = sys.argv[1]
    if mode == "--sets":
        labels = read_table(Path(sys.argv[2]))
        counts = select_counts(labels)
        print(counts, flush=True)
        expected = {"confirm_eval": 729, "confirm_risk": 427, "confirm_native": 315, "confirm_pool": 71}
        if counts != expected:
            raise SystemExit(f"set sizes {counts}")
        return
    if mode == "--compare":
        labels = read_table(Path(sys.argv[2]))
        pred = pd.read_parquet(sys.argv[3])
        for column in BOOL_COLS:
            if column in pred.columns and pred[column].dtype != bool:
                pred[column] = pred[column].astype(str).eq("True")
        if "library" in pred.columns:
            pred["library"] = pred["library"].fillna("")
        # Prediction file already carries the label columns. Join checks the keys.
        key_cols = ["cohort", "sample_id"]
        label_keys = labels[key_cols].astype(str)
        pred_keys = pred[key_cols].astype(str)
        merged_keys = label_keys.merge(pred_keys, on=key_cols, how="outer", indicator=True)
        if not (merged_keys["_merge"] == "both").all() or len(merged_keys) != len(labels):
            raise SystemExit(f"join mismatch labels {len(labels)} predictions {len(pred)} both {(merged_keys['_merge']=='both').sum()}")
        tables = Path(sys.argv[4])
        confirm = pred.loc[as_bool(pred, "include_confirm") & pred["layer"].eq("RNA-seq")].copy()
        primary, secondary = hypothesis(confirm)
        rows = []
        rows.extend(compare_hypotheses(primary, pd.read_csv(tables / "hypothesis_primary.tsv", sep="\t"), "primary"))
        rows.extend(compare_hypotheses(secondary, pd.read_csv(tables / "hypothesis_secondary.tsv", sep="\t"), "secondary"))
        # Layer metrics from the script's cohort-all rows.
        metrics = pd.read_csv(tables / "metrics.tsv", sep="\t")
        for layer, part in pred.groupby("layer", sort=False):
            indep = {row["method"]: row for row in layer_metrics(part)}
            script = metrics.loc[metrics["layer"].eq(layer) & metrics["cohort"].eq("all") & metrics["standard_site"].eq("all")]
            for record in script.to_dict("records"):
                got = indep[record["method"]]
                for field in ("n", "top1", "n_at_risk", "host_rate", "n_native_truth", "native_truth_top1"):
                    rows.append({
                        "family": "metrics", "hypothesis": f"{layer}:{record['method']}", "field": field,
                        "independent": got[field], "script": record[field], "match": close(got[field], record[field]),
                    })
        conf_n = {(row["layer"], row["level"]): row["n"] for row in conformal_n(pred)}
        conformal = pd.read_csv(tables / "conformal.tsv", sep="\t")
        for record in conformal.to_dict("records"):
            got = conf_n[(record["layer"], record["level"])]
            rows.append({
                "family": "conformal", "hypothesis": f"{record['layer']}:{record['method']}:{record['variant']}:{record['alpha']}",
                "field": record["level"], "independent": got, "script": int(record["n"]), "match": int(got) == int(record["n"]),
            })
        meta = pd.read_csv(tables / "meta_analysis.tsv", sep="\t")
        indep_meta = meta_records(pred)
        # Compare cohort b,c,n and DL estimate.
        for record in indep_meta:
            if record["kind"] == "cohort":
                script = meta.loc[meta["layer"].eq(record["layer"]) & meta["cohort"].eq(record["cohort"])]
            else:
                script = meta.loc[meta["layer"].eq(record["layer"]) & meta["status"].eq("계산")]
            if script.empty:
                rows.append({"family": "meta", "hypothesis": str(record), "field": "missing", "independent": record, "script": None, "match": False})
                continue
            script_row = script.iloc[0]
            if record["kind"] == "cohort":
                for field in ("n", "b", "c", "d", "v"):
                    rows.append({
                        "family": "meta", "hypothesis": f"{record['layer']}:{record['cohort']}", "field": field,
                        "independent": record[field], "script": script_row[field], "match": close(record[field], script_row[field]),
                    })
            else:
                rows.append({
                    "family": "meta", "hypothesis": record["layer"], "field": "estimate",
                    "independent": record["estimate"], "script": script_row["estimate"], "match": close(record["estimate"], script_row["estimate"]),
                })
        sens = pd.read_csv(tables / "sensitivity.tsv", sep="\t")
        for record in sensitivity_sizes(pred):
            script = sens.loc[sens["analysis"].eq(record["analysis"])]
            script_n = script.iloc[0]["n"] if len(script) else None
            rows.append({
                "family": "sensitivity", "hypothesis": record["analysis"], "field": "n",
                "independent": record["n"], "script": script_n, "match": close(record["n"], script_n),
            })
            if "n_favor" in record and "n_favor" in sens.columns:
                rows.append({
                    "family": "sensitivity", "hypothesis": record["analysis"], "field": "n_favor",
                    "independent": record["n_favor"], "script": script.iloc[0]["n_favor"],
                    "match": close(record["n_favor"], script.iloc[0]["n_favor"]),
                })
                rows.append({
                    "family": "sensitivity", "hypothesis": record["analysis"], "field": "n_against",
                    "independent": record["n_against"], "script": script.iloc[0]["n_against"],
                    "match": close(record["n_against"], script.iloc[0]["n_against"]),
                })
        out = pd.DataFrame(rows)
        dest = Path(sys.argv[5])
        dest.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(dest, sep="\t", index=False)
        bad = out.loc[~out["match"].astype(bool)]
        print("comparisons", len(out), "mismatches", len(bad), flush=True)
        if len(bad):
            print(bad.head(30).to_string(index=False), flush=True)
            raise SystemExit(1)
        return
    if mode == "--expect-planted":
        tables = Path(sys.argv[2])
        primary = pd.read_csv(tables / "hypothesis_primary.tsv", sep="\t")
        secondary = pd.read_csv(tables / "hypothesis_secondary.tsv", sep="\t")
        sens = pd.read_csv(tables / "sensitivity.tsv", sep="\t")
        expected = {
            "AH1": (30, 5, True),
            "AH2": (33, 13, True),
            "AS1": (8, 3, False),
            "AS2": (12, 1, True),
            "AS3": (8, 3, False),
        }
        lookup = {row.hypothesis: row for row in primary.itertuples(index=False)}
        lookup.update({row.hypothesis: row for row in secondary.itertuples(index=False)})
        for name, (favor, against, reject) in expected.items():
            row = lookup[name]
            p = one_sided(favor, against)
            if int(row.n_favor) != favor or int(row.n_against) != against:
                raise SystemExit(f"{name} counts {row.n_favor}/{row.n_against}")
            if abs(float(row.p) - p) > 1e-9:
                raise SystemExit(f"{name} p {row.p} != {p}")
            if bool(row.reject) != reject:
                raise SystemExit(f"{name} reject {row.reject}")
        ah3 = lookup["AH3"]
        if int(ah3.n_favor) != 3 or int(ah3.n_against) != 8 or not bool(ah3.tested):
            raise SystemExit(f"AH3 {ah3.n_favor} {ah3.n_against} tested {ah3.tested}")
        if abs(float(ah3.diff_sa_minus_comparator) - (-5 / int(ah3.n))) > 1e-9:
            raise SystemExit(f"AH3 diff {ah3.diff_sa_minus_comparator}")
        sens_row = sens.loc[sens["analysis"].eq("all_samples_cluster_AH1")].iloc[0]
        if int(sens_row.n_favor) != 40 or int(sens_row.n_against) != 5:
            raise SystemExit(f"sensitivity {sens_row.n_favor} {sens_row.n_against}")
        print("planted-expect ok", "AH3_n", int(ah3.n), flush=True)
        return
    raise SystemExit("usage")


if __name__ == "__main__":
    main()
