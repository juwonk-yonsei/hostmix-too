#!/usr/bin/env python3
"""Post hoc microarray and missing-gene ablation. The A15 plan must already be committed.

Does not set a thread cap, does not read n_threads, and does not import loader.py or genes.py.
Does not overwrite an existing stage-10 output.
"""
from __future__ import annotations

import hashlib
import importlib
import json
import subprocess
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from metrics import ORGANS  # noqa: E402
from stage2_common import load_gene_pack, load_tpm, sum_1e6  # noqa: E402
from stage2_fit import predict_bundle  # noqa: E402
from stage5_lib import rank_b0  # noqa: E402

ablation = importlib.import_module("45_ablation_A12")
mask = importlib.import_module("36_mask_A9")

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage10" / "ablation_microarray"
PLAN = ROOT / "config" / "analysis_plan_A15.yaml"
MODELS = ["BASE-Z", "PURE-Zs", "PURE-Zs-w", "MIX-Z0", "SA-Z"]
PLATFORM_ORDER = ["GPL96", "GPL20769", "GPL15659", "prad_fhcrc_agilent"]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(1 << 20)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def require_plan() -> str:
    if not PLAN.exists():
        raise SystemExit("config/analysis_plan_A15.yaml is missing")
    head = subprocess.run(["git", "cat-file", "-e", "HEAD:config/analysis_plan_A15.yaml"], cwd=ROOT)
    if head.returncode != 0:
        raise SystemExit("the plan is not in HEAD; scoring was not started")
    dirty = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", "config/analysis_plan_A15.yaml"], cwd=ROOT)
    if dirty.returncode != 0:
        raise SystemExit("the plan differs from the committed file; scoring was not started")
    line = subprocess.check_output(
        ["git", "log", "-1", "--format=%H %ci", "--", "config/analysis_plan_A15.yaml"], cwd=ROOT, text=True
    ).strip()
    print("plan", line, flush=True)
    return line


def as_bool(frame: pd.DataFrame, column: str) -> pd.Series:
    return frame[column].map(lambda value: str(value).strip().lower() in {"true", "1"})


def dedupe_patients(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty or frame["patient_id"].is_unique:
        return frame.reset_index(drop=True)
    part = frame.copy()
    library = part["library"].fillna("")
    part["_library_rank"] = np.where(library.eq("polyA"), 0, 1)
    part = part.sort_values(["_library_rank", "sample_id", "cohort"], kind="mergesort")
    return part.drop_duplicates("patient_id", keep="first").drop(columns="_library_rank").reset_index(drop=True)


def z_matrix(tpm: np.ndarray, genes: list[str], keep: set[str] | None, b0: list[str]) -> np.ndarray:
    if keep is None:
        block = tpm
        present = genes
    else:
        columns = [i for i, gene in enumerate(genes) if gene in keep]
        block = np.ascontiguousarray(tpm[:, columns])
        present = [genes[i] for i in columns]
    z_all = mask.rank_z(block)
    present_index = {gene: i for i, gene in enumerate(present)}
    z = np.zeros((tpm.shape[0], len(b0)), dtype=np.float32)
    for j, symbol in enumerate(b0):
        hit = present_index.get(symbol)
        if hit is not None:
            z[:, j] = z_all[:, hit]
    return z


def call_organs(bundle: dict, matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    organ, pred, *_rest = predict_bundle(bundle, matrix)
    return organ, np.asarray(pred, dtype=object)


def load_bundles() -> dict:
    stage2 = ROOT / "results" / "stage2" / "models"
    fitted = ROOT / "results" / "stage9" / "ablation" / "models"
    return {
        "BASE-Z": joblib.load(stage2 / "BASE_Z.joblib"),
        "SA-Z": joblib.load(stage2 / "SA_Z.joblib"),
        "PURE-Zs": joblib.load(fitted / "PURE-Zs.joblib"),
        "PURE-Zs-w": joblib.load(fitted / "PURE-Zs-w.joblib"),
        "MIX-Z0": joblib.load(fitted / "MIX-Z0.joblib"),
    }


def microarray_frame() -> pd.DataFrame:
    columns = [
        "cohort", "patient_id", "sample_id", "library", "organ", "native_organs",
        "layer", "layer_test", "selected_eval", "selected_risk", "selected_native",
        "BASE-Z__pred", "SA-Z__pred",
    ]
    columns += [f"{method}__p_{organ}" for method in ("BASE-Z", "SA-Z") for organ in ORGANS]
    pred = pd.read_parquet(ROOT / "results/stage7/confirm/predictions.parquet", columns=columns)
    part = pred.loc[pred["layer"].eq("마이크로어레이")].copy()
    for column in ("sample_id", "patient_id", "library", "organ", "native_organs"):
        part[column] = part[column].fillna("").astype(str)
    return part.reset_index(drop=True)


def attach_microarray_z(frame: pd.DataFrame) -> np.ndarray:
    rows = []
    cache: dict[str, tuple[np.ndarray, dict[str, int]]] = {}
    for cohort, sample in zip(frame["cohort"], frame["sample_id"]):
        if cohort not in cache:
            folder = ROOT / "data" / "processed" / "aux" / cohort
            ids = [line for line in (folder / "sample_ids.txt").read_text().splitlines() if line]
            cache[cohort] = (np.load(folder / "Z.npy"), {item: i for i, item in enumerate(ids)})
        matrix, index = cache[cohort]
        if sample not in index:
            raise SystemExit(f"{cohort} Z.npy has no {sample}")
        rows.append(matrix[index[sample]])
    return np.vstack(rows).astype(np.float32)


def identity_microarray(bundles: dict, frame: pd.DataFrame, matrix: np.ndarray) -> list[dict]:
    rows = []
    for method in ("BASE-Z", "SA-Z"):
        organ, pred = call_organs(bundles[method], matrix)
        stored = frame[f"{method}__pred"].astype(str).to_numpy()
        mismatch = int(np.sum(pred != stored))
        held = frame[[f"{method}__p_{organ_name}" for organ_name in ORGANS]].to_numpy(dtype=np.float64)
        gap = float(np.max(np.abs(organ - held)))
        print("microarray", method, "mismatch", mismatch, "proba_gap", gap, flush=True)
        rows.append({"check": "microarray_rescore", "method": method, "label_mismatch": mismatch, "max_abs_probability": gap, "pass": mismatch == 0 and gap <= 1e-12})
    return rows


def cohort_matrices(paths: dict):
    genes, b0_idx, _offsets, _indices, _names = load_gene_pack(paths)
    import importlib as _importlib
    score = _importlib.import_module("31_score_A8")
    met = score.met500()
    print("MET500 features", flush=True)
    met_samples, met_tpm = load_tpm(paths["data_processed"] / "met500" / "tpm_G.h5")
    met_row = {sample: i for i, sample in enumerate(met_samples)}
    met_z = rank_b0(sum_1e6(met_tpm[[met_row[sample] for sample in met["sample_id"]]]), b0_idx)
    del met_tpm
    print("POG features", flush=True)
    expr = pd.read_csv(paths["data_raw"] / "POG570" / "POG570_TPM_expression.txt.gz", sep="\t", index_col=0)
    pog_matrix, pog_ids = ablation.collapse_to_g(expr, genes, ablation.ensembl_base_to_symbol(paths))
    del expr
    pog_z = rank_b0(sum_1e6(pog_matrix), b0_idx)
    del pog_matrix
    pog = score.pog570()
    pog_pos = {sample: i for i, sample in enumerate(pog_ids)}
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
        folder = ablation.feature_cohort(str(cohort), str(library))
        if folder not in cache:
            ids = [line for line in (paths["data_processed"] / "aux" / folder / "sample_ids.txt").read_text().splitlines() if line]
            cache[folder] = (np.load(paths["data_processed"] / "aux" / folder / "Z.npy"), {item: i for i, item in enumerate(ids)})
        matrix, index = cache[folder]
        if sample not in index:
            raise SystemExit(f"{folder} missing {sample}")
        rows_z.append(matrix[index[sample]])
    aux_z = np.vstack(rows_z).astype(np.float32)
    return {"MET500": (met, met_z), "POG570": (pog, pog_z[pog_take]), "aux_rnaseq": (aux, aux_z)}


def identity_stored_cohorts(bundles: dict, paths: dict) -> list[dict]:
    stored = pd.read_csv(ROOT / "results/stage9/ablation/predictions.tsv", sep="\t", dtype=str)
    rows = []
    built = cohort_matrices(paths)
    for cohort, (frame, matrix) in built.items():
        for name in ("PURE-Zs", "PURE-Zs-w", "MIX-Z0"):
            _organ, pred = call_organs(bundles[name], matrix)
            held = stored.loc[stored["cohort"].eq(cohort) & stored["method"].eq(name), ["sample_id", "pred"]]
            joined = frame[["sample_id"]].copy()
            joined["sample_id"] = joined["sample_id"].astype(str)
            joined["pred"] = pred
            merged = joined.merge(held, on="sample_id", how="left", suffixes=("_new", "_stored"))
            if merged["pred_stored"].isna().any() or len(merged) != len(frame):
                raise SystemExit(f"stored predictions do not cover {cohort} {name}")
            mismatch = int(np.sum(merged["pred_new"].to_numpy(dtype=object) != merged["pred_stored"].to_numpy(dtype=object)))
            print("stored cohort", cohort, name, "mismatch", mismatch, flush=True)
            rows.append({"check": "ablation_cohort_rescore", "cohort": cohort, "method": name, "label_mismatch": mismatch, "pass": mismatch == 0})
    del built
    return rows


def top1(pred: np.ndarray, truth: np.ndarray) -> float:
    return float(np.mean(pred == truth))


def identity_and_restriction(bundles: dict) -> tuple[list[dict], pd.DataFrame]:
    genes = pd.read_csv(ROOT / "data/processed/genes/G_symbols.txt", header=None)[0].astype(str).tolist()
    b0 = pd.read_csv(ROOT / "results/B0_genes.txt")["symbol"].astype(str).tolist()
    platforms = mask.platform_symbols(genes)
    expr, _gene_index = mask.load_eval(genes)
    unmasked = {}
    original_top1 = {}
    for cohort, (tpm, meta) in expr.items():
        print("unmasked", cohort, flush=True)
        z = z_matrix(tpm, genes, None, b0)
        unmasked[cohort] = z
        truth = meta["truth"].astype(str).to_numpy()
        for method in ("BASE-Z", "SA-Z"):
            _organ, pred = call_organs(bundles[method], z)
            stored = meta[f"pred__{method}"].astype(str).to_numpy()
            mismatch = int(np.sum(pred != stored))
            if mismatch:
                raise SystemExit(f"unmasked {cohort} {method} mismatch {mismatch}")
            original_top1[(cohort, method)] = top1(pred, truth)
    stored_metrics = pd.read_csv(ROOT / "results/stage9/mask_metrics.tsv", sep="\t")
    identity = []
    locked = {}

    def one_drop(kind: str, name: str, cohort: str, meta: pd.DataFrame, z: np.ndarray, method: str, rep: str) -> dict:
        truth = meta["truth"].astype(str).to_numpy()
        if (cohort, method) not in original_top1:
            _base_organ, base_pred = call_organs(bundles[method], unmasked[cohort])
            original_top1[(cohort, method)] = top1(base_pred, truth)
        _organ, pred = call_organs(bundles[method], z)
        change = top1(pred, truth) - original_top1[(cohort, method)]
        row = {"model": method, "platform": name, "cohort": cohort, "kind": kind, "rep": rep, "top1": top1(pred, truth), "drop": -change}
        if method in {"BASE-Z", "SA-Z"}:
            held = stored_metrics.loc[
                stored_metrics["cohort"].eq(cohort) & stored_metrics["mask"].eq(name)
                & stored_metrics["rep"].eq(rep) & stored_metrics["method"].eq(method)
            ].iloc[0]
            gap = abs(float(change) - float(held["top1_change"]))
            identity.append({"check": f"{kind}_drop", "cohort": cohort, "mask": name, "rep": rep, "method": method, "max_abs": gap, "pass": gap <= 1e-12})
            if gap > 1e-12:
                print("drop mismatch", kind, cohort, name, rep, method, gap, flush=True)
        return row

    def walk(methods: list[str]) -> list[dict]:
        rows = []
        rng = np.random.default_rng(20261023)
        b0_index = np.arange(len(b0))
        for name in PLATFORM_ORDER:
            keep = platforms[name]
            for cohort, (tpm, meta) in expr.items():
                key = ("platform", name, cohort)
                if key not in locked:
                    print("platform", name, cohort, flush=True)
                    locked[key] = z_matrix(tpm, genes, keep, b0)
                for method in methods:
                    rows.append(one_drop("platform", name, cohort, meta, locked[key], method, "platform"))
        for name in PLATFORM_ORDER:
            n_drop = len([gene for gene in b0 if gene not in platforms[name]])
            for rep in range(10):
                chosen = {b0[int(i)] for i in rng.choice(b0_index, size=n_drop, replace=False)}
                keep = set(genes) - chosen
                for cohort, (tpm, meta) in expr.items():
                    key = ("random", name, cohort, rep)
                    if key not in locked:
                        print("random", name, rep, cohort, flush=True)
                        locked[key] = z_matrix(tpm, genes, keep, b0)
                    for method in methods:
                        rows.append(one_drop("random", name, cohort, meta, locked[key], method, str(rep)))
        return rows

    base_rows = walk(["BASE-Z", "SA-Z"])
    if not all(row["pass"] for row in identity):
        return identity, pd.DataFrame(base_rows)
    extra = walk(["PURE-Zs", "PURE-Zs-w", "MIX-Z0"])
    return identity, pd.DataFrame(base_rows + extra)


def host_rate(frame: pd.DataFrame, pred: np.ndarray) -> float:
    truth = frame["organ"].to_numpy(dtype=object)
    natives = [set(filter(None, str(text).split("|"))) for text in frame["native_organs"]]
    risk = np.array([bool(natives[i]) and truth[i] not in natives[i] for i in range(len(frame))])
    if not risk.any():
        return float("nan")
    pulled = risk & np.array([pred[i] in natives[i] for i in range(len(pred))])
    return float(pulled[risk].mean())


def score_microarray(bundles: dict, frame: pd.DataFrame, matrix: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    tested = frame.loc[as_bool(frame, "layer_test")].copy()
    tested["_row"] = np.arange(len(frame))[as_bool(frame, "layer_test").to_numpy()]
    sets = {
        "evaluation": dedupe_patients(tested.loc[as_bool(tested, "selected_eval")]),
        "at_risk": dedupe_patients(tested.loc[as_bool(tested, "selected_risk")]),
        "native_truth": dedupe_patients(tested.loc[as_bool(tested, "selected_native")]),
    }
    expected = {"evaluation": 536, "at_risk": 207, "native_truth": 341}
    for name, part in sets.items():
        if len(part) != expected[name]:
            raise SystemExit(f"microarray {name} n {len(part)}")
    metric_rows = []
    esophagus_rows = []
    calls = {}
    for method in MODELS:
        _organ, pred = call_organs(bundles[method], matrix)
        calls[method] = pred
    for method, pred in calls.items():
        for set_name, part in sets.items():
            taken = pred[part["_row"].to_numpy()]
            grouped = [("all", part, taken)]
            for cohort, sub in part.groupby("cohort", sort=False):
                grouped.append((str(cohort), sub, pred[sub["_row"].to_numpy()]))
            for cohort, sub, values in grouped:
                if set_name == "evaluation":
                    metric_rows.append({"model": method, "set": set_name, "cohort": cohort, "metric": "top1", "value": float(np.mean(values == sub["organ"].to_numpy(dtype=object))), "n": int(len(sub))})
                    n_eso = int(np.sum(values == "Esophagus"))
                    metric_rows.append({"model": method, "set": set_name, "cohort": cohort, "metric": "n_esophagus", "value": n_eso, "n": int(len(sub))})
                    esophagus_rows.append({"model": method, "cohort": cohort, "n_esophagus": n_eso, "n": int(len(sub))})
                elif set_name == "at_risk":
                    metric_rows.append({"model": method, "set": set_name, "cohort": cohort, "metric": "host_rate", "value": host_rate(sub, values), "n": int(len(sub))})
                else:
                    metric_rows.append({"model": method, "set": set_name, "cohort": cohort, "metric": "native_truth_top1", "value": float(np.mean(values == sub["organ"].to_numpy(dtype=object))), "n": int(len(sub))})
    return pd.DataFrame(metric_rows), pd.DataFrame(esophagus_rows)


def reading(metrics: pd.DataFrame, restriction: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    def value(model: str, set_name: str, metric: str) -> float:
        hit = metrics.loc[metrics["model"].eq(model) & metrics["set"].eq(set_name) & metrics["cohort"].eq("all") & metrics["metric"].eq(metric)]
        return float(hit.iloc[0]["value"])

    accuracy = {model: value(model, "evaluation", "top1") for model in MODELS}
    g = accuracy["BASE-Z"] - accuracy["SA-Z"]
    shares = {
        "s_mix": (accuracy["BASE-Z"] - accuracy["MIX-Z0"]) / g,
        "s_std": (accuracy["BASE-Z"] - accuracy["PURE-Zs"]) / g,
        "s_stdw": (accuracy["BASE-Z"] - accuracy["PURE-Zs-w"]) / g,
    }
    random = restriction.loc[restriction["kind"].eq("random")]
    mean_drop = random.groupby(["model", "platform", "cohort"])["drop"].mean().reset_index()
    r_value = {model: float(part["drop"].mean()) for model, part in mean_drop.groupby("model")}
    g2 = r_value["SA-Z"] - r_value["BASE-Z"]
    if g2 < 0.02:
        second = "Z"
        s2 = {}
    else:
        s2 = {
            "s2_mix": (r_value["MIX-Z0"] - r_value["BASE-Z"]) / g2,
            "s2_std": (r_value["PURE-Zs"] - r_value["BASE-Z"]) / g2,
            "s2_stdw": (r_value["PURE-Zs-w"] - r_value["BASE-Z"]) / g2,
        }
        second = category(s2["s2_mix"], s2["s2_std"], s2["s2_stdw"])
    share_rows = [{"contrast": "accuracy", "name": name, "value": val, "gap": g} for name, val in shares.items()]
    share_rows += [{"contrast": "random_removal", "name": name, "value": val, "gap": g2} for name, val in s2.items()]
    categories = pd.DataFrame([
        {"contrast": "accuracy", "category": category(shares["s_mix"], shares["s_std"], shares["s_stdw"]), "gap": g, **shares},
        {"contrast": "random_removal", "category": second, "gap": g2, **{key: r_value[key] for key in MODELS}, **s2},
    ])
    return pd.DataFrame(share_rows), categories


def category(s_mix: float, s_std: float, s_stdw: float) -> str:
    low = min(s_std, s_stdw)
    high = max(s_std, s_stdw)
    if s_mix < 0.5 and low >= 0.5:
        return "S"
    if s_mix >= 0.5 and high < 0.5:
        return "M"
    if s_mix >= 0.5 and low >= 0.5:
        return "B"
    return "O"


def restriction_table(records: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (model, platform, cohort), part in records.groupby(["model", "platform", "cohort"]):
        platform_drop = float(part.loc[part["kind"].eq("platform"), "drop"].iloc[0])
        random = part.loc[part["kind"].eq("random"), "drop"]
        rows.append({
            "model": model, "platform": platform, "cohort": cohort,
            "restriction_drop": platform_drop,
            "random_mean": float(random.mean()), "random_min": float(random.min()), "random_max": float(random.max()),
        })
    return pd.DataFrame(rows)


def main() -> None:
    if (OUT / "metrics.tsv").exists():
        raise SystemExit(f"{OUT / 'metrics.tsv'} already exists; nothing was overwritten")
    plan_line = require_plan()
    OUT.mkdir(parents=True, exist_ok=True)
    bundles = load_bundles()
    frame = microarray_frame()
    matrix = attach_microarray_z(frame)
    identity = identity_microarray(bundles, frame, matrix)
    paths = ablation.project_paths()
    identity.extend(identity_stored_cohorts(bundles, paths))
    if not all(row["pass"] for row in identity):
        pd.DataFrame(identity).to_csv(OUT / "identity.tsv", sep="\t", index=False)
        raise SystemExit("identity failed before gene-restricted scoring of the ablation models")
    more, records = identity_and_restriction(bundles)
    identity.extend(more)
    pd.DataFrame(identity).to_csv(OUT / "identity.tsv", sep="\t", index=False)
    if not all(row["pass"] for row in identity):
        raise SystemExit("identity failed; metrics were not written")
    metrics, esophagus = score_microarray(bundles, frame, matrix)
    metrics.to_csv(OUT / "metrics.tsv", sep="\t", index=False)
    esophagus.to_csv(OUT / "esophagus_calls.tsv", sep="\t", index=False)
    restriction = restriction_table(records)
    restriction.to_csv(OUT / "restriction.tsv", sep="\t", index=False)
    shares, categories = reading(metrics, records)
    shares.to_csv(OUT / "shares.tsv", sep="\t", index=False)
    categories.to_csv(OUT / "category.tsv", sep="\t", index=False)
    hashed = {
        "config/analysis_plan_A15.yaml": sha256_file(PLAN),
        "results/stage2/models/BASE_Z.joblib": sha256_file(ROOT / "results/stage2/models/BASE_Z.joblib"),
        "results/stage2/models/SA_Z.joblib": sha256_file(ROOT / "results/stage2/models/SA_Z.joblib"),
        "results/stage9/ablation/models/PURE-Zs.joblib": sha256_file(ROOT / "results/stage9/ablation/models/PURE-Zs.joblib"),
        "results/stage9/ablation/models/PURE-Zs-w.joblib": sha256_file(ROOT / "results/stage9/ablation/models/PURE-Zs-w.joblib"),
        "results/stage9/ablation/models/MIX-Z0.joblib": sha256_file(ROOT / "results/stage9/ablation/models/MIX-Z0.joblib"),
        "results/stage9/ablation/predictions.tsv": sha256_file(ROOT / "results/stage9/ablation/predictions.tsv"),
        "results/stage9/mask_metrics.tsv": sha256_file(ROOT / "results/stage9/mask_metrics.tsv"),
    }
    (OUT / "run.json").write_text(json.dumps({"plan_commit": plan_line, "sha256": hashed, "status": "post hoc"}, indent=2))
    print(categories.to_string(index=False), flush=True)
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
