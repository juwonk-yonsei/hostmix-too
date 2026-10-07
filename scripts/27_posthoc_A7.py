#!/usr/bin/env python3
"""Post-hoc tables after the auxiliary confirmation commit. Does not rewrite frozen features."""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from stage2_common import load_gene_pack, rank_z, score_ks, sum_1e6
from stage2_fit import predict_bundle

OUT = ROOT / "results" / "stage7" / "posthoc"
PRED = ROOT / "results" / "stage7" / "confirm" / "predictions.parquet"
POOL_OUT = {
    "kidney", "pancreas", "bladder", "thyroid", "ovary", "breast", "stomach",
    "colon_rectum", "prostate", "esophagus", "head_neck", "cervix",
}
PROXY = {"lymph_node", "bone"}
REAL_HOST = {"liver", "lung", "soft_tissue", "skin", "subcutaneous", "brain", "adrenal", "omentum", "peritoneum", "bone_marrow"}
BOOLS = (
    "in_eval", "in_risk", "in_native", "in_pool_out_risk", "excluded", "include_confirm",
    "selected_eval", "selected_risk", "selected_native", "selected_pool_out",
)


def paths() -> dict:
    raw = yaml.safe_load((ROOT / "config" / "paths.yaml").read_text())
    out = {"root": ROOT}
    for key, value in raw.items():
        if key in ("seed", "n_threads"):
            continue
        out[key] = ROOT / value
    return out


def as_bool(frame: pd.DataFrame, column: str) -> pd.Series:
    if frame[column].dtype == bool:
        return frame[column]
    return frame[column].astype(str).eq("True")


def dedupe(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty or frame["patient_id"].is_unique:
        return frame.reset_index(drop=True)
    part = frame.copy()
    library = part["library"].fillna("")
    part["_r"] = np.where(library.eq("polyA"), 0, 1)
    part = part.sort_values(["_r", "sample_id", "cohort"], kind="mergesort")
    return part.drop_duplicates("patient_id").drop(columns="_r").reset_index(drop=True)


def native_sets(frame: pd.DataFrame) -> list[set[str]]:
    return [set(filter(None, str(text).split("|"))) for text in frame["native_organs"].fillna("")]


def host_error(frame: pd.DataFrame, method: str) -> np.ndarray:
    pred = frame[f"{method}__pred"].to_numpy(dtype=object)
    truth = frame["organ"].to_numpy(dtype=object)
    native = native_sets(frame)
    return np.array([bool(native[i]) and truth[i] not in native[i] and pred[i] in native[i] for i in range(len(frame))])


def correct(frame: pd.DataFrame, method: str) -> np.ndarray:
    return frame[f"{method}__pred"].to_numpy(dtype=object) == frame["organ"].to_numpy(dtype=object)


def bootstrap_diff(left: np.ndarray, right: np.ndarray, rng: np.random.Generator) -> tuple[float, float, float]:
    point = float(left.mean() - right.mean()) if len(left) else float("nan")
    if len(left) == 0:
        return point, float("nan"), float("nan")
    draws = rng.integers(0, len(left), size=(2000, len(left)))
    diffs = left[draws].mean(axis=1) - right[draws].mean(axis=1)
    return point, float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def cluster_diff(frame: pd.DataFrame, left: np.ndarray, right: np.ndarray, rng: np.random.Generator) -> tuple[float, float, float]:
    patients = frame["patient_id"].to_numpy()
    unique = pd.unique(patients)
    groups = [np.flatnonzero(patients == patient) for patient in unique]
    point = float(left.mean() - right.mean())
    draws = rng.integers(0, len(unique), size=(2000, len(unique)))
    diffs = np.empty(2000)
    for i, draw in enumerate(draws):
        idx = np.concatenate([groups[j] for j in draw])
        diffs[i] = left[idx].mean() - right[idx].mean()
    return point, float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def effect_row(name: str, frame: pd.DataFrame, favor: np.ndarray, against: np.ndarray, left: np.ndarray, right: np.ndarray, rng) -> dict:
    diff, low, high = bootstrap_diff(left, right, rng)
    return {
        "analysis": name, "n": int(len(frame)), "n_favor": int(favor.sum()) if len(frame) else 0,
        "n_against": int(against.sum()) if len(frame) else 0, "diff": diff, "ci_low": low, "ci_high": high,
    }


def load_predictions() -> pd.DataFrame:
    frame = pd.read_parquet(PRED)
    for column in BOOLS:
        if column in frame.columns and frame[column].dtype != bool:
            frame[column] = frame[column].astype(str).eq("True")
    for column in ("library", "standard_site", "organ", "native_organs", "raw_site", "patient_id", "sample_id"):
        if column in frame.columns:
            frame[column] = frame[column].fillna("").astype(str)
    return frame


def discordant(pred: pd.DataFrame) -> None:
    confirm = pred.loc[as_bool(pred, "include_confirm") & pred["layer"].eq("RNA-seq")]
    pieces = []
    risk = dedupe(confirm.loc[as_bool(confirm, "selected_risk")])
    base, sa = host_error(risk, "BASE-Z"), host_error(risk, "SA-Z")
    for kind, mask in (("AH1_favor", base & ~sa), ("AH1_against", sa & ~base)):
        part = risk.loc[mask].copy()
        part.insert(0, "case", kind)
        pieces.append(part)
    pool = dedupe(confirm.loc[as_bool(confirm, "selected_pool_out")])
    sa, pool22 = host_error(pool, "SA-Z"), host_error(pool, "SA-pool22")
    for kind, mask in (("AS2_favor", sa & ~pool22), ("AS2_against", pool22 & ~sa)):
        part = pool.loc[mask].copy()
        part.insert(0, "case", kind)
        pieces.append(part)
    cols = [
        "case", "cohort", "patient_id", "sample_id", "raw_site", "standard_site", "organ",
        "BASE-Z__pred", "BASE-Z__pmax", "SA-Z__pred", "SA-Z__pmax", "SA-pool22__pred", "SA-pool22__pmax", "beta",
    ]
    table = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=cols)
    for method in ("BASE-Z", "SA-Z", "SA-pool22"):
        pcols = [c for c in table.columns if c.startswith(f"{method}__p_")]
        table[f"{method}__pmax"] = table[pcols].to_numpy(dtype=float).max(axis=1) if len(table) and pcols else np.nan
    keep = [c for c in cols if c in table.columns]
    table[keep].to_csv(OUT / "discordant_AH1_AS2.tsv", sep="\t", index=False)


def site_groups(pred: pd.DataFrame) -> None:
    confirm = pred.loc[as_bool(pred, "include_confirm") & pred["layer"].eq("RNA-seq")]
    risk = dedupe(confirm.loc[as_bool(confirm, "selected_risk")])

    def group(site: str) -> str:
        if site in PROXY:
            return "proxy_host"
        if site in POOL_OUT:
            return "pool_out"
        if site in REAL_HOST:
            return "real_host_tissue"
        return "other"

    risk = risk.copy()
    risk["site_group"] = [group(site) for site in risk["standard_site"]]
    rows = []
    for name, part in risk.groupby("site_group", sort=False):
        base, sa = host_error(part, "BASE-Z"), host_error(part, "SA-Z")
        rows.append({
            "site_group": name, "n": int(len(part)),
            "base_host": int(base.sum()), "sa_host": int(sa.sum()),
            "base_host_rate": float(base.mean()) if len(part) else None,
            "sa_host_rate": float(sa.mean()) if len(part) else None,
            "diff_sa_minus_base": float(sa.mean() - base.mean()) if len(part) else None,
            "n_favor": int((base & ~sa).sum()), "n_against": int((sa & ~base).sum()),
            "sites": ",".join(sorted(set(part["standard_site"]))),
        })
    pd.DataFrame(rows).to_csv(OUT / "ah1_by_site_group.tsv", sep="\t", index=False)


def gse_prefix(pred: pd.DataFrame, rng: np.random.Generator) -> None:
    part = pred.loc[pred["cohort"].eq("GSE209998")].copy()
    part["prefix"] = part["sample_id"].str.extract(r"^AUR-([A-Z0-9]{4})-", expand=False)
    part["patient_id"] = part["prefix"]
    specs = {
        "AH1": ("in_risk", "SA-Z", "BASE-Z", "host"),
        "AH2": ("in_eval", "SA-Z", "BASE-Z", "correct"),
        "AH3": ("in_native", "SA-Z", "BASE-Z", "correct"),
        "AS1": ("in_native", "SA-G", "SA-Z", "correct"),
        "AS2": ("in_pool_out_risk", "SA-pool22", "SA-Z", "host"),
        "AS3": ("in_eval", "SA-G", "SA-Z", "correct"),
    }
    rows = []
    for name, (mask, left_m, right_m, kind) in specs.items():
        sub = part.loc[as_bool(part, mask) & part["prefix"].notna()].copy()
        one = dedupe(sub)
        flag = host_error if kind == "host" else correct
        left_one, right_one = flag(one, left_m), flag(one, right_m)
        left_all, right_all = flag(sub, left_m), flag(sub, right_m)
        if name in ("AH1", "AS2"):
            favor, against = right_one & ~left_one, left_one & ~right_one
        else:
            favor, against = left_one & ~right_one, right_one & ~left_one
        row = effect_row(name, one, favor, against, left_one.astype(float), right_one.astype(float), rng)
        cdiff, clow, chigh = cluster_diff(sub, left_all.astype(float), right_all.astype(float), rng) if len(sub) else (float("nan"), float("nan"), float("nan"))
        row.update({
            "one_sample_n": int(len(one)), "cluster_n_samples": int(len(sub)),
            "cluster_diff": cdiff, "cluster_ci_low": clow, "cluster_ci_high": chigh,
            "n_prefix_all_rows": int(part["prefix"].nunique(dropna=True)),
            "diff_definition": f"{left_m} minus {right_m}",
        })
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "gse209998_prefix.tsv", sep="\t", index=False)


def beta_table(pred: pd.DataFrame) -> None:
    rows = []
    missing = {}
    path = ROOT / "results" / "stage6" / "ld_inputs.json"
    if path.exists():
        for record in json.loads(path.read_text()):
            missing[record.get("cohort")] = record.get("n_ld_missing")
    for cohort, part in pred.groupby("cohort", sort=False):
        beta = pd.to_numeric(part["beta"], errors="coerce")
        n = len(part)
        rows.append({
            "cohort": cohort, "n": n,
            "n_na": int(beta.isna().sum()), "n_le_0.02": int((beta <= 0.02).sum()), "n_gt_0.02": int((beta > 0.02).sum()),
            "frac_na": float(beta.isna().mean()), "frac_le_0.02": float((beta <= 0.02).mean()),
            "frac_gt_0.02": float((beta > 0.02).mean()),
            "median": float(beta.median()) if beta.notna().any() else None,
            "n_ld_missing": missing.get(cohort),
        })
    pd.DataFrame(rows).to_csv(OUT / "beta_by_cohort.tsv", sep="\t", index=False)


def imvigor(pred: pd.DataFrame) -> None:
    part = pred.loc[pred["cohort"].eq("blca_iatlas_imvigor210_2017") & pred["standard_site"].eq("kidney")].copy()
    samples = pd.read_csv(ROOT / "results" / "stage5" / "aux" / "cohorts" / "blca_iatlas_imvigor210_2017" / "samples.tsv", sep="\t", dtype=str)
    keep = ["sample_id", "BIOPSY_SITE", "SAMPLE_TYPE", "METASTASIZED", "CANCER_TYPE", "CANCER_TYPE_DETAILED", "ONCOTREE_CODE", "raw_site"]
    joined = part[["sample_id", "raw_site", "standard_site"]].merge(samples[keep], on="sample_id", how="left", suffixes=("_pred", "_samples"))
    rows = []
    for column in ("raw_site_pred", "raw_site_samples", "BIOPSY_SITE", "SAMPLE_TYPE", "METASTASIZED", "CANCER_TYPE", "CANCER_TYPE_DETAILED", "ONCOTREE_CODE"):
        if column not in joined.columns:
            continue
        for value, n in joined[column].fillna("").value_counts().items():
            rows.append({"field": column, "value": value, "n": int(n)})
    pd.DataFrame(rows).to_csv(OUT / "imvigor_kidney_raw_site.tsv", sep="\t", index=False)
    meta = {
        "sample_columns": list(samples.columns),
        "n_kidney_in_predictions": int(len(part)),
        "n_samples_file": int(len(samples)),
        "fields_tabulated": ["raw_site", "BIOPSY_SITE", "SAMPLE_TYPE", "METASTASIZED", "CANCER_TYPE", "CANCER_TYPE_DETAILED", "ONCOTREE_CODE"],
    }
    (OUT / "imvigor_kidney_meta.json").write_text(json.dumps(meta, indent=2))


def three_cohorts(pred: pd.DataFrame) -> None:
    use = pred.loc[as_bool(pred, "include_confirm") & pred["layer"].eq("RNA-seq") & as_bool(pred, "selected_eval")]
    rows = []
    for method in ("BASE-Z", "SA-Z", "SA-G"):
        truth = use["organ"].to_numpy(dtype=object)
        got = use[f"{method}__pred"].to_numpy(dtype=object)
        host = host_error(use, method)
        risk = np.array([bool(native) and truth[i] not in native for i, native in enumerate(native_sets(use))])
        rows.append({
            "cohort": "aux_rnaseq_confirm", "method": method, "n": int(len(use)),
            "top1": float(np.mean(got == truth)), "n_at_risk": int(risk.sum()),
            "host_rate": float(host[risk].mean()) if risk.any() else None,
            "source": "results/stage7/confirm/predictions.parquet selected_eval",
        })
    met = pd.read_csv(ROOT / "results" / "stage5" / "posthoc" / "met500_overall.tsv", sep="\t")
    met_g = pd.read_csv(ROOT / "results" / "stage5" / "posthoc" / "met500_rules.tsv", sep="\t")
    pog = pd.read_csv(ROOT / "results" / "stage5" / "posthoc" / "pog_overall.tsv", sep="\t") if (ROOT / "results" / "stage5" / "posthoc" / "pog_overall.tsv").exists() else None
    pog_g = pd.read_csv(ROOT / "results" / "stage5" / "posthoc" / "pog_rules.tsv", sep="\t")

    def add(table, cohort, method, source):
        hit = table.loc[table["method"].eq(method)]
        if hit.empty:
            return
        row = hit.iloc[0]
        rows.append({"cohort": cohort, "method": method, "n": int(row["n"]), "top1": float(row["top1"]),
                     "n_at_risk": int(row["n_at_risk"]), "host_rate": float(row["host_rate"]), "source": source})

    for method in ("BASE-Z", "SA-Z"):
        add(met, "MET500", method, "results/stage5/posthoc/met500_overall.tsv")
    g = met_g.loc[met_g["rule"].eq("G-beta(0.02)")].iloc[0]
    rows.append({"cohort": "MET500", "method": "SA-G", "n": int(g["n"]), "top1": float(g["top1"]),
                 "n_at_risk": int(g["n_at_risk"]), "host_rate": float(g["host_rate"]),
                 "source": "results/stage5/posthoc/met500_rules.tsv G-beta(0.02)"})
    pog_path = ROOT / "results" / "stage4" / "confirm" / "tables"
    # POG overall lives next to the rules file.
    pog_file = ROOT / "results" / "stage5" / "posthoc" / "pog_rules.tsv"
    # BASE/SA from the comparison table's source primary metrics if a file exists.
    for folder in (
        ROOT / "results" / "stage5" / "posthoc" / "pog_overall.tsv",
        ROOT / "results" / "stage4" / "pog_overall.tsv",
    ):
        if folder.exists():
            pog = pd.read_csv(folder, sep="\t")
            for method in ("BASE-Z", "SA-Z"):
                add(pog, "POG570", method, str(folder.relative_to(ROOT)))
            break
    g = pog_g.loc[pog_g["rule"].eq("G-beta(0.02)")].iloc[0]
    rows.append({"cohort": "POG570", "method": "SA-G", "n": int(g["n"]), "top1": float(g["top1"]),
                 "n_at_risk": int(g["n_at_risk"]), "host_rate": float(g["host_rate"]),
                 "source": "results/stage5/posthoc/pog_rules.tsv G-beta(0.02)"})
    if not any(row["cohort"] == "POG570" and row["method"] == "BASE-Z" for row in rows):
        # Numbers recorded in A report 4 §4.7.
        rows.append({"cohort": "POG570", "method": "BASE-Z", "n": 512, "top1": 0.681640625, "n_at_risk": 378,
                     "host_rate": 0.2275132275132275, "source": "A_보고서_4.md §4.7"})
        rows.append({"cohort": "POG570", "method": "SA-Z", "n": 512, "top1": 0.755859375, "n_at_risk": 378,
                     "host_rate": 0.023809523809523808, "source": "A_보고서_4.md §4.7"})
    pd.DataFrame(rows).to_csv(OUT / "three_cohorts.tsv", sep="\t", index=False)


def read_cbio(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    symbol = "Hugo_Symbol" if "Hugo_Symbol" in frame.columns else frame.columns[0]
    frame = frame.loc[frame[symbol].notna() & ~frame[symbol].isin(["", "NA"])]
    cols = [col for col in frame.columns if col not in {symbol, "Entrez_Gene_Id"}]
    values = frame[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    values.index = frame[symbol].astype(str).to_numpy()
    return values.groupby(level=0).sum()


def gpl96(downloads: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(downloads / "GPL96.annot.gz", "rt", errors="replace") as handle:
        header = None
        for line in handle:
            if line.startswith("#") or line.startswith("!") or line.startswith("^"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            rec = dict(zip(header, parts))
            symbol = rec.get("Gene symbol", "")
            if not symbol or symbol == "---" or "///" in symbol:
                continue
            rows.append((rec["ID"], symbol))
    return pd.DataFrame(rows, columns=["probe", "symbol"]).drop_duplicates("probe")


def series_matrix(path: Path, annot: pd.DataFrame) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as handle:
        lines = handle.read().splitlines()
    start = lines.index("!series_matrix_table_begin") + 1
    end = lines.index("!series_matrix_table_end")
    header = [part.strip().strip('"') for part in lines[start].split("\t")]
    probes, data = [], []
    for line in lines[start + 1:end]:
        parts = [part.strip().strip('"') for part in line.split("\t")]
        probes.append(parts[0])
        data.append([float(x) if x not in {"", "null", "NA"} else np.nan for x in parts[1:]])
    values = pd.DataFrame(data, index=probes, columns=header[1:])
    merged = values.join(annot.set_index("probe"), how="inner")
    means = merged.drop(columns=["symbol"]).mean(axis=1)
    merged = merged.assign(_mean=means)
    chosen = merged.sort_values("_mean", ascending=False).groupby("symbol", sort=False).head(1)
    matrix = chosen.drop(columns=["_mean", "symbol"])
    matrix.index = chosen["symbol"].to_numpy()
    return matrix


def pack_offsets(member_lists):
    offsets = np.zeros(len(member_lists) + 1, dtype=np.int64)
    chunks = []
    for i, members in enumerate(member_lists):
        chunks.append(np.asarray(members, dtype=np.int32))
        offsets[i + 1] = offsets[i] + chunks[-1].size
    indices = np.concatenate(chunks) if chunks else np.empty(0, dtype=np.int32)
    return offsets, indices


def z_and_k(block, present, genes, b0, offsets, indices):
    present_index = {gene: i for i, gene in enumerate(present)}
    z_all = rank_z(block)
    z = np.zeros((block.shape[0], len(b0)), dtype=np.float32)
    for j, symbol in enumerate(b0):
        if symbol in present_index:
            z[:, j] = z_all[:, present_index[symbol]].astype(np.float32)
    lists, ok = [], []
    n_sets = int(offsets.shape[0] - 1)
    for g in range(n_sets):
        members = []
        for gene_i in indices[offsets[g]:offsets[g + 1]]:
            symbol = genes[int(gene_i)]
            if symbol in present_index:
                members.append(present_index[symbol])
        if len(members) < 5:
            ok.append(False)
            lists.append(np.array([0], dtype=np.int32))
        else:
            ok.append(True)
            lists.append(np.array(members, dtype=np.int32))
    valid = [members for members, flag in zip(lists, ok) if flag]
    packed_off, packed_idx = pack_offsets(valid)
    scored = score_ks(block.astype(np.float32), packed_off, packed_idx)
    k = np.zeros((block.shape[0], n_sets), dtype=np.float32)
    names = []
    cursor = 0
    for g, flag in enumerate(ok):
        names.append(g)
        if flag:
            k[:, g] = scored[:, cursor]
            cursor += 1
    return z, k


def align(matrix: pd.DataFrame, genes: list[str], sample_ids: list[str]):
    matrix = matrix.reindex(columns=sample_ids)
    present = [gene for gene in genes if gene in matrix.index]
    block = np.nan_to_num(matrix.loc[present].to_numpy(dtype=np.float64).T, nan=0.0)
    return block, present


def tie_extra(block: np.ndarray) -> int:
    """Genes that share a float32 value with another gene, summed over samples."""
    values = np.ascontiguousarray(block, dtype=np.float32)
    total = 0
    for row in values:
        _uniq, counts = np.unique(row, return_counts=True)
        total += int(counts[counts > 1].sum())
    return total


def feature_sensitivity(pred: pd.DataFrame, cfg: dict) -> None:
    genes, b0_idx, offsets, indices, set_names = load_gene_pack(cfg)
    b0 = [genes[int(i)] for i in b0_idx]
    downloads = cfg["data_processed"] / "aux" / "downloads"
    lengths = pd.read_csv(downloads / "gencode_v23_exon_union.tsv", sep="\t").set_index("symbol")["exon_union_bp"]
    jobs = {
        "blca_iatlas_imvigor210_2017": ("rna", read_cbio(downloads / "blca_iatlas_imvigor210_2017__data_mrna_seq_tpm.txt")),
        "brca_iatlas_anders_2022": ("rna", read_cbio(downloads / "brca_iatlas_anders_2022__data_mrna_seq_tpm.txt")),
        "paad_iatlas_prince_2022": ("paad", read_cbio(downloads / "paad_iatlas_prince_2022__data_mrna_seq_expression.txt")),
        "GSE14018": ("microarray", series_matrix(downloads / "GSE14018_series_matrix.txt.gz", gpl96(downloads))),
        "prad_fhcrc": ("microarray", read_cbio(downloads / "prad_fhcrc_data_mrna_agilent_microarray.txt")),
    }
    models = {
        "BASE-Z": joblib.load(cfg["results"] / "stage2" / "models" / "BASE_Z.joblib"),
        "SA-Z": joblib.load(cfg["results"] / "stage2" / "models" / "SA_Z.joblib"),
    }
    rows = []
    tie_rows = []
    set_rows = []
    for cohort, (kind, matrix) in jobs.items():
        ids = [line for line in (cfg["data_processed"] / "aux" / cohort / "sample_ids.txt").read_text().splitlines() if line]
        matrix = matrix.reindex(columns=[col for col in ids if col in matrix.columns])
        missing_ids = [sample for sample in ids if sample not in matrix.columns]
        if missing_ids:
            raise SystemExit(f"{cohort} missing {len(missing_ids)} sample columns")
        stage5_block, present5 = align(matrix, genes, ids)
        if kind == "paad":
            linear = np.maximum(np.exp2(matrix.to_numpy(dtype=np.float64)) - 1.0, 0.0)
            linear = pd.DataFrame(np.nan_to_num(linear, nan=0.0, posinf=0.0, neginf=0.0), index=matrix.index, columns=matrix.columns)
            common = linear.index.intersection(lengths.index)
            rate = linear.loc[common].div(lengths.loc[common], axis=0)
            new_block, present_new = align(rate, genes, ids)
        else:
            values = matrix.to_numpy(dtype=np.float64)
            linear = np.exp2(values) if kind == "microarray" else np.maximum(np.exp2(values) - 1.0, 0.0)
            linear = pd.DataFrame(np.nan_to_num(linear, nan=0.0, posinf=0.0, neginf=0.0), index=matrix.index, columns=matrix.columns)
            new_block, present_new = align(linear, genes, ids)
        stage5_sum = sum_1e6(stage5_block)
        new_sum = sum_1e6(new_block)
        z5, k5 = z_and_k(stage5_sum, present5, genes, b0, offsets, indices)
        z_new, k_new = z_and_k(new_sum, present_new, genes, b0, offsets, indices)
        disk_k = np.load(cfg["data_processed"] / "aux" / cohort / "K.npy")
        disk_z = np.load(cfg["data_processed"] / "aux" / cohort / "Z.npy")
        # predictions from the stage-5-style Z
        old_pred = {}
        for method, model in models.items():
            _p, labels, *_rest = predict_bundle(model, z5)
            old_pred[method] = labels
        sub = pred.loc[pred["cohort"].eq(cohort)].copy()
        have = [sample for sample in ids if sample in set(sub["sample_id"])]
        sub = sub.drop_duplicates("sample_id").set_index("sample_id").loc[have].reset_index()
        pos = {sample: i for i, sample in enumerate(ids)}
        take = np.array([pos[sample] for sample in sub["sample_id"]])
        for method in ("BASE-Z", "SA-Z"):
            disagree = int((old_pred[method][take] != sub[f"{method}__pred"].to_numpy(dtype=object)).sum())
            rows.append({
                "cohort": cohort, "method": method, "n": int(len(sub)), "top1_disagree_vs_confirm": disagree,
                "max_abs_z_recomputed_new_vs_disk": float(np.max(np.abs(z_new - disk_z))) if z_new.shape == disk_z.shape else None,
                "max_abs_k_stage5_vs_new": float(np.max(np.abs(k5 - k_new))) if k5.shape == k_new.shape else None,
                "max_abs_k_new_vs_disk": float(np.max(np.abs(k_new - disk_k))) if k_new.shape == disk_k.shape else None,
            })
        for label, getter in (
            ("confirm_predictions", lambda method: sub[f"{method}__pred"].to_numpy(dtype=object)),
            ("stage5_style_features", lambda method: old_pred[method][take]),
        ):
            work = sub.copy()
            for method in ("BASE-Z", "SA-Z"):
                work[f"{method}__pred"] = getter(method)
            risk = dedupe(work.loc[as_bool(work, "selected_risk")])
            ev = dedupe(work.loc[as_bool(work, "selected_eval")])
            if len(risk):
                b, s = host_error(risk, "BASE-Z"), host_error(risk, "SA-Z")
                rows.append({"cohort": cohort, "method": "AH1_" + label, "n": int(len(risk)),
                             "top1_disagree_vs_confirm": None, "n_favor": int((b & ~s).sum()), "n_against": int((s & ~b).sum()),
                             "diff_sa_minus_base": float(s.mean() - b.mean()),
                             "include_confirm": bool(as_bool(risk, "include_confirm").all())})
            if len(ev):
                s, b = correct(ev, "SA-Z"), correct(ev, "BASE-Z")
                rows.append({"cohort": cohort, "method": "AH2_" + label, "n": int(len(ev)),
                             "top1_disagree_vs_confirm": None, "n_favor": int((s & ~b).sum()), "n_against": int((b & ~s).sum()),
                             "diff_sa_minus_base": float(s.mean() - b.mean()),
                             "include_confirm": bool(as_bool(ev, "include_confirm").all())})
        if cohort == "blca_iatlas_imvigor210_2017" and k5.shape == k_new.shape:
            delta = np.abs(k5.astype(np.float64) - k_new.astype(np.float64))
            flat = np.unravel_index(np.argmax(delta), delta.shape)
            tie_rows.append({
                "stage5_sum_tie_gene_copies": tie_extra(stage5_sum),
                "linear_sum_tie_gene_copies": tie_extra(new_sum),
                "max_abs_k": float(delta[flat]),
                "max_sample_id": ids[int(flat[0])],
                "max_set_index": int(flat[1]),
                "max_set_name": set_names[int(flat[1])] if int(flat[1]) < len(set_names) else "",
                "n_samples": int(delta.shape[0]), "n_sets": int(delta.shape[1]),
            })
            set_max = delta.max(axis=0)
            order = np.argsort(-set_max)[:15]
            for j in order:
                set_rows.append({
                    "set_index": int(j), "set_name": set_names[int(j)] if int(j) < len(set_names) else "",
                    "max_abs_k": float(set_max[j]), "mean_abs_k": float(delta[:, j].mean()),
                    "max_sample_id": ids[int(np.argmax(delta[:, j]))],
                })
        feat_dir = OUT / "features" / cohort
        feat_dir.mkdir(parents=True, exist_ok=True)
        np.save(feat_dir / "Z_stage5_style.npy", z5)
        np.save(feat_dir / "K_stage5_style.npy", k5)
        print("features", cohort, flush=True)
    pd.DataFrame(rows).to_csv(OUT / "feature_sensitivity.tsv", sep="\t", index=False)
    pd.DataFrame(tie_rows).to_csv(OUT / "blca_k_ties.tsv", sep="\t", index=False)
    pd.DataFrame(set_rows).to_csv(OUT / "blca_k_sets.tsv", sep="\t", index=False)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pred = load_predictions()
    rng = np.random.default_rng(20261001)
    discordant(pred)
    site_groups(pred)
    gse_prefix(pred, rng)
    beta_table(pred)
    imvigor(pred)
    three_cohorts(pred)
    feature_sensitivity(pred, paths())
    print("posthoc", OUT, flush=True)


if __name__ == "__main__":
    main()
