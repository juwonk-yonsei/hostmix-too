"""Post-hoc IMvigor210 contributions and TCGA ESCA histology. Does not retrain."""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

def _release_config():
    cfg = Path(__file__).resolve().parents[1] / "config.yaml"
    out = {}
    for line in cfg.read_text().splitlines():
        if ":" not in line or line.strip().startswith("#"):
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        if key in {"project_root", "figure_dir", "disk_mount"}:
            out[key] = value.strip().strip('"').strip("'")
    return out


def project_root():
    return Path(_release_config()["project_root"])


import h5py
import joblib
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from scipy.stats import norm

ROOT = project_root()
OUT = ROOT / "results/stage9"
NOTE = "사후"


def rank_b0(values: np.ndarray, b0_idx: np.ndarray) -> np.ndarray:
    out = np.empty((values.shape[0], len(b0_idx)), dtype=np.float32)
    n_genes = values.shape[1]
    for start in range(0, values.shape[0], 200):
        block = values[start:start + 200]
        ranks = rankdata(block, method="average", axis=1)
        scored = norm.ppf((ranks - 0.5) / n_genes).astype(np.float32)
        out[start:start + block.shape[0]] = scored[:, b0_idx]
    return out


def sum_1e6(tpm: np.ndarray) -> np.ndarray:
    totals = tpm.sum(axis=1, keepdims=True).astype(np.float64)
    totals[totals == 0] = np.nan
    return np.nan_to_num(tpm.astype(np.float64) / totals * 1e6, nan=0.0).astype(np.float32)


def class_index(bundle, label: str) -> int:
    classes = [str(item) for item in bundle["clf"].classes_]
    if label not in classes:
        raise SystemExit(f"missing class {label}")
    return classes.index(label)


def histology_bin(text: str) -> str:
    low = str(text).lower()
    if "adenosquamous" in low:
        return "other"
    if "squamous" in low and "adenocarcinoma" in low:
        return "other"
    if "squamous" in low:
        return "squamous"
    if "adenocarcinoma" in low:
        return "adenocarcinoma"
    return "other"


def gdc_diagnoses(patients: list[str]) -> pd.DataFrame:
    payload = {
        "filters": {"op": "in", "content": {"field": "cases.submitter_id", "value": patients}},
        "fields": "submitter_id,diagnoses.primary_diagnosis",
        "size": max(len(patients), 1),
        "format": "JSON",
    }
    request = urllib.request.Request(
        "https://api.gdc.cancer.gov/cases",
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(request, timeout=120) as handle:
        body = json.loads(handle.read().decode())
    rows = []
    for hit in body["data"]["hits"]:
        diagnoses = hit.get("diagnoses") or []
        text = diagnoses[0].get("primary_diagnosis", "") if diagnoses else ""
        rows.append({"patient": hit["submitter_id"], "primary_diagnosis": text, "histology": histology_bin(text)})
    return pd.DataFrame(rows)


def read_tpm(sample_ids: list[str]) -> np.ndarray:
    with h5py.File(ROOT / "data/processed/toil/tpm_G.h5", "r") as handle:
        samples = [item.decode() if isinstance(item, bytes) else str(item) for item in handle["samples"][:]]
        index = {sample: i for i, sample in enumerate(samples)}
        missing = [sample for sample in sample_ids if sample not in index]
        if missing:
            raise SystemExit(f"missing tpm rows {missing[:3]}")
        order = np.argsort([index[sample] for sample in sample_ids])
        sorted_ids = [sample_ids[i] for i in order]
        rows = handle["tpm"][[index[sample] for sample in sorted_ids]]
    back = np.empty(len(sample_ids), dtype=int)
    back[order] = np.arange(len(sample_ids))
    return sum_1e6(rows[back])


def cohort_z(cohort: str, sample_id: str, cache: dict) -> np.ndarray | None:
    if cohort not in cache:
        paths = [ROOT / "data/processed/aux" / cohort]
        if cohort == "prad_su2c_2019":
            paths = [
                ROOT / "data/processed/aux/prad_su2c_2019__fpkm_polya",
                ROOT / "data/processed/aux/prad_su2c_2019__fpkm_capture",
            ]
        loaded = []
        for path in paths:
            ids = (path / "sample_ids.txt").read_text().splitlines()
            z = np.load(path / "Z.npy")
            loaded.append(( {sample: i for i, sample in enumerate(ids)}, z))
        cache[cohort] = loaded
    for index, z in cache[cohort]:
        if sample_id in index:
            return z[index[sample_id]]
    return None


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = [gene for gene in (ROOT / "results/B0_genes.txt").read_text().splitlines() if gene and gene != "symbol"]
    g_symbols = pd.read_csv(ROOT / "data/processed/genes/G_symbols.txt", header=None)[0].astype(str).tolist()
    b0_idx = np.array([g_symbols.index(gene) for gene in genes], dtype=np.int32)
    lengths = pd.read_csv(ROOT / "data/processed/aux/downloads/gencode_v23_exon_union.tsv", sep="\t", dtype=str)
    length_kb = {}
    for symbol, bp in zip(lengths["symbol"], lengths["exon_union_bp"]):
        value = float(bp)
        if value > 0 and symbol not in length_kb:
            length_kb[symbol] = value / 1000.0
    b0_lengths = np.array([length_kb[gene] for gene in genes if gene in length_kb], dtype=float)

    confirm = pd.read_parquet(
        ROOT / "results/stage7/confirm/predictions.parquet",
        columns=["cohort", "sample_id", "organ", "selected_native", "selected_eval", "BASE-Z__pred", "SA-Z__pred", "LD-Z__pred"],
    )
    for column in ("selected_native", "selected_eval"):
        confirm[column] = confirm[column].astype(str).str.lower().isin(["true", "1"])
    native = confirm.loc[confirm["cohort"].eq("blca_iatlas_imvigor210_2017") & confirm["selected_native"]].copy()
    if len(native) != 194:
        raise SystemExit(f"native {len(native)}")
    ids = (ROOT / "data/processed/aux/blca_iatlas_imvigor210_2017/sample_ids.txt").read_text().splitlines()
    z = np.load(ROOT / "data/processed/aux/blca_iatlas_imvigor210_2017/Z.npy")
    index = {sample: i for i, sample in enumerate(ids)}
    native["row"] = native["sample_id"].astype(str).map(index)
    if native["row"].isna().any():
        raise SystemExit("IMvigor Z id missing")
    split = pd.read_csv(ROOT / "config/split_tcga.tsv", sep="\t", dtype=str)
    blca = split.loc[split["split"].eq("TCGA-test") & split["learning_label"].eq("BLCA"), "sample"].astype(str).tolist()
    print("tcga blca z", len(blca), flush=True)
    blca_z = rank_b0(read_tpm(blca), b0_idx)
    contrib_rows = []
    compare_rows = []
    for method, path, n_expected in (("BASE-Z", "BASE_Z.joblib", 145), ("SA-Z", "SA_Z.joblib", 65)):
        bundle = joblib.load(ROOT / "results/stage2/models" / path)
        coef = bundle["clf"].coef_
        esca = class_index(bundle, "ESCA")
        blca_i = class_index(bundle, "BLCA")
        called = native.loc[native[f"{method}__pred"].eq("Esophagus")]
        if len(called) != n_expected:
            raise SystemExit(f"{method} esophagus {len(called)}")
        block = z[[int(i) for i in called["row"]]]
        delta = coef[esca] - coef[blca_i]
        mean_contrib = (delta * block).mean(axis=0)
        order = np.argsort(-np.abs(mean_contrib))[:20]
        top_genes = [genes[i] for i in order]
        for rank, gene_i in enumerate(order, start=1):
            contrib_rows.append({
                "note": NOTE, "method": method, "n_samples": int(len(called)), "rank": rank,
                "gene": genes[gene_i], "mean_contribution": float(mean_contrib[gene_i]),
                "length_kb": length_kb.get(genes[gene_i]),
                "length_definition": "GENCODE v23 exon-union kb, first positive row per symbol",
            })
        top_len = np.array([length_kb[gene] for gene in top_genes if gene in length_kb], dtype=float)
        compare_rows.append({
            "note": NOTE, "method": method, "n_top": int(len(top_genes)),
            "median_length_kb_top20": float(np.median(top_len)),
            "n_b0_with_length": int(len(b0_lengths)),
            "median_length_kb_B0": float(np.median(b0_lengths)),
        })
        for gene in top_genes:
            gene_i = genes.index(gene)
            compare_rows_gene = {
                "note": NOTE, "method": method, "gene": gene,
                "imvigor_n": int(len(called)), "imvigor_mean_Z": float(block[:, gene_i].mean()),
                "imvigor_median_Z": float(np.median(block[:, gene_i])),
                "tcga_test_blca_n": int(len(blca)), "tcga_test_blca_mean_Z": float(blca_z[:, gene_i].mean()),
                "tcga_test_blca_median_Z": float(np.median(blca_z[:, gene_i])),
                "tcga_Z_denominator": "rank-normal on all G genes, then B0 columns",
            }
            contrib_rows.append({"note": NOTE, "method": method, "rank": "", "gene": gene, "z_compare": True, **{k: v for k, v in compare_rows_gene.items() if k not in {"note", "method", "gene"}}})
    pd.DataFrame([row for row in contrib_rows if row.get("rank") != ""]).to_csv(OUT / "imvigor_contribution.tsv", sep="\t", index=False)
    pd.DataFrame([row for row in contrib_rows if row.get("z_compare")]).to_csv(OUT / "imvigor_gene_Z.tsv", sep="\t", index=False)
    pd.DataFrame(compare_rows).to_csv(OUT / "imvigor_length.tsv", sep="\t", index=False)

    train = split.loc[split["split"].eq("TCGA-train") & split["learning_label"].eq("ESCA") & ~split["exclude_from_eval"].astype(str).str.lower().eq("true")].copy()
    if len(train) != 145:
        raise SystemExit(f"ESCA train {len(train)}")
    train["patient"] = train["sample"].astype(str).str.slice(0, 12)
    print("gdc", train["patient"].nunique(), flush=True)
    clinical = gdc_diagnoses(sorted(train["patient"].unique()))
    clinical.to_csv(OUT / "esca_gdc_diagnosis.tsv", sep="\t", index=False)
    merged = train.merge(clinical, on="patient", how="left")
    if merged["histology"].isna().any():
        raise SystemExit(f"unmatched ESCA {int(merged['histology'].isna().sum())}")
    counts = merged["histology"].value_counts()
    raw_counts = merged["primary_diagnosis"].value_counts()
    pd.DataFrame({
        "note": NOTE, "histology": counts.index.astype(str), "n": counts.to_numpy(), "n_train": len(merged),
    }).to_csv(OUT / "esca_histology_counts.tsv", sep="\t", index=False)
    pd.DataFrame({
        "note": NOTE, "primary_diagnosis": raw_counts.index.astype(str), "n": raw_counts.to_numpy(),
    }).to_csv(OUT / "esca_diagnosis_strings.tsv", sep="\t", index=False)
    print("esca z", flush=True)
    esca_z = rank_b0(read_tpm(merged["sample"].astype(str).tolist()), b0_idx)
    centers = {}
    for label in ("adenocarcinoma", "squamous"):
        mask = merged["histology"].eq(label).to_numpy()
        if mask.sum() == 0:
            raise SystemExit(f"no {label}")
        centers[label] = esca_z[mask].mean(axis=0)
    use = confirm.loc[confirm["selected_eval"]].copy()
    cache: dict = {}
    rows = []
    for method in ("BASE-Z", "SA-Z", "LD-Z"):
        called = use.loc[use[f"{method}__pred"].eq("Esophagus")]
        for row in called.itertuples(index=False):
            vector = cohort_z(row.cohort, str(row.sample_id), cache)
            if vector is None:
                rows.append({"note": NOTE, "method": method, "cohort": row.cohort, "sample_id": row.sample_id, "status": "Z missing"})
                continue
            scores = {}
            for label, center in centers.items():
                coef, _p = spearmanr(vector, center)
                scores[label] = float(coef)
            if scores["adenocarcinoma"] > scores["squamous"]:
                winner = "adenocarcinoma"
            elif scores["squamous"] > scores["adenocarcinoma"]:
                winner = "squamous"
            else:
                winner = "tie"
            rows.append({
                "note": NOTE, "method": method, "cohort": row.cohort, "sample_id": row.sample_id,
                "status": "scored", "closer": winner,
                "spearman_adenocarcinoma": scores["adenocarcinoma"], "spearman_squamous": scores["squamous"],
            })
    detail = pd.DataFrame(rows)
    detail.to_csv(OUT / "esophagus_histology_similarity.tsv", sep="\t", index=False)
    summary = []
    scored = detail.loc[detail["status"].eq("scored")]
    for (method, cohort), sub in scored.groupby(["method", "cohort"]):
        summary.append({
            "note": NOTE, "method": method, "cohort": cohort, "n": int(len(sub)),
            "n_closer_adenocarcinoma": int(sub["closer"].eq("adenocarcinoma").sum()),
            "n_closer_squamous": int(sub["closer"].eq("squamous").sum()),
            "n_tie": int(sub["closer"].eq("tie").sum()),
        })
    pd.DataFrame(summary).to_csv(OUT / "esophagus_histology_summary.tsv", sep="\t", index=False)
    print("extra diagnostics", OUT, counts.to_dict(), flush=True)


if __name__ == "__main__":
    main()
