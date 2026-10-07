"""Post-hoc platform-gene mask on POG570 and MET500. Does not retrain."""
from __future__ import annotations

import gzip
import sys
import tarfile
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


import joblib
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from scipy.stats import norm

sys.path.insert(0, str(project_root() / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ks_fast import ks_matrix, right_and_tie  # noqa: E402
from mapping_rules import SITE_NATIVE  # noqa: E402
from stage2_fit import predict_bundle  # noqa: E402

ROOT = project_root()
OUT = ROOT / "results/stage9"
DOWNLOADS = ROOT / "data/processed/aux/downloads"
NOTE = "사후"
METHODS = ["BASE-Z", "SA-Z", "BASE-K", "SA-K"]
PLATFORM_ORDER = ["GPL96", "GPL20769", "GPL15659", "prad_fhcrc_agilent"]
EXPECTED_G = {"GPL96": 11389, "GPL20769": 10693, "GPL15659": 15156, "prad_fhcrc_agilent": 16323}


def native_list(site: str) -> list[str]:
    return list(SITE_NATIVE.get(str(site), ([], [], []))[2])


def single_symbol(text: str) -> str:
    symbol = str(text).strip().strip('"')
    if symbol in {"", "---", "NA", "null", "None", "nan"}:
        return ""
    if "///" in symbol or ";" in symbol or "|" in symbol or " " in symbol:
        return ""
    return symbol


def rank_z(values: np.ndarray) -> np.ndarray:
    out = np.empty(values.shape, dtype=np.float32)
    n_genes = values.shape[1]
    for start in range(0, values.shape[0], 200):
        block = values[start:start + 200]
        ranks = rankdata(block, method="average", axis=1)
        out[start:start + block.shape[0]] = norm.ppf((ranks - 0.5) / n_genes).astype(np.float32)
    return out


def score_ks(tpm: np.ndarray, offsets, indices) -> np.ndarray:
    parts = []
    for start in range(0, tpm.shape[0], 2000):
        block = np.ascontiguousarray(tpm[start:start + 2000], dtype=np.float32)
        right, tie = right_and_tie(block)
        parts.append(ks_matrix(right, tie, offsets, indices))
    return np.vstack(parts) if parts else np.empty((0, int(offsets.shape[0] - 1)), dtype=np.float32)


def sum_1e6(tpm: np.ndarray) -> np.ndarray:
    totals = tpm.sum(axis=1, keepdims=True).astype(np.float64)
    totals[totals == 0] = np.nan
    return np.nan_to_num(tpm.astype(np.float64) / totals * 1e6, nan=0.0).astype(np.float32)


def gpl96_symbols() -> set[str]:
    probes = set()
    with gzip.open(DOWNLOADS / "GSE41258_series_matrix.txt.gz", "rt", errors="replace") as handle:
        started = False
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                started = True
                next(handle)
                continue
            if not started:
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            probes.add(line.split("\t", 1)[0].strip().strip('"'))
    symbols = set()
    with gzip.open(DOWNLOADS / "GPL96.annot.gz", "rt", errors="replace") as handle:
        header = None
        for line in handle:
            if line.startswith("#") or line.startswith("!") or line.startswith("^"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            rec = dict(zip(header, parts))
            if rec.get("ID") not in probes:
                continue
            symbol = single_symbol(rec.get("Gene symbol", ""))
            if symbol:
                symbols.add(symbol)
    return symbols


def first_tar_symbols(tar_path: Path, kind: str) -> set[str]:
    with tarfile.open(tar_path) as tar:
        for member in tar:
            name = Path(member.name).name
            if not name.startswith("GSM") or not name.endswith(".txt.gz"):
                continue
            raw = gzip.decompress(tar.extractfile(member).read())
            lines = raw.decode("latin1").splitlines()
            symbols = set()
            if kind == "genepix":
                header = lines[0].split("\t")
                i_sym = header.index("Gene Symbol")
                for line in lines[1:]:
                    parts = line.split("\t")
                    if i_sym < len(parts):
                        symbol = single_symbol(parts[i_sym])
                        if symbol:
                            symbols.add(symbol)
            else:
                symbol_of = {}
                with open(DOWNLOADS / "GPL15659_data.txt", encoding="latin1") as handle:
                    header = None
                    for line in handle:
                        if line.startswith("#") or line.startswith("!") or line.startswith("^"):
                            continue
                        parts = line.rstrip("\n").split("\t")
                        if header is None:
                            header = parts
                            continue
                        rec = dict(zip(header, parts))
                        symbol = single_symbol(rec.get("GENE_SYMBOL", rec.get("Gene symbol", "")))
                        if symbol and rec.get("ID"):
                            symbol_of[rec["ID"]] = symbol
                start = next(i for i, line in enumerate(lines) if line.startswith("FEATURES\t"))
                header = lines[start].split("\t")
                i_probe = header.index("ProbeName")
                i_ctrl = header.index("ControlType")
                for line in lines[start + 1:]:
                    if not line.startswith("DATA\t"):
                        continue
                    parts = line.split("\t")
                    if i_ctrl < len(parts) and parts[i_ctrl] not in {"0", "FALSE"}:
                        continue
                    symbol = symbol_of.get(parts[i_probe], "")
                    if symbol:
                        symbols.add(symbol)
            return symbols
    raise SystemExit(f"no sample in {tar_path}")


def fhcrc_symbols() -> set[str]:
    frame = pd.read_csv(DOWNLOADS / "prad_fhcrc_data_mrna_agilent_microarray.txt", sep="\t", usecols=["Hugo_Symbol"], dtype=str)
    return {single_symbol(value) for value in frame["Hugo_Symbol"] if single_symbol(value)}


def platform_symbols(genes: list[str]) -> dict[str, set[str]]:
    raw = {
        "GPL96": gpl96_symbols(),
        "GPL20769": first_tar_symbols(DOWNLOADS / "GSE71729_RAW.tar", "genepix"),
        "GPL15659": first_tar_symbols(DOWNLOADS / "GSE74685_RAW.tar", "agilent"),
        "prad_fhcrc_agilent": fhcrc_symbols(),
    }
    gene_set = set(genes)
    out = {}
    for name, symbols in raw.items():
        present = symbols & gene_set
        print("platform", name, "symbols", len(symbols), "G", len(present), "expected", EXPECTED_G[name], flush=True)
        if len(present) != EXPECTED_G[name]:
            raise SystemExit(f"{name} G intersect {len(present)} != {EXPECTED_G[name]}")
        out[name] = present
    return out


def features(tpm: np.ndarray, genes: list[str], keep: set[str] | None, b0: list[str], offsets, indices):
    if keep is None:
        block = tpm
        present = genes
    else:
        columns = [i for i, gene in enumerate(genes) if gene in keep]
        block = np.ascontiguousarray(tpm[:, columns])
        present = [genes[i] for i in columns]
    z_all = rank_z(block)
    present_index = {gene: i for i, gene in enumerate(present)}
    z = np.zeros((tpm.shape[0], len(b0)), dtype=np.float32)
    n_b0 = 0
    for j, symbol in enumerate(b0):
        hit = present_index.get(symbol)
        if hit is not None:
            z[:, j] = z_all[:, hit]
            n_b0 += 1
    lists = []
    set_ok = []
    n_sets = int(offsets.shape[0] - 1)
    for g in range(n_sets):
        members = [present_index[genes[int(gene_i)]] for gene_i in indices[offsets[g]:offsets[g + 1]] if genes[int(gene_i)] in present_index]
        if len(members) < 5:
            set_ok.append(False)
        else:
            set_ok.append(True)
            lists.append(np.asarray(members, dtype=np.int32))
    valid = lists
    packed_off = np.zeros(len(valid) + 1, dtype=np.int64)
    if valid:
        packed_idx = np.concatenate(valid)
        for i, members in enumerate(valid):
            packed_off[i + 1] = packed_off[i] + members.size
    else:
        packed_idx = np.empty(0, dtype=np.int32)
    scored = score_ks(block, packed_off, packed_idx) if valid else np.zeros((tpm.shape[0], 0), dtype=np.float32)
    k = np.zeros((tpm.shape[0], n_sets), dtype=np.float32)
    cursor = 0
    for g, ok in enumerate(set_ok):
        if ok:
            k[:, g] = scored[:, cursor]
            cursor += 1
    return z, k, {"n_G_kept": len(present), "n_B0_kept": n_b0, "n_ks_zero": int(n_sets - cursor)}


def load_models():
    models = {}
    for name in ("BASE_Z", "SA_Z", "BASE_K", "SA_K"):
        models[name] = joblib.load(ROOT / "results/stage2/models" / f"{name}.joblib")
    return models


def predict(models, z, k) -> dict[str, np.ndarray]:
    out = {}
    for method, bundle, matrix in (
        ("BASE-Z", models["BASE_Z"], z),
        ("SA-Z", models["SA_Z"], z),
        ("BASE-K", models["BASE_K"], k),
        ("SA-K", models["SA_K"], k),
    ):
        _organ, pred, _raw, _proba, _classes = predict_bundle(bundle, matrix)
        out[method] = np.asarray(pred, dtype=object)
    return out


def ensembl_map() -> dict[str, str]:
    hgnc = pd.read_csv(ROOT / "data/raw/hgnc/hgnc_complete_set.txt", sep="\t", dtype=str, usecols=["symbol", "locus_group"])
    coding = set(hgnc.loc[hgnc["locus_group"].eq("protein-coding gene"), "symbol"].dropna().astype(str))
    probe = pd.read_csv(ROOT / "data/raw/toil/gencode.v23.annotation.gene.probemap", sep="\t", dtype=str, usecols=["id", "gene"])
    probe = probe.loc[probe["gene"].isin(coding), ["id", "gene"]].drop_duplicates()
    by_base: dict[str, set[str]] = {}
    for ens, sym in zip(probe["id"].astype(str), probe["gene"].astype(str)):
        by_base.setdefault(ens.split(".")[0], set()).add(sym)
    return {base: next(iter(symbols)) for base, symbols in by_base.items() if len(symbols) == 1}


def load_eval(genes: list[str]):
    gene_index = {gene: i for i, gene in enumerate(genes)}
    met = pd.read_csv(ROOT / "results/stage5/posthoc/met500_samples.tsv", sep="\t", dtype=str)
    import h5py
    with h5py.File(ROOT / "data/processed/met500/tpm_G.h5", "r") as handle:
        samples = [item.decode() if isinstance(item, bytes) else str(item) for item in handle["samples"][:]]
        index = {sample: i for i, sample in enumerate(samples)}
        order = np.argsort([index[sample] for sample in met["id"]])
        sorted_ids = [met["id"].iloc[i] for i in order]
        rows = handle["tpm"][[index[sample] for sample in sorted_ids]]
    back = np.empty(len(met), dtype=int)
    back[order] = np.arange(len(met))
    met_tpm = sum_1e6(rows[back])
    met_meta = met[["id", "truth", "site", "pred__BASE-Z", "pred__SA-Z"]].copy()
    met_meta["cohort"] = "MET500"
    print("load POG", flush=True)
    expr = pd.read_csv(ROOT / "data/raw/POG570/POG570_TPM_expression.txt.gz", sep="\t", index_col=0)
    mapping = ensembl_map()
    bases = pd.Index(expr.index.astype(str)).str.split(".").str[0]
    symbol = np.asarray(bases.map(mapping), dtype=object)
    keep = pd.notna(symbol)
    grouped = expr.loc[keep].astype(np.float32).groupby(symbol[keep], sort=False).sum()
    aligned = grouped.reindex(genes).fillna(0).astype(np.float32)
    pog_all = sum_1e6(np.ascontiguousarray(aligned.to_numpy(dtype=np.float32).T))
    patient_ids = [str(column) for column in expr.columns]
    del expr, grouped, aligned
    labels = pd.read_csv(ROOT / "config/pog570_eval_labels.tsv", sep="\t", dtype=str)
    labels = labels.loc[~labels["organ"].isin(["exclude", "NA"])].copy()
    demo = pd.read_excel(ROOT / "data/raw/POG570/Table_S1_Demographics.xlsx", dtype=str)
    demo["standard_site"] = demo["BIOPSY_SITE"]
    from mapping_rules import map_site
    demo["standard_site"] = [map_site("pog570_biopsy_site", value) for value in demo["BIOPSY_SITE"]]
    merged = labels.merge(demo[["PATIENT_ID", "standard_site"]], on="PATIENT_ID", how="left")
    parquet = pd.read_parquet(ROOT / "results/stage3/pog570_predictions.parquet", columns=["patient_id", "BASE-Z__pred", "SA-Z__pred", "BASE-K__pred", "SA-K__pred"])
    merged = merged.merge(parquet, left_on="PATIENT_ID", right_on="patient_id", how="left")
    position = {patient: i for i, patient in enumerate(patient_ids)}
    take = [position[patient] for patient in merged["PATIENT_ID"]]
    pog_tpm = pog_all[take]
    pog_meta = pd.DataFrame({
        "id": merged["PATIENT_ID"].astype(str),
        "truth": merged["organ"].astype(str),
        "site": merged["standard_site"].astype(str),
        "pred__BASE-Z": merged["BASE-Z__pred"].astype(str),
        "pred__SA-Z": merged["SA-Z__pred"].astype(str),
        "pred__BASE-K": merged["BASE-K__pred"].astype(str),
        "pred__SA-K": merged["SA-K__pred"].astype(str),
        "cohort": "POG570",
    })
    if len(pog_meta) != 512 or len(met_meta) != 437:
        raise SystemExit(f"eval n {len(met_meta)} {len(pog_meta)}")
    print("G filled check", int((pog_tpm == 0).all(axis=0).sum()), flush=True)
    return {"MET500": (met_tpm, met_meta), "POG570": (pog_tpm, pog_meta)}, gene_index


def metric_rows(meta: pd.DataFrame, preds: dict[str, np.ndarray], original: dict[str, np.ndarray], cohort: str, mask: str, rep: str) -> list[dict]:
    rows = []
    truth = meta["truth"].astype(str).to_numpy()
    sites = meta["site"].astype(str).tolist()
    risk = np.array([bool(native_list(site)) and truth_i not in native_list(site) for truth_i, site in zip(truth, sites)])
    for method in METHODS:
        call = preds[method]
        base = original[method]
        correct = call == truth
        base_correct = base == truth
        host = np.array([
            bool(native_list(site)) and call_i in native_list(site) and call_i != truth_i
            for call_i, truth_i, site in zip(call, truth, sites)
        ])
        base_host = np.array([
            bool(native_list(site)) and call_i in native_list(site) and call_i != truth_i
            for call_i, truth_i, site in zip(base, truth, sites)
        ])
        moved = (base != "Esophagus") & (call == "Esophagus")
        rows.append({
            "note": NOTE, "cohort": cohort, "mask": mask, "rep": rep, "method": method,
            "n": int(len(meta)), "top1": float(correct.mean()),
            "top1_change": float(correct.mean() - base_correct.mean()),
            "n_at_risk": int(risk.sum()),
            "host_rate": float(host[risk].mean()) if risk.any() else float("nan"),
            "host_rate_change": float(host[risk].mean() - base_host[risk].mean()) if risk.any() else float("nan"),
            "n_moved_to_esophagus": int(moved.sum()),
        })
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(ROOT / "data/processed/genes/G_symbols.txt", header=None)[0].astype(str).tolist()
    b0 = pd.read_csv(ROOT / "results/B0_genes.txt")["symbol"].astype(str).tolist()
    packed = np.load(ROOT / "data/processed/genes/gene_sets.npz", allow_pickle=True)
    offsets, indices = packed["offsets"], packed["indices"]
    platforms = platform_symbols(genes)
    expr, _gene_index = load_eval(genes)
    models = load_models()
    stored = {}
    original = {}
    for cohort, (tpm, meta) in expr.items():
        print("unmasked", cohort, flush=True)
        z, k, info = features(tpm, genes, None, b0, offsets, indices)
        print("unmasked info", cohort, info, flush=True)
        preds = predict(models, z, k)
        baseline = {method: preds[method].copy() for method in METHODS}
        checks = (("BASE-Z", "pred__BASE-Z"), ("SA-Z", "pred__SA-Z"))
        if cohort == "POG570":
            checks = checks + (("BASE-K", "pred__BASE-K"), ("SA-K", "pred__SA-K"))
        for method, column in checks:
            mismatch = int((preds[method] != meta[column].astype(str).to_numpy()).sum())
            print("match", cohort, method, mismatch, flush=True)
            if method.endswith("-Z") and mismatch:
                raise SystemExit(f"unmasked {cohort} {method} mismatch {mismatch}")
            if column in meta.columns:
                baseline[method] = meta[column].astype(str).to_numpy()
        original[cohort] = baseline
        stored[cohort] = meta
    rows = []
    info_rows = []
    for cohort, (tpm, meta) in expr.items():
        rows.extend(metric_rows(meta, original[cohort], original[cohort], cohort, "unmasked", "0"))
    for name in PLATFORM_ORDER:
        keep = platforms[name]
        for cohort, (tpm, meta) in expr.items():
            print("mask", name, cohort, flush=True)
            z, k, info = features(tpm, genes, keep, b0, offsets, indices)
            preds = predict(models, z, k)
            rows.extend(metric_rows(meta, preds, original[cohort], cohort, name, "platform"))
            info_rows.append({"note": NOTE, "mask": name, "cohort": cohort, "rep": "platform", **info})
    b0_index = np.arange(len(b0))
    rng = np.random.default_rng(20261023)
    for name in PLATFORM_ORDER:
        n_drop = len([gene for gene in b0 if gene not in platforms[name]])
        for rep in range(10):
            chosen = set(b0[i] for i in rng.choice(b0_index, size=n_drop, replace=False))
            keep = set(genes) - chosen
            for cohort, (tpm, meta) in expr.items():
                print("random", name, rep, cohort, "drop", n_drop, flush=True)
                z, k, info = features(tpm, genes, keep, b0, offsets, indices)
                preds = predict(models, z, k)
                rows.extend(metric_rows(meta, preds, original[cohort], cohort, name, str(rep)))
                info_rows.append({"note": NOTE, "mask": name, "cohort": cohort, "rep": str(rep), "n_B0_dropped": n_drop, **info})
    pd.DataFrame(rows).to_csv(OUT / "mask_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(info_rows).to_csv(OUT / "mask_feature_info.tsv", sep="\t", index=False)
    summary = []
    frame = pd.DataFrame(rows)
    random_rows = frame.loc[frame["rep"].isin([str(i) for i in range(10)])]
    for (cohort, mask, method), sub in random_rows.groupby(["cohort", "mask", "method"]):
        summary.append({
            "note": NOTE, "cohort": cohort, "mask": mask, "method": method,
            "top1_mean": float(sub["top1"].mean()), "top1_min": float(sub["top1"].min()), "top1_max": float(sub["top1"].max()),
            "top1_change_mean": float(sub["top1_change"].mean()),
            "host_rate_mean": float(sub["host_rate"].mean()), "host_rate_min": float(sub["host_rate"].min()), "host_rate_max": float(sub["host_rate"].max()),
            "n_moved_to_esophagus_mean": float(sub["n_moved_to_esophagus"].mean()),
            "n_moved_to_esophagus_min": int(sub["n_moved_to_esophagus"].min()),
            "n_moved_to_esophagus_max": int(sub["n_moved_to_esophagus"].max()),
        })
    pd.DataFrame(summary).to_csv(OUT / "mask_random_summary.tsv", sep="\t", index=False)
    print("mask tables", OUT, flush=True)


if __name__ == "__main__":
    main()
