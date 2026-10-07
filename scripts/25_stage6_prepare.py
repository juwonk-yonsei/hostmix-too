"""Stage-6 evaluation-set and LD-input preparation. Does not fit or apply a classifier and does not compute LD beta."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import binom

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from stage2_common import load_gene_pack, sum_1e6
from stage5_lib import load_references


def load_paths():
    """Project paths. Does not read n_threads."""
    root = ROOT.parent
    raw = yaml.safe_load((root / "config" / "paths.yaml").read_text())
    out = {"root": root}
    for key, value in raw.items():
        if key in ("seed", "n_threads"):
            continue
        out[key] = root / value
    return out


def load_mod(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


geo = load_mod("auxgeo", "20_aux_geo.py")
aux = geo.aux
two = load_mod("auxtwo", "21_aux_twochannel.py")

POOL_OUT = {
    "kidney", "pancreas", "bladder", "thyroid", "ovary", "breast", "stomach",
    "colon_rectum", "prostate", "esophagus", "head_neck", "cervix",
}
PREFIX = r"^AUR-([A-Z0-9]{4})-"
PAPER_PATIENTS = 55
PAPER_SOURCE = "Garcia-Recio et al., Nat Cancer 2023, 10.1038/s43018-022-00491-x, 55 individuals"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def quantiles(matrix: pd.DataFrame) -> dict:
    values = matrix.to_numpy(dtype=np.float64).ravel()
    values = values[np.isfinite(values)]
    q0, q25, q99, q100 = np.quantile(values, [0, 0.25, 0.99, 1])
    linear = bool(q99 > 100 or ((q100 - q0) > 50 and q25 > 0))
    return {
        "n_values": int(values.size),
        "q0": float(q0), "q25": float(q25), "q99": float(q99), "q100": float(q100),
        "judgment": "선형" if linear else "로그",
    }


def length_then_sum(matrix: pd.DataFrame, lengths: pd.Series, genes: list[str]) -> tuple[np.ndarray, list[str], list[str], int]:
    common = matrix.index.intersection(lengths.index)
    missing = int(len(matrix.index.difference(lengths.index)))
    rate = matrix.loc[common].div(lengths.loc[common], axis=0)
    block, present, ids = aux.align_g(rate, genes)
    return sum_1e6(block), present, ids, missing


def sum_genes(matrix: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray, list[str], list[str]]:
    block, present, ids = aux.align_g(matrix, genes)
    return sum_1e6(block), present, ids


def to_linear_frame(matrix: pd.DataFrame, mode: str) -> pd.DataFrame:
    values = matrix.to_numpy(dtype=np.float64)
    if mode == "microarray":
        lin = np.exp2(values)
    else:
        lin = np.maximum(np.exp2(values) - 1.0, 0.0)
    lin = np.nan_to_num(lin, nan=0.0, posinf=0.0, neginf=0.0)
    return pd.DataFrame(lin, index=matrix.index, columns=matrix.columns)


def save_parquet(path: Path, block: np.ndarray, present: list[str], ids: list[str]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(block, columns=present)
    frame.insert(0, "sample_id", ids)
    frame.to_parquet(path, index=False)
    return sha256_file(path)


def missing_ld(genes: list[str], keep: np.ndarray, present: list[str]) -> int:
    have = set(present)
    return int(sum(1 for gene, flag in zip(genes, keep) if flag and gene not in have))


def compare_zk(data_root: Path, z_new: np.ndarray, k_new: np.ndarray, ids: list[str]) -> dict:
    old_ids = [line for line in (data_root / "sample_ids.txt").read_text().splitlines() if line]
    if old_ids != list(ids):
        return {"status": "sample_id_order_differs", "max_abs_z": None, "max_abs_k": None}
    old_z = np.load(data_root / "Z.npy")
    old_k = np.load(data_root / "K.npy")
    if old_z.shape != z_new.shape or old_k.shape != k_new.shape:
        return {"status": "shape_differs", "old_z": list(old_z.shape), "new_z": list(z_new.shape), "max_abs_z": None, "max_abs_k": None}
    return {
        "status": "compared",
        "max_abs_z": float(np.max(np.abs(old_z.astype(np.float64) - z_new.astype(np.float64)))),
        "max_abs_k": float(np.max(np.abs(old_k.astype(np.float64) - k_new.astype(np.float64)))),
    }


def rebuild_zk(name, block, present, ids, genes, b0, offsets, indices, data_root: Path, replace: bool, stage5_block, stage5_present) -> dict:
    z, k, info = aux.z_and_k(block, present, genes, b0, offsets, indices)
    z5, k5, _info5 = aux.z_and_k(stage5_block, stage5_present, genes, b0, offsets, indices)
    compared = {
        "status": "recomputed_stage5_procedure",
        "max_abs_z": float(np.max(np.abs(z.astype(np.float64) - z5.astype(np.float64)))),
        "max_abs_k": float(np.max(np.abs(k.astype(np.float64) - k5.astype(np.float64)))),
        "replaced": False,
    }
    if replace or compared["max_abs_z"] > 1e-6 or compared["max_abs_k"] > 1e-6:
        np.save(data_root / "Z.npy", z)
        np.save(data_root / "K.npy", k)
        (data_root / "sample_ids.txt").write_text("\n".join(ids) + "\n")
        info_path = data_root / "feature_info.json"
        old = json.loads(info_path.read_text()) if info_path.exists() else {}
        old.update(info)
        info_path.write_text(json.dumps(old))
        compared["replaced"] = True
        compared["new_z_sha256"] = sha256_file(data_root / "Z.npy")
        compared["new_k_sha256"] = sha256_file(data_root / "K.npy")
    return compared


def ordered(matrix: pd.DataFrame, data_root: Path) -> pd.DataFrame:
    ids = [line for line in (data_root / "sample_ids.txt").read_text().splitlines() if line]
    missing = [sample for sample in ids if sample not in matrix.columns]
    if missing:
        raise SystemExit(f"{data_root.name} missing {len(missing)} expression columns, e.g. {missing[:3]}")
    return matrix.loc[:, ids]


def gse209998_decision(samples: pd.DataFrame) -> dict:
    import re
    pattern = re.compile(PREFIX)
    prefix = samples["sample_id"].map(lambda value: pattern.match(str(value)).group(1) if pattern.match(str(value)) else None)
    if prefix.isna().any():
        raise SystemExit("GSE209998 barcode did not match the prefix expression")
    disease = samples["disease"].astype(str)
    primary = set(prefix.loc[disease.str.contains("Primary", case=False)])
    meta = set(prefix.loc[disease.str.contains("Metast", case=False)])
    both = primary & meta
    counts = prefix.value_counts()
    match_paper = int(prefix.nunique()) == PAPER_PATIENTS
    return {
        "regex": PREFIX,
        "n_samples": int(len(samples)),
        "n_prefix": int(prefix.nunique()),
        "paper_patients": PAPER_PATIENTS,
        "paper_source": PAPER_SOURCE,
        "n_prefix_with_multiple_samples": int((counts > 1).sum()),
        "n_prefix_primary_and_metastasis": int(len(both)),
        "samples_per_prefix": {str(k): int(v) for k, v in counts.value_counts().sort_index().items()},
        "check_multiple_and_paired": bool((counts > 1).any() and len(both) > 0),
        "check_count_equals_paper": match_paper,
        "decision": "확인 검정 포함" if match_paper and (counts > 1).any() and len(both) > 0 else "서술용",
        "patient_rule": "prefix" if match_paper else "barcode",
        "note": "고유 접두부 53과 논문 환자 55가 같지 않고, 빠진 2명의 식별자가 이 시리즈 메타데이터에 없다.",
    }


def native_map(path: Path) -> dict[str, str]:
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return dict(zip(table["standard_site"], table["native_organs"]))


def build_labels(paths, decision: dict, root: Path) -> pd.DataFrame:
    natives = native_map(Path(paths["config"]) / "mappings" / "site_native.tsv")
    rows = []

    def add_frame(cohort, layer, frame, library):
        frame = frame.copy()
        flagged = frame["flag_candidate"].astype(str).eq("True") if "flag_candidate" in frame.columns else pd.Series(False, index=frame.index)
        for row in frame.itertuples(index=False):
            site = str(getattr(row, "standard_site"))
            in_eval = str(getattr(row, "in_eval")) == "True"
            in_risk = str(getattr(row, "in_risk")) == "True"
            in_native = str(getattr(row, "in_native")) == "True"
            flag = bool(flagged.loc[row.Index]) if False else str(getattr(row, "flag_candidate", "False")) == "True"
            rule = str(getattr(row, "organ_rule"))
            if flag:
                reason = "duplicate"
            elif not in_eval:
                reason = f"organ_rule:{rule}"
            else:
                reason = ""
            excluded = reason != ""
            rows.append({
                "cohort": cohort,
                "layer": layer,
                "sample_id": row.sample_id,
                "patient_id": row.patient_id,
                "library": library if library else str(getattr(row, "library", "")),
                "organ": row.organ,
                "organ_rule": rule,
                "raw_site": row.raw_site,
                "standard_site": site,
                "native_organs": natives.get(site, ""),
                "in_eval": in_eval and not excluded,
                "in_risk": in_risk and not excluded,
                "in_native": in_native and not excluded,
                "pool_out_site": site in POOL_OUT,
                "in_pool_out_risk": site in POOL_OUT and in_risk and not excluded,
                "excluded": excluded,
                "exclude_reason": reason,
            })

    cohort_root = root / "cohorts"
    simple = {
        "blca_iatlas_imvigor210_2017": "RNA-seq",
        "brca_iatlas_anders_2022": "RNA-seq",
        "mel_dfci_2019": "RNA-seq",
        "paad_iatlas_prince_2022": "RNA-seq",
        "GSE50760": "RNA-seq",
        "GSE209998": "RNA-seq",
        "GSE41258": "마이크로어레이",
        "GSE14018": "마이크로어레이",
        "GSE71729": "마이크로어레이",
        "GSE74685": "마이크로어레이",
        "prad_fhcrc": "마이크로어레이",
    }
    for cohort, layer in simple.items():
        frame = pd.read_csv(cohort_root / cohort / "samples.tsv", sep="\t", dtype=str).reset_index(drop=True)
        if cohort == "GSE209998" and decision["patient_rule"] == "barcode":
            frame["patient_id"] = frame["sample_id"]
        add_frame(cohort, layer, frame, "")
    poly = pd.read_csv(cohort_root / "prad_su2c_2019__fpkm_polya" / "samples.tsv", sep="\t", dtype=str)
    cap = pd.read_csv(cohort_root / "prad_su2c_2019__fpkm_capture" / "samples.tsv", sep="\t", dtype=str)
    cap_ids = set(cap["sample_id"])
    chosen = []
    for row in poly.itertuples(index=False):
        other = str(getattr(row, "OTHER_SAMPLE_ID", ""))
        partner = cap.loc[cap["sample_id"] == row.sample_id]
        if partner.empty and other in cap_ids:
            partner = cap.loc[cap["sample_id"] == other]
        flag = str(row.flag_candidate) == "True" or (not partner.empty and (partner["flag_candidate"] == "True").any())
        item = row._asdict()
        item["flag_candidate"] = str(flag)
        item["library"] = "polyA"
        chosen.append(item)
        if not partner.empty:
            cap_ids.discard(str(partner["sample_id"].iloc[0]))
    for row in cap.itertuples(index=False):
        if row.sample_id not in cap_ids:
            continue
        item = row._asdict()
        item["library"] = "capture"
        chosen.append(item)
    add_frame("prad_su2c_2019", "RNA-seq", pd.DataFrame(chosen), "")
    labels = pd.DataFrame(rows)
    risk_n = labels.loc[labels["in_risk"]].groupby("cohort").size()
    confirm = {}
    for cohort, layer in labels.groupby("cohort")["layer"].first().items():
        enough = int(risk_n.get(cohort, 0)) >= 10
        if cohort == "GSE209998" and decision["decision"] == "서술용":
            confirm[cohort] = False
        elif layer == "RNA-seq" and enough:
            confirm[cohort] = True
        else:
            confirm[cohort] = False
    labels["include_confirm"] = [bool(confirm[c]) and (not ex) for c, ex in zip(labels["cohort"], labels["excluded"])]
    labels["layer_test"] = [int(risk_n.get(c, 0)) >= 10 and (not ex) for c, ex in zip(labels["cohort"], labels["excluded"])]
    for column, mask_col in (
        ("selected_eval", "in_eval"),
        ("selected_risk", "in_risk"),
        ("selected_native", "in_native"),
        ("selected_pool_out", "in_pool_out_risk"),
    ):
        labels[column] = False
        part = labels.loc[labels[mask_col]].copy()
        part["library_rank"] = np.where(part["library"].eq("polyA"), 0, 1)
        part = part.sort_values(["cohort", "patient_id", "library_rank", "sample_id"])
        keep = part.drop_duplicates(["cohort", "patient_id"]).index
        labels.loc[keep, column] = True
    return labels


def power_table(labels: pd.DataFrame) -> pd.DataFrame:
    confirm = labels.loc[labels["include_confirm"] & labels["layer"].eq("RNA-seq")]
    n = {
        "AH1": int(confirm.loc[confirm["selected_risk"]].shape[0]),
        "AH2": int(confirm.loc[confirm["selected_eval"]].shape[0]),
        "AS1": int(confirm.loc[confirm["selected_native"]].shape[0]),
        "AS2": int(confirm.loc[confirm["selected_pool_out"]].shape[0]),
        "AS3": int(confirm.loc[confirm["selected_eval"]].shape[0]),
    }
    scenarios = {
        "AH1": [(77 / 378, 0 / 378), (34 / 361, 1 / 361)],
        "AH2": [(61 / 512, 23 / 512), (29 / 437, 12 / 437)],
        "AS1": [(7 / 91, 1 / 91), (5.5 / 91, 2.5 / 91)],
        "AS2": [(0.084, 0.002), (0.043, 0.002)],
        "AS3": [(19 / 512, 9 / 512), (8 / 437, 5 / 437)],
    }
    alpha = {"AH1": 0.05, "AH2": 0.05, "AS1": 0.05 / 3, "AS2": 0.05 / 3, "AS3": 0.05 / 3}
    rng = np.random.default_rng(20261012)
    rows = []
    for name in ("AH1", "AH2", "AS1", "AS2", "AS3"):
        for i, (p_fav, p_against) in enumerate(scenarios[name], start=1):
            p_tie = 1.0 - p_fav - p_against
            if p_tie < -1e-12:
                raise SystemExit(f"probabilities exceed 1 for {name}")
            draws = rng.multinomial(n[name], [p_fav, p_against, max(p_tie, 0.0)], size=10000)
            favor, against = draws[:, 0], draws[:, 1]
            m = favor + against
            pvalues = np.ones(len(m))
            ok = m > 0
            pvalues[ok] = binom.sf(favor[ok] - 1, m[ok], 0.5)
            rows.append({
                "hypothesis": name, "scenario": i, "n": n[name],
                "p_favor": p_fav, "p_against": p_against, "alpha": alpha[name],
                "power": float(np.mean(pvalues <= alpha[name])),
            })
    return pd.DataFrame(rows)


def recount_as1(paths) -> dict:
    samples = pd.read_csv(Path(paths["results"]) / "stage5" / "posthoc" / "pog_samples.tsv", sep="\t", dtype=str)
    natives = native_map(Path(paths["config"]) / "mappings" / "site_native.tsv")
    organs = samples["site"].map(lambda site: set(filter(None, natives.get(str(site), "").split("|"))))
    native = np.array([truth in native_set and bool(native_set) for truth, native_set in zip(samples["truth"], organs)])
    beta = pd.to_numeric(samples["beta"], errors="coerce").to_numpy()
    use_base = np.isnan(beta) | (beta <= 0.02)
    pred = np.where(use_base, samples["pred__BASE-Z"], samples["pred__SA-Z"])
    truth = samples["truth"].to_numpy()
    only_g = native & (pred == truth) & (samples["pred__SA-Z"].to_numpy() != truth)
    only_sa = native & (samples["pred__SA-Z"].to_numpy() == truth) & (pred != truth)
    return {
        "n_native": int(native.sum()),
        "sa_g_only": samples.loc[only_g, "id"].tolist(),
        "sa_z_only": samples.loc[only_sa, "id"].tolist(),
        "n_sa_g_correct": int((native & (pred == truth)).sum()),
        "n_sa_z_correct": int((native & (samples["pred__SA-Z"].to_numpy() == truth)).sum()),
    }


def handle_matrix(name, matrix, kind, judgment, lengths, genes, b0, offsets, indices, paths, keep, replace_zk, data_name):
    data_root = Path(paths["data_processed"]) / "aux" / data_name
    matrix = ordered(matrix, data_root)
    if kind == "counts":
        block, present, ids, missing_len = length_then_sum(matrix, lengths, genes)
        stage5_block, present5 = block, present
    elif kind == "paad":
        linear = to_linear_frame(matrix, "rna")
        block, present, ids, missing_len = length_then_sum(linear, lengths, genes)
        stage5_block, present5, ids5 = sum_genes(matrix, genes)
        if ids5 != ids:
            raise SystemExit(f"{name} sample order differs")
    elif judgment == "로그":
        linear = to_linear_frame(matrix, "microarray" if kind == "microarray" else "rna")
        block, present, ids = sum_genes(linear, genes)
        stage5_block, present5, ids5 = sum_genes(matrix, genes)
        missing_len = 0
        if ids5 != ids or present5 != present:
            raise SystemExit(f"{name} alignment differs after the log transform")
    else:
        block, present, ids = sum_genes(matrix, genes)
        stage5_block, present5 = block, present
        missing_len = 0
    out = Path(paths["results"]).parents[0] / "data" / "aux" / "ld_input" / f"{name}.parquet"
    # paths are project-relative; data/aux is next to data/processed
    out = Path(paths["data_processed"]).parents[0] / "aux" / "ld_input" / f"{name}.parquet"
    digest = save_parquet(out, block, present, ids)
    record = {
        "cohort": name, "kind": kind, "judgment": judgment, "sha256": digest, "path": str(out.relative_to(Path(paths["data_processed"]).parents[1] if False else Path(paths["config"]).parents[0])),
        "n_samples": len(ids), "n_genes": len(present), "n_ld_missing": missing_ld(genes, keep, present),
        "length_missing_symbols": missing_len,
    }
    record["path"] = str(out.relative_to(Path(paths["config"]).parents[0]))
    if kind == "paad" or judgment == "로그":
        record["zk"] = rebuild_zk(name, block, present, ids, genes, b0, offsets, indices, data_root, replace_zk, stage5_block, present5)
    if kind == "paad":
        patients = pd.read_csv(Path(paths["results"]) / "stage5" / "aux" / "cohorts" / name / "samples.tsv", sep="\t", dtype=str)
        patients = patients.set_index("sample_id").loc[ids]
        refs = aux.load_reference_tpm(paths, genes)
        dup_dir = Path(paths["results"]) / "stage6" / "paad_duplicate"
        dup_dir.mkdir(parents=True, exist_ok=True)
        dup = aux.duplicate_screen(block, present, patients["patient_id"].tolist(), ids, name, genes, refs, dup_dir)
        patients = patients.reset_index()
        merged = patients.drop(columns=[c for c in ("external_max", "internal_max", "flag_candidate") if c in patients.columns], errors="ignore")
        merged = merged.merge(dup[["sample_id", "external_max", "internal_max", "flag_candidate"]], on="sample_id", how="left")
        merged.to_csv(Path(paths["results"]) / "stage5" / "aux" / "cohorts" / name / "samples.tsv", sep="\t", index=False)
        record["n_flag"] = int(dup["flag_candidate"].sum())
    return record


def main() -> None:
    paths = load_paths()
    out_dir = Path(paths["results"]) / "stage6"
    out_dir.mkdir(parents=True, exist_ok=True)
    downloads = Path(paths["data_processed"]) / "aux" / "downloads"
    genes, b0_idx, offsets, indices, _names = load_gene_pack(paths)
    b0 = [genes[int(i)] for i in b0_idx]
    _mu, _h, keep, ref_genes = load_references(Path(paths["results"]) / "stage3" / "references.npz")
    if ref_genes != genes:
        raise SystemExit("reference gene order does not match G")
    lengths = pd.read_csv(downloads / "gencode_v23_exon_union.tsv", sep="\t").set_index("symbol")["exon_union_bp"]
    gse_samples = pd.read_csv(Path(paths["results"]) / "stage5" / "aux" / "cohorts" / "GSE209998" / "samples.tsv", sep="\t", dtype=str)
    decision = gse209998_decision(gse_samples)
    (out_dir / "gse209998_patient.json").write_text(json.dumps(decision, indent=2, ensure_ascii=False))
    print("GSE209998", decision["decision"], decision["n_prefix"], flush=True)
    as1 = recount_as1(paths)
    (out_dir / "as1_recount.json").write_text(json.dumps(as1, indent=2))
    print("AS1", as1["n_native"], len(as1["sa_g_only"]), len(as1["sa_z_only"]), flush=True)

    records = []
    jobs = []
    # cBio RNA-seq except PAAD, which has its own rule
    for study, filename, kind in (
        ("blca_iatlas_imvigor210_2017", "data_mrna_seq_tpm.txt", "tpm"),
        ("brca_iatlas_anders_2022", "data_mrna_seq_tpm.txt", "tpm"),
        ("mel_dfci_2019", "data_mrna_seq_tpm.txt", "tpm"),
        ("prad_su2c_2019", "data_mrna_seq_fpkm_polya.txt", "fpkm"),
        ("prad_su2c_2019", "data_mrna_seq_fpkm_capture.txt", "fpkm"),
    ):
        matrix = aux.read_cbio_expression(downloads / f"{study}__{filename}")
        data_name = study if "fpkm_" not in filename else f"{study}__{filename.replace('data_mrna_seq_', '').replace('.txt', '')}"
        jobs.append((data_name, matrix, kind, data_name, False))
    paad = aux.read_cbio_expression(downloads / "paad_iatlas_prince_2022__data_mrna_seq_expression.txt")
    jobs.append(("paad_iatlas_prince_2022", paad, "paad", "paad_iatlas_prince_2022", True))
    samples, matrix = geo.gse50760(downloads)
    jobs.append(("GSE50760", matrix, "fpkm", "GSE50760", False))
    _samples, matrix = geo.gse209998(downloads, lengths)
    jobs.append(("GSE209998", matrix, "counts", "GSE209998", False))
    annot = geo.gpl96_annot(downloads)
    for acc in ("GSE14018", "GSE41258"):
        matrix, _ids = geo.series_matrix(downloads / f"{acc}_series_matrix.txt.gz", annot)
        jobs.append((acc, matrix, "microarray", acc, False))
    fhcrc = aux.read_cbio_expression(downloads / "prad_fhcrc_data_mrna_agilent_microarray.txt")
    jobs.append(("prad_fhcrc", fhcrc, "microarray", "prad_fhcrc", False))

    judgments = []
    for name, matrix, kind, data_name, replace in jobs:
        data_root = Path(paths["data_processed"]) / "aux" / data_name
        matrix = ordered(matrix, data_root)
        if kind == "paad":
            stat = {"judgment": "PAAD 고정: log2(x+1)로 보고 2^v-1", "q0": None, "q25": None, "q99": None, "q100": None, "n_values": int(np.isfinite(matrix.to_numpy(dtype=float)).sum())}
        else:
            stat = quantiles(matrix)
        stat["cohort"] = name
        judgments.append(stat)
        print("judge", name, stat["judgment"], flush=True)
        record = handle_matrix(name, matrix, kind, stat["judgment"], lengths, genes, b0, offsets, indices, paths, keep, replace, data_name)
        record["quantiles"] = {key: stat[key] for key in ("q0", "q25", "q99", "q100", "n_values")}
        records.append(record)
        if "zk" in record:
            print("zk", name, record["zk"]["max_abs_z"], record["zk"]["max_abs_k"], record["zk"]["replaced"], flush=True)
        del matrix

    print("two-channel", flush=True)
    for acc, parser in (
        ("GSE71729", two.parse_genepix),
        ("GSE74685", None),
    ):
        header = geo.header_frame(Path(paths["results"]) / "stage5" / "aux" / "geo_headers" / f"{acc}_header.txt")
        wanted = set(header["geo_accession"].astype(str))
        if acc == "GSE74685":
            symbols_of = two.load_gpl_symbols(downloads / "GPL15659_data.txt")
            parser = lambda raw, table=symbols_of: two.parse_agilent(raw, table)
        _probes, symbols, values, found = two.read_tar(downloads / f"{acc}_RAW.tar", wanted, parser)
        matrix = two.collapse(symbols, values, found)
        stat = quantiles(matrix)
        stat["cohort"] = acc
        judgments.append(stat)
        print("judge", acc, stat["judgment"], flush=True)
        record = handle_matrix(acc, matrix, "microarray", stat["judgment"], lengths, genes, b0, offsets, indices, paths, keep, False, acc)
        record["quantiles"] = {key: stat[key] for key in ("q0", "q25", "q99", "q100", "n_values")}
        records.append(record)
        del values, matrix

    # SU2C one biological sample: polyA row when both exist. Capture-only rows stay.
    poly = next(row for row in records if row["cohort"] == "prad_su2c_2019__fpkm_polya")
    cap = next(row for row in records if row["cohort"] == "prad_su2c_2019__fpkm_capture")
    poly_df = pd.read_parquet(Path(paths["config"]).parents[0] / poly["path"])
    cap_df = pd.read_parquet(Path(paths["config"]).parents[0] / cap["path"])
    samples_poly = pd.read_csv(Path(paths["results"]) / "stage5" / "aux" / "cohorts" / "prad_su2c_2019__fpkm_polya" / "samples.tsv", sep="\t", dtype=str)
    samples_cap = pd.read_csv(Path(paths["results"]) / "stage5" / "aux" / "cohorts" / "prad_su2c_2019__fpkm_capture" / "samples.tsv", sep="\t", dtype=str)
    # labels are built from samples; LD file uses the same chosen sample ids after labels exist
    labels = build_labels(paths, decision, Path(paths["results"]) / "stage5" / "aux")
    chosen_ids = labels.loc[labels["cohort"].eq("prad_su2c_2019"), ["sample_id", "library"]]
    parts = []
    poly_ids = set(chosen_ids.loc[chosen_ids["library"].eq("polyA"), "sample_id"])
    cap_ids = set(chosen_ids.loc[chosen_ids["library"].eq("capture"), "sample_id"])
    if poly_ids:
        parts.append(poly_df.loc[poly_df["sample_id"].isin(poly_ids)])
    if cap_ids:
        parts.append(cap_df.loc[cap_df["sample_id"].isin(cap_ids)])
    merged = pd.concat(parts, ignore_index=True)
    merged_path = Path(paths["data_processed"]).parents[0] / "aux" / "ld_input" / "prad_su2c_2019.parquet"
    merged.to_parquet(merged_path, index=False)
    present = [col for col in merged.columns if col != "sample_id"]
    records.append({
        "cohort": "prad_su2c_2019", "kind": "fpkm", "judgment": "라이브러리별 판정 후 polyA 우선",
        "sha256": sha256_file(merged_path),
        "path": str(merged_path.relative_to(Path(paths["config"]).parents[0])),
        "n_samples": int(len(merged)), "n_genes": len(present),
        "n_ld_missing": missing_ld(genes, keep, present),
        "n_polya": len(poly_ids), "n_capture_only": len(cap_ids),
        "library_judgments": {"polyA": poly["judgment"], "capture": cap["judgment"]},
    })
    labels_path = Path(paths["config"]) / "aux_eval_labels.tsv"
    labels.to_csv(labels_path, sep="\t", index=False)
    power = power_table(labels)
    power.to_csv(out_dir / "power.tsv", sep="\t", index=False)
    pd.DataFrame(judgments).to_csv(out_dir / "log_judgment.tsv", sep="\t", index=False)
    (out_dir / "ld_inputs.json").write_text(json.dumps(records, indent=2, ensure_ascii=False))
    print("labels", len(labels), "power rows", len(power), flush=True)
    print(power.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
