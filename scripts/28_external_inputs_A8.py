"""Build SCOPE and CUP-AI-Dx input matrices. Does not run either tool."""
from __future__ import annotations

import gzip
import hashlib
import json
import re
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/kangjw/WORK/Cisplatin/2026_NewProj/proj_A_host_too")
OUT = ROOT / "results/stage8/external/inputs"
SCOPE_GENES = ROOT / "results/stage7/external_tools/cancerscope/cancerscope/resources/scope_features_genenames.txt"
CUP_FEATURES = ROOT / "results/stage7/external_tools/CUP-AI-Dx/data/features_791.csv"
GTF = ROOT / "data/processed/aux/downloads/gencode.v23.annotation.gtf.gz"
LENGTHS = ROOT / "data/processed/aux/downloads/gencode_v23_exon_union.tsv"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
RNA = [
    "blca_iatlas_imvigor210_2017",
    "brca_iatlas_anders_2022",
    "mel_dfci_2019",
    "paad_iatlas_prince_2022",
    "GSE50760",
    "prad_su2c_2019",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def scope_features() -> list[str]:
    frame = pd.read_csv(SCOPE_GENES, sep="\t", dtype=str)
    symbols = frame["HUGO"].tolist()
    if len(symbols) != 17688:
        raise SystemExit(f"SCOPE features {len(symbols)}")
    return symbols


def cup_features() -> list[str]:
    raw = pd.read_csv(CUP_FEATURES, header=None, dtype=str)[0].tolist()
    if len(raw) != 791:
        raise SystemExit(f"CUP features {len(raw)}")
    return [item[1:] if item.startswith("X") else item for item in raw]


def gencode_symbol() -> dict[str, str]:
    pattern = re.compile(r'gene_id "([^"]+)";.*gene_name "([^"]+)"')
    found: dict[str, str] = {}
    with gzip.open(GTF, "rt") as handle:
        for line in handle:
            if line.startswith("#") or "\tgene\t" not in line:
                continue
            match = pattern.search(line)
            if not match:
                continue
            gene_id = match.group(1).split(".")[0]
            name = match.group(2)
            previous = found.get(gene_id)
            if previous is not None and previous != name:
                found[gene_id] = ""
            elif previous is None:
                found[gene_id] = name
    return {gene_id: name for gene_id, name in found.items() if name}


def hgnc_maps(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    current: dict[str, str] = {}
    for symbol, entrez in zip(frame["symbol"], frame["entrez_id"]):
        if symbol and entrez:
            current[symbol] = entrez
    previous_hits: dict[str, set[str]] = {}
    for raw, entrez in zip(frame["prev_symbol"], frame["entrez_id"]):
        if not raw or not entrez:
            continue
        for symbol in raw.split("|"):
            symbol = symbol.strip()
            if symbol:
                previous_hits.setdefault(symbol, set()).add(entrez)
    previous = {symbol: next(iter(ids)) for symbol, ids in previous_hits.items() if len(ids) == 1}
    return current, previous


def download_hgnc() -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "hgnc_complete_set.txt"
    if not path.exists():
        urllib.request.urlretrieve(HGNC_URL, path)
    return path


def lengths_kb() -> dict[str, float]:
    frame = pd.read_csv(LENGTHS, sep="\t", dtype=str)
    out = {}
    for symbol, bp in zip(frame["symbol"], frame["exon_union_bp"]):
        value = float(bp)
        if value > 0 and symbol not in out:
            out[symbol] = value / 1000.0
    return out


def read_selected() -> dict[str, pd.DataFrame]:
    met = pd.read_csv(ROOT / "results/stage5/posthoc/met500_samples.tsv", sep="\t", dtype=str)
    if len(met) != 437:
        raise SystemExit(f"MET500 n={len(met)}")
    pog_labels = pd.read_csv(ROOT / "config/pog570_eval_labels.tsv", sep="\t", dtype=str)
    pog = pog_labels.loc[~pog_labels["organ"].isin(["exclude", "NA"])].copy()
    if len(pog) != 512:
        raise SystemExit(f"POG eval n={len(pog)}")
    confirm = pd.read_parquet(
        ROOT / "results/stage7/confirm/predictions.parquet",
        columns=["cohort", "sample_id", "patient_id", "library", "include_confirm", "selected_eval"],
    )
    confirm["include_confirm"] = confirm["include_confirm"].astype(str).str.lower().isin(["true", "1"])
    confirm["selected_eval"] = confirm["selected_eval"].astype(str).str.lower().isin(["true", "1"])
    aux = confirm.loc[confirm["include_confirm"] & confirm["selected_eval"] & confirm["cohort"].isin(RNA)].copy()
    if len(aux) != 729:
        raise SystemExit(f"aux eval n={len(aux)}")
    return {"MET500": met, "POG570": pog, "aux": aux}


def write_matrix(path: Path, frame: pd.DataFrame) -> str:
    with gzip.open(path, "wt") as handle:
        frame.to_csv(handle, sep="\t", float_format="%.8e")
    return sha256_file(path)


def align(linear: pd.DataFrame, symbols: list[str]) -> tuple[pd.DataFrame, int]:
    present = [gene for gene in symbols if gene in linear.index]
    missing = len(symbols) - len(present)
    out = linear.reindex(symbols).fillna(0.0).astype(np.float64)
    return out, missing


def cup_from_linear(linear: pd.DataFrame, features: list[str], current: dict[str, str], previous: dict[str, str]) -> tuple[pd.DataFrame, int]:
    buckets: dict[str, list[str]] = {}
    for symbol in linear.index:
        entrez = current.get(symbol) or previous.get(symbol)
        if entrez:
            buckets.setdefault(entrez, []).append(symbol)
    rows = {}
    for entrez, names in buckets.items():
        rows[entrez] = linear.loc[names].sum(axis=0)
    collapsed = pd.DataFrame(rows).T if rows else pd.DataFrame(index=[], columns=linear.columns)
    logged = np.log2(collapsed.clip(lower=0) + 1.0)
    missing = sum(feature not in logged.index for feature in features)
    out = logged.reindex(features).fillna(0.0).T
    out.columns = [f"X{feature}" for feature in features]
    return out, missing


def collapse_symbols(values: pd.DataFrame, symbols: pd.Series) -> pd.DataFrame:
    keep = symbols.notna() & symbols.ne("")
    use = values.loc[keep].copy()
    use.index = symbols.loc[keep].to_numpy()
    return use.groupby(level=0).sum()


def load_wide(path: Path, gene_col: str, samples: list[str]) -> pd.DataFrame:
    header = pd.read_csv(path, sep="\t", nrows=0).columns.astype(str).tolist()
    missing = [sample for sample in samples if sample not in header]
    if missing:
        raise SystemExit(f"{path.name} missing {len(missing)} samples, first {missing[:5]}")
    frame = pd.read_csv(path, sep="\t", usecols=[gene_col, *samples], dtype={gene_col: str})
    frame[gene_col] = frame[gene_col].astype(str)
    values = frame.drop(columns=[gene_col]).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    values.index = frame[gene_col]
    return values


def ensembl_symbols(index: pd.Index, gene_to_symbol: dict[str, str]) -> pd.Series:
    keys = [str(gene).split(".")[0] for gene in index]
    return pd.Series([gene_to_symbol.get(key, "") for key in keys], index=index)


def gse50760(samples: list[str]) -> pd.DataFrame:
    frames = []
    with tarfile.open(ROOT / "data/processed/aux/downloads/GSE50760_RAW.tar") as tar:
        for member in tar.getmembers():
            if not member.name.endswith(".gz"):
                continue
            raw = gzip.decompress(tar.extractfile(member).read()).decode().splitlines()
            sample = raw[0].split("\t")[1]
            if sample not in samples:
                continue
            symbols, values = [], []
            for line in raw[1:]:
                gene, value = line.split("\t")
                symbols.append(gene)
                values.append(float(value))
            frames.append(pd.Series(values, index=symbols, name=sample))
    matrix = pd.concat(frames, axis=1).fillna(0.0).groupby(level=0).sum()
    missing = [sample for sample in samples if sample not in matrix.columns]
    if missing:
        raise SystemExit(f"GSE50760 missing {missing[:5]}")
    return matrix.loc[:, samples]


def build_one(name: str, linear: pd.DataFrame, full_sum: pd.Series | None, mode: str, scope_genes: list[str], cup_genes: list[str], current, previous, length_kb) -> dict:
    if mode == "fpkm":
        scope_linear = linear
        cup_source = linear.div(full_sum, axis=1) * 1e6
    elif mode == "tpm":
        scope_linear = linear
        cup_source = linear
    elif mode == "uq":
        kept = linear.loc[linear.index.isin(length_kb)]
        divided = kept.div([length_kb[gene] for gene in kept.index], axis=0)
        scope_linear = divided
        denom = divided.sum(axis=0).replace(0, np.nan)
        cup_source = divided.div(denom, axis=1) * 1e6
        cup_source = cup_source.fillna(0.0)
    else:
        raise SystemExit(mode)
    scope, scope_missing = align(scope_linear, scope_genes)
    cup, cup_missing = cup_from_linear(cup_source, cup_genes, current, previous)
    scope_path = OUT / f"{name}.scope.tsv.gz"
    cup_path = OUT / f"{name}.cup.tsv.gz"
    scope.index.name = "HUGO"
    cup.index.name = "sample"
    return {
        "cohort": name,
        "n_samples": int(scope.shape[1]),
        "n_scope_matched": 17688 - scope_missing,
        "n_scope_missing": scope_missing,
        "n_cup_missing": cup_missing,
        "sha256_scope_matrix": write_matrix(scope_path, scope),
        "sha256_cup_matrix": write_matrix(cup_path, cup),
        "scope_path": str(scope_path.relative_to(ROOT)),
        "cup_path": str(cup_path.relative_to(ROOT)),
    }


def undo_log(values: pd.DataFrame) -> pd.DataFrame:
    return (np.power(2.0, values) - 1.0).clip(lower=0.0)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    hgnc_path = download_hgnc()
    current, previous = hgnc_maps(hgnc_path)
    scope_genes = scope_features()
    cup_genes = cup_features()
    gene_to_symbol = gencode_symbol()
    length_kb = lengths_kb()
    selected = read_selected()
    records = []

    met_samples = selected["MET500"]["id"].tolist()
    met = load_wide(ROOT / "data/raw/met500/M.mx.txt.gz", "sample", met_samples)
    met_sum = met.sum(axis=0)
    met_linear = collapse_symbols(met, ensembl_symbols(met.index, gene_to_symbol))
    records.append(build_one("MET500", met_linear, met_sum, "fpkm", scope_genes, cup_genes, current, previous, length_kb))
    del met
    print("MET500", records[-1]["n_scope_matched"], records[-1]["n_cup_missing"], flush=True)

    pog_samples = selected["POG570"]["PATIENT_ID"].astype(str).tolist()
    pog = load_wide(ROOT / "data/raw/POG570/POG570_TPM_expression.txt.gz", "genes", pog_samples)
    pog_linear = collapse_symbols(pog, ensembl_symbols(pog.index, gene_to_symbol))
    records.append(build_one("POG570", pog_linear, None, "tpm", scope_genes, cup_genes, current, previous, length_kb))
    del pog
    print("POG570", records[-1]["n_scope_matched"], records[-1]["n_cup_missing"], flush=True)

    aux = selected["aux"]
    sources = {
        "blca_iatlas_imvigor210_2017": (ROOT / "data/processed/aux/downloads/blca_iatlas_imvigor210_2017__data_mrna_seq_tpm.txt", "log_tpm"),
        "brca_iatlas_anders_2022": (ROOT / "data/processed/aux/downloads/brca_iatlas_anders_2022__data_mrna_seq_tpm.txt", "log_tpm"),
        "mel_dfci_2019": (ROOT / "data/processed/aux/downloads/mel_dfci_2019__data_mrna_seq_tpm.txt", "tpm"),
        "paad_iatlas_prince_2022": (ROOT / "data/processed/aux/downloads/paad_iatlas_prince_2022__data_mrna_seq_expression.txt", "uq"),
        "prad_su2c_2019": (None, "fpkm"),
    }
    for cohort, (path, mode) in sources.items():
        part = aux.loc[aux["cohort"] == cohort]
        if cohort == "prad_su2c_2019":
            pieces = []
            for library, filename in (("polyA", "prad_su2c_2019__data_mrna_seq_fpkm_polya.txt"), ("capture", "prad_su2c_2019__data_mrna_seq_fpkm_capture.txt")):
                ids = part.loc[part["library"] == library, "sample_id"].tolist()
                if not ids:
                    continue
                piece = load_wide(ROOT / "data/processed/aux/downloads" / filename, "Hugo_Symbol", ids)
                pieces.append(piece)
            values = pd.concat(pieces, axis=1)
            values = values.groupby(level=0).sum()
            full_sum = values.sum(axis=0)
            linear = values
            use_mode = "fpkm"
        elif cohort == "paad_iatlas_prince_2022":
            values = load_wide(path, "Hugo_Symbol", part["sample_id"].tolist())
            values = undo_log(values).groupby(level=0).sum()
            full_sum = None
            linear = values
            use_mode = "uq"
        else:
            values = load_wide(path, "Hugo_Symbol", part["sample_id"].tolist())
            if mode == "log_tpm":
                values = undo_log(values)
            values = values.groupby(level=0).sum()
            full_sum = None
            linear = values
            use_mode = "tpm"
        records.append(build_one(cohort, linear, full_sum, use_mode, scope_genes, cup_genes, current, previous, length_kb))
        print(cohort, records[-1]["n_samples"], records[-1]["n_scope_matched"], records[-1]["n_cup_missing"], flush=True)

    gse_ids = aux.loc[aux["cohort"] == "GSE50760", "sample_id"].tolist()
    gse = gse50760(gse_ids)
    records.append(build_one("GSE50760", gse, gse.sum(axis=0), "fpkm", scope_genes, cup_genes, current, previous, length_kb))
    print("GSE50760", records[-1]["n_samples"], records[-1]["n_scope_matched"], records[-1]["n_cup_missing"], flush=True)

    summary = {
        "hgnc_sha256": sha256_file(hgnc_path),
        "hgnc_path": str(hgnc_path.relative_to(ROOT)),
        "cohorts": records,
    }
    (OUT / "input_matrix_summary.json").write_text(json.dumps(summary, indent=2))
    print("wrote", OUT / "input_matrix_summary.json", flush=True)


if __name__ == "__main__":
    main()
