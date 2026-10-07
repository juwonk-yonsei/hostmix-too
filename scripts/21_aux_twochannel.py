#!/usr/bin/env python3
"""Build Z/K for the two-channel GEO starters from the sample channel only.

GSE71729: ch1 is a common reference, ch2 is the sample, column CH2I_MEAN.
GSE74685: ch1 is a cell-line reference pool, ch2 is the sample, column gProcessedSignal.
Does not apply a classifier.
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("auxgeo", ROOT / "20_aux_geo.py")
geo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(geo)
aux = geo.aux

sys.path.insert(0, str(ROOT))
from stage2_common import load_gene_pack
from stage5_lib import load_paths

SITE_71729 = {
    "Pancreas_Primary": "pancreas",
    "Liver_Metastasis": "liver",
    "Lung_Metastasis": "lung",
    "LymphNode_Metastasis": "lymph_node",
    "Peritoneal_Metastasis": "peritoneum",
    "AbWall_Metastasis": "soft_tissue",
    "Fat_Metastasis": "soft_tissue",
    "Colon_Metastasis": "colon_rectum",
    "Diaphragm_Metastasis": "other",
    "Duo_Metastasis": "other",
}
SITE_74685 = {
    "LN": "lymph_node",
    "LUNG": "lung",
    "LIVER": "liver",
    "ADRENAL": "adrenal",
    "BONE": "bone",
    "SKIN": "skin",
    "PERITONEAL": "peritoneum",
    "PERITONEUM": "peritoneum",
    "RENAL": "kidney",
    "KIDNEY": "kidney",
    "RETROPERITONEAL": "other",
    "SCROTUM": "other",
    "SPLEEN": "other",
    "APPENDIX": "other",
}


def single_symbol(text: str) -> str:
    symbol = str(text).strip().strip('"')
    if symbol in {"", "---", "NA", "null", "None", "nan"}:
        return ""
    if "///" in symbol or ";" in symbol or "|" in symbol or " " in symbol:
        return ""
    return symbol


def field_column(header: pd.DataFrame, label: str) -> str:
    hits = [col for col in header.columns if str(col).endswith("__" + label)]
    if not hits:
        raise SystemExit(f"missing characteristic {label}")
    return hits[0]


def after_colon(value: str) -> str:
    text = str(value)
    return text.split(":", 1)[1].strip() if ":" in text else text.strip()


def samples_71729(header: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    source = dict(zip(header["geo_accession"], header["source_name_ch2"].astype(str)))
    rows = []
    for sample in sample_ids:
        raw = source.get(sample, "")
        if raw == "CellLine":
            organ, site, rule = "exclude", "", "2"
        elif raw.endswith("_Normal"):
            organ, site, rule = "exclude", "", "8"
        elif raw in SITE_71729:
            organ, site, rule = "Pancreas", SITE_71729[raw], "cohort"
        elif raw == "":
            organ, site, rule = "exclude", "", "3"
        else:
            organ, site, rule = "Pancreas", "other", "cohort"
        rows.append({
            "sample_id": sample, "patient_id": sample, "raw_site": raw,
            "standard_site": site, "organ": organ, "organ_rule": rule,
            "library": "GPL20769 CH2I_MEAN", "lcm": False,
            "patient_note": "no patient field in the series header; patient_id is the GSM",
        })
    return pd.DataFrame(rows)


def characteristic(header: pd.DataFrame, row: pd.Series, label: str) -> str:
    """Match the label inside the sample, because GEO characteristic rows are not aligned."""
    prefix = label.casefold() + ":"
    for col in header.columns:
        if not str(col).startswith("characteristics_ch2"):
            continue
        text = str(row[col])
        if text.casefold().startswith(prefix):
            return text.split(":", 1)[1].strip()
    return ""


def samples_74685(header: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    indexed = header.set_index("geo_accession")
    site_of, ne_of, patient_of = {}, {}, {}
    for sample in sample_ids:
        row = indexed.loc[sample]
        site_of[sample] = characteristic(header, row, "tumor site")
        ne_of[sample] = characteristic(header, row, "ne group")
        patient_of[sample] = characteristic(header, row, "patient")
    rows = []
    for sample in sample_ids:
        raw = site_of.get(sample, "")
        ne = ne_of.get(sample, "")
        site = SITE_74685.get(raw, "other" if raw else "")
        if "chga positive" in ne.casefold():
            organ, rule = "exclude", "1"
        elif raw == "":
            organ, site, rule = "exclude", "", "3"
        else:
            organ, rule = "Prostate", "cohort"
        rows.append({
            "sample_id": sample, "patient_id": patient_of.get(sample, sample),
            "raw_site": raw, "standard_site": site, "organ": organ, "organ_rule": rule,
            "ne_group": ne, "library": "GPL15659 gProcessedSignal", "lcm": False,
        })
    return pd.DataFrame(rows)


def parse_genepix(raw: bytes):
    lines = raw.decode("latin1").splitlines()
    header = lines[0].split("\t")
    i_name = header.index("NAME")
    i_sym = header.index("Gene Symbol")
    i_val = header.index("CH2I_MEAN")
    probes, symbols, values = [], [], []
    for line in lines[1:]:
        if not line:
            continue
        parts = line.split("\t")
        probes.append(parts[i_name])
        symbols.append(parts[i_sym] if i_sym < len(parts) else "")
        token = parts[i_val] if i_val < len(parts) else ""
        try:
            values.append(float(token))
        except ValueError:
            values.append(np.nan)
    return probes, symbols, np.asarray(values, dtype=np.float32)


def parse_agilent(raw: bytes, symbol_of: dict[str, str]):
    lines = raw.decode("latin1").splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("FEATURES\t"))
    header = lines[start].split("\t")
    i_probe = header.index("ProbeName")
    i_val = header.index("gProcessedSignal")
    i_ctrl = header.index("ControlType")
    probes, symbols, values = [], [], []
    for line in lines[start + 1:]:
        if not line or line.startswith("*") or line.startswith("TYPE\t"):
            if probes:
                break
            continue
        if not line.startswith("DATA\t"):
            continue
        parts = line.split("\t")
        if i_ctrl < len(parts) and parts[i_ctrl] not in {"0", "FALSE"}:
            continue
        probe = parts[i_probe]
        probes.append(probe)
        symbols.append(symbol_of.get(probe, ""))
        token = parts[i_val] if i_val < len(parts) else ""
        try:
            values.append(float(token))
        except ValueError:
            values.append(np.nan)
    return probes, symbols, np.asarray(values, dtype=np.float32)


def load_gpl_symbols(path: Path) -> dict[str, str]:
    symbol = {}
    with path.open(errors="replace") as handle:
        header = None
        for line in handle:
            if line.startswith("ID\t"):
                header = line.rstrip("\n").split("\t")
                break
        if header is None:
            raise SystemExit(f"no table header in {path}")
        i_id = header.index("ID")
        i_sym = header.index("GENE_SYMBOL")
        for line in handle:
            if line.startswith("!platform_table_end"):
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_id, i_sym):
                continue
            symbol[parts[i_id]] = parts[i_sym]
    return symbol


def read_tar(tar_path: Path, wanted: set[str], parser) -> tuple[list[str], list[str], np.ndarray, list[str]]:
    matrix = None
    probe_ids = None
    symbols = None
    position = None
    found = []
    with tarfile.open(tar_path) as tar:
        for member in tar:
            name = Path(member.name).name
            if not name.startswith("GSM") or not name.endswith(".txt.gz"):
                continue
            gsm = name.split("_", 1)[0]
            if gsm not in wanted:
                continue
            raw = gzip.decompress(tar.extractfile(member).read())
            probes, syms, values = parser(raw)
            if matrix is None:
                probe_ids = probes
                symbols = syms
                position = {probe: i for i, probe in enumerate(probes)}
                matrix = np.full((len(probes), len(wanted)), np.nan, dtype=np.float32)
            column = len(found)
            if probes == probe_ids:
                matrix[:, column] = values
            else:
                for probe, value in zip(probes, values):
                    hit = position.get(probe)
                    if hit is not None:
                        matrix[hit, column] = value
            found.append(gsm)
            if len(found) % 40 == 0:
                print("read", tar_path.name, len(found), flush=True)
    if matrix is None:
        raise SystemExit(f"no sample files in {tar_path}")
    return probe_ids, symbols, matrix[:, :len(found)], found


def collapse(symbols: list[str], values: np.ndarray, sample_ids: list[str]) -> pd.DataFrame:
    finite = np.isfinite(values)
    counts = finite.sum(axis=1)
    totals = np.where(finite, values, 0).sum(axis=1)
    means = np.divide(totals, counts, out=np.full(len(symbols), np.nan), where=counts > 0)
    best: dict[str, tuple[float, int]] = {}
    dropped_multi = 0
    for i, raw in enumerate(symbols):
        symbol = single_symbol(raw)
        if not symbol:
            if str(raw).strip() not in {"", "---", "NA", "null", "None", "nan"}:
                dropped_multi += 1
            continue
        score = float(means[i])
        if not np.isfinite(score):
            continue
        prev = best.get(symbol)
        if prev is None or score > prev[0]:
            best[symbol] = (score, i)
    genes = list(best)
    chosen = np.vstack([values[best[gene][1]] for gene in genes])
    frame = pd.DataFrame(chosen, index=genes, columns=sample_ids)
    frame.attrs["n_multi_gene_probes_dropped"] = dropped_multi
    frame.attrs["n_symbols"] = len(genes)
    return frame


def append_sites(path: Path, source: str, mapping: dict[str, str]) -> None:
    frame = pd.read_csv(path, sep="\t", dtype=str) if path.exists() else pd.DataFrame(columns=["source", "raw_value", "standard_site"])
    rows = [{"source": source, "raw_value": raw, "standard_site": standard} for raw, standard in mapping.items()]
    out = pd.concat([frame, pd.DataFrame(rows)], ignore_index=True)
    out = out.drop_duplicates(["source", "raw_value"], keep="last")
    out.to_csv(path, sep="\t", index=False)


def append_log(log_path: Path, acc: str, n: int, risk: int, memo: str) -> None:
    frame = pd.read_csv(log_path, sep="\t", dtype=str)
    note = "최종 판정: " + memo
    if (frame["memo"] == note).any():
        return
    base = frame.loc[frame["id"] == acc]
    title = base.iloc[0]["title"] if len(base) else ""
    row = {
        "path": "GEO", "query": "starter cohort, not limited to the title search",
        "id": acc, "title": title, "n": str(n), "platform": "microarray two-channel",
        "decision": "포함", "failed_criterion": "", "memo": note,
    }
    pd.concat([frame, pd.DataFrame([row])], ignore_index=True).to_csv(log_path, sep="\t", index=False)


def main() -> None:
    paths = load_paths()
    downloads = Path(paths["data_processed"]) / "aux" / "downloads"
    genes, b0_idx, offsets, indices, _names = load_gene_pack(paths)
    b0 = [genes[int(i)] for i in b0_idx]
    lengths = pd.Series(dtype=float)
    refs = None
    feature_rows = []
    download_rows = []

    done_71729 = Path(paths["results"]) / "stage5" / "aux" / "cohorts" / "GSE71729" / "samples.tsv"
    if done_71729.exists():
        data_root = Path(paths["data_processed"]) / "aux" / "GSE71729"
        info = json.loads((data_root / "feature_info.json").read_text())
        feature_rows.append({
            "cohort": "GSE71729",
            "Z": aux.sha256_file(data_root / "Z.npy"),
            "K": aux.sha256_file(data_root / "K.npy"),
            "n_samples": sum(1 for _ in (data_root / "sample_ids.txt").read_text().splitlines() if _),
            **info,
        })
        print("GSE71729 already built", flush=True)
    else:
        header = geo.header_frame(Path("results/stage5/aux/geo_headers/GSE71729_header.txt"))
        wanted = set(header["geo_accession"].astype(str))
        print("GSE71729 files", len(wanted), flush=True)
        _probes, symbols, values, found = read_tar(
            downloads / "GSE71729_RAW.tar", wanted, parse_genepix,
        )
        matrix = collapse(symbols, values, found)
        print("GSE71729 symbols", matrix.attrs["n_symbols"], "multi", matrix.attrs["n_multi_gene_probes_dropped"], flush=True)
        samples = samples_71729(header, found)
        feature_rows.append(aux.finish_cohort(
            "GSE71729", samples, matrix, "microarray", lengths, genes, b0, offsets, indices, paths,
            "microarray",
            "two-channel; ch1 Human Reference unused; sample channel CH2I_MEAN; highest-mean probe per symbol; deposited intensities rescaled to sum 1e6",
            refs,
        ))
    append_sites(Path(paths["config"]) / "mappings" / "site_dictionary_aux.tsv", "GSE71729", {**SITE_71729, "CellLine": "", **{f"{name}_Normal": "" for name in ["Pancreas", "Liver", "Lung", "Spleen", "LymphNode", "Diaphragm", "Fat", "Peritoneal", "PelvicWall", "Vessel"]}})

    symbol_of = load_gpl_symbols(downloads / "GPL15659_data.txt")
    header = geo.header_frame(Path("results/stage5/aux/geo_headers/GSE74685_header.txt"))
    wanted = set(header["geo_accession"].astype(str))
    print("GSE74685 files", len(wanted), "gpl symbols", len(symbol_of), flush=True)

    def parser(raw, table=symbol_of):
        return parse_agilent(raw, table)

    _probes, symbols, values, found = read_tar(downloads / "GSE74685_RAW.tar", wanted, parser)
    matrix = collapse(symbols, values, found)
    print("GSE74685 symbols", matrix.attrs["n_symbols"], "multi", matrix.attrs["n_multi_gene_probes_dropped"], flush=True)
    samples = samples_74685(header, found)
    feature_rows.append(aux.finish_cohort(
        "GSE74685", samples, matrix, "microarray", lengths, genes, b0, offsets, indices, paths,
        "microarray",
        "two-channel; ch1 cell-line reference pool unused; sample channel gProcessedSignal; GPL15659 GENE_SYMBOL; highest-mean probe per symbol; deposited intensities rescaled to sum 1e6",
        refs,
    ))
    append_sites(Path(paths["config"]) / "mappings" / "site_dictionary_aux.tsv", "GSE74685", SITE_74685)

    log_path = Path(paths["results"]) / "stage5" / "aux" / "screening_log.tsv"
    for acc in ("GSE71729", "GSE74685"):
        table = pd.read_csv(Path(paths["results"]) / "stage5" / "aux" / "cohorts" / acc / "samples.tsv", sep="\t")
        append_log(log_path, acc, len(table), int(table["in_risk"].sum()), f"n={len(table)} n_at_risk={int(table['in_risk'].sum())} after sample-channel extraction")
        tar = downloads / f"{acc}_RAW.tar"
        url = {
            "GSE71729": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71729/suppl/GSE71729_RAW.tar",
            "GSE74685": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74685/suppl/GSE74685_RAW.tar",
        }[acc]
        download_rows.append({
            "cohort": acc, "file": tar.name,
            "url": url,
            "bytes": tar.stat().st_size, "sha256": aux.sha256_file(tar),
            "license": "",
        })
    download_rows.append({
        "cohort": "GSE74685", "file": "GPL15659_data.txt",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15659&targ=self&view=data&form=text",
        "bytes": (downloads / "GPL15659_data.txt").stat().st_size,
        "sha256": aux.sha256_file(downloads / "GPL15659_data.txt"),
        "license": "",
    })
    existing = pd.read_csv(Path(paths["config"]) / "aux_downloads_A5.tsv", sep="\t", dtype=str)
    pd.concat([existing, pd.DataFrame(download_rows)], ignore_index=True).drop_duplicates(["cohort", "file"], keep="last").to_csv(
        Path(paths["config"]) / "aux_downloads_A5.tsv", sep="\t", index=False,
    )
    features = pd.read_csv(Path(paths["config"]) / "aux_features_A5.tsv", sep="\t", dtype=str)
    pd.concat([features, pd.DataFrame(feature_rows).astype(str)], ignore_index=True).drop_duplicates("cohort", keep="last").to_csv(
        Path(paths["config"]) / "aux_features_A5.tsv", sep="\t", index=False,
    )
    Path(paths["results"]).joinpath("stage5/aux/twochannel_feature_rows.json").write_text(json.dumps(feature_rows, indent=2))
    print("two-channel cohorts", len(feature_rows), flush=True)


if __name__ == "__main__":
    main()
