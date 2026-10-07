#!/usr/bin/env python3
"""Build features for the GEO starter cohorts whose processed matrices are already downloaded.

Does not apply a classifier. Two-channel archives are handled only when the tar is complete.
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
spec = importlib.util.spec_from_file_location("auxbuild", ROOT / "19_aux_build.py")
aux = importlib.util.module_from_spec(spec)
spec.loader.exec_module(aux)

sys.path.insert(0, str(ROOT))
from stage2_common import load_gene_pack
from stage5_lib import load_paths


def header_frame(path: Path) -> pd.DataFrame:
    rows = {}
    char_i = 0
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("!Sample_"):
            continue
        parts = [part.strip().strip('"') for part in line.split("\t")]
        key = parts[0][len("!Sample_"):]
        values = parts[1:]
        if key in {"characteristics_ch1", "characteristics_ch2", "description"}:
            label = values[0].split(":", 1)[0].strip() if values and ":" in values[0] else key
            col = f"{key}__{char_i}__{label}"
            char_i += 1
        else:
            col = key
            if col in rows:
                col = f"{key}__dup"
        rows[col] = values
    n = max(len(v) for v in rows.values())
    frame = pd.DataFrame({col: vals for col, vals in rows.items() if len(vals) == n})
    return frame


def gse50760(downloads: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    with tarfile.open(downloads / "GSE50760_RAW.tar") as tar:
        for member in tar.getmembers():
            if not member.name.endswith(".gz"):
                continue
            raw = gzip.decompress(tar.extractfile(member).read()).decode().splitlines()
            header = raw[0].split("\t")
            sample = header[1]
            symbols, values = [], []
            for line in raw[1:]:
                gene, value = line.split("\t")
                symbols.append(gene)
                values.append(float(value))
            frames.append(pd.Series(values, index=symbols, name=sample))
    matrix = pd.concat(frames, axis=1).fillna(0.0).groupby(level=0).sum()
    header = header_frame(Path("results/stage5/aux/geo_headers/GSE50760_header.txt"))
    # title is like GSM... and characteristics tissue
    tissue_col = [c for c in header.columns if c.endswith("__tissue")][0]
    title = header["title"].astype(str)
    tissue_values = header[tissue_col].astype(str).str.split(":").str[-1].str.strip()
    samples = []
    for sample in matrix.columns:
        token = sample.replace("_FPKM", "").replace(".", "-")
        hit = title.str.contains(token, regex=False)
        tissue = tissue_values.loc[hit].iloc[0] if hit.any() else ""
        patient = sample.replace("_FPKM", "").rsplit(".", 1)[0]
        site = {
            "primary colorectal cancer": "colon_rectum",
            "metastatic colorectal cancer to the liver": "liver",
        }.get(tissue, "")
        organ = "exclude" if "normal" in tissue else "Colorectal"
        rule = "8" if organ == "exclude" else "cohort"
        samples.append({
            "sample_id": sample, "patient_id": patient, "raw_site": tissue,
            "standard_site": site, "organ": organ, "organ_rule": rule, "library": "FPKM",
        })
    return pd.DataFrame(samples), matrix


def gse209998(downloads: Path, lengths) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = downloads / "GSE209998_AUR_129_raw_counts.txt.gz"
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    frame = frame.groupby(level=0).sum()
    header = header_frame(Path("results/stage5/aux/geo_headers/GSE209998_header.txt"))
    desc = header[[c for c in header.columns if c.startswith("description__")][1]] if any(c.startswith("description__") for c in header.columns) else header["description"]
    # the barcode description line is the one whose values contain AUR-
    barcode = None
    for col in header.columns:
        if header[col].astype(str).str.startswith("AUR-").mean() > 0.5:
            barcode = header[col].astype(str)
            break
    tissue_col = [c for c in header.columns if c.endswith("__tissue")][0]
    disease_col = [c for c in header.columns if c.endswith("__disease")][0]
    meta = pd.DataFrame({
        "barcode": barcode.to_numpy(),
        "tissue": header[tissue_col].astype(str).str.split(":").str[-1].str.strip().to_numpy(),
        "disease": header[disease_col].astype(str).str.split(":").str[-1].str.strip().to_numpy(),
        "gsm": header["geo_accession"].astype(str).to_numpy(),
    })
    site_map = {
        "Breast": "breast", "Liver": "liver", "Lymph node": "lymph_node", "Brain": "brain",
        "Lung": "lung", "Chest": "other", "Soft tissue": "soft_tissue", "Adrenal": "adrenal",
        "Pleura": "pleura", "Skin": "skin", "Ovary": "ovary", "Bone": "bone",
    }
    samples = []
    keep = []
    for sample in frame.columns:
        hit = meta.loc[meta["barcode"] == sample]
        if hit.empty:
            continue
        row = hit.iloc[0]
        tissue = row["tissue"]
        organ = "exclude" if row["disease"] == "Normal tissue" else "Breast"
        samples.append({
            "sample_id": sample, "patient_id": sample,
            "raw_site": tissue, "standard_site": "" if organ == "exclude" else site_map.get(tissue, "other"),
            "organ": organ, "organ_rule": "8" if organ == "exclude" else "cohort",
            "library": "raw_counts", "gsm": row["gsm"], "disease": row["disease"],
        })
        keep.append(sample)
    return pd.DataFrame(samples), frame[keep]


def gpl96_annot(downloads: Path) -> pd.DataFrame:
    path = downloads / "GPL96.annot.gz"
    if not path.exists():
        aux.download_to("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz", path)
    rows = []
    with gzip.open(path, "rt", errors="replace") as handle:
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


def series_matrix(path: Path, annot: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
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
    # highest mean deposited value per symbol
    means = merged.drop(columns=["symbol"]).mean(axis=1)
    merged = merged.assign(_mean=means)
    chosen = merged.sort_values("_mean", ascending=False).groupby("symbol", sort=False).head(1)
    matrix = chosen.drop(columns=["_mean", "symbol"])
    matrix.index = chosen["symbol"].to_numpy()
    return matrix, header[1:]


def from_titles_gse14018(header: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    title = dict(zip(header["geo_accession"], header["title"]))
    source = dict(zip(header["geo_accession"], header["source_name_ch1"]))
    rows = []
    for sample in sample_ids:
        raw = str(source.get(sample, ""))
        site = {"Lung": "lung", "Liver": "liver", "Brain": "brain", "Bone": "bone"}.get(raw, "other" if raw else "")
        rows.append({
            "sample_id": sample, "patient_id": sample, "raw_site": raw, "standard_site": site,
            "organ": "Breast", "organ_rule": "cohort", "library": "GPL96", "title": title.get(sample, ""),
        })
    return pd.DataFrame(rows)


def gse41258_samples(header: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    tissue_col = [c for c in header.columns if c.endswith("__tissue")][0]
    patient_col = [c for c in header.columns if "patient id" in c][0]
    tissue = dict(zip(header["geo_accession"], header[tissue_col].astype(str)))
    patient = dict(zip(header["geo_accession"], header[patient_col].astype(str).str.split(":").str[-1].str.strip()))
    site_of = {
        "Primary Tumor": "colon_rectum",
        "Liver Metastasis": "liver",
        "Lung Metastasis": "lung",
    }
    exclude_tissue = {"Normal Colon", "Normal Liver", "Normal Lung", "Polyp", "Polyp, high grade", "Microadenoma"}
    rows = []
    for sample in sample_ids:
        full = tissue.get(sample, "")
        label, _, rest = full.partition(":")
        raw = rest.strip()
        if label.strip() == "cell line":
            organ, site, rule = "exclude", "", "2"
        elif raw in exclude_tissue or raw == "":
            organ, site, rule = "exclude", "", "8"
        else:
            organ, site, rule = "Colorectal", site_of.get(raw, "other"), "cohort"
        rows.append({
            "sample_id": sample, "patient_id": patient.get(sample, sample), "raw_site": raw,
            "standard_site": site, "organ": organ, "organ_rule": rule, "library": "GPL96",
        })
    return pd.DataFrame(rows)


def main() -> None:
    paths = load_paths()
    downloads = Path(paths["data_processed"]) / "aux" / "downloads"
    gtf = downloads / "gencode.v23.annotation.gtf.gz"
    if not gtf.exists():
        raise SystemExit("gencode gtf missing; run 19_aux_build.py first")
    lengths = aux.gene_lengths(gtf, downloads / "gencode_v23_exon_union.tsv")
    genes, b0_idx, offsets, indices, _names = load_gene_pack(paths)
    b0 = [genes[int(i)] for i in b0_idx]
    print("reference", flush=True)
    refs = aux.load_reference_tpm(paths, genes)
    download_rows = []
    feature_rows = []

    samples, matrix = gse50760(downloads)
    feature_rows.append(aux.finish_cohort(
        "GSE50760", samples, matrix, "fpkm", lengths, genes, b0, offsets, indices, paths,
        "RNA-seq unspecified", "per-sample FPKM", refs,
    ))
    samples, matrix = gse209998(downloads, lengths)
    feature_rows.append(aux.finish_cohort(
        "GSE209998", samples, matrix, "counts", lengths, genes, b0, offsets, indices, paths,
        "RNA-seq unspecified", "GSE209998_AUR_129_raw_counts.txt.gz fractional counts, exon-union length", refs,
    ))
    annot = gpl96_annot(downloads)
    for acc, builder, kind in (
        ("GSE14018", from_titles_gse14018, "microarray"),
        ("GSE41258", gse41258_samples, "microarray"),
    ):
        header = header_frame(Path(f"results/stage5/aux/geo_headers/{acc}_header.txt"))
        matrix, ids = series_matrix(downloads / f"{acc}_series_matrix.txt.gz", annot)
        samples = builder(header, ids)
        feature_rows.append(aux.finish_cohort(
            acc, samples, matrix, kind, lengths, genes, b0, offsets, indices, paths,
            "microarray", "GEO series matrix deposited values; probe with the highest mean kept", refs,
        ))
    out = Path(paths["results"]) / "stage5" / "aux" / "geo_feature_rows.json"
    out.write_text(json.dumps(feature_rows, indent=2))
    print("geo cohorts", len(feature_rows), flush=True)


if __name__ == "__main__":
    main()
