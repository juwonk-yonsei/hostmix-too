#!/usr/bin/env python3
"""Collect included auxiliary cohorts and build Z/K features. No model is applied."""
from __future__ import annotations

import gzip
import hashlib
import json
import re
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from stage2_common import load_gene_pack, rank_z, score_ks
from stage3_rules import native_organs
from stage5_lib import load_paths, load_tpm, sum_1e6

from loader import load_pog570_expression
from genes import probemap_symbol

UA = {"User-Agent": "projA-aux-build"}
DATAHUB = "https://raw.githubusercontent.com/cBioPortal/datahub/master/public"
LFS = "https://github.com/cBioPortal/datahub.git/info/lfs/objects/batch"

# Latest decision for studies already in the first screen. GEO starters are added below.
FINAL = {
    "aml_ohsu_2018": ("제외", "2", "every SAMPLE_SITE is Bone Marrow Aspirate, Peripheral Blood, or Leukapheresis"),
    "asclc_msk_2024": ("제외", "4", "cohort title is atypical small cell; rule 1 excludes the cohort"),
    "brain_cptac_2020": ("제외", "6", "collection sites are intracranial or spinal; truth organ Brain; n_at_risk 0"),
    "gbm_cptac_2021": ("제외", "6", "TUMOR_SITE_CURATED values are brain lobes; n_at_risk 0"),
    "gbm_iatlas_prins_2019": ("제외", "6", "BIOPSY_SITE is Brain for all 30 samples"),
    "lgg_ctf_synodos_2025": ("제외", "6", "BIOPSY_SITE values are brain regions"),
    "mel_iatlas_gide_2019": ("제외", "6", "biopsy site filled for 5 samples; n_at_risk below 10"),
    "mel_iatlas_hugo_ucla_2016": ("제외", "6", "non-missing non-native sites are fewer than 10"),
    "mel_iatlas_liu_2019": ("제외", "7", "same DFCI Nature Medicine 2019 melanoma study as mel_dfci_2019"),
    "mel_tsam_liang_2017": ("제외", "6", "metastases with a non-native site: lymph node 5 and lung 1"),
    "nbl_target_2018_pub": ("제외", "3", "TUMOR_TISSUE_SITE is empty for 1088 of 1089 samples; relapse site is not the sequenced sample"),
    "nepc_wcm_2016": ("제외", "4", "cohort is neuroendocrine prostate cancer; rule 1"),
    "schw_ctf_synodos_2025": ("제외", "4", "schwannoma is not an organ in the stage-3 table; TISSUE_SITE is Nerves"),
    "sclc_ucologne_2015": ("제외", "4", "small cell lung carcinoma; rule 1"),
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def get_bytes(url: str, timeout: int = 180, data: bytes | None = None, headers: dict | None = None) -> bytes:
    last = None
    for attempt in range(4):
        request = urllib.request.Request(url, data=data, headers={**UA, **(headers or {})})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as handle:
                return handle.read()
        except Exception as exc:
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise last


def lfs_download(pointer: bytes) -> bytes:
    text = pointer.decode()
    oid = re.search(r"oid sha256:([0-9a-f]+)", text).group(1)
    size = int(re.search(r"size (\d+)", text).group(1))
    body = json.dumps({"operation": "download", "transfers": ["basic"], "objects": [{"oid": oid, "size": size}]}).encode()
    raw = get_bytes(LFS, data=body, headers={"Accept": "application/vnd.git-lfs+json", "Content-Type": "application/vnd.git-lfs+json"})
    href = json.loads(raw)["objects"][0]["actions"]["download"]["href"]
    return get_bytes(href, timeout=600)


def datahub_file(study: str, name: str) -> bytes:
    raw = get_bytes(f"{DATAHUB}/{study}/{name}")
    if raw.startswith(b"version https://git-lfs"):
        raw = lfs_download(raw)
    return raw


def download_to(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    path.write_bytes(get_bytes(url, timeout=600))


def append_decisions(log_path: Path) -> None:
    frame = pd.read_csv(log_path, sep="\t", dtype=str)
    rows = []
    for study, (decision, failed, memo) in FINAL.items():
        note = "최종 판정: " + memo
        if (frame["memo"] == note).any():
            continue
        base = frame.loc[frame["id"] == study].iloc[0]
        rows.append({
            "path": base["path"], "query": base["query"], "id": study, "title": base["title"],
            "n": base["n"], "platform": base["platform"], "decision": decision,
            "failed_criterion": failed, "memo": note,
        })
    starters = [
        ("GSE209998", "포함", "", "processed raw-count file on GEO FTP; n=129 so salmon is not used"),
        ("GSE50760", "포함", "", "per-sample FPKM files in GSE50760_RAW.tar"),
        ("GSE41258", "포함", "", "GPL96 series matrix; polyp, normal, and cell-line samples removed at labeling"),
        ("GSE14018", "포함", "", "GPL96 series matrix; site taken from the sample title"),
        ("GSE71729", "포함", "", "two-channel; sample is ch2 Cy5; CH2I_MEAN is present in the raw file"),
        ("GSE74685", "포함", "", "two-channel; sample is ch2 Cy3; gProcessedSignal is present in the raw file"),
    ]
    for acc, decision, failed, memo in starters:
        note = "출발 후보 확인: " + memo
        if (frame["id"] == acc).any() and (frame["memo"] == note).any():
            continue
        if (frame["memo"] == note).any():
            continue
        rows.append({
            "path": "GEO", "query": "starter cohort, not limited to the title search",
            "id": acc, "title": "", "n": "", "platform": "",
            "decision": decision, "failed_criterion": failed, "memo": note,
        })
    if rows:
        pd.concat([frame, pd.DataFrame(rows)], ignore_index=True).to_csv(log_path, sep="\t", index=False)
    print("decision rows", len(rows), flush=True)


def read_cbio_clinical(path: Path) -> pd.DataFrame:
    lines = path.read_text(errors="replace").splitlines()
    header = next(i for i, line in enumerate(lines) if "SAMPLE_ID" in line and not line.startswith("#"))
    names = lines[header].split("\t")
    rows = [line.split("\t") for line in lines[header + 1:] if line]
    width = len(names)
    rows = [row + [""] * (width - len(row)) if len(row) < width else row[:width] for row in rows]
    return pd.DataFrame(rows, columns=names)


def read_cbio_expression(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    symbol_col = "Hugo_Symbol" if "Hugo_Symbol" in frame.columns else frame.columns[0]
    frame = frame.loc[frame[symbol_col].notna() & (frame[symbol_col] != "") & (frame[symbol_col] != "NA")]
    sample_cols = [col for col in frame.columns if col not in {symbol_col, "Entrez_Gene_Id"}]
    values = frame[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    values.index = frame[symbol_col].astype(str).to_numpy()
    grouped = values.groupby(level=0).sum()
    return grouped


def gene_lengths(gtf: Path, cache: Path) -> pd.Series:
    if cache.exists():
        table = pd.read_csv(cache, sep="\t")
        return table.set_index("symbol")["exon_union_bp"]
    intervals: dict[str, dict[str, list[tuple[int, int]]]] = {}
    names: dict[str, set[str]] = {}
    with gzip.open(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            attrs = dict(re.findall(r'(\S+) "([^"]+)"', parts[8]))
            gene = attrs.get("gene_id", "").split(".")[0]
            symbol = attrs.get("gene_name", "")
            if not gene or not symbol:
                continue
            intervals.setdefault(gene, {}).setdefault(parts[0], []).append((int(parts[3]), int(parts[4])))
            names.setdefault(symbol, set()).add(gene)
    length = {}
    for gene, chroms in intervals.items():
        total = 0
        for spans in chroms.values():
            spans = sorted(spans)
            cur_s, cur_e = spans[0]
            for start, end in spans[1:]:
                if start <= cur_e + 1:
                    cur_e = max(cur_e, end)
                else:
                    total += cur_e - cur_s + 1
                    cur_s, cur_e = start, end
            total += cur_e - cur_s + 1
        length[gene] = total
    rows = []
    for symbol, genes in names.items():
        if len(genes) != 1:
            continue
        gene = next(iter(genes))
        rows.append({"symbol": symbol, "gene_id": gene, "exon_union_bp": length[gene]})
    table = pd.DataFrame(rows)
    table.to_csv(cache, sep="\t", index=False)
    return table.set_index("symbol")["exon_union_bp"]


def to_tpm_from_counts(counts: pd.DataFrame, lengths: pd.Series) -> tuple[pd.DataFrame, int]:
    common = counts.index.intersection(lengths.index)
    missing = int(len(counts.index.difference(lengths.index)))
    rate = counts.loc[common].div(lengths.loc[common], axis=0)
    denom = rate.sum(axis=0).replace(0, np.nan)
    return rate.div(denom, axis=1) * 1e6, missing


def align_g(matrix: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray, list[str], list[str]]:
    """Return samples x present-genes, present symbols, and the sample ids."""
    present = [gene for gene in genes if gene in matrix.index]
    block = matrix.loc[present].to_numpy(dtype=np.float64).T
    block = np.nan_to_num(block, nan=0.0)
    return block, present, list(matrix.columns)


def z_and_k(block: np.ndarray, present: list[str], genes: list[str], b0: list[str], offsets, indices) -> tuple[np.ndarray, np.ndarray, dict]:
    gene_index = {gene: i for i, gene in enumerate(genes)}
    present_index = {gene: i for i, gene in enumerate(present)}
    z_all = rank_z(block)
    b0_cols = []
    z = np.zeros((block.shape[0], len(b0)), dtype=np.float32)
    for j, symbol in enumerate(b0):
        if symbol in present_index:
            z[:, j] = z_all[:, present_index[symbol]].astype(np.float32)
            b0_cols.append(symbol)
    lists = []
    set_ok = []
    n_short = 0
    n_sets = int(offsets.shape[0] - 1)
    for g in range(n_sets):
        members = []
        for gene_i in indices[offsets[g]:offsets[g + 1]]:
            symbol = genes[int(gene_i)]
            if symbol in present_index:
                members.append(present_index[symbol])
        if len(members) < 5:
            n_short += 1
            set_ok.append(False)
            lists.append(np.array([0], dtype=np.int32))
        else:
            set_ok.append(True)
            lists.append(np.array(members, dtype=np.int32))
    # score only the valid sets, then scatter
    valid = [members for members, ok in zip(lists, set_ok) if ok]
    packed_off, packed_idx = pack_offsets(valid)
    scored = score_ks(block.astype(np.float32), packed_off, packed_idx)
    k = np.zeros((block.shape[0], n_sets), dtype=np.float32)
    cursor = 0
    for g, ok in enumerate(set_ok):
        if ok:
            k[:, g] = scored[:, cursor]
            cursor += 1
    info = {
        "n_G_intersect": len(present),
        "n_B0_intersect": len(b0_cols),
        "n_ks_sets_zero": int(n_short),
        "n_ks_sets": n_sets,
    }
    del gene_index
    return z, k, info


def pack_offsets(member_lists):
    offsets = np.zeros(len(member_lists) + 1, dtype=np.int64)
    chunks = []
    for i, members in enumerate(member_lists):
        chunks.append(np.asarray(members, dtype=np.int32))
        offsets[i + 1] = offsets[i] + chunks[-1].size
    indices = np.concatenate(chunks) if chunks else np.empty(0, dtype=np.int32)
    return offsets, indices


def set_counts(samples: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for site, sub in samples.groupby("standard_site", dropna=False):
        for kind, mask in (
            ("eval", sub["in_eval"]),
            ("risk", sub["in_risk"]),
            ("native", sub["in_native"]),
        ):
            part = sub.loc[mask]
            rows.append({
                "standard_site": site, "set": kind,
                "n_samples": int(len(part)),
                "n_patients": int(part["patient_id"].nunique()) if len(part) else 0,
            })
    return pd.DataFrame(rows)


def mark_sets(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    natives = [native_organs(site) for site in out["standard_site"]]
    out["in_eval"] = out["organ"].ne("exclude") & out["organ"].ne("NA") & out["standard_site"].ne("")
    out["in_native"] = [bool(out["in_eval"].iloc[i] and out["organ"].iloc[i] in natives[i]) for i in range(len(out))]
    out["in_risk"] = [bool(out["in_eval"].iloc[i] and natives[i] and out["organ"].iloc[i] not in natives[i]) for i in range(len(out))]
    return out


def spearman_max(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    def norm_ranks(block):
        ranks = rankdata(block, method="average", axis=1)
        ranks = ranks - ranks.mean(axis=1, keepdims=True)
        denom = np.linalg.norm(ranks, axis=1, keepdims=True)
        denom[denom == 0] = 1
        return ranks / denom
    return norm_ranks(left) @ norm_ranks(right).T


def load_reference_tpm(paths, genes: list[str]):
    met_samples, met = load_tpm(paths["data_processed"] / "met500/tpm_G.h5")
    pog = load_pog570_expression(unlock=True)
    mapping = probemap_symbol()
    by_base: dict[str, set[str]] = {}
    for ens, sym in zip(mapping["id"].astype(str), mapping["gene"].astype(str)):
        by_base.setdefault(ens.split(".")[0], set()).add(sym)
    base_to_symbol = {base: next(iter(symbols)) for base, symbols in by_base.items() if len(symbols) == 1}
    bases = pd.Index(pog.index.astype(str)).str.split(".").str[0]
    symbol = np.asarray(bases.map(base_to_symbol), dtype=object)
    keep = pd.notna(symbol)
    grouped = pog.loc[keep].astype(np.float32).groupby(symbol[keep], sort=False).sum()
    del pog
    pog_ids = [str(col) for col in grouped.columns]
    pog_block = grouped.reindex(genes).fillna(0).to_numpy(dtype=np.float32).T
    del grouped
    return list(met_samples), sum_1e6(met), pog_ids, sum_1e6(pog_block)


def duplicate_screen(block, present, patients, sample_ids, cohort, genes, refs, out_dir: Path) -> pd.DataFrame:
    """RNA-seq only. Returns one row per auxiliary sample with the external maximum."""
    gene_index = {gene: i for i, gene in enumerate(genes)}
    take = [gene_index[gene] for gene in present]
    met_samples, met, pog_ids, pog_block = refs
    met = met[:, take]
    pog_block = pog_block[:, take]
    aux = np.log2(block + 1.0)
    met_l = np.log2(met + 1.0)
    pog_l = np.log2(pog_block + 1.0)
    ext_met = spearman_max(aux, met_l)
    ext_pog = spearman_max(aux, pog_l)
    internal = spearman_max(aux, aux)
    patient = np.asarray(patients)
    for i in range(len(patient)):
        internal[i, patient == patient[i]] = -np.inf
    ext_max = np.maximum(ext_met.max(axis=1), ext_pog.max(axis=1))
    int_max = np.where(np.isfinite(internal).any(axis=1), np.max(np.where(np.isfinite(internal), internal, -np.inf), axis=1), np.nan)
    which_met = ext_met.max(axis=1) >= ext_pog.max(axis=1)
    partner = []
    partner_cohort = []
    partner_r = []
    for i in range(len(sample_ids)):
        if which_met[i]:
            j = int(np.argmax(ext_met[i]))
            partner.append(met_samples[j])
            partner_cohort.append("MET500")
            partner_r.append(float(ext_met[i, j]))
        else:
            j = int(np.argmax(ext_pog[i]))
            partner.append(pog_ids[j])
            partner_cohort.append("POG570")
            partner_r.append(float(ext_pog[i, j]))
    flag = (ext_max >= 0.98) | ((ext_max >= 0.95) & (ext_max > int_max))
    table = pd.DataFrame({
        "cohort": cohort, "sample_id": sample_ids, "patient_id": patients,
        "external_max": ext_max, "internal_max": int_max,
        "flag_candidate": flag,
    })
    pairs = table.loc[flag, ["cohort", "sample_id", "patient_id"]].copy()
    pairs["partner_cohort"] = [partner_cohort[i] for i in np.flatnonzero(flag)]
    pairs["partner_id"] = [partner[i] for i in np.flatnonzero(flag)]
    pairs["correlation"] = [partner_r[i] for i in np.flatnonzero(flag)]
    pairs.to_csv(out_dir / "duplicate_pairs.tsv", sep="\t", index=False)
    quant = []
    for column in ("external_max", "internal_max"):
        values = table[column].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if len(values) == 0:
            continue
        qs = np.percentile(values, [0, 1, 50, 99, 100])
        quant.append({"cohort": cohort, "which": column, "min": qs[0], "p01": qs[1], "p50": qs[2], "p99": qs[3], "max": qs[4], "n": len(values), "n_flag_098": int((ext_max >= 0.98).sum()), "n_flag_095": int(((ext_max >= 0.95) & (ext_max > int_max)).sum())})
    pd.DataFrame(quant).to_csv(out_dir / "duplicate_quantiles.tsv", sep="\t", index=False)
    return table


def save_features(cohort_dir: Path, sample_ids, z, k, info) -> dict:
    cohort_dir.mkdir(parents=True, exist_ok=True)
    z_path = cohort_dir / "Z.npy"
    k_path = cohort_dir / "K.npy"
    id_path = cohort_dir / "sample_ids.txt"
    np.save(z_path, z)
    np.save(k_path, k)
    id_path.write_text("\n".join(sample_ids) + "\n")
    info_path = cohort_dir / "feature_info.json"
    info_path.write_text(json.dumps(info))
    return {
        "Z": sha256_file(z_path),
        "K": sha256_file(k_path),
        "n_samples": len(sample_ids),
        **info,
    }


def series_table(header_path: Path) -> pd.DataFrame:
    buckets: dict[str, list[str]] = {}
    for line in header_path.read_text(errors="replace").splitlines():
        if not line.startswith("!Sample_"):
            continue
        parts = line.split("\t")
        key = parts[0][len("!Sample_"):]
        values = [part.strip().strip('"') for part in parts[1:]]
        buckets.setdefault(key, [])
        # repeated keys (characteristics, description) are stored with a suffix
        name = key if key not in buckets or len(buckets[key]) == 1 and buckets[key][0] == [] else key
        if key in ("characteristics_ch1", "characteristics_ch2", "description"):
            label = values[0].split(":", 1)[0].strip() if values and ":" in values[0] else key
            col = f"{key}::{label}::{len([c for c in buckets if c.startswith(key + '::' + label)])}"
            # simpler: stack by occurrence index
        buckets[key].append(values)
    n = max(len(vals) for rows in buckets.values() for vals in rows)
    frame = pd.DataFrame({"i": np.arange(n)})
    for key, rows in buckets.items():
        for j, vals in enumerate(rows):
            col = key if j == 0 else f"{key}__{j}"
            frame[col] = vals
    return frame


def finish_cohort(name, samples, matrix, kind, lengths, genes, b0, offsets, indices, paths, platform, library_note, refs) -> dict:
    """matrix is genes x samples, already on the scale named by kind (tpm or counts)."""
    out_root = Path(paths["results"]) / "stage5" / "aux" / "cohorts" / name
    data_root = Path(paths["data_processed"]) / "aux" / name
    out_root.mkdir(parents=True, exist_ok=True)
    samples = mark_sets(samples)
    length_missing = 0
    if kind == "counts":
        tpm, length_missing = to_tpm_from_counts(matrix, lengths)
    else:
        tpm = matrix.copy()
        # FPKM/RPKM/TPM are rescaled to sum 1e6. Log-valued matrices are not.
        values = tpm.to_numpy(dtype=float)
        if values.size == 0:
            raise SystemExit(f"{name} expression matrix is empty")
        if float(np.nanmin(values)) < 0:
            raise SystemExit(f"{name} has negative expression; not a linear abundance")
        totals = tpm.sum(axis=0).replace(0, np.nan)
        tpm = tpm.div(totals, axis=1) * 1e6
    block, present, sample_ids = align_g(tpm, genes)
    order = samples.set_index("sample_id").loc[sample_ids].reset_index()
    if len(order) != len(sample_ids):
        raise SystemExit(f"{name} sample alignment failed")
    z, k, info = z_and_k(block, present, genes, b0, offsets, indices)
    info["length_missing_symbols"] = length_missing
    info["platform"] = platform
    info["library_note"] = library_note
    info["kind"] = kind
    hashes = save_features(data_root, sample_ids, z, k, info)
    order.to_csv(out_root / "samples.tsv", sep="\t", index=False)
    set_counts(order).to_csv(out_root / "set_counts.tsv", sep="\t", index=False)
    if kind in {"tpm", "fpkm", "rpkm", "counts"} and platform.startswith("RNA-seq"):
        dup = duplicate_screen(block, present, order["patient_id"].tolist(), sample_ids, name, genes, refs, out_root)
        order = order.merge(dup[["sample_id", "external_max", "internal_max", "flag_candidate"]], on="sample_id", how="left")
        order.to_csv(out_root / "samples.tsv", sep="\t", index=False)
        kept = order.loc[~order["flag_candidate"].fillna(False)]
        set_counts(kept).to_csv(out_root / "set_counts_after_flag.tsv", sep="\t", index=False)
    else:
        order.to_csv(out_root / "samples.tsv", sep="\t", index=False)
    print("cohort", name, "n", len(order), "risk", int(order["in_risk"].sum()), "intersect", info["n_G_intersect"], flush=True)
    return {"cohort": name, **hashes}


MISSING_SITE = {"", "nan", "na", "n/a", "unknown", "not available", "not applicable", "none"}


def process_cbio(paths, lengths, genes, b0, offsets, indices, download_rows, refs) -> list[dict]:
    clinical_dir = Path(paths["results"]) / "stage5" / "aux" / "clinical"
    dest = Path(paths["data_processed"]) / "aux" / "downloads"
    specs = [
        ("blca_iatlas_imvigor210_2017", "data_mrna_seq_tpm.txt", "tpm", "RNA-seq unspecified", "BIOPSY_SITE", "Bladder"),
        ("brca_iatlas_anders_2022", "data_mrna_seq_tpm.txt", "tpm", "RNA-seq unspecified", "BIOPSY_SITE", "Breast"),
        ("mel_dfci_2019", "data_mrna_seq_tpm.txt", "tpm", "RNA-seq unspecified", "BIOPSY_SITE", "Melanoma"),
        ("paad_iatlas_prince_2022", "data_mrna_seq_expression.txt", "tpm", "RNA-seq unspecified", "BIOPSY_SITE", "Pancreas"),
        ("skcm_mskcc_2014", "data_mrna_seq_rpkm.txt", "rpkm", "RNA-seq unspecified", "SIMPLIFIED_TUMOR_SITE", "Melanoma"),
        ("prad_su2c_2019", "data_mrna_seq_fpkm_polya.txt", "fpkm", "RNA-seq polyA", "TISSUE_SITE", "Prostate"),
        ("prad_su2c_2019", "data_mrna_seq_fpkm_capture.txt", "fpkm", "RNA-seq capture", "TISSUE_SITE", "Prostate"),
    ]
    site_map = json.loads(Path(__file__).with_name("aux_site_map.json").read_text())
    made = []
    for study, filename, kind, platform, site_col, organ in specs:
        raw_path = dest / f"{study}__{filename}"
        if not raw_path.exists():
            print("download", study, filename, flush=True)
            raw_path.write_bytes(datahub_file(study, filename))
        license_path = dest / f"{study}__LICENSE"
        if not license_path.exists():
            license_path.write_bytes(datahub_file(study, "LICENSE"))
        download_rows.append({
            "cohort": study, "file": filename, "url": f"{DATAHUB}/{study}/{filename}",
            "bytes": raw_path.stat().st_size, "sha256": sha256_file(raw_path),
            "license": license_path.read_text(errors="replace").strip(),
        })
        clinical = read_cbio_clinical(clinical_dir / f"{study}_sample.txt")
        expr = read_cbio_expression(raw_path)
        name = study if "fpkm_" not in filename else f"{study}__{filename.replace('data_mrna_seq_', '').replace('.txt', '')}"
        samples = clinical.copy()
        samples["sample_id"] = samples["SAMPLE_ID"]
        samples["patient_id"] = samples["PATIENT_ID"] if "PATIENT_ID" in samples.columns else samples["SAMPLE_ID"]
        raw_site = samples[site_col].fillna("").astype(str)
        samples["raw_site"] = raw_site
        mapped = []
        for value in raw_site:
            token = value.strip()
            if token.casefold() in MISSING_SITE:
                mapped.append("")
            else:
                key = f"{study}|{token}"
                if key not in site_map:
                    print("unmapped site", key, flush=True)
                    mapped.append("other")
                else:
                    mapped.append(site_map[key])
        samples["standard_site"] = mapped
        samples["organ"] = organ
        samples["organ_rule"] = "cohort"
        if study == "prad_su2c_2019":
            ne = samples.get("NEUROENDOCRINE_FEATURES", pd.Series("", index=samples.index)).fillna("")
            pathol = samples.get("PATHOLOGY_CLASSIFICATION", pd.Series("", index=samples.index)).fillna("")
            drop = ne.eq("Yes") | pathol.str.contains("Small cell|NE features", case=False, na=False)
            samples.loc[drop, "organ"] = "exclude"
            samples.loc[drop, "organ_rule"] = "1"
        keep_ids = [col for col in expr.columns if col in set(samples["sample_id"])]
        expr = expr[keep_ids]
        samples = samples.loc[samples["sample_id"].isin(keep_ids)].drop_duplicates("sample_id")
        samples["library"] = filename
        made.append(finish_cohort(name, samples, expr, kind, lengths, genes, b0, offsets, indices, paths, platform, filename, refs))
    return made


def main() -> None:
    paths = load_paths()
    root = Path(paths["results"]) / "stage5" / "aux"
    append_decisions(root / "screening_log.tsv")
    downloads = Path(paths["data_processed"]) / "aux" / "downloads"
    downloads.mkdir(parents=True, exist_ok=True)
    gtf = downloads / "gencode.v23.annotation.gtf.gz"
    download_to("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_23/gencode.v23.annotation.gtf.gz", gtf)
    lengths = gene_lengths(gtf, downloads / "gencode_v23_exon_union.tsv")
    print("lengths", len(lengths), flush=True)
    genes, b0_idx, offsets, indices, _names = load_gene_pack(paths)
    b0 = [genes[int(i)] for i in b0_idx]
    print("reference TPM", flush=True)
    refs = load_reference_tpm(paths, genes)
    download_rows = []
    feature_rows = process_cbio(paths, lengths, genes, b0, offsets, indices, download_rows, refs)
    pd.DataFrame(download_rows).to_csv(paths["config"] / "aux_downloads_A5.tsv", sep="\t", index=False)
    pd.DataFrame(feature_rows).to_csv(paths["config"] / "aux_features_A5.tsv", sep="\t", index=False)
    print("cBio cohorts", len(feature_rows), flush=True)


if __name__ == "__main__":
    main()
