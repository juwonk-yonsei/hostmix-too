#!/usr/bin/env python3
"""Download auxiliary clinical metadata and GEO series headers. No expression, no models."""
from __future__ import annotations

import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd

from stage5_lib import load_paths

DATAHUB = "https://raw.githubusercontent.com/cBioPortal/datahub/master/public"
LFS = "https://github.com/cBioPortal/datahub.git/info/lfs/objects/batch"
UA = {"User-Agent": "projA-aux-meta"}

# Studies whose only matched attribute was "Primary Tumor Site" do not have a
# collection-site field. Normal, organoid, and preclinical titles fail criterion 2.
SECOND_PASS = {
    "brain_cptac_gdc": ("제외", "3", "matched attribute display is Primary Tumor Site"),
    "hcc_clca_2024": ("제외", "3", "matched attribute display is Primary Tumor Site"),
    "luad_cptac_gdc": ("제외", "3", "matched attribute display is Primary Tumor Site"),
    "lusc_cptac_gdc": ("제외", "3", "matched attribute display is Primary Tumor Site"),
    "mel_ucla_2016": ("제외", "3", "matched attribute display is Primary Tumor Site"),
    "ohnca_cptac_gdc": ("제외", "3", "matched attribute display is Primary Tumor Site"),
    "pancreas_cptac_gdc": ("제외", "3", "matched attribute display is Primary Tumor Site"),
    "rcc_cptac_gdc": ("제외", "3", "matched attribute display is Primary Tumor Site"),
    "skcm_dfci_2015": ("제외", "3", "matched attribute display is Primary Tumor Site"),
    "uec_cptac_gdc": ("제외", "3", "matched attribute display is Primary Tumor Site"),
    "normal_skin_fibroblast_2024": ("제외", "2", "title says Normal Fibroblasts"),
    "normal_skin_keratinocytes_2024": ("제외", "2", "title says Normal Keratinocytes"),
    "normal_skin_melanocytes_2024": ("제외", "2", "title says Normal Melanocytes"),
    "prad_organoids_msk_2022": ("제외", "2", "study id contains organoids"),
    "pptc_2019": ("제외", "2", "title says Preclinical Testing Consortium"),
}

CLINICAL = [
    "aml_ohsu_2018",
    "asclc_msk_2024",
    "blca_iatlas_imvigor210_2017",
    "brain_cptac_2020",
    "brca_iatlas_anders_2022",
    "gbm_cptac_2021",
    "gbm_iatlas_prins_2019",
    "lgg_ctf_synodos_2025",
    "mel_dfci_2019",
    "mel_iatlas_gide_2019",
    "mel_iatlas_hugo_ucla_2016",
    "mel_iatlas_liu_2019",
    "mel_tsam_liang_2017",
    "nbl_target_2018_pub",
    "nepc_wcm_2016",
    "paad_iatlas_prince_2022",
    "prad_fhcrc",
    "prad_su2c_2019",
    "schw_ctf_synodos_2025",
    "sclc_ucologne_2015",
    "skcm_mskcc_2014",
]

GEO_STARTERS = ["GSE209998", "GSE50760", "GSE41258", "GSE14018", "GSE71729", "GSE74685"]


def get_bytes(url: str, timeout: int = 120, data: bytes | None = None, headers: dict | None = None) -> bytes:
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
    body = json.dumps({
        "operation": "download",
        "transfers": ["basic"],
        "objects": [{"oid": oid, "size": size}],
    }).encode()
    raw = get_bytes(
        LFS, timeout=120, data=body,
        headers={"Accept": "application/vnd.git-lfs+json", "Content-Type": "application/vnd.git-lfs+json"},
    )
    payload = json.loads(raw)
    href = payload["objects"][0]["actions"]["download"]["href"]
    return get_bytes(href, timeout=300)


def datahub_file(study: str, name: str) -> bytes | None:
    url = f"{DATAHUB}/{study}/{name}"
    try:
        raw = get_bytes(url)
    except Exception as exc:
        print("miss", study, name, type(exc).__name__, flush=True)
        return None
    if raw.startswith(b"version https://git-lfs"):
        raw = lfs_download(raw)
    return raw


def attribute_text(study: str) -> str:
    url = f"https://www.cbioportal.org/api/studies/{urllib.parse.quote(study)}/clinical-attributes?projection=SUMMARY"
    request = urllib.request.Request(url, headers={"Accept": "application/json", **UA})
    with urllib.request.urlopen(request, timeout=60) as handle:
        attrs = json.load(handle)
    sample_level = [
        f"{row.get('clinicalAttributeId')}={row.get('displayName')}"
        for row in attrs if row.get("patientAttribute") is not True
    ]
    return "; ".join(sample_level)


def append_second_pass(log_path: Path) -> None:
    frame = pd.read_csv(log_path, sep="\t", dtype=str)
    existing = set(frame["memo"].fillna(""))
    rows = []
    for study, (decision, failed, memo) in SECOND_PASS.items():
        extra = ""
        if failed == "3":
            extra = " | sample attributes: " + attribute_text(study)
            time.sleep(0.05)
        note = "2차 판정: " + memo + extra
        if note in existing:
            continue
        base = frame.loc[frame["id"] == study].iloc[0]
        rows.append({
            "path": base["path"], "query": base["query"], "id": study, "title": base["title"],
            "n": base["n"], "platform": base["platform"], "decision": decision,
            "failed_criterion": failed, "memo": note,
        })
        print("second", study, decision, failed, flush=True)
    if rows:
        frame = pd.concat([frame, pd.DataFrame(rows)], ignore_index=True)
        frame.to_csv(log_path, sep="\t", index=False)
    print("second pass rows", len(rows), flush=True)


def clinical_value_counts(root: Path) -> None:
    out = root / "clinical"
    out.mkdir(parents=True, exist_ok=True)
    summaries = []
    for study in CLINICAL:
        path = out / f"{study}_sample.txt"
        if path.exists() and path.stat().st_size > 200:
            raw = path.read_bytes()
        else:
            raw = datahub_file(study, "data_clinical_sample.txt")
        if raw is None:
            summaries.append({"study": study, "status": "missing"})
            continue
        path = out / f"{study}_sample.txt"
        path.write_bytes(raw)
        text = raw.decode("utf-8", errors="replace").splitlines()
        header_i = next(i for i, line in enumerate(text) if line.startswith("#") is False and "SAMPLE_ID" in line.upper())
        # cBio files have a # display row then a column-name row.
        names = text[header_i].split("\t")
        records = [line.split("\t") for line in text[header_i + 1:] if line and not line.startswith("#")]
        frame = pd.DataFrame(records, columns=names[:len(records[0])] if records else names)
        if len(frame.columns) != len(names):
            frame = pd.DataFrame(records)
            frame.columns = names[:frame.shape[1]]
        site_cols = [col for col in frame.columns if re.search(r"site|biopsy|tissue|organ|sample_type|metast", str(col), re.I)]
        for col in site_cols:
            counts = frame[col].fillna("").value_counts()
            for value, n in counts.items():
                summaries.append({"study": study, "column": col, "value": value, "n": int(n), "n_rows": len(frame)})
        print("clinical", study, len(frame), site_cols, flush=True)
        time.sleep(0.2)
    pd.DataFrame(summaries).to_csv(root / "clinical_value_counts.tsv", sep="\t", index=False)


def stream_geo_header(url: str, limit_bytes: int = 8_000_000) -> tuple[str, int]:
    """Read only the series-matrix header. Stop once the data table starts."""
    import zlib
    request = urllib.request.Request(url, headers=UA)
    decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
    text = ""
    got = 0
    with urllib.request.urlopen(request, timeout=180) as handle:
        while "!series_matrix_table_begin" not in text and got < limit_bytes:
            chunk = handle.read(65536)
            if not chunk:
                break
            got += len(chunk)
            text += decoder.decompress(chunk).decode("utf-8", errors="replace")
    head = []
    for line in text.splitlines():
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_") or line.startswith("!Series_"):
            head.append(line)
    return "\n".join(head), got


def geo_headers(root: Path) -> None:
    rows = []
    dest = root / "geo_headers"
    dest.mkdir(parents=True, exist_ok=True)
    for acc in GEO_STARTERS:
        prefix = acc[:-3] + "nnn"
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/matrix/{acc}_series_matrix.txt.gz"
        print("geo header", acc, flush=True)
        try:
            head, got = stream_geo_header(url)
        except Exception as exc:
            rows.append({"id": acc, "status": type(exc).__name__, "url": url})
            print("geo fail", acc, exc, flush=True)
            continue
        (dest / f"{acc}_header.txt").write_text(head)
        rows.append({"id": acc, "status": "header", "url": url, "bytes_read": got, "header_lines": head.count("\n") + 1})
        time.sleep(0.4)
    pd.DataFrame(rows).to_csv(root / "geo_header_status.tsv", sep="\t", index=False)


def main() -> None:
    paths = load_paths()
    root = Path(paths["results"]) / "stage5" / "aux"
    log_path = root / "screening_log.tsv"
    append_second_pass(log_path)
    clinical_value_counts(root)
    geo_headers(root)


if __name__ == "__main__":
    main()
